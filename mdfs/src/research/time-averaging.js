// Research · D.16: time averaging and correlations. fix ave/time (one, running, window) and fix ave/correlate (auto, running)
// reproduced bit for bit from per-step thermo output; the trap() function; a Green-Kubo viscosity and how noisy it is.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const SETUP = `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 4 0 4 0 4
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
run             1000
reset_timestep  0`;
const NR = 50;
const printRows = (fid, n) => Array.from({ length: n }, (_, k) => `print "C ${k + 1} $(f_${fid}[${k + 1}][1]) $(f_${fid}[${k + 1}][2]) $(f_${fid}[${k + 1}][3]:%.17g) $(f_${fid}[${k + 1}][4]:%.17g) $(f_${fid}[${k + 1}][5]:%.17g)"`).join('\n');
const SHOWN = `fix             2 all ave/time 2 5 20 c_thermo_temp c_thermo_press file one.out format " %.17g"
fix             4 all ave/time 2 5 20 c_thermo_temp c_thermo_press ave running file run.out format " %.17g"
fix             5 all ave/time 2 5 20 c_thermo_temp c_thermo_press ave window 3 file win.out format " %.17g"
fix             3 all ave/correlate 1 ${NR} 200 c_thermo_press[4] c_thermo_press[5] c_thermo_press[6] type auto ave running
thermo_style    custom step temp press pxy pxz pyz
thermo_modify   format float %.17g
thermo          1
run             1000
print           "C 1 $(f_3[1][1]) $(f_3[1][2]) $(f_3[1][3]:%.17g) …"   # one line per row, 50 rows`;
const MAIN = `${SETUP}\n${SHOWN.split('\nprint')[0]}\n${printRows(3, NR)}`;
const GKN = 100;
const GKSHOWN = `${SETUP}
fix             3 all ave/correlate 5 ${GKN} 20000 c_thermo_press[4] c_thermo_press[5] c_thermo_press[6] type auto ave running
variable        T equal temp
fix             6 all ave/time 100 200 20000 v_T
thermo          5000
run             20000
print           "TRAP $(trap(f_3[3]):%.17g) $(trap(f_3[4]):%.17g) $(trap(f_3[5]):%.17g) T $(f_6:%.17g) V $(vol:%.17g)"`;
const GK = `${GKSHOWN}\n${printRows(3, GKN)}`;

const aveFile = t => t.split('\n').filter(l => l.trim() && !l.startsWith('#')).map(l => l.trim().split(/\s+/).map(Number));
const corrRows = lines => lines.filter(l => l.startsWith('C ')).map(l => l.split(/\s+/).slice(1).map(Number));

// FixAveTime::invoke_scalar, two values, nevery 2, nrepeat 5, nfreq 20
function aveTime(th, mode, nwindow = 3) {
  const out = [], hist = []; let tot = [0.0, 0.0], norm = 0;
  for (let nf = 20; nf <= 1000; nf += 20) {
    const v = [0.0, 0.0];
    for (let s = nf - 8; s <= nf; s += 2) { v[0] += th[s][0]; v[1] += th[s][1]; }
    v[0] /= 5.0; v[1] /= 5.0;
    if (mode === 'one') { tot = v.slice(); norm = 1; }
    else if (mode === 'running') { tot = [tot[0] + v[0], tot[1] + v[1]]; norm++; }
    else { tot = [tot[0] + v[0], tot[1] + v[1]]; if (hist.length >= nwindow) { const o = hist[hist.length - nwindow]; tot = [tot[0] - o[0], tot[1] - o[1]]; } hist.push(v); norm = Math.min(hist.length, nwindow); }
    out.push([nf, tot[0] / norm, tot[1] / norm]);
  }
  return out;
}
// FixAveCorrelate::end_of_step + accumulate, type auto, ave running, prefactor 1
function aveCorrelate(vals, nrepeat) {
  const corr = Array.from({ length: nrepeat }, () => [0.0, 0.0, 0.0]), count = new Array(nrepeat).fill(0), buf = [];
  for (const v of vals) {
    buf.push(v); const ns = Math.min(buf.length, nrepeat);
    for (let k = 0; k < ns; k++) count[k]++;
    for (let k = 0; k < ns; k++) for (let c = 0; c < 3; c++) corr[k][c] += buf[buf.length - 1 - k][c] * buf[buf.length - 1][c];
  }
  return corr.map((c, k) => ({ count: count[k], c: c.map(x => 1.0 * x / count[k]) }));
}
// Variable special function trap(): sum everything, then remove half of the two ends
const trap = col => { let v = 0.0; for (const x of col) v += x; v -= 0.5 * col[0] + 0.5 * col[col.length - 1]; return v; };

export default {
  id: 'time-averaging',
  title: 'Time averaging and correlations',
  minutes: 20,
  summary: 'fix ave/time and fix ave/correlate turn per-step numbers into averages and time correlation functions. You work out by hand which steps each one samples, reproduce the three averaging modes and a 50-lag autocorrelation of the shear stress bit for bit from per-step thermo output, check LAMMPS\'s trap() function, and then use it for a Green–Kubo viscosity: the three stress components agree within a few per cent when the integral is stopped where the correlation has died, and disagree by a factor of two when it runs on through the noisy tail.',
  learn: [
    'nevery, nrepeat, nfreq: exactly which steps are sampled',
    'ave one, running and window: what is summed and what is divided',
    'fix ave/correlate: the circular buffer, counts per lag, ave running',
    'Green–Kubo with trap(), and where to stop the integral',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: samples, blocks and lags</h2>
    <div class="prose"><p><code>fix ave/time nevery nrepeat nfreq</code> averages <i>nrepeat</i> samples taken <i>nevery</i> steps apart, ending on every multiple of <i>nfreq</i> (${gh('src/fix_ave_time.cpp', 1129)}). For <code>2 5 20</code>: steps 12, 14, 16, 18, 20 make the block written at step 20, then 32 … 40, and so on. nfreq must be a multiple of nevery and nrepeat·nevery ≤ nfreq (line 131), so blocks never overlap.</p>
<pre class="lmp-mini">vector[i] += value;                         // each sample (invoke_scalar)
… after nrepeat samples:  vector[i] /= nrepeat;
ave one:     total = vector;               norm = 1
ave running: total += vector;              norm++
ave window M: total += vector − vector from M blocks ago;  norm = min(blocks, M)
output = total / norm</pre>
<p><code>fix ave/correlate nevery nrepeat nfreq</code> (${gh('src/fix_ave_correlate.cpp', 389)}) keeps the last <i>nrepeat</i> samples in a circular buffer. Every nevery steps it adds the new sample times each stored one into lag bins 0 … nrepeat−1 and counts how many products each lag has; every nfreq steps it publishes corr/count as a global array, row k = lag k·nevery, columns: lag, count, then one column per pair. <code>ave running</code> keeps accumulating over the whole run; <code>ave one</code> starts again after each output. The file it writes uses %g (6 digits), so below the array is read through variables with 17 digits.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: three averages and a correlation, by hand</h2>
    <div class="prose"><p>256 LJ atoms, melted for 1000 steps, then 1000 more with thermo every step (17 digits). Everything the fixes see is in that thermo output: thermo and the fixes ask the same compute (<code>thermo_press</code>), which runs once per step.</p></div>
    <pre class="lmp-mini">${esc(SETUP + '\n' + SHOWN)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run and compare</button><span class="lmp-status muted ta-sta"></span></div>
    <div class="tbl ta-a"></div>
    <div class="prose ta-asum"></div>
    <div class="tbl ta-c"></div>
    <div class="prose ta-csum"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: a Green–Kubo viscosity</h2>
    <div class="prose"><p>The shear viscosity is the time integral of the stress autocorrelation, η = V/(k<sub>B</sub>T) ∫₀^∞ ⟨P<sub>xy</sub>(0) P<sub>xy</sub>(t)⟩ dt; xz and yz give two more estimates. 20 000 steps (100 τ), samples every 5 steps, 100 lags (2.5 τ). LAMMPS's <code>trap()</code> does the integral in units of the sample spacing; multiply by 5·dt.</p></div>
    <pre class="lmp-mini">${esc(GKSHOWN)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run 20 000 steps</button><span class="lmp-status muted ta-stb"></span></div>
    <div class="tbl ta-b"></div>
    <div class="prose ta-bsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> <code>fix ave/time 10 5 40</code> defined at step 0: which steps go into the first output?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>nextvalid: 40 − (5 − 1)·10 = 0. Steps 0, 10, 20, 30, 40, written at step 40. Step 0 is in it: the initial configuration counts as a sample.</div></div>
    <div class="ex"><p><b>2.</b> And <code>fix ave/time 10 5 30</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>An error: 5 × 10 = 50 &gt; 30 ("Illegal fix ave/time nfreq value", line 131). Blocks would overlap.</div></div>
    <div class="ex"><p><b>3.</b> In section 2 the lag-0 count is 1001 and lag 49 has 952. Why?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>nevery = 1 from step 0 to 1000 is 1001 samples, each a lag-0 product. A lag-49 product needs a sample 49 steps earlier: 1001 − 49 = 952.</div></div>
    <div class="ex"><p><b>4.</b> With <code>ave one</code> instead of running, what are the counts in the block written at step 400?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>After each output LAMMPS clears the sums, sets nsample = 1 and calls accumulate() again, so the step-200 sample starts the new block. Steps 200 … 400 are 201 samples: count 201 − k at lag k. Neighbouring blocks share one sample.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.ta-sta'); b.disabled = true; st.textContent = 'Running 2000 steps…';
      const res = await runLammps(MAIN, { read: ['one.out', 'run.out', 'win.out'] }); b.disabled = false;
      if (res.error || !res.files['one.out']) return err(st, res.error); st.textContent = '';
      const th = {}, i0 = res.lines.findIndex(l => /^\s*Step\s+Temp\s+Press\s+Pxy/.test(l));
      for (let k = i0 + 1; k < res.lines.length; k++) { const w = res.lines[k].trim().split(/\s+/); if (w.length !== 6 || !/^\d+$/.test(w[0])) break; th[+w[0]] = w.slice(1).map(Number); }
      const rows = []; let tot = 0, n = 0;
      for (const [mode, f] of [['one', 'one.out'], ['running', 'run.out'], ['window', 'win.out']]) {
        const L = aveFile(res.files[f]), H = aveTime(th, mode); let ok = 0;
        L.forEach((r, k) => { if (r[0] === H[k][0] && r[1] === H[k][1] && r[2] === H[k][2]) ok++; }); tot += ok; n += L.length;
        for (const k of [0, 1, 2, 3, 49]) rows.push(row([k === 0 ? `<b>${mode}</b>` : '', L[k][0], `<span class="mono">${s17(L[k][1])}</span>`, `<span class="mono">${s17(L[k][2])}</span>`, L[k][1] === H[k][1] && L[k][2] === H[k][2] ? OK : BAD]));
      }
      $('.ta-a').innerHTML = `<table>${head(['ave', 'step', 'temp', 'press', 'hand'])}${rows.join('')}</table>`;
      const b0 = [12, 14, 16, 18, 20].map(s => th[s][0]);
      $('.ta-asum').innerHTML = `<p>${tot} of ${n} output lines identical to the last bit. The first block by hand: T at steps 12 … 20 = ${b0.map(x => x.toFixed(5)).join(', ')}, summed in that order and divided by 5.0 → ${s17(b0.reduce((s, x) => s + x, 0) / 5)}. The window average is a running sum with the oldest block subtracted, so after 50 blocks it carries the round-off of 50 additions and 47 subtractions; the hand code does the same operations in the same order and matches anyway.</p>`;
      const C = corrRows(res.lines), vals = []; for (let s = 0; s <= 1000; s++) vals.push(th[s].slice(2, 5));
      const H = aveCorrelate(vals, NR); let ok = 0;
      C.forEach((r, k) => { if (r[1] === k && r[2] === H[k].count && r[3] === H[k].c[0] && r[4] === H[k].c[1] && r[5] === H[k].c[2]) ok++; });
      $('.ta-c').innerHTML = `<table>${head(['lag (steps)', 'count', '⟨Pxy Pxy⟩', '⟨Pxz Pxz⟩', '⟨Pyz Pyz⟩', 'hand'])}${[0, 1, 2, 5, 10, 20, 30, 40, 49].map(k => { const r = C[k]; return row([r[1], r[2], `<span class="mono">${s17(r[3])}</span>`, `<span class="mono">${s17(r[4])}</span>`, `<span class="mono">${s17(r[5])}</span>`, r[2] === H[k].count && r[3] === H[k].c[0] && r[4] === H[k].c[1] && r[5] === H[k].c[2] ? OK : BAD]); }).join('')}</table>`;
      $('.ta-csum').innerHTML = `<p>${ok} of ${NR} rows identical: lag, count and the three correlations, for the array as it stands after step 1000. The products are added newest-times-older into lag k (accumulate, line 545), and each is divided by its own count. Over 49 steps (0.25 τ) the stress correlation falls from ${C[0][3].toFixed(4)} to ${C[49][3].toFixed(4)}: it has not decayed yet, and 50 lags of one step are far too short for a viscosity.</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.ta-stb'); b.disabled = true; st.textContent = 'Running 20 000 steps (a few seconds)…';
      const res = await runLammps(GK); b.disabled = false;
      if (res.error) return err(st, res.error); st.textContent = '';
      const C = corrRows(res.lines), w = res.lines.find(l => l.startsWith('TRAP')).split(/\s+/), lt = [+w[1], +w[2], +w[3]], T = +w[5], V = +w[7];
      const dt = 0.005 * 5, eta = [], mid = [], rows = [], names = ['xy', 'xz', 'yz'];
      for (let c = 0; c < 3; c++) {
        const col = C.map(r => r[3 + c]), ht = trap(col), e3 = V / T * dt * ht; eta.push(e3);
        const run = k => { const cc = col.slice(0, k + 1); return V / T * dt * trap(cc); };
        mid.push(run(50)); rows.push(row([names[c], col[0].toFixed(5), (col[20] / col[0]).toFixed(3), (col[99] / col[0]).toFixed(3), run(20).toFixed(2), run(50).toFixed(2), e3.toFixed(2), ht === lt[c] ? OK : BAD]));
      }
      const stat = a => { const m = (a[0] + a[1] + a[2]) / 3; return [m, Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / 2)]; }, [m, sd] = stat(eta), [m2, sd2] = stat(mid);
      $('.ta-b').innerHTML = `<table>${head(['P', 'C(0)', 'C(0.5τ)/C(0)', 'C(2.5τ)/C(0)', 'η to 0.5τ', 'η to 1.25τ', 'η to 2.5τ', 'trap() = hand'])}${rows.join('')}</table>`;
      $('.ta-bsum').innerHTML = `<p>T = ${T.toFixed(4)} (fix ave/time of temp over the run), V = ${V.toFixed(3)}σ³. trap() is sum minus half of each end (variable.cpp:4777); the hand sum agrees to the last bit for all three columns. The correlation drops to a few per cent of C(0) within 0.5τ, and the long tail is noise of either sign, so the running integral wanders instead of settling.</p><p>The three components are equally valid measurements of the same number, so their spread is an honest error bar. Integrated to 2.5τ: ${eta.map(x => x.toFixed(2)).join(', ')} (mean ${m.toFixed(2)} ± ${sd.toFixed(2)}). Stopped at 1.25τ, where the correlation has died: ${mid.map(x => x.toFixed(2)).join(', ')} (mean ${m2.toFixed(2)} ± ${sd2.toFixed(2)}). Past the decay the integral only collects noise, and the noise grows with the length of the tail. In practice: run much longer (the error falls only as 1/√time), average the off-diagonal components, plot the running integral and read it where it levels off, not at the last lag.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
