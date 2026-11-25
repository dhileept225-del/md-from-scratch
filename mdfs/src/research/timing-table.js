// Research · C.8: the timing table at the end of every run. Each line traced to finish.cpp, the derived numbers recomputed, the neighbour count predicted from FCC shells, and %varavg decoded from a real 2-process run.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

const base = (skin, nm) => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 3.0 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        ${skin} bin
neigh_modify    ${nm}
fix             1 all nve
`;
const MELT = 'every 20 delay 0 check no', SAFE = 'every 1 delay 0 check yes';

// measured natively (LAMMPS 22 Jul 2025, mpirun, 2 cores), melt, run 250
const NATIVE = {
  1: { loop: 0.334123, rows: [['Pair', 0.26807, 0.26807, 0.26807, 0.0, 80.23], ['Neigh', 0.054413, 0.054413, 0.054413, 0.0, 16.29], ['Comm', 0.0051831, 0.0051831, 0.0051831, 0.0, 1.55], ['Output', 6.1717e-05, 6.1717e-05, 6.1717e-05, 0.0, 0.02], ['Modify', 0.0051297, 0.0051297, 0.0051297, 0.0, 1.54]], other: 0.001263 },
  2: { loop: 0.181168, rows: [['Pair', 0.13253, 0.1326, 0.13268, 0.0, 73.19], ['Neigh', 0.02787, 0.028142, 0.028414, 0.2, 15.53], ['Comm', 0.015851, 0.016181, 0.016511, 0.3, 8.93], ['Output', 4.1995e-05, 0.00010224, 0.00016249, 0.0, 0.06], ['Modify', 0.0022632, 0.0023001, 0.0023371, 0.1, 1.27]], other: 0.00184 },
};

function parseTiming(lines) {
  const t = lines.join('\n'), g = re => { const m = t.match(re); return m ? m : null; };
  const lp = g(/Loop time of ([\d.e+-]+) on (\d+) procs for (\d+) steps with (\d+) atoms/);
  const pf = g(/Performance: ([\d.]+) tau\/day, ([\d.]+) timesteps\/s, ([\d.]+) (\w*atom-step\/s)/);
  const rows = [];
  for (const l of lines) { const m = l.match(/^(\w+)\s*\|\s*([\d.e+-]*)\s*\|\s*([\d.e+-]+)\s*\|\s*([\d.e+-]*)\s*\|\s*([\d.]*)\s*\|\s*([\d.]+)\s*$/); if (m) rows.push({ name: m[1], min: m[2], avg: +m[3], max: m[4], va: m[5], pct: +m[6] }); }
  const nb = g(/Total # of neighbors = (\d+)/), an = g(/Ave neighs\/atom = ([\d.]+)/), bu = g(/Neighbor list builds = (\d+)/);
  return lp && { loop: +lp[1], loopS: lp[1], procs: +lp[2], steps: +lp[3], atoms: +lp[4], perf: pf && { tau: +pf[1], ts: +pf[2], as: +pf[3], unit: pf[4] }, rows, neigh: nb && +nb[1], ave: an && an[1], builds: bu && +bu[1] };
}
// neighbours of one FCC site within rc
function shellCount(rc) {
  const a = Math.cbrt(4 / 0.8442), B = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]]; let n = 0;
  for (let i = -4; i <= 4; i++) for (let j = -4; j <= 4; j++) for (let k = -4; k <= 4; k++) for (const b of B) { const d = a * Math.hypot(i + b[0], j + b[1], k + b[2]); if (d > 0 && d < rc) n++; }
  return n;
}
const close = (h, l, rel) => Math.abs(h - l) <= rel * Math.abs(l) + 1e-12;

export default {
  id: 'timing-table',
  title: 'The timing table',
  minutes: 15,
  summary: 'Every run ends with a block of timings and counts that most people skip. You trace each line to finish.cpp, recompute the Performance line and the %total column from the loop time, predict the total number of neighbours from FCC shells before running, watch the skin trade Pair time against Neigh time, and decode %varavg from a real 2-process run: it is not the percentage spread it looks like.',
  learn: [
    'Timer::stamp: how the time is split into Pair, Neigh, Comm, Output, Modify and Other',
    'Performance: tau/day, timesteps/s and atom-step/s from the loop time',
    'Total # of neighbors = N × (neighbours within rc + skin) / 2, by hand',
    '%varavg = √(variance / mean) × 100, from min, avg and max',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a stopwatch with categories</h2>
    <div class="prose"><p>The loop calls <code>timer->stamp()</code> before a piece of work and <code>timer->stamp(Timer::PAIR)</code> after it (${gh('src/verlet.cpp', 308)}). ${gh('src/timer.cpp', 53)} adds the wall time since the previous stamp to that category and to ALL:</p>
<pre class="lmp-mini">wall_array[which] += delta_wall;     // 65
wall_array[ALL]   += delta_wall;     // 67</pre>
<p>At the end, ${gh('src/finish.cpp', 117)} prints the loop time, and for each category ${gh('src/finish.cpp', 600)}:</p>
<pre class="lmp-mini">%total  = time / time_loop * 100                     // 631
%varavg = sqrt(⟨t²⟩/⟨t⟩ − ⟨t⟩) * 100   (if ⟨t⟩ &gt; 0.001)   // 626, ⟨⟩ over MPI ranks
Other   = TOTAL − ALL                                 // 201: time no stamp claimed
t_step = time_loop / nsteps;  tau/day = 86400 / t_step * dt     // 142, 158</pre></div>
  </section>

  <section class="block">
    <h2>2 · Match: the numbers derived from the loop time</h2>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run melt, 250 steps</button><span class="lmp-status muted tt-sta"></span></div>
    <pre class="lmp-mini tt-raw" hidden></pre>
    <div class="calc tt-hand" hidden></div>
  </section>

  <section class="block">
    <h2>3 · Match: count the neighbours before LAMMPS does</h2>
    <div class="prose"><p>At step 0 the atoms sit on the FCC lattice (a = 1.6796). The list holds every pair closer than rc + skin = 2.5 + skin, each pair once (a half list, page A.7). So Total # of neighbors = 4000 × n / 2, with n the FCC neighbours within 2.5 + skin.</p></div>
    <div class="row"><label class="lc-in">skin<select class="lc-u tt-sk"><option>0.1</option><option selected>0.3</option><option>0.5</option><option>0.7</option><option>1.0</option><option>1.5</option></select></label><button class="btn sm primary" data-k="b">▶ run 0</button><span class="lmp-status muted tt-stb"></span></div>
    <div class="calc tt-pred"></div>
    <div class="prose tt-outb"></div>
  </section>

  <section class="block">
    <h2>4 · Modify: the skin trades Pair against Neigh</h2>
    <div class="prose"><p>Same melt, safe lists (<code>every 1 delay 0 check yes</code>), four skins. A small skin means short lists (cheap Pair) but frequent rebuilds (expensive Neigh); a big skin the opposite.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run 4 × 250 steps</button><span class="lmp-status muted tt-stc"></span></div>
    <div class="tbl tt-outc"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: two processes, and what %varavg means</h2>
    <div class="prose"><p>The browser runs on one process, so every min = avg = max and %varavg = 0. Here is the same melt run natively with <code>mpirun -np 1</code> and <code>-np 2</code>:</p></div>
    <div class="row"><label class="lc-in">processes<select class="lc-u tt-np"><option>1</option><option selected>2</option></select></label></div>
    <div class="tbl tt-nat"></div>
    <div class="prose"><p>With two ranks, ⟨t⟩ = (min + max)/2 and the variance is ((max − min)/2)², so LAMMPS's formula becomes <b>%varavg = (max − min) / 2 / √avg × 100</b>. That has units of √seconds: it is not the relative spread, which for Neigh is (max − min)/2/avg = 0.97 %, not 0.2. Use min and max columns to judge load imbalance; Comm growing from 1.6 % to 8.9 % is the real cost of splitting the box.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> 32 000 LJ atoms, dt = 0.005, 100 000 steps in a loop time of 10 s. What does the Performance line say?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>t_step = 10⁻⁴ s; tau/day = 86400 / 10⁻⁴ × 0.005 = <b>4 320 000</b>; timesteps/s = <b>10 000</b>; atom-step/s = 3.2 × 10⁸ → <b>320 Matom-step/s</b>.</div></div>
    <div class="ex"><p><b>2.</b> In the browser the "CPU use" line says 0.0 %. Why?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The loop's CPU time is the difference of two <code>platform::cputime()</code> calls at the start and end of the run (${gh('src/timer.cpp', 105)} and 124), which asks the operating system through <code>getrusage</code> (${gh('src/platform.cpp', 156)}). WebAssembly in a browser has no such clock, so it returns 0 both times: 0 / time_loop × 100 = 0.</div></div>
    <div class="ex"><p><b>3.</b> A run shows Neigh at 40 % with <code>every 1 delay 0 check yes</code>. Two things to try?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A larger skin (fewer rebuilds; section 4), and <code>every 2</code> or more so it does not even check every step. Then confirm "Dangerous builds = 0".</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.tt-sta'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(`${base(0.3, MELT)}run             250`); b.disabled = false;
      const T = !res.error && parseTiming(res.lines); if (!T) return err(st, res); st.textContent = '';
      const i0 = res.lines.findIndex(l => /^Loop time/.test(l)), i1 = res.lines.findIndex(l => /^Dangerous/.test(l));
      $('.tt-raw').hidden = false; $('.tt-raw').textContent = res.lines.slice(i0, i1 + 1).join('\n');
      const ts = T.loop / T.steps, tau = 86400 / ts * 0.005, stp = 1 / ts, as = T.atoms * stp / 1e6;
      const sum = T.rows.filter(r => r.name !== 'Other').reduce((s, r) => s + r.avg, 0), oth = T.rows.find(r => r.name === 'Other');
      const L = [`from "Loop time of ${T.loopS}" for ${T.steps} steps, ${T.atoms} atoms:`,
        `  t_step = ${T.loopS} / ${T.steps} = ${ts.toPrecision(6)} s`,
        `  tau/day = 86400 / t_step × 0.005 = ${tau.toFixed(3)}    LAMMPS ${T.perf.tau}  ${close(tau, T.perf.tau, 1e-5) ? '✓' : '✗'}`,
        `  timesteps/s = 1 / t_step = ${stp.toFixed(3)}    LAMMPS ${T.perf.ts}  ${close(stp, T.perf.ts, 1e-5) ? '✓' : '✗'}`,
        `  atom-step/s = 4000 / t_step = ${as.toFixed(3)} M    LAMMPS ${T.perf.as} ${T.perf.unit}  ${close(as, T.perf.as, 1e-3) ? '✓' : '✗'}`, '',
        ...T.rows.map(r => `  ${r.name.padEnd(7)} ${String(r.avg).padEnd(11)} / ${T.loopS} × 100 = ${(r.avg / T.loop * 100).toFixed(2).padStart(6)}   LAMMPS ${r.pct.toFixed(2).padStart(6)}  ${Math.abs(r.avg / T.loop * 100 - r.pct) <= 0.011 ? '✓' : '✗'}`), '',
        `  Other = loop − (sum of the rest) = ${T.loopS} − ${sum.toPrecision(6)} = ${(T.loop - sum).toPrecision(4)}    LAMMPS ${oth ? oth.avg : '?'}  ${oth && Math.abs(T.loop - sum - oth.avg) < 2e-5 ? '✓' : '✗'}`,
        `  (the printed times have 5 digits, so the hand values agree to that precision)`];
      const coarse = T.rows.every(r => Math.abs(r.avg * 1e4 - Math.round(r.avg * 1e4)) < 1e-6);
      if (coarse) L.push('', '  Every time above is a whole number of 0.1 ms: your browser rounds its clock (performance.now)', '  to blunt timing attacks, so each stamp is rounded too. A native build times to the microsecond.');
      $('.tt-hand').hidden = false; $('.tt-hand').innerHTML = L.map(esc).join('\n').replace(/✓/g, OK).replace(/✗/g, BAD);
    };

    const pred = () => {
      const sk = +$('.tt-sk').value, rc = 2.5 + sk, n = shellCount(rc), a = Math.cbrt(4 / 0.8442);
      const sh = []; const seen = new Map();
      const B = [[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]];
      for (let i = -4; i <= 4; i++) for (let j = -4; j <= 4; j++) for (let k = -4; k <= 4; k++) for (const b of B) { const q = (i + b[0]) ** 2 + (j + b[1]) ** 2 + (k + b[2]) ** 2; if (q > 0 && a * Math.sqrt(q) < rc) seen.set(q, (seen.get(q) || 0) + 1); }
      [...seen.keys()].sort((x, y) => x - y).forEach(q => sh.push(`${(a * Math.sqrt(q)).toFixed(3)}×${seen.get(q)}`));
      $('.tt-pred').textContent = `rc + skin = ${rc}\nshells inside: ${sh.join('  ')}\nn = ${n} neighbours per atom\nTotal # of neighbors = 4000 × ${n} / 2 = ${4000 * n / 2}\nAve neighs/atom = ${n / 2}`;
      $('.tt-outb').innerHTML = ''; return { n, total: 4000 * n / 2 };
    };
    $('.tt-sk').onchange = pred; pred();
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.tt-stb'), sk = +$('.tt-sk').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(`${base(sk, MELT)}run             0`); b.disabled = false;
      const T = !res.error && parseTiming(res.lines); if (!T || T.neigh == null) return err(st, res); st.textContent = '';
      const P = pred();
      $('.tt-outb').innerHTML = `<p>LAMMPS: Total # of neighbors = <b>${T.neigh}</b> ${T.neigh === P.total ? OK : BAD}, Ave neighs/atom = ${T.ave} ${+T.ave === P.n / 2 ? OK : BAD}, Neighbor list builds = ${T.builds} (the build in setup is not counted).</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.tt-stc'); b.disabled = true; const rs = [];
      for (const sk of [0.1, 0.3, 0.6, 1.0]) {
        st.textContent = `Running skin ${sk}…`;
        const res = await runLammps(`${base(sk, SAFE)}run             250`);
        const T = !res.error && parseTiming(res.lines); if (!T) { b.disabled = false; return err(st, res); }
        const g = n => T.rows.find(r => r.name === n) || { avg: 0, pct: 0 };
        rs.push(row([sk, T.builds, (T.neigh / 4000).toFixed(1), `${g('Pair').avg.toFixed(3)} (${g('Pair').pct.toFixed(1)} %)`, `${g('Neigh').avg.toFixed(3)} (${g('Neigh').pct.toFixed(1)} %)`, T.loop.toFixed(3)]));
        $('.tt-outc').innerHTML = `<table>${head(['skin', 'builds', 'neighbours per atom', 'Pair s', 'Neigh s', 'loop s'])}${rs.join('')}</table>`;
      }
      b.disabled = false; st.textContent = 'Times depend on your device; the builds and neighbour counts do not.';
    };

    const nat = () => {
      const D = NATIVE[$('.tt-np').value];
      $('.tt-nat').innerHTML = `<table>${head(['section', 'min', 'avg', 'max', '%varavg', 'by hand: (max−min)/2/√avg×100', '%total'])}${D.rows.map(([n, mi, av, ma, va, pc]) => { const h = av > 0.001 ? (ma - mi) / 2 / Math.sqrt(av) * 100 : 0; return row([n, mi, av, ma, va.toFixed(1), `${h.toFixed(3)} ${Math.abs(h - va) <= 0.05 + 1e-9 ? OK : BAD}${av <= 0.001 ? ' (avg ≤ 0.001: set to 0)' : ''}`, pc.toFixed(2)]); }).join('')}${row(['Other', '', D.other, '', '', '', (D.other / D.loop * 100).toFixed(2)])}</table><p class="muted">Loop time ${D.loop} s.</p>`;
    };
    $('.tt-np').onchange = nat; nat();
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
