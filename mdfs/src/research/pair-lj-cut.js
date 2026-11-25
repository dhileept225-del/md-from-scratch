// Research · D.5: pair style lj/cut as a whole class. settings, coeff, init_one, compute, single; a pair_write table matched to 15 digits; a re-issued pair_style that resets explicit cutoffs; the cost of the cutoff jump in NVE.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

// C's %.15g
function g15(v) {
  if (v === 0) return '0';
  const ex = v.toExponential(14), e = +ex.split('e')[1];
  if (e < -4 || e >= 15) { let [m, x] = ex.split('e'); m = m.replace(/\.?0+$/, ''); const xs = Math.abs(+x); return `${m}e${+x < 0 ? '-' : '+'}${xs < 10 ? '0' : ''}${xs}`; }
  return v.toFixed(14 - e).replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
}
// pow(x, n) correctly rounded (what glibc returns; JavaScript's Math.pow can be 1 ulp off)
function powr(x, n) {
  if (x === 0) return 0;
  const dv = new DataView(new ArrayBuffer(8)); dv.setFloat64(0, x);
  const hi = dv.getUint32(0), lo = dv.getUint32(4), ex = (hi >>> 20) & 0x7ff;
  let m = (BigInt(hi & 0xfffff) << 32n) | BigInt(lo), e;
  if (ex === 0) e = -1074; else { m |= 1n << 52n; e = ex - 1075; }
  const big = m ** BigInt(n), bits = big.toString(2).length, sh = Math.max(0, bits - 60);
  // keep 60 bits plus a sticky bit, so Number() rounds exactly once
  let top = big >> BigInt(sh); if (sh && (big & ((1n << BigInt(sh)) - 1n))) top |= 1n;
  let v = Number(top), k = e * n + sh;
  while (k > 0) { const s = Math.min(k, 1000); v *= 2 ** s; k -= s; }
  while (k < 0) { const s = Math.min(-k, 1000); v /= 2 ** s; k += s; }
  return x < 0 && n % 2 ? -v : v;
}
// PairLJCut::init_one and single(), in the C++ order
function coeffs(eps, sig, cut, shift, pw = powr) {
  const lj1 = 48.0 * eps * pw(sig, 12), lj2 = 24.0 * eps * pw(sig, 6), lj3 = 4.0 * eps * pw(sig, 12), lj4 = 4.0 * eps * pw(sig, 6);
  let offset = 0.0; if (shift && cut > 0.0) { const ratio = sig / cut; offset = 4.0 * eps * (pw(ratio, 12) - pw(ratio, 6)); }
  return { lj1, lj2, lj3, lj4, offset, cutsq: cut * cut };
}
function single(c, rsq) {
  const r2inv = 1.0 / rsq, r6inv = r2inv * r2inv * r2inv, forcelj = r6inv * (c.lj1 * r6inv - c.lj2), fforce = 1.0 * forcelj * r2inv;
  const philj = r6inv * (c.lj3 * r6inv - c.lj4) - c.offset;
  return { e: 1.0 * philj, fforce };
}
// Pair::write_file, style r
function table(c, n, inner, outer) {
  const out = [];
  for (let i = 0; i < n; i++) { const r = inner + (outer - inner) * i / (n - 1), rsq = r * r; let e = 0, f = 0; if (rsq < c.cutsq) { const s = single(c, rsq); e = s.e; f = s.fforce * r; } out.push([i + 1, r, e, f]); }
  return out;
}
const parseTable = (txt, key) => { const ls = txt.split('\n'), i = ls.findIndex(l => l.trim() === key); const o = []; for (let k = i + 3; k < ls.length; k++) { const w = ls[k].trim().split(/\s+/); if (w.length !== 4 || !/^\d+$/.test(w[0])) break; o.push(w); } return o; };

const scriptA = p => `units           real
atom_style      atomic
region          box block 0 30 0 30 0 30 units box
create_box      1 box
mass            1 39.948
pair_style      lj/cut ${p.cut}
pair_coeff      1 1 ${p.eps} ${p.sig}
${p.shift ? 'pair_modify     shift yes\n' : ''}pair_write      1 1 ${p.n} r ${p.lo} ${p.hi} w.table LJ`;
const SCRIPT_B = `units           real
atom_style      atomic
region          box block 0 30 0 30 0 30 units box
create_box      2 box
mass            * 39.948
pair_style      lj/cut 8.5
pair_coeff      1 1 0.2381 3.405
pair_coeff      2 2 0.2 3.0
pair_coeff      1 2 0.25 3.2 10.0
pair_modify     shift yes
pair_write      1 2 12 r 3.0 11.0 w.table BEFORE
pair_style      lj/cut 8.5
pair_write      1 2 12 r 3.0 11.0 w.table AFTER`;
const scriptC = rc => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 3.0 87287 loop geom
pair_style      lj/cut ${rc}
pair_coeff      1 1 1.0 1.0 ${rc}
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
fix             1 all nve
thermo_style    custom step temp etotal
thermo_modify   norm yes format float %.10g
thermo          50
run             1000`;

export default {
  id: 'pair-lj-cut',
  title: 'Pair style lj/cut',
  minutes: 20,
  summary: 'pair_style lj/cut is about 700 lines, and five functions do the work: settings and coeff read your input, init_one turns ε and σ into four numbers and an offset, compute runs every step, single answers "what is E and F at this r?" for other commands. You reproduce a pair_write table to all 15 printed digits, discover that typing pair_style a second time silently resets cutoffs you set by hand, and measure what the jump at the cutoff costs in energy conservation.',
  learn: [
    'The life of a pair style: settings, coeff, init_one, compute, single',
    'lj1 … lj4 and the offset, as init_one computes them',
    'pair_write reproduced by hand to 15 digits',
    'The cutoff is a step in E and F, and what it does to NVE',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: five functions</h2>
    <div class="tbl"><table>${head(['when', 'function', 'what it does'])}
      ${row(['<code>pair_style lj/cut 8.5</code>', gh('src/pair_lj_cut.cpp', 422), 'stores the global cutoff; <b>resets every cutoff already set</b> (line 430)'])}
      ${row(['<code>pair_coeff 1 2 ε σ [rc]</code>', gh('src/pair_lj_cut.cpp', 442), 'fills ε, σ, rc for the type pairs, marks setflag'])}
      ${row(['each <code>run</code>, every i ≤ j', gh('src/pair_lj_cut.cpp', 501), 'mixes unset pairs (page D.6), lj1…lj4, offset, tail; returns rc → cutsq'])}
      ${row(['every step', gh('src/pair_lj_cut.cpp', 69), 'the neighbour loop (page C.4)'])}
      ${row(['pair_write, compute pair/local, fix bond/swap …', gh('src/pair_lj_cut.cpp', 660), 'E and F/r for one pair at one rsq'])}</table></div>
    <div class="prose"><pre class="lmp-mini">lj1 = 48.0 * epsilon * pow(sigma,12.0);      // 509
lj2 = 24.0 * epsilon * pow(sigma,6.0);
lj3 =  4.0 * epsilon * pow(sigma,12.0);
lj4 =  4.0 * epsilon * pow(sigma,6.0);
offset = 4.0 * epsilon * (pow(σ/rc,12.0) - pow(σ/rc,6.0));   // only with pair_modify shift yes

double PairLJCut::single(..., double rsq, ..., double &fforce)
{
  r2inv = 1.0 / rsq;   r6inv = r2inv * r2inv * r2inv;
  forcelj = r6inv * (lj1 * r6inv - lj2);
  fforce = factor_lj * forcelj * r2inv;              // F/r, like fpair
  philj = r6inv * (lj3 * r6inv - lj4) - offset;
  return factor_lj * philj;
}</pre>
<p>No square root anywhere: everything is in powers of 1/r², which is why LJ is cheap.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: a pair_write table, by hand</h2>
    <div class="prose"><p><code>pair_write 1 1 N r lo hi file key</code> (${gh('src/pair.cpp', 1806)}) calls single() at N evenly spaced r and writes r, E and F = fforce · r with <code>%.15g</code>. Argon in real units by default (ε in kcal/mol, σ in Å):</p></div>
    <div class="row">
      <label class="lc-in">ε<input type="number" class="lmp-in sm pl-eps" value="0.2381" step="0.01"></label>
      <label class="lc-in">σ<input type="number" class="lmp-in sm pl-sig" value="3.405" step="0.05"></label>
      <label class="lc-in">rc<input type="number" class="lmp-in sm pl-cut" value="8.5" step="0.5"></label>
      <label class="lc-in">shift<select class="lc-u pl-sh"><option value="0">no</option><option value="1">yes</option></select></label>
      <label class="lc-in">N<select class="lc-u pl-n"><option>8</option><option selected>12</option><option>25</option></select></label>
      <button class="btn sm primary" data-k="a">▶ pair_write</button><span class="lmp-status muted pl-sta"></span></div>
    <div class="calc pl-co"></div>
    <div class="tbl pl-outa"></div>
    <div class="prose pl-suma"></div>
  </section>

  <section class="block">
    <h2>3 · Modify: type pair_style twice</h2>
    <div class="prose"><p>Two types; the 1–2 pair gets its own cutoff 10.0. Then the script repeats <code>pair_style lj/cut 8.5</code>, perhaps meaning to change nothing. ${gh('src/pair_lj_cut.cpp', 430)}:</p>
<pre class="lmp-mini">if (allocated) {
  for (i = 1; i &lt;= atom-&gt;ntypes; i++)
    for (j = i; j &lt;= atom-&gt;ntypes; j++)
      if (setflag[i][j]) cut[i][j] = cut_global;     // your 10.0 becomes 8.5
}</pre>
<p>Predict the second table: zero beyond 8.5, and below 8.5 every energy moves by offset(10.0) − offset(8.5), forces unchanged.</p></div>
    <pre class="lmp-mini">${esc(SCRIPT_B)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run</button><span class="lmp-status muted pl-stb"></span></div>
    <div class="tbl pl-outb"></div>
    <div class="prose pl-sumb"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: the step at the cutoff</h2>
    <div class="prose"><p>Plain lj/cut jumps from E(rc) to 0 and from F(rc) to 0. Shift removes the energy step but not the force step, so the trajectory is the same either way. Every time a pair crosses rc, the integrator sees a force that changed abruptly, and the total energy takes a small kick.</p></div>
    <div class="tbl pl-jump"></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Melt, NVE, 1000 steps, rc = 2.5 and 4.0</button><span class="lmp-status muted pl-stc"></span></div>
    <div class="tbl pl-outc"></div>
    <div class="prose pl-sumc"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> lj1 for argon in real units?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>48 × 0.2381 × 3.405¹² = 48 × 0.2381 × 2.4389 × 10⁶ ≈ <b>2.787 × 10⁷</b> kcal/mol·Å¹².</div></div>
    <div class="ex"><p><b>2.</b> At what r is F = 0 in the table of section 2, and what is E there?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>r = 2^(1/6) σ = 3.822 Å, E = −ε = −0.2381 kcal/mol (−ε − offset with shift).</div></div>
    <div class="ex"><p><b>3.</b> You want 1–2 at cutoff 10 and everything else at 8.5, and a later part of the script changes the global cutoff. What is the safe order?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Issue <code>pair_style lj/cut 8.5</code> again first, then repeat <code>pair_coeff 1 2 … 10.0</code> after it. Any later pair_style line resets explicit cutoffs; pair_coeff lines after it set them again.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pA = () => { const g = (s, d, lo, hi) => { const v = +$(s).value; return isFinite(v) && v >= lo && v <= hi ? v : d; }; const sig = g('.pl-sig', 3.405, 0.5, 10), cut = g('.pl-cut', 8.5, sig * 1.2, 20); return { eps: g('.pl-eps', 0.2381, 0.001, 100), sig, cut, shift: $('.pl-sh').value === '1', n: +$('.pl-n').value, lo: +(0.88 * sig).toFixed(3), hi: +(cut + 0.5).toFixed(3) }; };
    const showCo = () => { const p = pA(), c = coeffs(p.eps, p.sig, p.cut, p.shift); $('.pl-co').textContent = `lj1 = ${g15(c.lj1)}   lj2 = ${g15(c.lj2)}\nlj3 = ${g15(c.lj3)}   lj4 = ${g15(c.lj4)}\noffset = ${g15(c.offset)}   table from r = ${p.lo} to ${p.hi}`; $('.pl-outa').innerHTML = ''; $('.pl-suma').innerHTML = ''; };
    root.querySelectorAll('.pl-eps,.pl-sig,.pl-cut,.pl-sh,.pl-n').forEach(e => e.oninput = e.onchange = showCo); showCo();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.pl-sta'), p = pA(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptA(p), { read: ['w.table'] }); b.disabled = false;
      const t = res.files?.['w.table']; if (res.error || !t) return err(st, res); st.textContent = '';
      const H = table(coeffs(p.eps, p.sig, p.cut, p.shift), p.n, p.lo, p.hi), Lr = parseTable(t, 'LJ'); let n = 0, ok = 0, near = 0;
      const d14 = (a, b) => Math.abs(+a - +b) <= 6e-15 * Math.max(Math.abs(+a), Math.abs(+b));
      $('.pl-outa').innerHTML = `<table>${head(['i', 'r', 'E by hand', 'F by hand', 'LAMMPS'])}${H.map((h, i) => { const l = Lr[i] || [], s = [g15(h[1]), g15(h[2]), g15(h[3])], same = s.map((q, k) => q === l[k + 1]), cl = s.map((q, k) => !same[k] && d14(q, l[k + 1])); same.forEach((x, k) => { n++; if (x) ok++; else if (cl[k]) near++; }); return row([h[0], s[0], s[1], s[2], same.every(Boolean) ? `${OK} 3 of 3` : same.every((x, k) => x || cl[k]) ? `◐ LAMMPS ${l.slice(1).filter((_, k) => !same[k]).join(' ')}` : `${BAD} ${(l.slice(1) || []).join(' ')}`]); }).join('')}</table>`;
      $('.pl-suma').innerHTML = `<p>${ok} of ${n} printed numbers identical, as strings of 15 significant digits${near ? `; the other ${near} (◐) differ only in the 15th digit` : ''}.</p>${near ? '<p class="muted">The last digit depends on <code>pow(σ, 12.0)</code>. Your hand code uses the correctly rounded value, which is what native LAMMPS gets from the Linux C library: the same five tables run natively (σ = 1, 2.64, 3.166, 3.405 with and without shift) match the hand values in every one of 15 digits. The C library inside this browser build returns pow off by an ulp or so for some σ, and that bit shows up here. Same physics, a different last bit.</p>' : ''}`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.pl-stb'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(SCRIPT_B, { read: ['w.table'] }); b.disabled = false;
      const t = res.files?.['w.table']; if (res.error || !t) return err(st, res); st.textContent = '';
      const B = parseTable(t, 'BEFORE'), Af = parseTable(t, 'AFTER');
      const c10 = coeffs(0.25, 3.2, 10.0, true), c85 = coeffs(0.25, 3.2, 8.5, true), HB = table(c10, 12, 3.0, 11.0), HA = table(c85, 12, 3.0, 11.0), dOff = c10.offset - c85.offset;
      let ok = 0, n = 0;
      $('.pl-outb').innerHTML = `<table>${head(['r', 'E before', 'E after', 'ΔE', 'F before = after?', 'by hand'])}${B.map((w, i) => { const a = Af[i], de = +a[2] - +w[2], cl = (h, l) => g15(h) === l || Math.abs(h - +l) <= 6e-15 * Math.abs(h), same = cl(HB[i][2], w[2]) && cl(HA[i][2], a[2]) && cl(HA[i][3], a[3]); n++; if (same) ok++; return row([w[1], w[2], a[2], +w[2] && +a[2] ? de.toPrecision(6) : '', w[3] === a[3] ? 'same' : `${w[3]} → ${a[3]}`, same ? OK : BAD]); }).join('')}</table>`;
      $('.pl-sumb').innerHTML = `<p>${ok} of ${n} rows match your prediction (to the 15th digit, as in section 2). offset(10.0) − offset(8.5) = ${dOff.toPrecision(6)} kcal/mol: the shift moved with the cutoff. Pairs from 8.5 to 10 Å now feel nothing. LAMMPS prints no warning.</p>`;
    };

    $('.pl-jump').innerHTML = `<table>${head(['rc (σ)', 'E(rc) / ε', 'F(rc) / (ε/σ)', 'F(rc) / F at 1.5σ'])}${[2.5, 3.0, 4.0].map(rc => { const c = coeffs(1, 1, rc, false), s = single(c, rc * rc), s15 = single(c, 2.25); return row([rc, s.e.toPrecision(4), (s.fforce * rc).toPrecision(4), ((s.fforce * rc) / (s15.fforce * 1.5)).toPrecision(3)]); }).join('')}</table>`;
    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.pl-stc'); b.disabled = true; const rs = [];
      for (const rc of [2.5, 4.0]) {
        st.textContent = `Running rc = ${rc}…`;
        const res = await runLammps(scriptC(rc)); if (res.error || !res.thermo.rows.length) { b.disabled = false; return err(st, res); }
        const c = res.thermo.cols, et = res.thermo.rows.map(r => +r[c.indexOf('TotEng')]).filter((_, i) => i >= 4);   // from step 200
        const mean = et.reduce((s, x) => s + x, 0) / et.length, sd = Math.sqrt(et.reduce((s, x) => s + (x - mean) ** 2, 0) / et.length), drift = et[et.length - 1] - et[0];
        rs.push(row([rc, et.length, mean.toFixed(5), sd.toExponential(2), drift.toExponential(2)]));
      }
      b.disabled = false; st.textContent = '';
      $('.pl-outc').innerHTML = `<table>${head(['rc', 'samples (steps 200–1000)', 'mean TotEng per atom', 'std', 'drift 200 → 1000'])}${rs.join('')}</table>`;
      $('.pl-sumc').innerHTML = '<p>The mean differs because rc = 4 counts more attraction (page C.5 tail). The spread and drift come from dt = 0.005 and from the force step; the force step at 4.0 is 27 times smaller. A smooth cutoff (lj/smooth/linear, lj/gromacs, a switching function) removes it; this browser build does not include those styles.</p>';
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
