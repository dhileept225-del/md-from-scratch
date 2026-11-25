// Research · C.7: the other loops. Respa::recurse (two time steps) and MinSD + backtracking line search, both reproduced by hand bit for bit.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const OK = '<span class="m-ok">✓</span>';

// ---------- RESPA: a bonded pair and a lone atom ----------
const DATA = `three atoms: a bonded pair and a lone atom

3 atoms
1 bonds
1 atom types
1 bond types

0 10 xlo xhi
0 10 ylo yhi
0 10 zlo zhi

Masses

1 1.0

Atoms # bond

1 1 1 5.0 5.0 5.0
2 1 1 6.05 5.0 5.0
3 2 1 5.5 6.2 5.0

Velocities

1 0.3 0.0 0.0
2 -0.3 0.1 0.0
3 0.0 -0.4 0.2

Bonds

1 1 1 2
`;
const X0 = [[5.0, 5.0, 5.0], [6.05, 5.0, 5.0], [5.5, 6.2, 5.0]], V0 = [[0.3, 0, 0], [-0.3, 0.1, 0], [0, -0.4, 0.2]];
const scriptR = (rs, dt, n, nin) => `units           lj
atom_style      bond
read_data       tri.data
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
bond_style      harmonic
bond_coeff      1 2000.0 1.0
special_bonds   lj 0.0 0.0 0.0
neighbor        0.3 bin
${rs ? `run_style       respa 2 ${nin} bond 1 pair 2\n` : ''}timestep        ${dt}
fix             1 all nve
thermo_style    custom step pe ke etotal
thermo_modify   norm no format float %.17g
thermo          1
run             ${n}`;
const z3 = () => [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
function fBond(x) {
  const f = z3(), d = [0, 1, 2].map(k => x[0][k] - x[1][k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2], r = Math.sqrt(rsq), dr = r - 1.0, rk = 2000.0 * dr;
  const fb = r > 0 ? -2.0 * rk / r : 0.0;
  for (let k = 0; k < 3; k++) { f[0][k] += d[k] * fb; f[1][k] -= d[k] * fb; }
  return { f, e: rk * dr };
}
function fPair(x) {
  const f = z3(); let e = 0;
  for (const [i, j] of [[0, 2], [1, 2]]) {
    const d = [0, 1, 2].map(k => x[i][k] - x[j][k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (rsq < 6.25) { const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48.0 * r6 - 24.0) * r2; for (let k = 0; k < 3; k++) { f[i][k] += d[k] * fp; f[j][k] -= d[k] * fp; } e += r6 * (4.0 * r6 - 4.0); }
  }
  return { f, e };
}
function thermoR(v, B, P) {
  let t = 0; for (let i = 0; i < 3; i++) t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * 1.0;
  const T = t * (1.0 / (6.0 * 1.0)); let ke = T; ke *= 0.5 * 6.0 * 1.0;
  let pe = 0.0 + (P.e + 0.0); pe += B.e;
  return { pe, ke, et: pe + ke };
}
function handR(rs, dt, n, nin) {
  const x = X0.map(a => [...a]), v = V0.map(a => [...a]);
  let B = fBond(x), P = fPair(x); const out = [thermoR(v, B, P)]; let nb = 1, np = 1;
  const kick = (F, c) => { for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) v[i][k] += c * F[i][k]; };
  for (let s = 1; s <= n; s++) {
    if (rs) {
      const so = dt, si = so / nin, dfo = 0.5 * so * 1.0 / 1.0, dfi = 0.5 * si * 1.0 / 1.0;
      kick(P.f, dfo);
      for (let it = 0; it < nin; it++) {
        for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { v[i][k] += dfi * B.f[i][k]; x[i][k] += si * v[i][k]; }
        B = fBond(x); nb++; kick(B.f, dfi);
      }
      P = fPair(x); np++; kick(P.f, dfo);
    } else {
      const df = 0.5 * dt * 1.0 / 1.0, F = P.f.map((a, i) => a.map((p, k) => p + B.f[i][k]));
      for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { v[i][k] += df * F[i][k]; x[i][k] += dt * v[i][k]; }
      P = fPair(x); B = fBond(x); np++; nb++;
      kick(P.f.map((a, i) => a.map((p, k) => p + B.f[i][k])), df);
    }
    out.push(thermoR(v, B, P));
  }
  return { out, nb, np };
}

// ---------- minimize: two atoms, steepest descent, backtracking ----------
const scriptM = (x2, dmax) => `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10 units box
create_box      1 box
create_atoms    1 single 5.0 5.0 5.0 units box
create_atoms    1 single ${x2} 5.0 5.0 units box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
min_style       sd
min_modify      line backtrack dmax ${dmax}
thermo_style    custom step pe fnorm
thermo_modify   norm no format float %.17g
thermo          1
minimize        0.0 1.0e-8 100 1000`;
function ef(x) {
  const d = [0, 1, 2].map(k => x[k] - x[3 + k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2], f = [0, 0, 0, 0, 0, 0]; let e = 0;
  if (rsq < 6.25) { const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48.0 * r6 - 24.0) * r2; for (let k = 0; k < 3; k++) { f[k] += d[k] * fp; f[3 + k] -= d[k] * fp; } e = r6 * (4.0 * r6 - 4.0) - 0.0; }
  return { e: 0.0 + (e + 0.0), f };
}
function handM(x2, dmax) {
  let x = [5.0, 5.0, 5.0, x2, 5.0, 5.0], { e, f } = ef(x), h = [...f], neval = 0;
  const fn = g => { let s = 0; for (const q of g) s += q * q; return Math.sqrt(s); };
  const it = [{ e, fn: fn(f) }];
  for (let k = 1; k <= 100; k++) {
    const eo = e; let fdh = 0; for (let i = 0; i < 6; i++) fdh += f[i] * h[i];
    if (fdh <= 0) return { it, stop: 'search direction is not downhill', neval };
    let hm = 0; for (let i = 0; i < 6; i++) hm = Math.max(hm, Math.abs(h[i]));
    if (hm === 0) return { it, stop: 'forces are zero', neval };
    let a = Math.min(1.0, dmax / hm), nb = 0, stop = null, di, de; const x0 = [...x];
    for (;;) {
      x = x0.map((q, i) => q + a * h[i]); ({ e, f } = ef(x)); neval++;
      di = -0.4 * a * fdh; de = e - eo;
      if (de <= di) break;
      a *= 0.5; nb++;
      if (a <= 0 || di >= -1e-8) { x = [...x0]; ({ e, f } = ef(x)); neval++; stop = de < 0 ? 'energy tolerance' : 'linesearch alpha is zero'; break; }
    }
    it.push({ e, fn: fn(f), a, nb, di, stop });
    if (stop) return { it, stop, neval };
    let s2 = 0; for (const q of f) s2 += q * q;
    if (s2 < 1e-8 * 1e-8) return { it, stop: 'force tolerance', neval };
    h = [...f];
  }
  return { it, stop: 'max iterations', neval };
}

export default {
  id: 'other-loops',
  title: 'The other loops: RESPA and minimize',
  minutes: 25,
  summary: 'Verlet::run is not the only loop. Respa::recurse nests a fast loop for stiff bonds inside a slow loop for pair forces; Min::run takes steps along the force until the energy stops falling. You reproduce both by hand for 2 and 3 atoms, every printed number bit for bit, and find out why a minimization with ftol 1e-8 stops at a force of 1e-3 and reports "energy tolerance" though you set etol to 0.',
  learn: [
    'Respa::recurse: kick with slow forces, n inner Verlet steps with fast forces, kick again',
    'Why RESPA is nearly as accurate as a small dt at the cost of a large one',
    'MinSD and the backtracking line search: alpha, Armijo test, EMACH',
    'Why minimize stops early, and what the stopping message means',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: RESPA, a loop inside the loop</h2>
    <div class="prose"><p><code>run_style respa 2 4 bond 1 pair 2</code>: two levels, the inner one 4 times shorter (${gh('src/respa.cpp', 330)}: <code>step[0] = step[1] / loop[0]</code>), bonds at level 1 (inner), pair at level 2 (outer). ${gh('src/respa.cpp', 577)}, trimmed:</p>
<pre class="lmp-mini">void Respa::recurse(int ilevel)
{
  copy_flevel_f(ilevel);                  // f = the stored force of this level
  for (int iloop = 0; iloop &lt; loop[ilevel]; iloop++) {
    modify-&gt;initial_integrate_respa(vflag, ilevel, iloop);   // 584
    ...
    if (ilevel) recurse(ilevel - 1);      // 642: the faster level, loop[ilevel-1] times
    force_clear();                        // 649
    if (level_pair == ilevel) force-&gt;pair-&gt;compute(eflag, vflag);   // 663
    if (level_bond == ilevel) force-&gt;bond-&gt;compute(eflag, vflag);   // 679
    modify-&gt;final_integrate_respa(ilevel, iloop);            // 710
  }
  copy_f_flevel(ilevel);                  // store this level's force
}</pre>
<p>And fix nve (${gh('src/fix_nve.cpp', 150)}) does a half kick only at outer levels, kick + drift at the innermost one:</p>
<pre class="lmp-mini">dtf = 0.5 * step_respa[ilevel] * force-&gt;ftm2v;
if (ilevel == 0) initial_integrate(vflag);   // v += dtf F/m; x += dt v
else final_integrate();                      // v += dtf F/m only</pre>
<p>So one outer step is: v += ½Δt·F<sub>pair</sub>; then 4 × (velocity Verlet with F<sub>bond</sub> and Δt/4); then v += ½Δt·F<sub>pair</sub> with the new pair force.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: three atoms, by hand</h2>
    <div class="prose"><p>Atoms 1–2: a stiff bond, E = 2000 (r − 1)², vibration period 0.070. Atom 3: LJ with both (the bonded pair has no LJ: <code>special_bonds lj 0 0 0</code>).</p></div>
    <div class="row">
      <label class="lc-in">loop<select class="lc-u ol-rs"><option value="1" selected>run_style respa</option><option value="0">run_style verlet</option></select></label>
      <label class="lc-in">dt<select class="lc-u ol-dt"><option>0.00125</option><option selected>0.005</option><option>0.01</option></select></label>
      <label class="lc-in">inner loops<select class="lc-u ol-in"><option>2</option><option selected>4</option><option>8</option></select></label>
      <button class="btn sm primary" data-k="r">▶ Run 5 steps</button><span class="lmp-status muted ol-str"></span></div>
    <pre class="lmp-mini ol-scr"></pre>
    <div class="tbl ol-outr"></div>
    <div class="prose"><p><b>Same time t = 5, by hand</b> (your code is now LAMMPS, bit for bit): largest change of TotEng, and how many times each force was computed.</p></div>
    <div class="tbl ol-scan"></div>
  </section>

  <section class="block">
    <h2>3 · Understand: minimize, steepest descent with backtracking</h2>
    <div class="prose"><p>${gh('src/min_sd.cpp', 56)}: each iteration takes a line search along h = F, checks the stopping rules, and sets h = F again. The line search, ${gh('src/min_linesearch.cpp', 176)}:</p>
<pre class="lmp-mini">fdothall = Σ f·h                                  // = |F|² for sd
alpha = MIN(ALPHA_MAX, dmax/hmaxall);              // 212: no atom moves more than dmax
while (true) {
  ecurrent = alpha_step(alpha,1);                  // x = x0 + alpha·h, new E and F
  de_ideal = -BACKTRACK_SLOPE*alpha*fdothall;      // 260: want at least 0.4 of the linear drop
  if (ecurrent - eoriginal &lt;= de_ideal) return 0;
  alpha *= ALPHA_REDUCE;                           // 272: halve
  if (alpha &lt;= 0.0 || de_ideal &gt;= -EMACH) {        // 279: EMACH = 1e-8
    ecurrent = alpha_step(0.0,0);                  // go back
    if (de &lt; 0.0) return ETOL; else return ZEROALPHA;
  }
}</pre></div>
  </section>

  <section class="block">
    <h2>4 · Match: minimize two atoms, by hand</h2>
    <div class="prose"><p><code>min_style sd</code>, <code>min_modify line backtrack</code>, <code>minimize 0.0 1.0e-8 100 1000</code>: etol 0, ftol 10⁻⁸. Predict every iteration:</p></div>
    <div class="row">
      <label class="lc-in">x₂<input type="number" class="lmp-in sm ol-x2" value="6.4" step="0.05" min="5.9" max="7.4"></label>
      <label class="lc-in">dmax<select class="lc-u ol-dm"><option>0.02</option><option selected>0.1</option><option>0.5</option></select></label>
      <button class="btn sm primary" data-k="m">▶ Minimize</button><span class="lmp-status muted ol-stm"></span></div>
    <div class="tbl ol-handm"></div>
    <div class="prose ol-summ"></div>
  </section>

  <section class="block prose">
    <h2>5 · Observe</h2>
    <ul>
      <li><b>RESPA buys accuracy cheaply.</b> At dt = 0.005 plain Verlet takes 14 steps per bond vibration and loses about 0.25 in TotEng; RESPA with 4 inner loops does almost as well as Verlet at dt = 0.00125 but computes the pair force 4 times less often. With thousands of atoms the pair force is the expensive part.</li>
      <li><b>"Energy tolerance" with etol = 0.</b> Near the minimum, |F|² is about 10⁻⁶, so the hoped-for drop 0.4·alpha·|F|² falls below EMACH = 10⁻⁸ and the line search gives up and reports ETOL. The force is then about 10⁻³, far above ftol. Real force tolerances need <code>min_style cg</code> with the default quadratic line search, or <code>fire</code>. Try x₂ = 6.1 with dmax 0.02: there the last trial raised the energy (de &gt; 0), so the same exit reports "linesearch alpha is zero".</li>
      <li><b>Force evaluations count line-search trials.</b> LAMMPS's "Iterations, force evaluations" adds every alpha tried, plus the step back on failure.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> With <code>run_style respa 3 2 2 bond 1 angle 2 pair 3</code> and timestep 4 fs, what is the bond time step?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>step[2] = 4, step[1] = 4/2 = 2, step[0] = 2/2 = <b>1 fs</b>. Angles at 2 fs, pair at 4 fs.</div></div>
    <div class="ex"><p><b>2.</b> First iteration of section 4 with x₂ = 6.4, dmax = 0.1: what alpha is tried first?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>hmax = |F| on one atom = 1.672, so alpha = min(1, 0.1/1.672) = <b>0.0598</b>, and each atom moves exactly 0.1. Energy falls enough, so no backtrack.</div></div>
    <div class="ex"><p><b>3.</b> Why are minimizations of large systems in LJ units sensitive to <code>thermo_modify norm</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>With norm yes, both the energy and fdothall are divided by N (${gh('src/min_linesearch.cpp', 198)} and ${gh('src/min.cpp', 599)}), so etol and EMACH compare against per-atom changes: the line search gives up at a different point than with norm no.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pR = () => ({ rs: $('.ol-rs').value === '1', dt: +$('.ol-dt').value, nin: +$('.ol-in').value });
    const scan = () => {
      const rs = [];
      for (const [r, dt, nin] of [[0, 0.01, 1], [0, 0.005, 1], [0, 0.00125, 1], [1, 0.01, 8], [1, 0.005, 4], [1, 0.005, 2]]) {
        const n = Math.round(5 / dt), H = handR(!!r, dt, n, nin), e0 = H.out[0].et; let m = 0; for (const o of H.out) m = Math.max(m, Math.abs(o.et - e0));
        rs.push(row([r ? `respa, ${nin} inner` : 'verlet', dt, n, H.np, H.nb, m.toExponential(2)]));
      }
      $('.ol-scan').innerHTML = `<table>${head(['loop', 'dt', 'steps', 'pair calls', 'bond calls', 'largest |ΔTotEng|'])}${rs.join('')}</table>`;
    };
    const showR = () => { const p = pR(); $('.ol-scr').textContent = scriptR(p.rs, p.dt, 5, p.nin).split('\n').slice(9, 12).join('\n'); $('.ol-in').disabled = !p.rs; $('.ol-outr').innerHTML = ''; };
    ['.ol-rs', '.ol-dt', '.ol-in'].forEach(s => $(s).onchange = showR); showR(); scan();
    $('[data-k=r]').onclick = async e => {
      const b = e.currentTarget, st = $('.ol-str'), p = pR(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptR(p.rs, p.dt, 5, p.nin), { files: { 'tri.data': DATA } }); b.disabled = false;
      if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = handR(p.rs, p.dt, 5, p.nin), c = res.thermo.cols; let n = 0, ok = 0;
      $('.ol-outr').innerHTML = `<table>${head(['step', 'PotEng by hand', 'KinEng by hand', 'TotEng by hand', 'LAMMPS'])}${res.thermo.rows.map(r => { const h = H.out[+r[0]], l = ['PotEng', 'KinEng', 'TotEng'].map(k => parseFloat(r[c.indexOf(k)])); const same = [h.pe, h.ke, h.et].map((q, i) => q === l[i]); same.forEach(s => { n++; if (s) ok++; }); return row([r[0], ...[h.pe, h.ke, h.et].map(q => `<span class="mono">${s17(q)}</span>`), same.every(Boolean) ? `${OK} all 3 identical` : '<span class="m-bad">✗</span>']); }).join('')}</table><p>${ok} of ${n} identical.</p>`;
    };

    const showM = () => {
      let x2 = +$('.ol-x2').value; if (!(x2 >= 5.9 && x2 <= 7.4)) x2 = 6.4; const H = handM(x2, +$('.ol-dm').value);
      $('.ol-handm').innerHTML = `<table>${head(['iter', 'alpha', 'halvings', 'PotEng by hand', 'Fnorm by hand', 'LAMMPS'])}${H.it.map((q, i) => row([i, i ? s17(q.a) : '', i ? q.nb : '', `<span class="mono">${s17(q.e)}</span>`, `<span class="mono">${s17(q.fn)}</span>`, `<span class="ol-l" data-i="${i}"></span>`])).join('')}</table>`;
      $('.ol-summ').innerHTML = `<p>Your prediction: stops with <b>${esc(H.stop)}</b> after ${H.it.length - 1} iterations and ${H.neval} force evaluations.</p>`;
      return H;
    };
    $('.ol-x2').oninput = showM; $('.ol-dm').onchange = showM; showM();
    $('[data-k=m]').onclick = async e => {
      const b = e.currentTarget, st = $('.ol-stm'); let x2 = +$('.ol-x2').value; if (!(x2 >= 5.9 && x2 <= 7.4)) x2 = 6.4; const dm = +$('.ol-dm').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptM(x2, dm)); b.disabled = false;
      if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = showM(); let n = 0, ok = 0;
      for (const r of res.thermo.rows) { const i = +r[0], h = H.it[i], cel = root.querySelector(`.ol-l[data-i="${i}"]`); if (!h || !cel) continue; const s = parseFloat(r[1]) === h.e && parseFloat(r[2]) === h.fn; n += 2; if (s) ok += 2; cel.innerHTML = s ? `${OK} both identical` : `<span class="m-bad">✗ ${r[1]} ${r[2]}</span>`; }
      const txt = res.lines.join('\n'), crit = (txt.match(/Stopping criterion = (.*)/) || [])[1] || '?', ev = (txt.match(/Iterations, force evaluations = (\d+) (\d+)/) || []);
      $('.ol-summ').innerHTML += `<p>LAMMPS: stopping criterion = <b>${esc(crit.trim())}</b> ${crit.trim() === H.stop ? OK : '<span class="m-bad">✗</span>'}; iterations, force evaluations = ${ev[1]} ${ev[2]} ${+ev[1] === H.it.length - 1 && +ev[2] === H.neval ? OK : '<span class="m-bad">✗</span>'}. ${ok} of ${n} printed numbers identical.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
