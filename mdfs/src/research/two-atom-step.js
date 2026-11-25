// Research · C.4: one timestep for two atoms, entirely by hand. Force, energy, virial, kicks, drift, temperature, pressure — each in the same floating-point order as the C++, matched bit for bit.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const script = p => `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10 units box
create_box      1 box
create_atoms    1 single 5.0 5.0 5.0 units box
create_atoms    1 single ${p.r} 5.0 5.0 units box
mass            1 1.0
set             atom 1 vx ${p.v1x}
set             atom 2 vx ${p.v2x} vy ${p.v2y}
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
timestep        ${p.dt}
fix             1 all nve
thermo_style    custom step temp pe ke etotal press
thermo_modify   norm no format float %.17g
thermo          1
dump            1 all custom 1 t.dump id x y z vx vy vz fx fy fz
dump_modify     1 format float %.17g sort id
run             ${p.n}`;

// the C++ in the same order, in doubles
function force(x) {
  const f = [[0, 0, 0], [0, 0, 0]], d = [0, 1, 2].map(k => x[0][k] - x[1][k]);
  const rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
  let e = 0, r2inv = 0, r6inv = 0, forcelj = 0, fpair = 0;
  if (rsq < 6.25) {
    r2inv = 1.0 / rsq; r6inv = r2inv * r2inv * r2inv;
    forcelj = r6inv * (48.0 * r6inv - 24.0); fpair = forcelj * r2inv;
    for (let k = 0; k < 3; k++) { f[0][k] += d[k] * fpair; f[1][k] -= d[k] * fpair; }
    e = r6inv * (4.0 * r6inv - 4.0) - 0.0;
  }
  const vir = [0, 0, 0];
  for (let i = 0; i < 2; i++) for (let k = 0; k < 3; k++) vir[k] += f[i][k] * x[i][k];
  return { f, e, vir, d, rsq, r2inv, r6inv, forcelj, fpair };
}
function thermo(v, F) {
  let t = 0; for (let i = 0; i < 2; i++) t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * 1.0;
  const dof = 3.0, T = t * (1.0 / (dof * 1.0));
  let ke = T; ke *= 0.5 * dof * 1.0;
  const pe = 0.0 + (F.e + 0.0);
  const P = (dof * 1.0 * T + F.vir[0] + F.vir[1] + F.vir[2]) / 3.0 * (1.0 / (10.0 * 10.0 * 10.0)) * 1.0;
  return { t, T, pe, ke, et: pe + ke, P };
}
function hand(p) {
  const x = [[5.0, 5.0, 5.0], [p.r, 5.0, 5.0]], v = [[p.v1x, 0, 0], [p.v2x, p.v2y, 0]];
  const dt = p.dt, dtfm = 0.5 * dt * 1.0 / 1.0;
  let F = force(x); const out = [{ x: x.map(a => [...a]), v: v.map(a => [...a]), F, th: thermo(v, F) }];
  for (let s = 1; s <= p.n; s++) {
    const vh = [];
    for (let i = 0; i < 2; i++) { for (let k = 0; k < 3; k++) { v[i][k] += dtfm * F.f[i][k]; x[i][k] += dt * v[i][k]; } vh.push([...v[i]]); }
    F = force(x);
    for (let i = 0; i < 2; i++) for (let k = 0; k < 3; k++) v[i][k] += dtfm * F.f[i][k];
    out.push({ x: x.map(a => [...a]), v: v.map(a => [...a]), vh, F, th: thermo(v, F) });
  }
  return out;
}

export default {
  id: 'two-atom-step',
  title: 'One timestep for two atoms, by hand',
  minutes: 25,
  summary: 'Page C.2 borrowed the new force from LAMMPS. Here nothing is borrowed: two LJ atoms, and you compute the force, the energy, the virial, both half-kicks, the drift, and then Temp, PotEng, KinEng, TotEng and Press exactly the way the C++ does, operation by operation. Every number LAMMPS prints is matched to all 17 digits. Then change dt and watch TotEng drift grow as dt².',
  learn: [
    'pair lj/cut compute(): rsq, r2inv, r6inv, forcelj, fpair, in that order',
    'Pressure from the virial Σ x·f, temperature with dof = 3N − 3',
    'One full step, every number matched bit for bit',
    'TotEng is not conserved exactly; the error scales as dt²',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the five formulas LAMMPS uses</h2>
    <div class="prose"><p><b>Force and energy</b>, ${gh('src/pair_lj_cut.cpp', 107)}, with lj1 = 48εσ¹², lj2 = 24εσ⁶, lj3 = 4εσ¹², lj4 = 4εσ⁶ (${gh('src/pair_lj_cut.cpp', 509)}):</p>
<pre class="lmp-mini">delx = xtmp - x[j][0];                       // 107
rsq = delx*delx + dely*dely + delz*delz;     // 110
r2inv = 1.0 / rsq;                           // 114
r6inv = r2inv * r2inv * r2inv;
forcelj = r6inv * (lj1*r6inv - lj2);         // = 48/r¹² − 24/r⁶
fpair = factor_lj * forcelj * r2inv;         // F/r
f[i][0] += delx * fpair;   f[j][0] -= delx * fpair;
evdwl = r6inv * (lj3*r6inv - lj4) - offset;  // 128, offset = 0 without pair_modify shift</pre>
<p><b>Virial</b>, ${gh('src/pair.cpp', 1762)}: <code>virial[0] += f[i][0]*x[i][0]</code> over all atoms (and ghosts). <b>Pressure</b>, ${gh('src/compute_pressure.cpp', 256)}:</p>
<pre class="lmp-mini">P = (dof·kB·T + W_xx + W_yy + W_zz) / 3 / V</pre>
<p><b>Temperature</b>, ${gh('src/compute_temp.cpp', 92)} and 65: <code>T = Σ m v² / (dof · kB)</code> with dof = 3N − 3 = 3 for two atoms (the 3 is <code>extra_dof</code>: the centre of mass is not counted). <b>KinEng</b> is not recomputed from v: thermo takes ½ · dof · kB · T (${gh('src/thermo.cpp', 1958)}), and rounding follows that path.</p></div>
  </section>

  <section class="block">
    <h2>2 · Modify: set up the two atoms</h2>
    <div class="prose"><p>Atom 1 sits at (5, 5, 5); atom 2 at (x₂, 5, 5) in a 10 × 10 × 10 periodic box. 2^(1/6) = 1.1225 is the bottom of the well: closer pushes apart, farther pulls together; beyond 2.5 nothing.</p></div>
    <div class="row">
      <label class="lc-in">x₂<input type="number" class="lmp-in sm" data-p="r" value="6.1" step="0.05" min="5.8" max="8"></label>
      <label class="lc-in">v₁ₓ<input type="number" class="lmp-in sm" data-p="v1x" value="0.5" step="0.1"></label>
      <label class="lc-in">v₂ₓ<input type="number" class="lmp-in sm" data-p="v2x" value="-0.5" step="0.1"></label>
      <label class="lc-in">v₂ᵧ<input type="number" class="lmp-in sm" data-p="v2y" value="0.2" step="0.1"></label>
      <label class="lc-in">dt<select class="lc-u" data-p="dt"><option>0.001</option><option>0.0025</option><option selected>0.005</option><option>0.01</option><option>0.02</option></select></label>
      <label class="lc-in">steps<select class="lc-u" data-p="n"><option>1</option><option selected>2</option><option>10</option><option>20</option><option>40</option><option>100</option></select></label>
    </div>
    <div class="calc ta-hand"></div>
  </section>

  <section class="block">
    <h2>3 · Match: run the real LAMMPS</h2>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run</button><span class="lmp-status muted ta-st"></span></div>
    <pre class="lmp-mini ta-scr" hidden></pre>
    <div class="tbl ta-thermo"></div>
    <div class="tbl ta-dump"></div>
    <div class="prose ta-sum"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe</h2>
    <ul>
      <li><b>6.1 is not 6.1.</b> The dump prints x₂ = 6.0999999999999996 at step 0: 6.1 has no exact binary form. Your hand calculation starts from the same double, which is why it can match.</li>
      <li><b>Newton's third law is exact.</b> f₂ = −f₁ to the last bit, because both come from the same <code>delx * fpair</code>, once added and once subtracted. So the total momentum stays exactly what you set.</li>
      <li><b>TotEng moves.</b> With dt = 0.005, TotEng changes in the 6th digit within 2 steps. Compare at the same simulated time t = n·dt, not the same number of steps. Your hand code (now proven identical to LAMMPS) does it instantly for your atoms:<div class="tbl ta-dt"></div>Halving dt cuts the error 4 times: the dt² error of velocity Verlet (lesson 1.5), seen in the real code.</li>
      <li><b>Temp is tiny, Press is tiny.</b> Two atoms in a volume of 1000: the pressure is dominated by the box, not the pair. dof = 3 also makes the temperature look odd: with v₁ₓ = 0.5 alone, T = 0.25/3 ≠ ½mv².</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Put x₂ = 5 + 2^(1/6) = 6.122462048309373. What is fx on atom 1, roughly, and why not exactly 0?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>LAMMPS gives f₁ₓ = 3.2 × 10⁻¹⁵, not 0: the double nearest 2^(1/6) is not 2^(1/6), and 48·r6inv − 24 subtracts two numbers near 24, keeping only the round-off. PotEng prints as exactly −1.</div></div>
    <div class="ex"><p><b>2.</b> Both atoms at rest, x₂ = 6.0. After one step with dt = 0.005, how far has atom 1 moved?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>x₁ changes by dt · (½ dt · f) = ½ dt² f. At r = 1, f = 48 − 24 = 24 (repulsive), so Δx₁ = −½ · 0.005² · 24 = <b>−0.0003</b>.</div></div>
    <div class="ex"><p><b>3.</b> With <code>pair_modify shift yes</code>, which printed numbers change?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only PotEng and TotEng: each pair inside the cutoff gets −offset = −4(2.5⁻¹² − 2.5⁻⁶) = +0.016316891136. Forces, x, v, Temp and Press are untouched, because offset is not in fpair.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const P = () => { const o = {}; root.querySelectorAll('[data-p]').forEach(e => o[e.dataset.p] = +e.value); if (!(o.r > 5.5 && o.r <= 8)) o.r = 6.1; for (const k of ['v1x', 'v2x', 'v2y']) if (!isFinite(o[k])) o[k] = 0; return o; };
    const showHand = () => {
      const p = P(), H = hand(p), a = H[0], b = H[1], F0 = a.F, F1 = b.F, dtfm = 0.5 * p.dt;
      const L = [`STEP 0 → 1   (dt = ${p.dt}, dtfm = ½ dt / m = ${dtfm})`, '',
        `force at t = 0:  delx = 5 − ${s17(p.r)} = ${s17(F0.d[0])}`,
        `  rsq = ${s17(F0.rsq)}   r2inv = ${s17(F0.r2inv)}   r6inv = ${s17(F0.r6inv)}`,
        `  forcelj = r6inv(48 r6inv − 24) = ${s17(F0.forcelj)}   fpair = ${s17(F0.fpair)}`,
        `  f₁ₓ = delx · fpair = ${s17(F0.f[0][0])}     f₂ₓ = ${s17(F0.f[1][0])}`,
        `  evdwl = r6inv(4 r6inv − 4) = ${s17(F0.e)}`, '',
        `half kick:  v₁ₓ = ${s17(a.v[0][0])} + ${dtfm} × ${s17(F0.f[0][0])} = ${s17(b.vh[0][0])}`,
        `            v₂ₓ = ${s17(a.v[1][0])} + ${dtfm} × ${s17(F0.f[1][0])} = ${s17(b.vh[1][0])}`,
        `drift:      x₁ = 5 + ${p.dt} × v₁ₓ = ${s17(b.x[0][0])}`,
        `            x₂ = ${s17(p.r)} + ${p.dt} × v₂ₓ = ${s17(b.x[1][0])}     y₂ = ${s17(b.x[1][1])}`, '',
        `new force:  rsq = ${s17(F1.rsq)}   fpair = ${s17(F1.fpair)}`,
        `            f₁ = (${s17(F1.f[0][0])}, ${s17(F1.f[0][1])}, 0)`,
        `half kick:  v₁ = (${s17(b.v[0][0])}, ${s17(b.v[0][1])}, 0)`,
        `            v₂ = (${s17(b.v[1][0])}, ${s17(b.v[1][1])}, 0)`, '',
        `thermo at step 1:`,
        `  Σ m v² = ${s17(b.th.t)}     Temp = Σ m v² / 3 = ${s17(b.th.T)}`,
        `  PotEng = ${s17(b.th.pe)}     KinEng = ½ · 3 · T = ${s17(b.th.ke)}`,
        `  TotEng = ${s17(b.th.et)}`,
        `  W = Σ x·f = ${s17(F1.vir[0])} + ${s17(F1.vir[1])} + 0`,
        `  Press = (3 T + W) / 3 / 1000 = ${s17(b.th.P)}`];
      $('.ta-hand').textContent = L.join('\n');
      $('.ta-thermo').innerHTML = ''; $('.ta-dump').innerHTML = ''; $('.ta-sum').innerHTML = '';
    };
    const dtScan = () => {
      const p = P(), T = p.dt * p.n, rs = []; let prev = null;
      for (const dt of [0.01, 0.005, 0.0025, 0.00125]) { const n = Math.round(T / dt); if (n < 1 || n > 4000) continue; const et = hand({ ...p, dt, n }).map(h => h.th.et); const m = Math.max(...et.map(x => Math.abs(x - et[0]))); rs.push(row([dt, n, m.toExponential(3), prev && m ? (prev / m).toFixed(2) + '×' : ''])); prev = m; }
      $('.ta-dt').innerHTML = `<table>${head(['dt', `steps for t = ${+T.toFixed(6)}`, 'largest |ΔTotEng|', 'smaller by'])}${rs.join('')}</table>`;
    };
    root.querySelectorAll('[data-p]').forEach(e => e.oninput = e.onchange = () => { showHand(); dtScan(); }); showHand(); dtScan();
    const cell = (h, l) => { const ok = h === l; return [`<span class="mono">${s17(h)}</span>`, `<span class="mono">${ok ? '✓' : s17(l)}</span>`, ok]; };
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.ta-st'), p = P(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(script(p), { read: ['t.dump'] }); b.disabled = false;
      const txt = res.files?.['t.dump'];
      if (res.error || !txt || !res.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; return; }
      st.textContent = ''; $('.ta-scr').hidden = false; $('.ta-scr').textContent = script(p);
      const H = hand(p); let n = 0, ok = 0;
      const tc = res.thermo.cols, keys = ['Temp', 'PotEng', 'KinEng', 'TotEng', 'Press'], hk = ['T', 'pe', 'ke', 'et', 'P'];
      const trs = res.thermo.rows.map(r => { const s = +r[0], h = H[s].th; const cs = keys.map((k, i) => cell(h[hk[i]], parseFloat(r[tc.indexOf(k)]))); cs.forEach(c => { n++; if (c[2]) ok++; }); return row([s, ...cs.map(c => c[0] + (c[2] ? ' <span class="m-ok">✓</span>' : ` <span class="m-bad">✗ ${c[1]}</span>`))]); });
      $('.ta-thermo').innerHTML = `<table>${head(['step', 'Temp', 'PotEng', 'KinEng', 'TotEng', 'Press'])}${trs.join('')}</table>`;
      const D = {};
      for (const f of txt.split('ITEM: TIMESTEP').slice(1)) { const ls = f.trim().split('\n'), s = +ls[0], m = {}; for (const l of ls.slice(8)) { const w = l.trim().split(/\s+/).map(Number); m[w[0]] = w.slice(1); } D[s] = m; }
      let dr = '';
      for (const s of Object.keys(D).map(Number)) for (const id of [1, 2]) {
        const L = D[s][id], h = [...H[s].x[id - 1], ...H[s].v[id - 1], ...H[s].F.f[id - 1]];
        let all = true; for (let k = 0; k < 9; k++) { n++; if (h[k] === L[k]) ok++; else all = false; }
        if (s <= 2) dr += row([s, id, ...[0, 1, 3, 4, 6, 7].map(k => `<span class="mono">${s17(L[k])}</span>`), all ? '<span class="m-ok">✓ all 9</span>' : '<span class="m-bad">✗</span>']);
      }
      $('.ta-dump').innerHTML = `<table>${head(['step', 'atom', 'x', 'y', 'vx', 'vy', 'fx', 'fy', 'by hand'])}${dr}</table>`;
      const et = H.map(h => h.th.et), dE = Math.max(...et.map(x => Math.abs(x - et[0])));
      $('.ta-sum').innerHTML = `<p><b>${ok} of ${n}</b> printed numbers are identical to your hand calculation. ${ok === n ? '<span class="m-ok">✓</span>' : ''}</p><p>Largest change of TotEng from step 0 over ${p.n} step${p.n > 1 ? 's' : ''}: <b>${dE.toExponential(3)}</b> at dt = ${p.dt}.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
