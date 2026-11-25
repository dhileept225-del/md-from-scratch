// Research · D.13: constraints, SHAKE and RATTLE. fix shake's prediction-and-quadratic step reproduced bit for bit for a diatomic next to an argon atom (native runs: the browser build has no RIGID package), RATTLE's velocity correction, and a live run showing why constraints allow a 2 fs step.
import { runLammps } from './lammps-runner.js';
import { NATIVE } from './shake-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const FTM2V = 1.0 / 48.88821291 / 48.88821291, MVV2E = 48.88821291 * 48.88821291, KB = 0.0019872067;
const M = [14.007, 14.007, 39.948], B0 = 1.098, DT = 2.0;
const X0 = [[9.45, 10.0, 10.0], [10.55, 10.02, 10.0], [10.1, 13.4, 10.3]], V0 = [[0.004, 0.002, -0.001], [-0.003, -0.004, 0.002], [0.0, -0.002, 0.0005]];
const DATA = `a diatomic and an argon

3 atoms
1 bonds
2 atom types
1 bond types

0 20 xlo xhi
0 20 ylo yhi
0 20 zlo zhi

Masses

1 14.007
2 39.948

Atoms # full

${X0.map((x, i) => `${i + 1} 1 ${i < 2 ? 1 : 2} 0.0 ${x.join(' ')}`).join('\n')}

Velocities

${V0.map((v, i) => `${i + 1} ${v.join(' ')}`).join('\n')}

Bonds

1 1 1 2
`;
const SCRIPT = fx => `units           real
atom_style      full
read_data       dia.data
pair_style      lj/cut 10.0
pair_coeff      1 1 0.07 3.3
pair_coeff      2 2 0.238 3.4
pair_coeff      1 2 0.13 3.35
bond_style      harmonic
bond_coeff      1 1000.0 1.098
timestep        ${fx === 'none' ? 'DT' : '2.0'}
fix             1 all nve
${fx === 'none' ? '' : `fix             2 all ${fx} 1.0e-10 100 0 b 1\n`}thermo_style    custom step temp pe ke etotal
thermo_modify   format float %.17g`;

function forces(x) {
  const f = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]; let e = 0; const ep = 0.13, s = 3.35;
  const lj1 = 48.0 * ep * Math.pow(s, 12.0), lj2 = 24.0 * ep * Math.pow(s, 6.0), lj3 = 4.0 * ep * Math.pow(s, 12.0), lj4 = 4.0 * ep * Math.pow(s, 6.0);
  for (const [i, j] of [[0, 2], [1, 2]]) {   // 1-2 is a bond: excluded by special_bonds 0 0 0
    const d = [0, 1, 2].map(k => x[i][k] - x[j][k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2], r2 = 1.0 / rsq, r6 = r2 * r2 * r2, fp = r6 * (lj1 * r6 - lj2) * r2;
    for (let k = 0; k < 3; k++) { f[i][k] += d[k] * fp; f[j][k] -= d[k] * fp; } e += r6 * (lj3 * r6 - lj4);
  }
  return { f, e };
}
// FixShake::unconstrained_update + FixShake::shake (2-atom cluster)
function shake(x, v, f, dtfsq) {
  const xs = [0, 1].map(i => [0, 1, 2].map(k => x[i][k] + DT * v[i][k] + (dtfsq / M[i]) * f[i][k]));
  const r = [0, 1, 2].map(k => x[0][k] - x[1][k]), s = [0, 1, 2].map(k => xs[0][k] - xs[1][k]);
  const rsq = r[0] * r[0] + r[1] * r[1] + r[2] * r[2], ssq = s[0] * s[0] + s[1] * s[1] + s[2] * s[2], i0 = 1.0 / M[0], i1 = 1.0 / M[1];
  const a = (i0 + i1) * (i0 + i1) * rsq, b = 2.0 * (i0 + i1) * (s[0] * r[0] + s[1] * r[1] + s[2] * r[2]), c = ssq - B0 * B0;
  let det = b * b - 4.0 * a * c; if (det < 0) det = 0;
  const l1 = (-b + Math.sqrt(det)) / (2.0 * a), l2 = (-b - Math.sqrt(det)) / (2.0 * a); let lam = Math.abs(l1) <= Math.abs(l2) ? l1 : l2; lam /= dtfsq;
  for (let k = 0; k < 3; k++) { f[0][k] += lam * r[k]; f[1][k] -= lam * r[k]; }
  return lam;
}
const temp = v => { let t = 0; for (let i = 0; i < 3; i++) t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * M[i]; return t * (MVV2E / (5.0 * KB)); };
function hand() {
  const x = X0.map(q => [...q]), v = V0.map(q => [...q]); let { f, e } = forces(x); const out = [];
  const dth = 0.5 * DT * DT * FTM2V, dts = DT * DT * FTM2V, dtf = 0.5 * DT * FTM2V;
  // setup: correct_coordinates (v = f = 0, half dt²), then shake_end_of_step (half dt²)
  const z = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], fc = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], lam0 = shake(x, z, fc, dth);
  for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) x[i][k] = x[i][k] + (dth / M[i]) * fc[i][k];
  shake(x, v, f, dth);
  out.push({ x: x.map(q => [...q]), v: v.map(q => [...q]), f: f.map(q => [...q]), T: temp(v), e, lam0 });
  for (let s = 1; s <= 10; s++) {
    for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) { v[i][k] += (dtf / M[i]) * f[i][k]; x[i][k] += DT * v[i][k]; }
    ({ f, e } = forces(x)); shake(x, v, f, dts);
    for (let i = 0; i < 3; i++) for (let k = 0; k < 3; k++) v[i][k] += (dtf / M[i]) * f[i][k];
    out.push({ x: x.map(q => [...q]), v: v.map(q => [...q]), f: f.map(q => [...q]), T: temp(v), e });
  }
  return out;
}
const bondInfo = d => { const a = d[0].map(Number), b = d[1].map(Number), r = [0, 1, 2].map(k => a[k] - b[k]), len = Math.hypot(...r), vr = [0, 1, 2].map(k => a[3 + k] - b[3 + k]); return { len, vbond: (r[0] * vr[0] + r[1] * vr[1] + r[2] * vr[2]) / len }; };

export default {
  id: 'constraints',
  title: 'Constraints: SHAKE and RATTLE',
  minutes: 20,
  summary: 'fix shake does not move atoms back after the fact: in post_force it predicts where each atom will be one step later and adds the force that makes the bond come out exactly right. For a two-atom bond that force is the root of a quadratic. You reproduce it, and ten 2-fs steps of a diatomic next to an argon atom, bit for bit against native LAMMPS (the browser build has no RIGID package); see why the first SHAKE step throws away energy that RATTLE removes at setup; and run the same molecule live without constraints to see what 2 fs does to a stiff bond.',
  learn: [
    'Prediction: x(t+dt) = x + dt·v(t−½dt) + dt²·f/m, and why dt² not ½dt²',
    'The 2-atom SHAKE quadratic and the constraint force λ r₀₁',
    'Setup: correct_coordinates, and why SHAKE leaves bond velocities but RATTLE does not',
    'Degrees of freedom lost and time step gained',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a force that keeps the bond</h2>
    <div class="prose"><p>In post_force (stage 5, page C.2) the velocities are still at v(t − ½dt), so the next position will be x + dt·v + dt²·f/m. ${gh('src/RIGID/fix_shake.cpp', 1946)} predicts it, and ${gh('src/RIGID/fix_shake.cpp', 2032)} finds the λ with |s₀₁ + (1/m₀ + 1/m₁) λ r₀₁|² = b₀²:</p>
<pre class="lmp-mini">xshake = x + dtv*v + dtfsq/m * f;          // dtfsq = dt²·ftm2v (½dt² at setup and for RATTLE, line 507)
a = (1/m0 + 1/m1)² |r01|²;   b = 2 (1/m0 + 1/m1) (s01·r01);   c = |s01|² − bond²;
lamda = the smaller root of a λ² + b λ + c = 0;   lamda /= dtfsq;
f[i0] += lamda * r01;   f[i1] −= lamda * r01;</pre>
<p>Clusters of 3, 4 atoms or water (shake3, shake4, shake3angle) solve the same conditions together by Newton iteration up to the tolerance. The bond itself is removed from the bond list, so E_bond is 0, and each constraint removes one degree of freedom from compute temp (page D.12): here 3·3 − 3 − 1 = 5.</p>
<p>At setup (${gh('src/RIGID/fix_shake.cpp', 3466)}) correct_coordinates moves the atoms onto the constraint with v = f = 0. correct_velocities is empty for SHAKE; RATTLE overrides it and removes the velocity along each bond (${gh('src/RIGID/fix_rattle.cpp', 313)}: l₀₁ = −r₀₁·v₀₁ / (|r₀₁|²(1/m₀ + 1/m₁))).</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: ten SHAKE steps, by hand</h2>
    <div class="prose"><p>N₂-like diatomic (m = 14.007, b₀ = 1.098 Å, started at 1.100 Å with the atoms moving towards each other along the bond) and an argon atom 3.4 Å away, dt = 2 fs, <code>fix shake 1.0e-10 100 0 b 1</code>. The native run's dump (x, v, f of all three atoms every step) and thermo are on this page; your browser redoes the whole thing.</p></div>
    <pre class="lmp-mini">${esc(SCRIPT('shake'))}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Compute by hand and compare</button></div>
    <div class="calc cs-setup" hidden></div>
    <div class="tbl cs-out"></div>
    <div class="prose cs-sum"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: SHAKE vs RATTLE, and dt without constraints</h2>
    <div class="tbl cs-sr"></div>
    <div class="prose cs-srsum"></div>
    <div class="prose"><p>Now the same molecule <b>without</b> constraints, live: the bond is a spring, 1000 (r − 1.098)² kcal/mol, period about 18 fs. Run 1000 fs at three time steps and compare with native SHAKE at 2 fs (TotEng from step 10 on).</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run dt = 0.5, 1, 2 fs</button><span class="lmp-status muted cs-stb"></span></div>
    <div class="tbl cs-outb"></div>
    <div class="prose cs-sumb"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why dt² and not ½dt² in the prediction during the run?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>At post_force, v is v(t − ½dt). The next position is x + dt·v(t + ½dt) = x + dt·(v(t − ½dt) + dt·f/m): a whole dt of force. At setup, v is v(0), so only half a kick is ahead: ½dt².</div></div>
    <div class="ex"><p><b>2.</b> 1000 rigid water molecules with fix shake on bonds and angle: how many degrees of freedom?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>9000 − 3000 − 3 = <b>5997</b> (page D.12).</div></div>
    <div class="ex"><p><b>3.</b> SHAKE's first step changed TotEng from 0.902 to 0.480 kcal/mol here. Where did the energy go?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The atoms were moving towards each other along the bond (relative 0.007 Å/fs). The constraint force cancels that motion in one step: ½ μ v_bond² ≈ ½ × 7.0 × 0.0070² × 2390 = 0.41 kcal/mol vanishes. RATTLE removes it at setup, before step 0 is printed.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    $('[data-k=a]').onclick = () => {
      const H = hand(); let n = 0, ok = 0;
      $('.cs-setup').hidden = false;
      $('.cs-setup').textContent = `setup: start bond ${Math.hypot(...[0, 1, 2].map(k => X0[0][k] - X0[1][k])).toFixed(6)} Å → after correct_coordinates ${Math.hypot(...[0, 1, 2].map(k => H[0].x[0][k] - H[0].x[1][k])).toPrecision(17)} Å`;
      $('.cs-out').innerHTML = `<table>${head(['step', 'x₁ by hand', 'fx₁ (incl. constraint)', 'T', 'PotEng', 'bond length', 'all 27 dump values + T + PE'])}${H.map((h, s) => {
        const d = NATIVE.sd[s], th = NATIVE.st[s]; let same = 0;
        for (let i = 0; i < 3; i++) for (let k = 0; k < 9; k++) { const hv = k < 3 ? h.x[i][k] : k < 6 ? h.v[i][k - 3] : h.f[i][k - 6]; if (hv === parseFloat(d[i][k])) same++; }
        if (h.T === parseFloat(th[0])) same++; if (h.e === parseFloat(th[1])) same++; n += 29; ok += same;
        return row([s, `<span class="mono">${s17(h.x[0][0])}</span>`, `<span class="mono">${s17(h.f[0][0])}</span>`, h.T.toFixed(6), h.e.toFixed(9), Math.hypot(...[0, 1, 2].map(k => h.x[0][k] - h.x[1][k])).toPrecision(16), same === 29 ? `${OK} 29 of 29` : `${BAD} ${same} of 29`]); }).join('')}</table>`;
      $('.cs-sum').innerHTML = `<p>${ok} of ${n} numbers identical to the native run, to the last bit. The bond length is 1.098 to 15 digits at every step: for a single bond the quadratic is solved exactly, no iteration.</p>`;
    };

    { const rows = []; for (let s = 0; s <= 10; s += 2) { const a = bondInfo(NATIVE.sd[s]), b = bondInfo(NATIVE.rd[s]); rows.push(row([s, (+NATIVE.st[s][0]).toFixed(4), a.vbond.toExponential(2), (+NATIVE.st[s][3]).toFixed(6), (+NATIVE.rt[s][0]).toFixed(4), b.vbond.toExponential(2), (+NATIVE.rt[s][3]).toFixed(6)])); }
      $('.cs-sr').innerHTML = `<table>${head(['step', 'T (SHAKE)', 'v along bond', 'TotEng', 'T (RATTLE)', 'v along bond', 'TotEng'])}${rows.join('')}</table>`;
      $('.cs-srsum').innerHTML = '<p>Native runs. SHAKE keeps the bond length but not the bond-direction velocity at the printed moment: at step 0 the atoms still approach at 0.007 Å/fs and T counts that motion. The first step removes it, and TotEng drops. RATTLE removes it during setup and keeps r₀₁·v₀₁ at round-off every step, so its kinetic energy and temperature are consistent with the constraint at every printed step.</p>'; }

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.cs-stb'); b.disabled = true; const rows = [];
      for (const dt of [0.5, 1.0, 2.0]) {
        st.textContent = `Running dt = ${dt} fs…`;
        const n = Math.round(1000 / dt), sc = SCRIPT('none').replace('DT', dt) + `\nthermo          ${Math.round(20 / dt)}\nrun             ${n}`;
        const res = await runLammps(sc, { files: { 'dia.data': DATA } }); if (res.error || !res.thermo.rows.length) { b.disabled = false; return err(st, res); }
        const c = res.thermo.cols, et = res.thermo.rows.filter(r => +r[0] * dt >= 20).map(r => +r[c.indexOf('TotEng')]), m = et.reduce((s, x) => s + x, 0) / et.length, sd = Math.sqrt(et.reduce((s, x) => s + (x - m) ** 2, 0) / et.length);
        rows.push(row([`harmonic bond, dt = ${dt} fs`, n, sd.toExponential(2), (Math.max(...et) - Math.min(...et)).toExponential(2)]));
      }
      const ls = NATIVE.ls.map(Number).slice(1), m = ls.reduce((s, x) => s + x, 0) / ls.length, sd = Math.sqrt(ls.reduce((s, x) => s + (x - m) ** 2, 0) / ls.length);
      rows.push(row(['SHAKE, dt = 2 fs (native)', 500, sd.toExponential(2), (Math.max(...ls) - Math.min(...ls)).toExponential(2)]));
      b.disabled = false; st.textContent = '';
      $('.cs-outb').innerHTML = `<table>${head(['', 'steps for 1 ps', 'σ(TotEng) kcal/mol', 'max − min'])}${rows.join('')}</table>`;
      $('.cs-sumb').innerHTML = '<p>With the bond as a spring, 2 fs is only 9 steps per vibration and TotEng swings by 0.06 kcal/mol; halving dt cuts the swing about 4 times (the dt² of Verlet), and even 0.5 fs, four times the work, stays far noisier than the constrained run at 2 fs, which only has to follow the slow motions. That is the whole case for SHAKE and RATTLE: freeze the fastest motion (X–H bonds, here N≡N) and take the step the slower motions allow.</p>';
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
