// Research · A.8: read_data files. One water molecule written by hand, read by LAMMPS, every energy and the dipole checked.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

const DATA = `One water molecule, written by hand

3 atoms
2 bonds
1 angles

2 atom types
1 bond types
1 angle types

-10.0 10.0 xlo xhi
-10.0 10.0 ylo yhi
-10.0 10.0 zlo zhi

Masses

1 15.9994
2 1.008

Atoms # full

1 1 1 -0.8476  0.000000  0.000000  0.000000
2 1 2  0.4238  1.020000  0.000000  0.000000
3 1 2  0.4238 -0.255732  0.987421  0.000000

Bonds

1 1 1 2
2 1 1 3

Angles

1 1 2 1 3
`;
const KB = 450.0, R0 = 1.0, KA = 55.0, TH0 = 109.47, QQR2E = 332.06371;
const SCRIPT = special => `units           real
atom_style      full
${special ? 'special_bonds   lj/coul 0.0 1.0 0.0\n' : ''}read_data       water.data
bond_style      harmonic
bond_coeff      1 ${KB.toFixed(1)} ${R0.toFixed(1)}
angle_style     harmonic
angle_coeff     1 ${KA.toFixed(1)} ${TH0}
pair_style      lj/cut/coul/cut 10.0
pair_coeff      1 1 0.1553 3.166
pair_coeff      2 2 0.0 0.0
variable        mux equal q[1]*x[1]+q[2]*x[2]+q[3]*x[3]
variable        muy equal q[1]*y[1]+q[2]*y[2]+q[3]*y[3]
thermo_style    custom step ebond eangle evdwl ecoul v_mux v_muy
thermo_modify   format float %.8g
run             0
`;

function parse(text) {
  const atoms = {}; const L = text.split('\n'); const k = L.findIndex(l => /^Atoms/.test(l.trim()));
  if (k < 0) return null;
  for (let i = k + 1; i < L.length; i++) { const w = L[i].trim().split(/\s+/); if (!L[i].trim()) { if (Object.keys(atoms).length) break; continue; } if (w.length < 7 || isNaN(+w[0])) break; atoms[+w[0]] = { type: +w[2], q: +w[3], r: [+w[4], +w[5], +w[6]] }; }
  return atoms;
}
function hand(atoms, special) {
  const O = atoms[1], H1 = atoms[2], H2 = atoms[3]; if (!O || !H1 || !H2) return null;
  const d = (a, b) => a.r.map((x, i) => x - b.r[i]), len = v => Math.hypot(...v);
  const v1 = d(H1, O), v2 = d(H2, O), r1 = len(v1), r2 = len(v2);
  const th = Math.acos((v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2]) / (r1 * r2)), dth = th - TH0 * Math.PI / 180;
  const Eb = KB * ((r1 - R0) ** 2 + (r2 - R0) ** 2), Ea = KA * dth * dth;
  const rHH = len(d(H1, H2)), Ec = special ? QQR2E * H1.q * H2.q / rHH : 0;
  const mu = [0, 1].map(k => O.q * O.r[k] + H1.q * H1.r[k] + H2.q * H2.r[k]);
  return { r1, r2, th: th * 180 / Math.PI, dth, Eb, Ea, rHH, Ec, mu };
}

export default {
  id: 'read-data',
  title: 'read_data files',
  minutes: 20,
  summary: 'Real systems are rarely built with lattice commands: they are read from a data file that lists every atom, bond and angle. You write one by hand for a single water molecule, let LAMMPS read it, and check the bond and angle energies, the excluded Coulomb energy and the dipole moment by hand, then edit the file and check again.',
  learn: [
    'The layout of a data file: header, Masses, Atoms, Bonds, Angles, and the columns of atom_style full',
    'Predict bond, angle and Coulomb energies and the dipole of a water molecule by hand',
    'Why the Coulomb energy inside a molecule is zero by default (special_bonds)',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: a data file is the system, written out</h2>
    <p>Instead of <code>lattice</code>, <code>region</code>, <code>create_box</code> and <code>create_atoms</code>, block 2 of the script (page A.3) can be one line: <code>read_data water.data</code>. The file has a header with counts and the box, then sections:</p>
    <div class="tbl wide"><table>${head(['Part', 'Contains', 'Water'])}
${row(['header', 'numbers of atoms, bonds, angles, and of their types; the box', '3 atoms, 2 bonds, 1 angle; 2 atom types'])}
${row(['<code>Masses</code>', 'one mass per atom type', 'O 15.9994, H 1.008 g/mol'])}
${row(['<code>Atoms # full</code>', 'one line per atom; the columns depend on atom_style', 'id, molecule, type, charge, x, y, z'])}
${row(['<code>Bonds</code>', 'bond id, bond type, the two atom ids', 'O–H1, O–H2'])}
${row(['<code>Angles</code>', 'angle id, angle type, three atom ids (the middle one is the vertex)', 'H1–O–H2'])}</table></div>
    <div class="tbl wide"><table>${head(['atom_style', 'columns of an Atoms line'])}
${row(['atomic', 'id type x y z'])}
${row(['charge', 'id type q x y z'])}
${row(['full', 'id molecule type q x y z'])}</table></div>
    <p>Three optional integers at the end of an Atoms line are image flags (page D.2).</p>
  </section>

  <section class="block">
    <h2>2 · Modify: the data file and the script</h2>
    <p class="muted">Edit any number in the data file (for example a coordinate). The hand calculation below updates at once; then run LAMMPS.</p>
    <div class="rd-grid">
      <div><label class="lab">water.data</label><textarea class="lmp-in rd-data" spellcheck="false" wrap="off"></textarea></div>
      <div><label class="lab">the input script</label><pre class="lmp-mini rd-script"></pre>
        <label class="rd-sp-lab"><input type="checkbox" class="rd-sp"> count the H–H Coulomb energy (<code>special_bonds lj/coul 0.0 1.0 0.0</code>)</label></div>
    </div>
  </section>

  <section class="block prose">
    <h2>3 · Predict by hand</h2>
    <p>LAMMPS's harmonic bond and angle are E = K(r − r₀)² and E = K(θ − θ₀)²: there is <b>no ½</b> in front, unlike the textbook ½k x². Force fields that publish k must be divided by 2 before use.</p>
    <div class="calc rd-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Read it with LAMMPS</button><span class="lmp-status muted rd-st"></span></div>
    <div class="lmp-compare rd-cmp"></div>
  </section>

  <section class="block prose">
    <h2>What the molecule shows</h2>
    <ul>
      <li>Every energy is a few lines of arithmetic once the coordinates are known; the file is the whole system.</li>
      <li>The Coulomb energy is <b>zero</b> by default even though the atoms carry charges: atoms joined by a bond (1–2) or by an angle (1–3) are excluded (<code>special_bonds</code> 0 0 0). The two H atoms, 1.61 Å apart, would otherwise add about +37 kcal/mol, far more than the bond and angle energies. Each force field states its own exclusion rule; the data file and the script must follow it.</li>
      <li>The dipole μ = Σ q r comes out as 0.529 e·Å = 2.54 D for this geometry (1 e·Å = 4.8032 D).</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You change the header to <code>4 atoms</code> but leave three Atoms lines. What happens?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>LAMMPS reads the counts first and reads exactly that many Atoms lines, so the count and the file disagree. This version stops with <code>ERROR: Did not assign all atoms correctly</code>. Counts in the header must match the sections.</div></div>
    <div class="ex"><p><b>2.</b> A paper gives the O–H bond as k = 1059.2 kcal/mol/Å² in E = ½k(r − r₀)². What do you write in bond_coeff?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>K = k/2 = <b>529.6</b>, because LAMMPS uses E = K(r − r₀)² with no ½. Writing 1059.2 makes the bond twice as stiff.</div></div>
    <div class="ex"><p><b>3.</b> Move H2 so that the angle is exactly 109.47°. What is E_angle, and does the dipole get larger or smaller?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>E_angle = 0. Opening the angle moves the two H atoms apart, so their charges point less in the same direction: the dipole gets smaller. Try x = −0.340, y = 0.962 for H2 (r ≈ 1.02, θ ≈ 109.47°) in the editor.</div></div>
  </section>`,

  mount(root) {
    const ta = root.querySelector('.rd-data'), sp = root.querySelector('.rd-sp'), sc = root.querySelector('.rd-script');
    ta.value = DATA;
    let H = null;
    function draw() {
      sc.textContent = SCRIPT(sp.checked);
      const a = parse(ta.value); H = a && hand(a, sp.checked);
      if (!H) { root.querySelector('.rd-hand').textContent = 'The page needs atoms 1 (O), 2 and 3 (H) in the Atoms section.'; return; }
      root.querySelector('.rd-hand').innerHTML = esc(`bond lengths   r₁ = |H1 − O| = ${g8(H.r1)} Å      r₂ = |H2 − O| = ${g8(H.r2)} Å
angle          θ = ${g8(H.th)}°,  θ − θ₀ = ${g8(H.dth)} rad
E_bond  = 450 × [(r₁ − 1)² + (r₂ − 1)²]          = `) + `<b>${g8(H.Eb)}</b>` + esc(`
E_angle = 55 × (θ − θ₀)²                          = `) + `<b>${g8(H.Ea)}</b>` + esc(`
E_vdwl  = 0   (only O has LJ, and there is one O)
E_coul  = `) + (sp.checked ? esc(`332.06371 × qH × qH / r_HH = 332.06371 × 0.4238² / ${g8(H.rHH)} = `) + `<b>${g8(H.Ec)}</b>` : `<b>0</b>` + esc('   (all pairs are bonded or share an angle: excluded)')) + esc(`
dipole  μx = Σ q x = `) + `<b>${g8(H.mu[0])}</b>` + esc(`  μy = Σ q y = `) + `<b>${g8(H.mu[1])}</b>` + esc(` e·Å
        |μ| = ${g8(Math.hypot(...H.mu))} e·Å = ${g8(Math.hypot(...H.mu) * 4.8032047)} D`);
      root.querySelector('.rd-cmp').innerHTML = '';
    }
    ta.oninput = draw; sp.onchange = draw; draw();
    root.querySelector('[data-k=run]').onclick = async e => {
      if (!H) return; const b = e.currentTarget, st = root.querySelector('.rd-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(SCRIPT(sp.checked), { files: { 'water.data': ta.value } }); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = res.lines.filter(l => /^\s+\d+ (atoms|bonds|angles)$/.test(l)).map(l => l.trim()).join(' · ') + ' read';
      const r = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols;
      const want = [['E_bond', H.Eb], ['E_angle', H.Ea], ['E_vdwl', 0], ['E_coul', H.Ec], ['v_mux', H.mu[0]], ['v_muy', H.mu[1]]];
      root.querySelector('.rd-cmp').innerHTML = `<div class="tbl"><table><tr><th></th><th>by hand</th><th>LAMMPS</th><th></th></tr>${want.map(([k, h]) => { const lv = r[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k.replace('v_mu', 'dipole ')}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
