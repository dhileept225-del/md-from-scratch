// ML potentials from scratch · M0.3 What a potential must obey.
// Water model: SPC/Fw (Wu, Tepper, Voth, J. Chem. Phys. 124, 024503, 2006), LAMMPS form
// E = Kb (r - r0)^2 per bond + Ka (θ - θ0)^2, Kb = 529.581 kcal/mol/Å², r0 = 1.012 Å, Ka = 37.95 kcal/mol/rad², θ0 = 113.24°.
import { tex } from '../tex.js';

const Kb = 529.581, r0 = 1.012, Ka = 37.95, T0 = 113.24 * Math.PI / 180;
const COL = { ink: '#F2F5FA', muted: '#8FA3C7', cyan: '#38E1FF', amber: '#FFC24B', coral: '#FF6B6B', green: '#7EE787' };

// energy and forces, 2D. P = [O, H1, H2] as [x, y]. naive: bond 1 read from the raw x of H1.
function model(P, naive) {
  const [O, H1, H2] = P, F = [[0, 0], [0, 0], [0, 0]];
  const u = [H1[0] - O[0], H1[1] - O[1]], v = [H2[0] - O[0], H2[1] - O[1]];
  const ru = Math.hypot(...u), rv = Math.hypot(...v);
  const c = Math.max(-1, Math.min(1, (u[0] * v[0] + u[1] * v[1]) / (ru * rv))), th = Math.acos(c), s = Math.sqrt(1 - c * c) || 1e-12;
  let E = Kb * (rv - r0) ** 2 + Ka * (th - T0) ** 2;
  // bond 2
  const g2 = 2 * Kb * (rv - r0) / rv; F[2][0] -= g2 * v[0]; F[2][1] -= g2 * v[1]; F[0][0] += g2 * v[0]; F[0][1] += g2 * v[1];
  // bond 1
  if (naive) { E += Kb * (H1[0] - r0) ** 2; F[1][0] -= 2 * Kb * (H1[0] - r0); }
  else { E += Kb * (ru - r0) ** 2; const g1 = 2 * Kb * (ru - r0) / ru; F[1][0] -= g1 * u[0]; F[1][1] -= g1 * u[1]; F[0][0] += g1 * u[0]; F[0][1] += g1 * u[1]; }
  // angle: dθ/du = (c û - v̂) / (|u| s)
  const dE = 2 * Ka * (th - T0);
  const du = [(c * u[0] / ru - v[0] / rv) / (ru * s), (c * u[1] / ru - v[1] / rv) / (ru * s)];
  const dv = [(c * v[0] / rv - u[0] / ru) / (rv * s), (c * v[1] / rv - u[1] / ru) / (rv * s)];
  for (let k = 0; k < 2; k++) { F[1][k] -= dE * du[k]; F[2][k] -= dE * dv[k]; F[0][k] += dE * (du[k] + dv[k]); }
  return { E, F };
}

export default {
  id: 'invariances',
  minutes: 12,
  summary: 'Any potential, learned or not, has to give the same energy when you rotate the system, move it, or swap two identical atoms, and its forces must turn with the system. These rules are not optional extras: they are what makes momentum conserved. A model that reads raw coordinates breaks all three, which is why every ML potential in this course is built on distances and angles.',
  learn: [
    'The three symmetries: rotation, translation, swapping identical atoms',
    'Energy is invariant; forces are equivariant (they rotate with the system)',
    'Why the forces must add up to zero, and what breaks if they do not',
    'How to test any potential, including MACE, in a few lines',
  ],
  simHint: 'Drag the molecule, rotate it, or swap the two H labels: with the physical model the energy and force lengths do not move, and the arrows turn with the molecule. Stretch the bond to change the shape: now E changes, as it should. Then switch to the naive model and do the same.',

  idea: `
  <p>The energy of a molecule depends on its <b>shape</b>, nothing else. Three things leave the shape unchanged, so they must leave the energy unchanged:</p>
  <ol class="steps">
    <li><b>Translation:</b> move every atom by the same vector. Water in one corner of the box is the same as water in the other.</li>
    <li><b>Rotation:</b> turn the whole molecule. Space has no preferred direction.</li>
    <li><b>Permutation:</b> swap two atoms of the same element. The two hydrogens of water carry labels in your file, but nature does not know which is "H1".</li>
  </ol>
  <p>The energy is <b>invariant</b>: it does not change at all. The forces are different: they are arrows, so when the molecule turns, the arrows must turn with it by the same angle. That is called <b>equivariant</b>. Their lengths stay the same.</p>
  <p>These rules have physical consequences. Because the energy does not change when everything moves together, the forces on an isolated molecule add up to <b>zero</b>, so the total momentum is conserved. Because it does not change on rotation, the total torque is zero, so angular momentum is conserved. A model that breaks the symmetries makes an isolated molecule accelerate or spin up on its own.</p>
  <p>A classical force field obeys all three by construction: it is written in bond lengths, angles and distances, which do not care where the molecule is or how it is turned, and every hydrogen of water gets the same parameters. An ML model has no such guarantee. Feed a neural network the raw x, y, z of each atom, in file order, and it will learn a function that changes when you rotate the input. Parts M2 to M4 are about the inputs and architectures that make the rules hold exactly.</p>
  <div class="callout"><b>Key idea.</b> Same shape, same energy; turned shape, turned forces. Test it on every potential you use: it takes three lines (see Code).</div>`,

  hand: `
  <p>Flexible water (SPC/Fw): each O–H bond costs ${tex('K_b(r-r_0)^2')} and the angle costs ${tex('K_a(\\theta-\\theta_0)^2')}, with ${tex('K_b = 529.581')} kcal/mol/Å², ${tex('r_0 = 1.012')} Å, ${tex('K_a = 37.95')} kcal/mol/rad², ${tex('\\theta_0 = 113.24^\\circ')}.</p>
  <p>Put O at the origin, stretch H1 to 1.100 Å along x, and leave H2 at 1.012 Å and 113.24°. Only bond 1 is off its rest value:</p>
  <div class="calc">r₁ − r₀ = 1.100 − 1.012 = 0.088 Å
E = 529.581 × 0.088² = 529.581 × 0.007744 = 4.1011 kcal/mol
force on H1 = −2 K_b (r₁ − r₀) = −2 × 529.581 × 0.088 = −93.206 kcal/mol/Å, pointing at O
force on O  = +93.206 kcal/mol/Å, pointing at H1;  force on H2 = 0
sum of forces = −93.206 + 93.206 + 0 = 0</div>
  <p><b>Rotate by 90°.</b> H1 moves to (0, 1.100). Every distance and the angle are unchanged, so E = 4.1011 kcal/mol again. The force on H1 still has length 93.206 and still points at O, which is now straight down: the arrow turned by 90° with the molecule.</p>
  <p><b>Move by 1 Å, or swap the H labels.</b> Same distances, same angle, same 4.1011 kcal/mol.</p>
  <p><b>Now a naive model</b> that read H1's bond length as its raw x coordinate. It agrees at the start, because H1 sits on the x axis with O at the origin:</p>
  <div class="calc">start:      x_H1 = 1.100            E = 529.581 × (1.100 − 1.012)² = 4.101
rotate 90°: x_H1 = 0                E = 529.581 × (0 − 1.012)² = 529.581 × 1.024144 = 542.37
move 1 Å:   x_H1 = 2.100            E = 529.581 × (2.100 − 1.012)² = 529.581 × 1.183744 = 626.89
swap H:     x_H1 = 1.012 cos 113.24° = −0.39932
            E = 529.581 × (−1.41132)² + 4.1011 (the real bond, now "H2") = 1054.83 + 4.10 = 1058.93</div>
  <p>Same molecule, four energies. Press the button to compute all eight numbers with the page's own code.</p>
  <p><button class="btn sm" data-k="hand">Compute the table</button></p>
  <div class="ml-hand-out"></div>`,

  algorithm: `
  <p>An invariance test for any potential, given a structure:</p>
  <ol>
    <li>Compute ${tex('E_0')} and forces ${tex('F_0')}.</li>
    <li>Make a random rotation ${tex('R')}: QR-decompose a random 3×3 matrix, fix the signs, and flip one column if ${tex('\\det R = -1')} (that would be a mirror).</li>
    <li>Make a random translation ${tex('t')} and a permutation ${tex('P')} that swaps atoms of the same element.</li>
    <li>Transform: new positions ${tex('x\' = P(xR^{\\mathsf T} + t)')}. Compute ${tex('E_1, F_1')}.</li>
    <li>Check ${tex('|E_1 - E_0|')} and ${tex('\\max|F_1 - P F_0 R^{\\mathsf T}|')}: both should be at rounding level, about ${tex('10^{-15}')} in double precision, ${tex('10^{-6}')} in single.</li>
    <li>Check ${tex('|\\sum_i F_i|')}: zero for an isolated system. With periodic boundaries it is also zero; a nonzero sum means momentum drifts.</li>
  </ol>`,

  maths() {
    return `
  <p>For any rotation ${tex('R')}, translation ${tex('t')} and permutation ${tex('P')} of identical atoms, an exact potential satisfies</p>
  ${tex('E(\\{R\\,r_i + t\\}) = E(\\{r_i\\}), \\qquad E(\\{r_{P(i)}\\}) = E(\\{r_i\\}).', true)}
  <p><b>Forces rotate.</b> Write ${tex('r_i\' = R\\,r_i + t')}. Differentiate ${tex('E(\\{r\'\\}) = E(\\{r\\})')} with respect to ${tex('r_i')} using the chain rule, ${tex('\\partial r_i\' / \\partial r_i = R')}:</p>
  ${tex('R^{\\mathsf T}\\,\\nabla_{i}E(\\{r\'\\}) = \\nabla_i E(\\{r\\}) \\;\\Rightarrow\\; F_i(\\{r\'\\}) = R\\,F_i(\\{r\\}),', true)}
  <p>using ${tex('R R^{\\mathsf T} = 1')}. That is equivariance.</p>
  <p><b>Forces add to zero.</b> Take ${tex('R = 1')} and a small shift ${tex('t = \\varepsilon\\,\\hat n')}. Since ${tex('E')} does not change, ${tex('0 = \\frac{d}{d\\varepsilon}E = \\sum_i \\nabla_i E \\cdot \\hat n = -\\hat n \\cdot \\sum_i F_i')} for every direction ${tex('\\hat n')}, so ${tex('\\sum_i F_i = 0')}: total momentum is conserved.</p>
  <p><b>Torques add to zero.</b> A small rotation about ${tex('\\hat n')} moves each atom by ${tex('\\varepsilon\\,\\hat n \\times r_i')}. The same argument gives ${tex('\\hat n \\cdot \\sum_i r_i \\times F_i = 0')}, so the total torque vanishes and angular momentum is conserved (for an isolated system; periodic boxes break full rotation symmetry).</p>
  <p>Mirror images (${tex('\\det R = -1')}) are a separate question: most potentials are mirror-invariant too, but chiral molecules show that this is a choice, not a law.</p>`;
  },

  code: {
    note: 'Test a real ML potential: MACE-MP-0 on the stretched water above. Run here (Nov 2026): in float64, |ΔE| = 1.8×10⁻¹⁵ eV, force mismatch 7.1×10⁻¹⁵ eV/Å, |ΣF| = 4.4×10⁻¹⁶ eV/Å. In float32 the force mismatch is 5.0×10⁻⁶ eV/Å: rounding, not a broken model.',
    src: `# Invariance check: same energy, and rotated forces, after a rotation,
# a translation and a swap of the two hydrogens?
import numpy as np
from ase import Atoms
from mace.calculators import mace_mp

def random_rotation(rng):
    q, r = np.linalg.qr(rng.normal(size=(3, 3)))
    q *= np.sign(np.diag(r))                  # unique QR
    if np.linalg.det(q) < 0: q[:, 0] *= -1    # proper rotation, det = +1
    return q

def check(atoms, calc, rng):
    atoms.calc = calc
    E0, F0 = atoms.get_potential_energy(), atoms.get_forces()
    R, t = random_rotation(rng), rng.normal(size=3) * 5.0
    perm = [0, 2, 1]                          # swap H1 and H2
    b = atoms.copy(); b.positions = (atoms.positions @ R.T + t)[perm]; b.calc = calc
    E1, F1 = b.get_potential_energy(), b.get_forces()
    print(f'|dE| = {abs(E1 - E0):.1e} eV')
    print(f'max |F_after - R F_before| = {np.abs(F1 - (F0 @ R.T)[perm]).max():.1e} eV/A')
    print(f'|sum of forces| = {np.linalg.norm(F0.sum(0)):.1e} eV/A')

# water: one O-H bond stretched to 1.100 A, H-O-H = 113.24 degrees
th = np.radians(113.24)
water = Atoms('OH2', positions=[[0, 0, 0], [1.100, 0, 0],
                                [1.012 * np.cos(th), 1.012 * np.sin(th), 0]])
calc = mace_mp(model='small', default_dtype='float64', device='cpu')
check(water, calc, np.random.default_rng(0))`,
  },

  quiz: [
    { q: 'You rotate a water molecule by 90°. What happens?', options: ['Energy and forces both stay exactly the same', 'Energy stays the same; each force keeps its length and turns by 90°', 'Energy changes; forces stay the same', 'Both change unless the rotation is about the oxygen'], answer: 1, why: 'Energy is invariant; forces are equivariant: F(Rr) = R F(r). The rotation point does not matter, because a rotation about any other point is a rotation plus a translation.' },
    { q: 'Why must the forces on an isolated molecule add up to zero?', options: ['Because the atoms have equal masses', 'Because moving every atom by the same vector does not change the energy', 'Because of the cutoff', 'They do not have to'], answer: 1, why: 'Translation invariance gives ΣF = 0, so total momentum is conserved. A model that breaks it makes a molecule accelerate on its own.' },
    { q: 'A model uses the raw x coordinate of the atom labelled H1. Which operations change its energy?', options: ['Only rotation', 'Rotation and translation', 'Rotation, translation and swapping the hydrogens', 'None, if it was trained well'], answer: 2, why: 'At 90° it gives 542 kcal/mol, moved 1 Å 627, swapped 1059, instead of 4.10 each time. Training data cannot fix an input that carries the wrong information.' },
    { q: 'MACE in float32 gives a force mismatch of 5×10⁻⁶ eV/Å after a rotation. What does that mean?', options: ['The model is not equivariant', 'Single-precision rounding; float64 gives about 10⁻¹⁵', 'The rotation matrix was wrong', 'The structure needs relaxing first'], answer: 1, why: 'Float32 has about 7 significant digits. A real symmetry error is orders of magnitude larger, like the naive model above.' },
  ],

  mount(root) {
    const el = root.querySelector('#sim');
    el.classList.add('sim');
    el.innerHTML = `
    <div class="sim-view"><canvas style="touch-action:none;cursor:grab"></canvas></div>
    <div class="sim-ctrl">
      <div class="seg mini" style="margin-bottom:10px"><button data-m="0" class="on">Physical (distances, angle)</button><button data-m="1">Naive (raw x of H1)</button></div>
      <label>Rotate <output data-o="rot"></output><input type="range" data-k="rot" min="-180" max="180" step="1" value="0"></label>
      <label>H1–O bond <output data-o="r1"></output><input type="range" data-k="r1" min="0.95" max="1.15" step="0.001" value="1.1"></label>
      <div style="display:flex;gap:8px;flex-wrap:wrap;margin:6px 0 10px"><button class="btn sm ghost" data-b="swap">Swap H labels</button><button class="btn sm ghost" data-b="reset">Reset position</button></div>
      <div class="readouts">
        <div><span>Energy</span><b data-o="E"></b></div>
        <div><span>Same shape, start pose</span><b data-o="E0"></b></div>
        <div><span>|ΣF|</span><b data-o="sum"></b></div>
        <div><span>Force check</span><b data-o="fchk"></b></div>
      </div>
      <p class="muted small">Drag the molecule to move it. kcal/mol and kcal/mol/Å. Arrows: forces (cyan on O, amber on H), same scale; very long arrows are shortened and labelled.</p>
    </div>`;
    const cv = el.querySelector('canvas'), cx = cv.getContext('2d'), view = el.querySelector('.sim-view'), o = {};
    el.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
    const st = { naive: false, rot: 0, r1: 1.1, swap: false, sx: 0, sy: 0 };
    let W = 600, H = 420, sc = 110;
    // start pose: O at origin, H1 on +x, H2 at 113.24°
    const base = () => [[0, 0], [st.r1, 0], [r0 * Math.cos(T0), r0 * Math.sin(T0)]];
    const cen = P => [(P[0][0] + P[1][0] + P[2][0]) / 3, (P[0][1] + P[1][1] + P[2][1]) / 3];
    function pose() {
      const B = base(), c = cen(B), a = st.rot * Math.PI / 180, ca = Math.cos(a), sa = Math.sin(a);
      const P = B.map(([x, y]) => [ca * (x - c[0]) - sa * (y - c[1]) + c[0] + st.sx, sa * (x - c[0]) + ca * (y - c[1]) + c[1] + st.sy]);
      const lab = st.swap ? [P[0], P[2], P[1]] : P;     // order fed to the model: O, "H1", "H2"
      return { B, P, lab, ca, sa };
    }
    const f3 = x => x.toFixed(3);
    function draw() {
      const { B, P, lab, ca, sa } = pose();
      const m = model(lab, st.naive), ref = model(B, st.naive);
      // expected forces: start-pose forces, rotated, in the same atom order as the picture
      const rf = ref.F.map(([x, y]) => [ca * x - sa * y, sa * x + ca * y]);
      const got = st.swap ? [m.F[0], m.F[2], m.F[1]] : m.F;      // back to picture order
      const mis = Math.max(...got.map((f, i) => Math.hypot(f[0] - rf[i][0], f[1] - rf[i][1])));
      const sum = Math.hypot(m.F[0][0] + m.F[1][0] + m.F[2][0], m.F[0][1] + m.F[1][1] + m.F[2][1]);
      o.rot.textContent = `${st.rot}°`; o.r1.textContent = `${st.r1.toFixed(3)} Å`;
      o.E.textContent = f3(m.E); o.E0.textContent = f3(ref.E);
      const same = Math.abs(m.E - ref.E) < 1e-6 * Math.max(1, ref.E);
      o.E.style.color = same ? COL.green : COL.coral;
      o.sum.textContent = f3(sum); o.sum.style.color = sum < 1e-6 ? COL.green : COL.coral;
      o.fchk.textContent = mis < 1e-6 ? 'rotated ✓' : `off by ${mis.toFixed(1)}`; o.fchk.style.color = mis < 1e-6 ? COL.green : COL.coral;
      // picture (y up)
      cx.clearRect(0, 0, W, H);
      const X = p => W / 2 + (p[0] - 0.3) * sc, Y = p => H / 2 - (p[1] - 0.35) * sc;
      // axes cross at the world origin (the naive model's reference)
      cx.strokeStyle = 'rgba(143,163,199,.25)'; cx.lineWidth = 1; cx.setLineDash([4, 6]);
      cx.beginPath(); cx.moveTo(0, Y([0, 0])); cx.lineTo(W, Y([0, 0])); cx.moveTo(X([0, 0]), 0); cx.lineTo(X([0, 0]), H); cx.stroke(); cx.setLineDash([]);
      cx.fillStyle = COL.muted; cx.font = '12px system-ui'; cx.fillText('x', W - 14, Y([0, 0]) - 6); cx.fillText('origin', X([0, 0]) + 6, Y([0, 0]) + 16);
      cx.strokeStyle = '#8FA3C7'; cx.lineWidth = 6; cx.lineCap = 'round';
      [1, 2].forEach(i => { cx.beginPath(); cx.moveTo(X(P[0]), Y(P[0])); cx.lineTo(X(P[i]), Y(P[i])); cx.stroke(); });
      const fs = 0.9;   // px per kcal/mol/Å
      const cap = 0.32 * Math.min(W, H);
      const arrow = (p, f, c) => { const L0 = Math.hypot(...f) * fs; if (L0 < 2) return; const L = Math.min(L0, cap), k = L / L0, ox = -f[1] / L0 * fs * 4, oy = -f[0] / L0 * fs * 4, x0 = X(p) + ox, y0 = Y(p) + oy, x1 = x0 + f[0] * fs * k, y1 = y0 - f[1] * fs * k, a = Math.atan2(y1 - y0, x1 - x0);
        if (L0 > cap) { cx.fillStyle = c; cx.font = '600 12px system-ui'; cx.fillText(`${(L0 / fs).toFixed(0)}`, x1 + 6, y1 - 6); }
        cx.strokeStyle = c; cx.fillStyle = c; cx.lineWidth = 3; cx.beginPath(); cx.moveTo(x0, y0); cx.lineTo(x1, y1); cx.stroke();
        cx.beginPath(); cx.moveTo(x1, y1); cx.lineTo(x1 - 11 * Math.cos(a - 0.4), y1 - 11 * Math.sin(a - 0.4)); cx.lineTo(x1 - 11 * Math.cos(a + 0.4), y1 - 11 * Math.sin(a + 0.4)); cx.fill(); };
      const names = st.swap ? ['O', 'H2', 'H1'] : ['O', 'H1', 'H2'];
      P.forEach((p, i) => {
        cx.fillStyle = i ? COL.ink : COL.coral; cx.beginPath(); cx.arc(X(p), Y(p), i ? 15 : 22, 0, 7); cx.fill();
        cx.fillStyle = '#06101F'; cx.font = '600 12px system-ui'; cx.textAlign = 'center'; cx.fillText(names[i], X(p), Y(p) + 4); cx.textAlign = 'left';
      });
      P.forEach((p, i) => arrow(p, got[i], i ? COL.amber : COL.cyan));
      cx.fillStyle = COL.muted; cx.font = '12px system-ui';
      cx.fillText(st.naive ? 'naive: bond 1 = x of the atom labelled H1' : 'physical: SPC/Fw distances and angle', 10, 18);
    }
    const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; sc = Math.min(W, H) / 4.2; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
    const ro = new ResizeObserver(resize); ro.observe(view);
    el.querySelector('[data-k=rot]').oninput = e => { st.rot = +e.target.value; draw(); };
    el.querySelector('[data-k=r1]').oninput = e => { st.r1 = +e.target.value; draw(); };
    el.querySelectorAll('[data-m]').forEach(b => b.onclick = () => { st.naive = b.dataset.m === '1'; el.querySelectorAll('[data-m]').forEach(x => x.classList.toggle('on', x === b)); draw(); });
    el.querySelector('[data-b=swap]').onclick = () => { st.swap = !st.swap; draw(); };
    el.querySelector('[data-b=reset]').onclick = () => { st.sx = st.sy = 0; st.rot = 0; el.querySelector('[data-k=rot]').value = 0; draw(); };
    let drag = null;
    cv.onpointerdown = e => { drag = [e.clientX, e.clientY, st.sx, st.sy]; cv.setPointerCapture(e.pointerId); cv.style.cursor = 'grabbing'; };
    cv.onpointermove = e => { if (!drag) return; st.sx = drag[2] + (e.clientX - drag[0]) / sc; st.sy = drag[3] - (e.clientY - drag[1]) / sc; draw(); };
    cv.onpointerup = cv.onpointercancel = () => { drag = null; cv.style.cursor = 'grab'; };

    const hb = root.querySelector('[data-k=hand]'), out = root.querySelector('.ml-hand-out');
    if (hb) hb.onclick = () => {
      const A = [[0, 0], [1.1, 0], [r0 * Math.cos(T0), r0 * Math.sin(T0)]];
      const cases = [['start', A], ['rotate 90°', A.map(([x, y]) => [-y, x])], ['move 1 Å', A.map(([x, y]) => [x + 1, y])], ['swap H', [A[0], A[2], A[1]]]];
      out.innerHTML = `<div class="tbl"><table><tr><th>pose</th><th>physical E</th><th>naive E</th><th>physical |ΣF|</th></tr>${cases.map(([n, P]) => { const p = model(P, false), q = model(P, true);
        return `<tr><td>${n}</td><td>${p.E.toFixed(4)}</td><td>${q.E.toFixed(2)}</td><td>${Math.hypot(p.F[0][0] + p.F[1][0] + p.F[2][0], p.F[0][1] + p.F[1][1] + p.F[2][1]).toExponential(1)}</td></tr>`; }).join('')}</table></div>
        <p>Force on H1 at the start: ${model(A, false).F[1].map(v => v.toFixed(3)).join(', ')} kcal/mol/Å. Same as by hand: 4.1011, 542.37, 626.89, 1058.93.</p>`;
    };
    return () => ro.disconnect();
  },
};
