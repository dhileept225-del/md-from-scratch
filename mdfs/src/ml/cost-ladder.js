// ML potentials from scratch · M0.1 The cost ladder: force field, DFT, ML.
// Timings measured on one 2-core cloud machine (Nov 2026) for 32 water molecules (96 atoms) in a 9.865 Å box.
import { tex } from '../tex.js';

const COL = { cyan: '#38E1FF', coral: '#FF6B6B', amber: '#FFC24B', ink: '#F2F5FA', muted: '#8FA3C7', green: '#7EE787', violet: '#B79CFF' };
// seconds per MD step (energy + forces), 2 CPU cores
const T = {
  ff: 0.79008e-3,        // LAMMPS, flexible TIP3P-like water, PPPM 1e-5, 10 000 steps: 7.9009 s
  mace: 0.70391,         // MACE-MP-0 "small" (3.85 M parameters), PyTorch on CPU, float64, best of 5 calls
  maceM: 2.53687,        // MACE-MP-0 "medium" (4.69 M parameters)
  xtb: 2.71825,          // GFN2-xTB (tblite), periodic box, one energy + gradient
};
// DFT: PBE/def2-SVP energy + gradient (PySCF), gas-phase clusters of the first n molecules, SCF from scratch
const DFT = [[1, 3, 1.7715], [2, 6, 1.3801], [4, 12, 4.6277], [8, 24, 29.894], [16, 48, 224.83]];   // [waters, atoms, s]
const kDFT = Math.log(DFT[4][2] / DFT[3][2]) / Math.log(2);          // exponent from the two largest clusters
const tDFT96 = DFT[4][2] * 2 ** kDFT;                                // extrapolated to 32 waters
const STEPS_NS = 1e-9 / 1e-15;                                       // 1 ns at a 1 fs time step
const fmtT = s => s < 1e-3 ? `${(s * 1e6).toFixed(0)} µs` : s < 1 ? `${(s * 1e3).toFixed(s < 0.01 ? 2 : 0)} ms` : s < 120 ? `${s.toFixed(s < 10 ? 2 : 0)} s` : s < 7200 ? `${(s / 60).toFixed(0)} min` : s < 2 * 86400 ? `${(s / 3600).toFixed(1)} h` : s < 2 * 365.25 * 86400 ? `${(s / 86400).toFixed(s < 30 * 86400 ? 1 : 0)} days` : `${(s / (365.25 * 86400)).toFixed(0)} years`;
const N0 = 96;
const SUP = '⁰¹²³⁴⁵⁶⁷⁸⁹', sci = (x, d = 3) => { const [m, e] = x.toExponential(d - 1).split('e'); const n = +e; return n === 0 ? m : `${m} × 10${n < 0 ? '⁻' : ''}${String(Math.abs(n)).split('').map(c => SUP[c]).join('')}`; };
const METHODS = [
  { k: 'ff', name: 'Force field', col: COL.green, t: T.ff, p: 1, note: 'LAMMPS, fixed charges, PPPM' },
  { k: 'mace', name: 'ML (MACE-MP-0 small)', col: COL.cyan, t: T.mace, p: 1, note: 'PyTorch, CPU' },
  { k: 'xtb', name: 'Tight binding (GFN2-xTB)', col: COL.amber, t: T.xtb, p: 2, note: 'tblite' },
  { k: 'dft', name: 'DFT (PBE/def2-SVP)', col: COL.coral, t: tDFT96, p: kDFT, note: 'PySCF, extrapolated' },
];
const cost = (m, n) => m.t * (n / N0) ** m.p;

export default {
  id: 'cost-ladder',
  minutes: 12,
  summary: 'Why would anyone fit a potential with machine learning? Because the methods that are accurate enough are too slow, and the methods that are fast enough are not accurate. This lesson puts numbers on that ladder: one MD step for 96 atoms of water, timed on the same two-core machine with a force field, an ML potential, tight binding and DFT. Then it turns them into the wall time for one nanosecond, by hand.',
  learn: [
    'What one MD step costs with a force field, an ML potential, xTB and DFT',
    'Turning seconds per step into days per nanosecond, by hand',
    'How cost grows with the number of atoms: N for force fields and ML, about N³ for DFT',
    'Where an ML potential sits on the ladder, and what it must be trained to reproduce',
  ],
  simHint: 'Move the sliders: the number of atoms and the length of the simulation. The bars show wall time on this two-core machine; the dots on the curve are the measured points.',

  idea: `
<p>Every MD step needs the forces on all atoms. How you get them decides everything else.</p>
<ul>
<li><b>Force field</b> (lessons 4.3, 5.1): a fixed formula — springs, Lennard-Jones, point charges. Fast, but it only knows the physics someone wrote into it. It cannot break bonds, and its charges do not respond.</li>
<li><b>DFT</b> (density functional theory): solve for the electrons at every step. It breaks bonds, polarizes, transfers charge, but each step costs minutes to hours.</li>
<li><b>Tight binding</b> (xTB): a cheap, approximate quantum method in between.</li>
<li><b>ML potential</b>: a flexible function, trained on DFT energies and forces, that maps atomic positions to the same numbers. Once trained it costs far less than DFT, and it has no electrons.</li>
</ul>
<p>Measured here, for 32 water molecules (96 atoms) on the same two cores: force field <b>${fmtT(T.ff)}</b> per step, MACE-MP-0 small <b>${fmtT(T.mace)}</b>, GFN2-xTB <b>${fmtT(T.xtb)}</b>, DFT about <b>${fmtT(tDFT96)}</b>. Each level of the ladder costs roughly a thousand times more than the one below, and ML lands near the bottom of the quantum part of the ladder. Its accuracy is whatever its training data had, and that is what the rest of this section is about.</p>
<p>Two honest caveats. A GPU makes the ML potential 10–100× faster, and production DFT-MD codes (CP2K, VASP) on a few hundred cores are far faster than PySCF on two. The ratios between the levels survive both.</p>`,

  hand: `
<p>One nanosecond at a 1 fs time step is</p>
${tex(String.raw`N_\text{steps} = \frac{1\ \text{ns}}{1\ \text{fs}} = \frac{10^{-9}}{10^{-15}} = 10^{6}\ \text{steps}`, true)}
<p>Wall time = steps × seconds per step. With the measured seconds per step for 96 atoms:</p>
<div class="tbl"><table><tr><th>method</th><th>s per step</th><th>× 10⁶ steps</th><th>wall time for 1 ns</th></tr>
<tr><td>Force field</td><td>7.90 × 10⁻⁴</td><td>790 s</td><td>13 min</td></tr>
<tr><td>MACE-MP-0 small</td><td>0.704</td><td>7.04 × 10⁵ s</td><td>8.1 days</td></tr>
<tr><td>GFN2-xTB</td><td>2.72</td><td>2.72 × 10⁶ s</td><td>31 days</td></tr>
<tr><td>DFT</td><td>≈ ${tDFT96.toFixed(0)}</td><td>≈ ${sci(tDFT96 * 1e6)} s</td><td>≈ ${(tDFT96 * 1e6 / (365.25 * 86400)).toFixed(0)} years</td></tr></table></div>
<p>Divide by 86 400 s per day: 7.04 × 10⁵ / 86 400 = 8.15 days. By 3.156 × 10⁷ s per year: DFT ≈ ${(tDFT96 * 1e6 / 3.15576e7).toFixed(1)} years.</p>
<p><b>Where the DFT number comes from.</b> 32 waters were too slow to time here, so the smaller clusters were timed and the trend extended. From 8 to 16 waters the time went from 29.9 s to 224.8 s: a factor 7.52 for twice the atoms, so</p>
${tex(String.raw`k = \frac{\ln(224.8/29.9)}{\ln 2} = \frac{\ln 7.52}{0.693} = ${kDFT.toFixed(2)}, \qquad t_{32} = 224.8 \times 2^{${kDFT.toFixed(2)}} = ${tDFT96.toFixed(0)}\ \text{s}`, true)}
<div class="row"><button class="btn sm primary" data-k="hand">▶ Redo the table from the measured seconds</button></div>
<div class="ml-hand-out"></div>`,

  algorithm: `
<pre class="lmp-mini">for step in 1 … N_steps:                 # 10⁶ for 1 ns at 1 fs
    v ← v + ½ Δt F/m                   # lesson 1.5, the same for every method
    x ← x + Δt v
    F ← forces(x)                      # ← the only line that differs, and nearly all the cost
    v ← v + ½ Δt F/m

forces(x):
    force field : sum of pair, bond and Coulomb terms           ~ N  (with neighbour lists, PPPM ~ N log N)
    ML          : descriptors of each atom's neighbourhood → network → E; F = −∂E/∂x   ~ N
    xTB, DFT    : solve for the electrons self-consistently, then Hellmann–Feynman forces  ~ N² … N³</pre>
<p>The integrator is identical; only <code>forces(x)</code> changes. An ML potential is a drop-in replacement for that one function, which is why every MD code (LAMMPS, GROMACS, ASE, OpenMM) can use one.</p>`,

  maths: () => `
<h4>1. Wall time</h4>
${tex(String.raw`t_\text{wall} = \frac{t_\text{sim}}{\Delta t}\; t_\text{step}(N)`, true)}
<h4>2. Scaling with system size</h4>
<p>Write the cost per step as a power law, ${tex(String.raw`t_\text{step}(N) = t_0 (N/N_0)^k`)}. Short-ranged methods (force fields with neighbour lists, ML potentials with a cutoff) have k = 1: double the atoms, double the work. Conventional DFT diagonalizes a matrix whose size grows with the number of electrons, k ≈ 3; the measured exponent between 24 and 48 atoms was ${kDFT.toFixed(2)}. From two measurements,</p>
${tex(String.raw`k = \frac{\ln\left(t_2/t_1\right)}{\ln\left(N_2/N_1\right)}`, true)}
<h4>3. Why the gap widens</h4>
<p>The ratio of DFT to ML cost grows as ${tex(String.raw`N^{k_\text{DFT}-1} \approx N^{2}`)}: for 96 atoms it is about ${(tDFT96 / T.mace).toFixed(0)}, for 960 atoms a hundred times more. This is the regime ML potentials were invented for: DFT-like forces for systems and times DFT cannot reach.</p>`,

  code: {
    note: 'The timing script (abridged). It needs LAMMPS, PyTorch, mace-torch, PySCF and tblite; on another machine your seconds will differ, the ratios much less.',
    src: `import time
from ase.io import read
from mace.calculators import mace_mp
from tblite.ase import TBLite
from pyscf import gto, dft

water = read("water32.xyz")            # 32 H2O in a 9.865 A box, 96 atoms

def per_step(atoms, calc, n=5):
    atoms = atoms.copy(); atoms.calc = calc; atoms.get_forces()   # warm-up
    best = 1e9
    for _ in range(n):
        atoms.positions += 0.001                                    # force a new evaluation
        t = time.time(); atoms.get_forces(); best = min(best, time.time() - t)
    return best

t_mace = per_step(water, mace_mp(model="small", device="cpu", default_dtype="float64"))
t_xtb = per_step(water, TBLite(method="GFN2-xTB", verbosity=0), n=1)

def t_dft(nwat):                       # gas-phase cluster of the first nwat molecules
    at = water[: 3 * nwat]
    mol = gto.M(atom=[[s, tuple(x)] for s, x in zip(at.get_chemical_symbols(), at.positions)],
                basis="def2-svp", verbose=0)
    mf = dft.RKS(mol); mf.xc = "pbe"
    t = time.time(); mf.kernel(); mf.nuc_grad_method().kernel(); return time.time() - t

t8, t16 = t_dft(8), t_dft(16)
import math
k = math.log(t16 / t8) / math.log(2)
t_dft32 = t16 * 2 ** k

t_ff = 7.9009 / 10000                  # LAMMPS "Loop time" for 10 000 steps, same 96 atoms
for name, t in [("force field", t_ff), ("MACE small", t_mace), ("GFN2-xTB", t_xtb), ("DFT", t_dft32)]:
    wall = t * 1e6                     # 1 ns at 1 fs
    print(f"{name:12s} {t:10.4g} s/step   1 ns: {wall/86400:10.4g} days")`,
  },

  quiz: [
    { q: 'An ML potential takes 0.5 s per step on your laptop. How long is 100 ps at a 0.5 fs time step?', options: ['50 s', 'About 1 hour', 'About 28 hours', 'About 10 days'], answer: 2, why: '100 ps / 0.5 fs = 2 × 10⁵ steps; × 0.5 s = 10⁵ s = 27.8 hours.' },
    { q: 'You double the number of atoms. Roughly how much longer does a DFT step take, and an ML step?', options: ['Both twice as long', 'DFT about 8×, ML about 2×', 'DFT about 2×, ML about 8×', 'Both about 8×'], answer: 1, why: 'DFT cost grows about as N³ (2³ = 8; measured here 7.5×), a short-ranged ML potential as N.' },
    { q: 'What does an ML potential replace in an MD code?', options: ['The integrator', 'The thermostat', 'The force (and energy) evaluation', 'The neighbour list'], answer: 2, why: 'The integrator, thermostats and analysis are unchanged; only forces(x) is a trained model instead of a formula or a DFT calculation.' },
    { q: 'Why is an ML potential not "free DFT"?', options: ['It is slower than DFT', 'It is only as good as its training data, and can fail where it has none', 'It cannot compute forces', 'It needs a GPU'], answer: 1, why: 'It reproduces DFT where it was trained. Outside that region it can be badly wrong without warning — parts M5 and the labs show how to detect that.' },
  ],

  mount(root) {
    const el = root.querySelector('#sim');
    el.classList.add('sim');
    el.innerHTML = `
    <div class="sim-view"><canvas></canvas></div>
    <div class="sim-ctrl">
      <label>Atoms <output data-o="n"></output><input type="range" data-k="n" min="1" max="5" step="0.01" value="1.982"></label>
      <label>Simulated time <output data-o="t"></output><input type="range" data-k="t" min="-3" max="3" step="1" value="0"></label>
      <div class="readouts">${METHODS.map(m => `<div><span style="color:${m.col}">${m.name}</span><b data-o="${m.k}"></b></div>`).join('')}
        <div class="wide"><span>DFT ÷ ML</span><b data-o="ratio"></b></div></div>
      <p class="muted small">Time step 1 fs. Force field and ML scale as N, xTB as N², DFT as N^${kDFT.toFixed(2)} (measured). Two CPU cores.</p>
    </div>`;
    const cv = el.querySelector('canvas'), cx = cv.getContext('2d'), view = el.querySelector('.sim-view'), o = {};
    el.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
    const st = { n: 96, t: 1e-9 };
    el.querySelector('[data-k=n]').oninput = e => { st.n = Math.round(10 ** +e.target.value); draw(); };
    el.querySelector('[data-k=t]').oninput = e => { st.t = 1e-9 * 10 ** +e.target.value; draw(); };
    let W = 600, H = 420;
    const resize = () => { const d = Math.min(devicePixelRatio || 1, 2); W = view.clientWidth; H = view.clientHeight; cv.width = W * d; cv.height = H * d; cx.setTransform(d, 0, 0, d, 0, 0); draw(); };
    const ro = new ResizeObserver(resize); ro.observe(view);
    const tl = t => t >= 1e-6 ? `${t * 1e6} µs` : t >= 1e-9 ? `${t * 1e9} ns` : `${t * 1e12} ps`;
    function draw() {
      const steps = st.t / 1e-15;
      o.n.textContent = st.n.toLocaleString('en'); o.t.textContent = tl(st.t);
      METHODS.forEach(m => { o[m.k].textContent = fmtT(cost(m, st.n) * steps); });
      o.ratio.textContent = sci(cost(METHODS[3], st.n) / cost(METHODS[1], st.n), 2);
      // log-log plot: seconds per step vs atoms
      cx.clearRect(0, 0, W, H);
      const x0 = 56, x1 = W - 14, y0 = H - 40, y1 = 16, lx = n => x0 + (Math.log10(n) - 0) / 5 * (x1 - x0), ly = s => y0 - (Math.log10(s) + 5) / 12 * (y0 - y1);
      cx.strokeStyle = 'rgba(143,163,199,.25)'; cx.fillStyle = COL.muted; cx.font = '12px system-ui'; cx.lineWidth = 1;
      for (let e = 0; e <= 5; e++) { const X = lx(10 ** e); cx.beginPath(); cx.moveTo(X, y0); cx.lineTo(X, y1); cx.stroke(); cx.fillText(`10${'⁰¹²³⁴⁵'[e]}`, X - 8, y0 + 16); }
      for (let e = -5; e <= 7; e += 2) { const Y = ly(10 ** e); cx.beginPath(); cx.moveTo(x0, Y); cx.lineTo(x1, Y); cx.stroke(); cx.fillText(fmtT(10 ** e), 4, Y + 4); }
      cx.fillText('atoms', x1 - 36, y0 + 30); cx.fillText('wall time per step', x0 + 6, y1 + 12);
      for (const m of METHODS) {
        cx.strokeStyle = m.col; cx.lineWidth = 2; cx.beginPath();
        for (let e = 0; e <= 5; e += 0.05) { const s = cost(m, 10 ** e), X = lx(10 ** e), Y = ly(s); if (Y < y1 - 2) break; e ? cx.lineTo(X, Y) : cx.moveTo(X, Y); }
        cx.stroke();
      }
      cx.fillStyle = COL.coral; DFT.forEach(([, n, s]) => { cx.beginPath(); cx.arc(lx(n), ly(s), 4, 0, 7); cx.fill(); });
      [[METHODS[0], T.ff], [METHODS[1], T.mace], [METHODS[2], T.xtb]].forEach(([m, s]) => { cx.fillStyle = m.col; cx.beginPath(); cx.arc(lx(N0), ly(s), 4, 0, 7); cx.fill(); });
      cx.strokeStyle = COL.ink; cx.setLineDash([4, 4]); cx.beginPath(); cx.moveTo(lx(st.n), y0); cx.lineTo(lx(st.n), y1); cx.stroke(); cx.setLineDash([]);
    }
    resize();

    const b = root.querySelector('[data-k=hand]'), out = root.querySelector('.ml-hand-out');
    if (b) b.onclick = () => {
      const rows = METHODS.map(m => { const w = m.t * STEPS_NS; return `<tr><td>${m.name}</td><td>${sci(m.t, 4)}</td><td>${sci(w, 4)} s</td><td>${fmtT(w)}</td></tr>`; }).join('');
      out.innerHTML = `<div class="tbl"><table><tr><th>method</th><th>s per step</th><th>s for 10⁶ steps</th><th>wall time</th></tr>${rows}</table></div><p>Same numbers as the hand table, from the same measured seconds. MACE-MP-0 "medium" (${T.maceM.toFixed(2)} s per step) would take ${fmtT(T.maceM * STEPS_NS)}: the model size is part of the price.</p>`;
    };
    return () => ro.disconnect();
  },
};
