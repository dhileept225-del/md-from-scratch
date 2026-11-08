import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-beyond',
  lab: true,
  summary: 'Three last checks with numbers small enough to follow. First, a real MD code: LAMMPS run on the two atoms of Lab 1 must reproduce our hand calculation digit for digit, in any unit system. Second, a machine-learned potential in miniature: fit an energy curve to three "DFT" points by least squares, then see why ML potentials need flexible models and much more data. Third, metadynamics: add three hills by hand, then watch a long run fill a free-energy well and measure the barrier.',
  minutes: 35,
  sim: 'lab7',
  learn: [
    'Check a real code against a hand calculation, and convert between LAMMPS unit systems without mistakes',
    'Fit a potential by linear least squares by hand, and see interpolation, extrapolation and the role of data',
    'Add metadynamics hills by hand, and compare standard and well-tempered runs on a known barrier',
  ],
  simHint: 'ML fit: choose a model and the number of training points (amber dots); tap the plot to add a point where the error is large, like active learning. Metadynamics: press Run and watch the hills (grey) fill the left well until the walker escapes; the amber line is the free energy recovered from the bias. Compare standard and well-tempered: one keeps overshooting, the other settles at 5.',

  setup: `
<p>Three independent parts, each with its own tiny system:</p>
<div class="calc">steps 1–2  the two argon atoms of Lab 1 (r₀ = 1.05σ, at rest,
           dt = 0.01τ), now in LAMMPS
steps 3–5  a diatomic energy curve E(r), known exactly
           (a Morse curve, D = 1, a = 2, r_e = 1.12),
           pretending we only see a few points from "DFT"
steps 6–7  one coordinate s in a double well
           U(s) = 5 (s² − 1)²: minima at ±1, barrier 5 k_B T</div>
<p>Everything is checkable: the hand numbers of Lab 1, the exact Morse curve, and the exact barrier of 5.</p>`,

  steps: [
    {
      title: 'Lab 1 in LAMMPS: the same numbers',
      why: 'A real MD code does exactly what we did by hand, only for millions of atoms. Before trusting it with a research system, it is worth checking that it reproduces a calculation we understand completely. This also shows how a code reports its results.',
      work: `<div class="calc">units lj · atom_style atomic · 2 atoms at x = 0 and 1.05
pair_style lj/cut 2.5 · fix nve · timestep 0.01

total energy after 0, 1, 2, 3 steps:
by hand (Lab 1)   −0.757512  −0.757582  −0.757771  −0.758028
LAMMPS            −0.757512  −0.757582  −0.757771  −0.758028

LAMMPS default output in lj units:
                  −0.378756  −0.378791  −0.378885  −0.379014
                  (per atom! thermo_modify norm yes)</div>`,
      see: 'All six digits agree, because LAMMPS <code>fix nve</code> is the velocity Verlet we did by hand in Lab 1. The one trap: in <code>units lj</code> LAMMPS divides energies by the number of atoms in its thermo output by default, which gives exactly half of our numbers for 2 atoms. Use <code>thermo_modify norm no</code> when comparing totals. Checking a code on a tiny system with a known answer is the cheapest bug test there is, and it is the first step in the Research section on modifying LAMMPS.',
    },
    {
      title: 'Units: the same run in real and metal units',
      why: 'Research inputs use physical units: <code>real</code> (kcal/mol, Å, fs; for biomolecules and liquids) or <code>metal</code> (eV, Å, ps; for materials and most ML potentials). A unit slip in ε or in the time step silently gives a different simulation, so let us convert once by hand and check that LAMMPS agrees.',
      work: `<div class="calc">argon: ε/k_B = 119.8 K, σ = 3.405 Å, m = 39.948 g/mol
time unit τ = σ √(m/ε) = 2.1564 ps

real:   ε = 119.8 × 0.0019872 = 0.23807 kcal/mol
        dt = 0.01 τ = <b>21.56 fs</b>,  start at 1.05σ = 3.57525 Å
        etotal (kcal/mol): −0.180338 −0.180355 −0.180400 −0.180461
        divided by ε:      −0.757512 −0.757582 −0.757771 −0.758028

metal:  ε = 119.8 × 8.617333×10⁻⁵ = 0.010324 eV
        dt = <b>0.02156</b> (ps!)
        etotal / ε:        −0.757512 −0.757582 −0.757771 −0.758028</div>`,
      see: 'Three unit systems, one trajectory: dividing by ε gives back the reduced numbers exactly. The classic mistakes are all visible here: the time step is in fs for <code>real</code> but in ps for <code>metal</code> (writing 1 instead of 0.001 gives a 1000 times larger step), energies are in kcal/mol vs eV (a factor of 23.06), and k_B differs between them. Always check the <code>units</code> line of an input before reading any number from it.',
    },
    {
      title: 'A machine-learned potential in miniature: least squares by hand',
      why: 'An ML potential is a flexible function fitted to energies (and forces) from quantum chemistry. The simplest version is a linear model: energy = a sum of fixed functions of the geometry, each with an unknown coefficient. The coefficients come from least squares, a 2 × 2 system here.',
      work: `<div class="calc">three "DFT" energies (really the Morse curve):
r = 1.0  E = −0.9264
r = 1.2  E = −0.9781
r = 1.5  E = −0.7166

model  E(r) = A / r¹² − B / r⁶   (linear in A and B)
rows (1/r¹², −1/r⁶):
r = 1.0  (1.00000, −1.00000)
r = 1.2  (0.11216, −0.33490)
r = 1.5  (0.00771, −0.08779)

normal equations  XᵀX (A, B) = Xᵀ E:
[ 1.01264  −1.03824 ] [A]   [−1.04165]
[−1.03824   1.11986 ] [B] = [ 1.31691]
determinant = 0.056080
A = <b>3.5799</b>,  B = <b>4.4949</b>  (Cramer's rule)

fit at the training points:
r = 1.0  −0.9150 vs −0.9264
r = 1.2  −1.1038 vs −0.9781
r = 1.5  −0.3670 vs −0.7166   (error 0.35)</div>`,
      see: 'The whole "training" is two sums of squares and one 2 × 2 solve; real ML potentials do the same with thousands of coefficients and millions of energies and forces (MACE, NequIP and ACE are built on this idea, with smarter functions of the geometry). But this model has only the Lennard-Jones shape: with two coefficients it cannot follow the Morse curve, and it misses the third point by 0.35. A model that is too rigid cannot be fixed by fitting.',
    },
    {
      title: 'Interpolation, extrapolation and active learning',
      why: 'A fitted potential is only tested where it was trained. During a simulation, atoms visit geometries outside the training set, for example two atoms pushed closer than any training point. What happens there?',
      work: `<div class="calc">training range r = 1.0 to 1.5
r       true       fit       error
0.90   −0.6945    4.2174   <b>+4.91</b>   outside
0.95   −0.8360    0.5103   +1.35    outside
1.10   −0.9983   −1.3966   −0.40
1.30   −0.9086   −0.7776   +0.13
2.00   −0.3145   −0.0694   +0.25    outside
2.50   −0.1226   −0.0184   +0.10    outside

add the two worst points (r = 0.92, 2.2), refit:
A = 1.2768,  B = 2.5058
largest error on [0.92, 3.0]:  3.08  ->  <b>0.51</b>
but at r = 1.30 the error grows: +0.13 -> +0.44</div>`,
      see: 'Outside the training data the error explodes: at r = 0.90 the fit says +4.2 where the truth is −0.69, a fake wall or, with another model, a fake hole the atoms fall into. Adding the geometries where the model is worst ("active learning": run, find uncertain configurations, compute them with DFT, refit) reduces the worst error six-fold. The remaining error is the rigid model again: improving one region worsens another. Limitation for ML potentials: they are reliable only inside the space covered by their training set, so the training data must include the high-temperature, compressed and reactive geometries the simulation will visit.',
    },
    {
      title: 'Flexible models need data',
      why: 'The cure for a rigid model is more flexibility. But a flexible model can pass through every training point and still be wrong in between. We test both effects with two flexible models.',
      work: `<div class="calc">parabola c₀ + c₁r + c₂r² through the same 3 points:
training error 3×10⁻¹⁵ (exact)
r = 1.1   true −0.9983   parabola −0.9749
r = 1.3   true −0.9086   parabola −0.9362
r = 2.0   true −0.3145   parabola +0.6235
r = 3.0   true −0.0460   parabola +6.6946</div>
<div class="tbl"><table>${head(['K Gaussian basis functions', 'training points', 'max training error', 'max test error (0.95–2.0)'])}
${row(['3', '3', '0.0000', '0.8703'])}
${row(['8', '3', '0.0000', '<b>0.2642</b>'])}
${row(['8', '20', '0.0044', '<b>0.0044</b>'])}
${row(['15', '20', '0.0008', '0.0017'])}
${row(['15', '6', '0.0000', '0.0143'])}</table></div>`,
      see: 'Zero training error means nothing by itself: the parabola fits all three points exactly and then rises to +6.7 at r = 3, where the true energy goes to 0. With only 3 points, even 8 basis functions leave a test error of 0.26; with 20 points the same model is accurate to 0.004 everywhere in the range. The lesson behind every ML potential: always judge a model on <b>test</b> data it has not seen, the amount and spread of training data matters more than the model, and physics built into the model (the energy must go to zero at large r, smoothness, symmetry) is what makes it safe outside the data.',
    },
    {
      title: 'Metadynamics: three hills by hand',
      why: 'A barrier of 5 k_B T is crossed only rarely, so a plain simulation spends almost all its time in one well. Metadynamics adds small Gaussian hills at the places the system visits, along a chosen collective variable s. The hills fill the current well until the system has to leave it, and their sum, the bias, records the free energy.',
      work: `<div class="calc">U(s) = 5 (s² − 1)²,  k_B T = 1
hill: w exp(−(s − sᵢ)² / 2σ²),  w = 0.5,  σ = 0.2
three hills at sᵢ = −1.00, −0.95, −1.05

s = −1.00:  0.5000 + 0.4846 + 0.4846 = V 1.4692
            U + V = 0 + 1.4692 = 1.4692
s = −0.80:  0.3033 + 0.3774 + 0.2289 = V 0.9096
            U + V = 0.6480 + 0.9096 = 1.5576
s = −0.50:  V = 0.0731,  U + V = 2.8856
s =  0.00:  V = 0,       U + V = 5.0000

barrier seen from the left well: 5.0000 − 1.4692
                               = <b>3.5308</b>  (was 5.0)
escape factor  e^(−5) = 0.0067  ->  e^(−3.53) = 0.029</div>`,
      see: 'Three hills already lower the effective barrier from 5.0 to 3.53, which makes escape about 4 times more likely (0.029 vs 0.0067). Each new hill is added where the system is, so the well fills from the bottom up; once it is full, the walker spills over and starts filling the other well. At the end, the bias is a mirror image of the free energy: F(s) ≈ −V(s) + constant. The hills are only as good as the variable s: if s misses a slow motion (a "hidden" variable), filling along s does not help, and choosing the collective variable is the real difficulty of metadynamics.',
    },
    {
      title: 'A long run: standard vs well-tempered',
      why: 'Does the bias really converge to the free energy, and how fast? We run a Langevin walker in the double well for 10⁶ steps, adding a hill every 200 steps, and read the barrier from the bias as the run goes on. Well-tempered metadynamics makes each new hill smaller where the bias is already high, w × exp(−V/ΔT).',
      work: `<div class="tbl"><table>${head(['steps (×10⁵)', '1', '2', '3', '4', '6', '8', '9', '10'])}
${row(['standard: barrier', '7.59', '6.07', '5.99', '2.66', '4.20', '1.82', '5.49', '3.44'])}
${row(['well-tempered ΔT = 4: barrier', '4.85', '4.71', '4.81', '4.76', '4.73', '5.07', '5.22', '<b>5.00</b>'])}</table></div>
<div class="calc">exact barrier 5.0
well-tempered: F = −(1 + k_B T/ΔT) V = −1.25 V
crossings in 10⁶ steps:  standard 79,  well-tempered 54,
                         no bias <b>5</b></div>`,
      see: 'Both runs cross the barrier 10–16 times more often than the unbiased walker (5 crossings), which is the whole point. But standard metadynamics never settles: with full-height hills forever, the bias keeps overshooting, and the barrier estimate swings between 1.8 and 7.6. Well-tempered hills shrink as the bias grows, so the estimate settles at 5.00 ± 0.3. Limitations: hill height, width, stride and ΔT must be chosen (too big gives a rough, overfilled surface; too small converges slowly); the dynamics are biased, so rates need extra theory; and everything depends on a good collective variable. In practice this is done with PLUMED, coupled to LAMMPS or GROMACS.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['tool', 'merit (number)', 'limitation (number)'])}
${row(['LAMMPS vs hand', 'identical to 6 digits, in lj, real and metal', 'default per-atom output in lj units (factor 2 here)'])}
${row(['unit systems', 'one trajectory, three unit sets', 'dt in fs (real) vs ps (metal): ×1000 traps'])}
${row(['rigid fitted model', '2 coefficients, one 2 × 2 solve', 'training error 0.35; +4.9 outside the data'])}
${row(['active learning', 'worst error 3.08 -> 0.51', 'cannot fix a rigid model form'])}
${row(['flexible ML model', 'test error 0.004 with 20 points', '0.26 with 3 points; parabola +6.7 at r = 3'])}
${row(['metadynamics', '54–79 crossings vs 5 unbiased', 'standard: barrier swings 1.8–7.6'])}
${row(['well-tempered', 'barrier 5.00 (exact 5.0)', 'needs a good collective variable'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the hand Verlet steps and the same run in LAMMPS (lj, real and metal units), the least-squares fit and its tests, the flexible-model table, the three hills and the long metadynamics runs. Run <code>pip install numpy lammps</code>, then <code>python lab_beyond.py</code> (about 10 s). Without LAMMPS installed, steps 1–2 are skipped.',
    src: `# Lab 7.4: beyond by hand. (1) the two atoms of Lab 1 in LAMMPS, (2) a tiny machine-learned potential fitted by least
# squares, (3) metadynamics: three hills by hand, then a long run. Needs: pip install numpy lammps
import numpy as np, math

# ---------- Steps 1-2: the same two atoms in LAMMPS ----------
def hand_verlet(n=3, r0=1.05, dt=0.01):
    F = lambda r: 24 * (2 * r**-13 - r**-7); U = lambda r: 4 * (r**-12 - r**-6)
    x = np.array([0.0, r0]); v = np.zeros(2); f = F(r0) * np.array([-1, 1]); E = [U(r0)]
    for _ in range(n):
        v += 0.5 * dt * f; x += dt * v; r = x[1] - x[0]; f = F(r) * np.array([-1, 1]); v += 0.5 * dt * f; E.append(U(r) + 0.5 * (v**2).sum())
    return E
def lammps_run(units, sig, eps, mass, dt, n=3, norm='no'):
    from lammps import lammps
    l = lammps(cmdargs=['-log', 'none', '-screen', 'none']); L = 10 * sig
    l.commands_string(f"""
units {units}
atom_style atomic
boundary f f f
region box block {-L} {L} {-L} {L} {-L} {L}
create_box 1 box
create_atoms 1 single 0 0 0
create_atoms 1 single {1.05 * sig} 0 0
mass 1 {mass}
pair_style lj/cut {2.5 * sig}
pair_coeff 1 1 {eps} {sig}
timestep {dt}
fix 1 all nve
thermo_style custom step pe ke etotal
thermo_modify norm {norm}
run 0""")
    out = [l.get_thermo('etotal')]
    for _ in range(n): l.command('run 1'); out.append(l.get_thermo('etotal'))
    l.close(); return out
print("Step 1: total energy after 0, 1, 2, 3 steps (r0 = 1.05, dt = 0.01, reduced units)")
print("   by hand (Lab 1):      " + "  ".join(f"{e:.6f}" for e in hand_verlet()))
try:
    E = lammps_run('lj', 1.0, 1.0, 1.0, 0.01); print("   LAMMPS, units lj:     " + "  ".join(f"{e:.6f}" for e in E))
    En = lammps_run('lj', 1.0, 1.0, 1.0, 0.01, norm='yes'); print("   LAMMPS default output (thermo_modify norm yes in lj units): " + "  ".join(f"{e:.6f}" for e in En) + "   <- per atom!")
    kB = 0.0019872; eps = 119.8 * kB; tau = 3.405e-10 * math.sqrt(39.948 * 1.66054e-27 / (119.8 * 1.380649e-23))
    print(f"\\nStep 2: argon in LAMMPS 'real' units (kcal/mol, Angstrom, fs): eps = 119.8 K x {kB} = {eps:.5f} kcal/mol, sigma = 3.405 A")
    print(f"   time unit tau = sigma sqrt(m/eps) = {tau * 1e12:.4f} ps, so dt = 0.01 tau = {0.01 * tau * 1e15:.2f} fs;  start at 1.05 sigma = {1.05 * 3.405:.5f} A")
    Er = lammps_run('real', 3.405, eps, 39.948, 0.01 * tau * 1e15)
    print("   LAMMPS etotal (kcal/mol): " + "  ".join(f"{e:.6f}" for e in Er))
    print("   divided by eps:           " + "  ".join(f"{e / eps:.6f}" for e in Er))
    kBe = 8.617333e-5; Em = lammps_run('metal', 3.405, 119.8 * kBe, 39.948, 0.01 * tau * 1e12)
    print(f"   'metal' units (eV, ps): eps = {119.8 * kBe:.6f} eV, dt = {0.01 * tau * 1e12:.5f} ps;  etotal / eps: " + "  ".join(f"{e / (119.8 * kBe):.6f}" for e in Em))
except ImportError:
    print("   (install the LAMMPS Python module to compare: pip install lammps)")

# ---------- Steps 3-5: a tiny machine-learned potential ----------
D, a, re = 1.0, 2.0, 1.12
true = lambda r: D * (1 - np.exp(-a * (r - re)))**2 - D                     # the "quantum chemistry" we pretend not to know
rt = np.array([1.0, 1.2, 1.5]); Et = true(rt)
print("\\nStep 3: three reference energies (pretend they came from DFT):")
for r_, e_ in zip(rt, Et): print(f"   r = {r_:.1f}: E = {e_:.4f}")
X = np.c_[rt**-12, -rt**-6]                                                  # model E = A r^-12 - B r^-6, linear in A and B
print("   model E(r) = A / r^12 - B / r^6;  design matrix rows (1/r^12, -1/r^6):")
for r_, row in zip(rt, X): print(f"   r = {r_:.1f}: ({row[0]:.5f}, {row[1]:.5f})")
XtX = X.T @ X; Xty = X.T @ Et
print(f"   normal equations  [{XtX[0, 0]:.5f} {XtX[0, 1]:.5f}] [A]   [{Xty[0]:.5f}]")
print(f"                     [{XtX[1, 0]:.5f} {XtX[1, 1]:.5f}] [B] = [{Xty[1]:.5f}]")
det = XtX[0, 0] * XtX[1, 1] - XtX[0, 1] * XtX[1, 0]
A = (Xty[0] * XtX[1, 1] - XtX[0, 1] * Xty[1]) / det; B = (XtX[0, 0] * Xty[1] - XtX[1, 0] * Xty[0]) / det
print(f"   determinant = {det:.6f};  A = {A:.4f},  B = {B:.4f}  (Cramer's rule)")
fit = lambda r, A=A, B=B: A * r**-12 - B * r**-6
print("   check at the training points: " + ", ".join(f"r {r_}: fit {fit(r_):.4f} vs {e_:.4f}" for r_, e_ in zip(rt, Et)))
print("\\nStep 4: test inside and outside the training range [1.0, 1.5]")
for r_ in (0.90, 0.95, 1.1, 1.3, 2.0, 2.5, 3.0):
    print(f"   r = {r_:.2f}: true {true(r_):8.4f}   fit {fit(r_):8.4f}   error {fit(r_) - true(r_):+8.4f}" + ("   (outside)" if r_ < 1.0 or r_ > 1.5 else ""))
rt2 = np.r_[rt, 0.92, 2.2]; X2 = np.c_[rt2**-12, -rt2**-6]; A2, B2 = np.linalg.lstsq(X2, true(rt2), rcond=None)[0]
print(f"   add the two worst points (r = 0.92, 2.2) and refit: A = {A2:.4f}, B = {B2:.4f}")
for r_ in (0.90, 0.95, 1.3, 2.5):
    print(f"   r = {r_:.2f}: true {true(r_):8.4f}   new fit {fit(r_, A2, B2):8.4f}   error {fit(r_, A2, B2) - true(r_):+8.4f}")
grid = np.linspace(0.92, 3.0, 300)
print(f"   largest error on [0.92, 3.0]: first fit {np.abs(fit(grid) - true(grid)).max():.3f}, after adding points {np.abs(fit(grid, A2, B2) - true(grid)).max():.3f}  (the model form itself limits it)")
c = np.polyfit(rt, Et, 2)
print(f"\\nStep 5: a 'flexible' model, a parabola c0 + c1 r + c2 r^2 through the same 3 points: c = {c[2]:.3f}, {c[1]:.3f}, {c[0]:.3f}")
for r_ in (1.1, 1.3, 2.0, 3.0): print(f"   r = {r_:.1f}: true {true(r_):8.4f}   parabola {np.polyval(c, r_):8.4f}")
print(f"   training error: parabola {np.abs(np.polyval(c, rt) - Et).max():.1e}, LJ-form {np.abs(fit(rt) - Et).max():.3f}")
print("   a model like real ML potentials: E = sum of K Gaussian basis functions of r (centres spread over 0.9-2.2), linear least squares")
test = np.linspace(0.95, 2.0, 200)
def rbf(K, rtrain, lam=1e-8):
    mu = np.linspace(0.9, 2.2, K); wd = (2.2 - 0.9) / K
    Phi = lambda r: np.exp(-(r[:, None] - mu[None])**2 / (2 * wd**2))
    P = Phi(rtrain); cc = np.linalg.solve(P.T @ P + lam * np.eye(K), P.T @ true(rtrain))
    return np.abs(P @ cc - true(rtrain)).max(), np.abs(Phi(test) @ cc - true(test)).max()
print("   K basis   training points   max training error   max test error (0.95-2.0)")
for K, n in ((3, 3), (8, 3), (8, 20), (15, 20), (15, 6)):
    tr = np.linspace(0.95, 2.0, n) if n > 3 else rt
    e1, e2 = rbf(K, tr); print(f"   {K:7d}   {n:15d}   {e1:18.4f}   {e2:14.4f}")

# ---------- Steps 6-7: metadynamics ----------
U = lambda s: (s**2 - 1)**2 * 5.0                                            # double well, barrier 5 (energy units), minima at +-1
dU = lambda s: 20.0 * s * (s**2 - 1)
kT = 1.0; w, sg = 0.5, 0.2
hills = [-1.0, -0.95, -1.05]
V = lambda s, hl=hills, ws=None: sum((ws[i] if ws else w) * np.exp(-(s - h)**2 / (2 * sg**2)) for i, h in enumerate(hl))
print(f"\\nStep 6: double well U(s) = 5 (s^2 - 1)^2, barrier {U(0.0):.1f} = {U(0.0) / kT:.0f} kT. Hills: height w = {w}, width sigma = {sg}, at s = {hills}")
for s_ in (-1.0, -0.8, -0.5, 0.0):
    terms = [w * np.exp(-(s_ - h)**2 / (2 * sg**2)) for h in hills]
    print(f"   s = {s_:5.2f}: hills " + " + ".join(f"{t_:.4f}" for t_ in terms) + f" = V {sum(terms):.4f};  U + V = {U(s_):.4f} + {sum(terms):.4f} = {U(s_) + sum(terms):.4f}")
print(f"   barrier seen from the left well: {U(0.0) + V(0.0) - (U(-1.0) + V(-1.0)):.4f} (was {U(0.0) - U(-1.0):.1f})")
print(f"   unbiased escape estimate: rate ~ exp(-5) = {math.exp(-5):.4f} per attempt; each hill lowers the barrier by about w = {w}")

def metad(nsteps=1000000, dt=0.001, gamma=5.0, stride=200, wt=None, seed=0):
    r = np.random.default_rng(seed); s = -1.0; v = 0.0; c = np.exp(-gamma * dt); sn = np.sqrt((1 - c * c) * kT); noise = r.normal(size=nsteps)
    grid = np.linspace(-2.5, 2.5, 1001); dg = grid[1] - grid[0]; Vg = np.zeros_like(grid); dVg = np.zeros_like(grid)
    cross = 0; side = -1; out = []; nh = 0
    for k in range(nsteps):
        i = min(max(int((s + 2.5) / dg), 0), 999); dV = dVg[i] + (dVg[i + 1] - dVg[i]) * ((s + 2.5) / dg - i)
        v += -0.5 * dt * (dU(s) + dV); s += dt * v
        i = min(max(int((s + 2.5) / dg), 0), 999); dV = dVg[i] + (dVg[i + 1] - dVg[i]) * ((s + 2.5) / dg - i)
        v += -0.5 * dt * (dU(s) + dV); v = c * v + sn * noise[k]
        if s * side < 0 and abs(s) > 0.8: cross += 1; side = -side
        if (k + 1) % stride == 0:
            ww = w if wt is None else w * np.exp(-float(np.interp(s, grid, Vg)) / wt); g = np.exp(-(grid - s)**2 / (2 * sg**2)); nh += 1
            Vg += ww * g; dVg += -ww * (grid - s) / sg**2 * g
        if (k + 1) % 100000 == 0:
            F = -Vg if wt is None else -(1 + kT / wt) * Vg
            out.append((k + 1, nh, F[np.argmin(np.abs(grid))] - F[(grid > -1.5) & (grid < -0.5)].min(), cross))
    return out
print("\\nStep 7: long metadynamics run (Langevin, kT = 1, dt = 0.001, a hill every 200 steps). Barrier F(0) - F(left minimum) estimated from the bias (exact 5.0):")
for name, wt in (("standard", None), ("well-tempered, Delta T = 4 kT", 4.0)):
    print(f"   {name}:")
    for k, nh, bar, cr in metad(wt=wt):
        print(f"     after {k:6d} steps, {nh:4d} hills: barrier {bar:6.2f}, crossings so far {cr}")
r = np.random.default_rng(1); s = -1.0; v = 0.0; dt = 0.001; c = np.exp(-5 * dt); sn = np.sqrt((1 - c * c) * kT); cross = 0; side = -1
for k in range(1000000):
    v += -0.5 * dt * dU(s); s += dt * v; v += -0.5 * dt * dU(s); v = c * v + sn * r.normal()
    if s * side < 0 and abs(s) > 0.8: cross += 1; side = -side
print(f"   without bias, same length: {cross} crossings")
`,
  },

  exercises: [
    { q: 'An input uses units metal with timestep 0.002. What is the time step in fs? What would the same number mean in units real?', a: 'metal: time in ps, so 0.002 ps = 2 fs, a normal step. real: time in fs, so 0.002 fs, 1000 times smaller: the run would be correct but would move about 1000 times too slowly through time.' },
    { q: 'Fit E(r) = −B / r⁶ (a single coefficient) to the three "DFT" points by least squares. What is B? (Hint: B = Σ(xE) / Σ(x²) with x = −1/r⁶.)', a: 'x = −1, −0.33490, −0.08779. Σ xE = 0.9264 + 0.32757 + 0.06291 = 1.31688. Σ x² = 1 + 0.11216 + 0.00771 = 1.11987. B = 1.1759. At r = 1.0 the model gives −1.176 vs −0.926: an attraction-only model cannot describe the repulsive wall.' },
    { q: 'After many well-tempered hills with ΔT = 4 k_B T, the bias at the barrier top is 4.0 higher than at the bottom of the left well. What barrier do you report?', a: 'F = −(1 + k_B T/ΔT) V = −1.25 V, so the barrier is 1.25 × 4.0 = 5.0 k_B T. Reporting 4.0 would underestimate it by 20 %.' },
  ],
};
