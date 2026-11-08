import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-dielectric',
  lab: true,
  summary: 'The dielectric constant ε says how strongly a liquid screens charges, and it controls ion pairing, solvation and conductivity. In a simulation it comes from the fluctuations of the total dipole of the box, and how the molecules rotate decides how fast those fluctuations relax. We compute the box dipole and ε by hand, measure rotational correlation functions for small-step and jump-like rotation, and see why ε is one of the slowest quantities to converge.',
  minutes: 35,
  sim: 'lab14',
  learn: [
    'Compute the box dipole and ε from its fluctuations, and read the Kirkwood factor G_K',
    'Compute P1 and P2 rotational correlation functions, and tell diffusive rotation from jumps by τ₁/τ₂',
    'See why the box dipole relaxes on a different time than one molecule, and why ε converges slowly',
  ],
  simHint: 'Dipoles on a lattice rotate by small random steps; neighbours prefer to be parallel (J > 0) or antiparallel (J < 0). Watch the box dipole M wander and the Kirkwood factor G_K settle: around 1 for J = 0, much larger for J > 0, smaller for J < 0. The lower plot compares how fast one dipole and the whole box forget their direction.',

  setup: `
<p>Every polar molecule carries a dipole μᵢ. Their sum is the box dipole M = Σ μᵢ. With Ewald sums and conducting ("tin-foil") boundaries, the default in most codes:</p>
<div class="calc">ε = 1 + (⟨M²⟩ − ⟨M⟩²) / (3 ε₀ V k_B T)

Kirkwood factor  G_K = ⟨M²⟩ / (N μ²)
                 = 1 for independent dipoles
                 > 1 if neighbours point the same way
rotation:  C₁(t) = ⟨u(0)·u(t)⟩,   C₂(t) = ⟨P₂(u(0)·u(t))⟩,  P₂(x) = (3x² − 1)/2</div>
<p>To keep the physics visible we use a model: 216 unit dipoles on a 6 × 6 × 6 lattice (periodic), each rotating by small random steps (rotational diffusion D_r = 1), with a coupling J between nearest neighbours (k_BT = 1). Real dipolar liquids also need the long-range dipole–dipole interaction, which Ewald handles.</p>`,

  steps: [
    {
      title: 'The box dipole and ε by hand',
      why: 'ε is built from M, so first compute M for a few molecules, then put real numbers for water into the fluctuation formula. Comparing with independent dipoles shows how much the molecules\' cooperation matters.',
      work: `<div class="calc">three dipoles (D): (1.2, 0.8, −1.5), (−0.4, 2.1, 0.9), (0.7, −1.0, 1.9)
M = (1.5, 1.9, 1.3) D;  |M|² = 7.55 D²;  Σ|μ|² = 14.81 D²

water-like: μ = 2.35 D = 7.839×10⁻³⁰ C m, n = 33.4 nm⁻³, T = 300 K
independent dipoles (G_K = 1):
ε − 1 = n μ² / (3 ε₀ k_B T)
      = 3.340×10²⁸ × 6.145×10⁻⁵⁹ / 1.1002×10⁻³¹ = <b>18.65</b>
ε = <b>19.7</b>

water models give ε ≈ 70  ->  G_K = 69 / 18.65 = <b>3.70</b></div>`,
      see: 'Three dipoles pointing in different directions partly cancel (|M|² = 7.55 against Σ|μ|² = 14.81). If water molecules pointed independently, ε would be only about 20; the real value near 70–80 needs G_K ≈ 3.7: hydrogen bonds make neighbouring dipoles point the same way, so their fluctuations add up. This is why ε is a sensitive test of a water model: a small change in how molecules orient together changes ε a lot. Limitation: the formula assumes conducting boundaries; with a reaction field or vacuum boundary it changes form, and for ionic systems the translational current also contributes to M and has to be separated.',
    },
    {
      title: 'The Kirkwood factor in the lattice model',
      why: 'With the model we can switch the neighbour coupling on and off and watch G_K respond. For J = 0 the answer must be exactly 1.',
      work: `<div class="tbl"><table>${head(['coupling J', 'G_K = ⟨M²⟩/N', 'water-like ε = 1 + 18.65 G_K'])}
${row(['0 (independent)', '0.968 ± 0.040', '19.1 ± 0.7'])}
${row(['+0.4 (parallel)', '<b>3.104 ± 0.155</b>', '58.9 ± 2.9'])}
${row(['−0.4 (antiparallel)', '0.514 ± 0.015', '10.6 ± 0.3'])}</table></div>
<div class="calc">200 000 steps each, errors from 10 blocks;  ⟨M⟩²/N ≤ 0.016</div>`,
      see: 'Independent dipoles give G_K = 0.97 ± 0.04, consistent with 1. A modest preference for parallel neighbours (J = 0.4 k_BT) triples G_K and would raise ε from 19 to 59; antiparallel neighbours halve it. ε measures correlations between molecules, not the dipole of one molecule. That is also why ε is hard to fit: small errors in the orientational correlations are multiplied by 18.65 (for water) in ε.',
    },
    {
      title: 'Rotational correlation functions: P1, P2, diffusion or jumps',
      why: 'How does one molecule change direction? C₁ (from the first Legendre polynomial) is what dielectric spectroscopy sees; C₂ is what NMR, Raman and fluorescence anisotropy see. For small random rotations theory says τ₁/τ₂ = 3; for big jumps to new directions, τ₁/τ₂ = 1. So the ratio tells the mechanism.',
      work: `<div class="calc">one molecule's axis, t = 0..3 ps:
t   u(0)·u(t) = cos θ    P₂ = (3cos²θ − 1)/2
0   1.000                1.000
1   0.933                0.806
2   0.640                0.114
3   0.800                0.460</div>
<div class="tbl"><table>${head(['rotation', 'C₁ (t = 0.1, 0.2, 0.5)', 'C₂ (t = 0.1, 0.2, 0.5)', 'τ₁', 'τ₂', 'τ₁/τ₂'])}
${row(['small steps, D_r = 1', '0.820 0.673 0.371', '0.552 0.304 0.050', '0.508', '0.169', '<b>3.01</b>'])}
${row(['jumps, 2 per time unit', '0.817 0.667 0.365', '0.816 0.666 0.361', '0.500', '0.486', '<b>1.03</b>'])}</table></div>
<div class="calc">theory: diffusion τ₁ = 1/(2D_r) = 0.500, τ₂ = 1/(6D_r) = 0.167
        jumps τ₁ = τ₂ = 1/rate = 0.500</div>`,
      see: 'P₂ falls faster than cos θ for small steps, because squaring makes it sensitive to small tilts (0.933 → 0.806 already at t = 1). Both mechanisms give the same τ₁ = 0.5, so C₁ alone cannot tell them apart; the ratio τ₁/τ₂ can: 3.01 for diffusion, 1.03 for jumps, exactly as theory predicts. In water the ratio is below 3, a signature of the large-angle jumps by which molecules swap hydrogen-bond partners. Limitation: real molecules are not spheres; a non-spherical molecule has different rates about different axes, and a single τ hides that.',
    },
    {
      title: 'One molecule vs the whole box',
      why: 'Dielectric spectroscopy measures how fast the box dipole M relaxes, not one molecule. Are they the same? Compare the single-dipole τ₁ with τ_M, the area of ⟨M(0)·M(t)⟩/⟨M²⟩.',
      work: `<div class="tbl"><table>${head(['coupling J', 'single τ₁', 'box τ_M', 'τ_M / τ₁', 'G_K'])}
${row(['0', '0.488', '0.440', '0.90', '0.97'])}
${row(['+0.4', '0.573', '<b>1.539</b>', '<b>2.68</b>', '3.10'])}
${row(['−0.4', '0.558', '0.259', '0.46', '0.51'])}</table></div>`,
      see: 'For independent dipoles the box forgets at the same rate as one molecule (0.44 vs 0.49, equal within noise). With parallel coupling the box dipole lives 2.7 times longer than any single molecule\'s orientation: a fluctuation of M involves many aligned molecules and decays slowly. With antiparallel coupling it decays faster. τ_M / τ₁ follows G_K (2.68 vs 3.10, 0.46 vs 0.51). So a dielectric relaxation time is not a molecular rotation time; comparing the two needs a theory that includes these correlations. For ionic liquids the translational current adds its own contribution to the measured spectrum.',
    },
    {
      title: 'Convergence: why ε needs very long runs',
      why: 'G_K is built from M, one vector for the whole box, relaxing on the slow time τ_M. Like the conductivity (Lab 6.7) and the viscosity (Lab 6.10) it gets no help from the number of molecules. How precise is it after a given run length?',
      work: `<div class="tbl"><table>${head(['run length (steps)', 'time units', 'in units of τ_M', 'spread of G_K', 'spread of water-like ε'])}
${row(['10 000', '20', '13', '0.62 (20 %)', '± 11.6'])}
${row(['20 000', '40', '26', '0.49 (16 %)', '± 9.1'])}
${row(['50 000', '100', '65', '0.39 (12 %)', '± 7.2'])}
${row(['200 000', '400', '260', '0.16 (5 %, blocks)', '± 2.9'])}</table></div>`,
      see: 'Even 65 relaxation times of M give G_K only to ±12 %, which would be ±7 in a water-like ε. Getting ε to a few percent takes hundreds of τ_M; for water that means tens of nanoseconds, and for viscous ionic liquids, where M relaxes in nanoseconds, microseconds or many independent runs. A quoted ε without a run length and an error bar is hard to interpret. Bigger boxes do not help, because M is one number per frame whatever N is.',
    },
    {
      title: 'Practical notes',
      why: 'A few choices in real codes change the formula or the number, and are easy to overlook.',
      work: `<div class="calc">boundary condition       formula
tin-foil Ewald (default)  ε = 1 + ⟨δM²⟩/(3ε₀Vk_BT)
reaction field ε_RF       (ε − 1)(2ε_RF + 1)/(2ε_RF + ε) = ⟨δM²⟩/(3ε₀Vk_BT)
vacuum boundary           (ε − 1)/(ε + 2) = ⟨δM²⟩/(9ε₀Vk_BT)

δM² = M² − ⟨M⟩²: subtract the mean (it should be ≈ 0)
molecules split across the box wall: use unwrapped molecules
ions: M also contains the ionic current; separate the
rotational part or use the frequency-dependent formalism</div>`,
      see: 'With the same ⟨δM²⟩ the three boundary conditions give different ε, so the formula must match the electrostatics of the run. A nonzero ⟨M⟩ usually signals an error (molecules broken across the periodic boundary, or a slowly drifting system) and must be subtracted. GROMACS: <code>gmx dipoles</code> (which also prints G_K); LAMMPS: <code>compute dipole</code> with your own post-processing.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['quantity', 'merit (number)', 'limitation (number)'])}
${row(['ε from ⟨δM²⟩', 'water: G_K = 3.7 explains ε ≈ 70 vs 19.7', 'formula depends on the boundary condition'])}
${row(['Kirkwood G_K', 'measures orientational correlation (0.51–3.10)', 'errors × 18.65 in a water-like ε'])}
${row(['C₁ and C₂', 'τ₁/τ₂ = 3.01 (diffusion) vs 1.03 (jumps)', 'non-spherical molecules need more than one τ'])}
${row(['box dipole relaxation', 'τ_M/τ₁ = 2.68 when neighbours align', 'not a single-molecule rotation time'])}
${row(['convergence', '±5 % after ~260 τ_M', '±12 % after 65 τ_M; box size does not help'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the box dipole and the water-like ε by hand, the rotational correlation functions for small steps and jumps, the Kirkwood factor and relaxation times for three couplings, and the convergence table. Run <code>pip install numpy</code>, then <code>python lab_dielectric.py</code> (about a minute).',
    src: `# Lab: dielectric constant and rotation by hand. Box dipole and the fluctuation formula, the Kirkwood factor, rotational
# correlation functions P1 and P2, single vs collective relaxation, convergence. Model: unit dipoles on a cubic lattice
# (6 x 6 x 6, periodic) with nearest-neighbour coupling J (kT = 1), rotational Brownian dynamics.
import numpy as np, time

# ---------- Step 1: the box dipole and the fluctuation formula ----------
mu = np.array([[1.2, 0.8, -1.5], [-0.4, 2.1, 0.9], [0.7, -1.0, 1.9]])        # three molecular dipoles, debye
M = mu.sum(0)
print("Step 1: three molecular dipoles (D):", mu.tolist())
print(f"   box dipole M = sum = ({M[0]:.1f}, {M[1]:.1f}, {M[2]:.1f}) D;  |M|^2 = {M @ M:.2f} D^2;  sum of |mu|^2 = {(mu**2).sum():.2f} D^2")
e0, kB, D2Cm = 8.8541878e-12, 1.380649e-23, 3.33564e-30
T, n, mw = 300.0, 33.4e27, 2.35 * D2Cm
y3 = n * mw**2 / (3 * e0 * kB * T)
print(f"   water-like: mu = 2.35 D = {mw:.3e} C m, n = 33.4 /nm^3, T = 300 K")
print(f"   independent dipoles: eps - 1 = n mu^2 / (3 eps0 kB T) = {n:.3e} x {mw**2:.3e} / ({3 * e0 * kB * T:.4e}) = {y3:.2f}  ->  eps = {1 + y3:.1f}")
print(f"   simulated SPC/E-like water eps ~ 70  ->  G_K = <M^2>/(N mu^2) = (70 - 1)/{y3:.2f} = {69 / y3:.2f}")

# ---------- model ----------
nl = 6; N = nl**3
idx = np.arange(N).reshape(nl, nl, nl)
nbr = np.stack([np.roll(idx, s, a).ravel() for a in range(3) for s in (1, -1)], 1)      # 6 neighbours each
def run(J, nsteps, Dr=1.0, dt=0.002, seed=0, jump=0.0, every=10):
    r = np.random.default_rng(seed); u = r.normal(size=(N, 3)); u /= np.linalg.norm(u, axis=1)[:, None]
    Ms, U = [], []
    for s in range(-5000, nsteps):
        if jump:                                                # jump model: with probability jump*dt a dipole picks a new random direction
            m = r.random(N) < jump * dt; k = m.sum()
            if k: w = r.normal(size=(k, 3)); u[m] = w / np.linalg.norm(w, axis=1)[:, None]
        else:
            h = J * u[nbr].sum(1); tq = h - (h * u).sum(1)[:, None] * u
            xi = r.normal(size=(N, 3)); xi -= (xi * u).sum(1)[:, None] * u
            u = u + Dr * tq * dt + np.sqrt(2 * Dr * dt) * xi; u /= np.linalg.norm(u, axis=1)[:, None]
        if s >= 0 and s % every == 0: Ms.append(u.sum(0)); U.append(u.copy())
    return np.array(Ms), np.array(U)
def corr(U, kmax):                                             # single-dipole P1 and P2 correlation, all origins
    c1, c2 = [], []
    for k in range(kmax):
        c = (U[:len(U) - k] * U[k:]).sum(-1); c1.append(c.mean()); c2.append((1.5 * c**2 - 0.5).mean())
    return np.array(c1), np.array(c2)
def acfM(Ms, kmax):
    return np.array([(Ms[:len(Ms) - k] * Ms[k:]).sum(-1).mean() for k in range(kmax)])

# ---------- Step 3: rotation by hand, then free rotational diffusion ----------
u0 = np.array([0.0, 0.0, 1.0]); us = [np.array([0.0, 0.0, 1.0]), np.array([0.36, 0.0, 0.933]), np.array([0.48, 0.6, 0.64]), np.array([-0.6, 0.0, 0.8])]
print("\\nStep 3: one molecule's axis at t = 0, 1, 2, 3 ps:", [np.round(v, 3).tolist() for v in us])
for t_, v in enumerate(us):
    c = u0 @ v; print(f"   t = {t_}: cos(theta) = u(0).u(t) = {c:.3f}   P2 = (3 cos^2 - 1)/2 = {1.5 * c * c - 0.5:.3f}")
dt, every = 0.002, 10; fdt = dt * every
t0 = time.time()
for name, kw in (("small random rotations, D_r = 1", dict(J=0.0)), ("large jumps, rate 2 per time unit", dict(J=0.0, jump=2.0))):
    Ms, U = run(nsteps=50000, seed=1, **kw); c1, c2 = corr(U, 150)
    t1, t2 = fdt * (c1.sum() - 0.5), fdt * (c2.sum() - 0.5)
    print(f"   {name}: C1 at t = 0.1, 0.2, 0.5: {c1[5]:.3f}, {c1[10]:.3f}, {c1[25]:.3f};  C2: {c2[5]:.3f}, {c2[10]:.3f}, {c2[25]:.3f}")
    print(f"      tau1 = area of C1 = {t1:.3f},  tau2 = area of C2 = {t2:.3f},  tau1/tau2 = {t1 / t2:.2f}")
print(f"   theory, rotational diffusion: tau1 = 1/(2 D_r) = 0.500, tau2 = 1/(6 D_r) = 0.167, ratio 3;  jumps: tau1 = tau2 = 1/rate = 0.500, ratio 1  ({time.time() - t0:.0f} s)")

# ---------- Step 2 and 4: Kirkwood factor and collective relaxation ----------
print(f"\\nStep 2: Kirkwood factor G_K = <M^2>/(N mu^2) for the lattice model (N = {N}), 200 000 steps each")
res = {}
t0 = time.time()
for J in (0.0, 0.4, -0.4):
    Ms, U = run(J, 200000, seed=3); GK = (Ms**2).sum(1).mean() / N; m2 = (Ms.mean(0)**2).sum() / N
    blocks = [(b**2).sum(1).mean() / N for b in np.array_split(Ms, 10)]
    c1, _ = corr(U[:5000], 200); aM = acfM(Ms, 600); t1 = fdt * (c1.sum() - 0.5); tM = fdt * (aM / aM[0]).sum() - 0.5 * fdt
    res[J] = (GK, np.std(blocks, ddof=1) / np.sqrt(10), t1, tM, Ms)
    print(f"   J = {J:+.1f}: G_K = {GK:.3f} +- {res[J][1]:.3f} (<M>^2/N = {m2:.3f});  water-like eps = 1 + {y3:.2f} G_K = {1 + y3 * GK:.1f} +- {y3 * res[J][1]:.1f}")
print(f"   ({time.time() - t0:.0f} s)")
print("\\nStep 4: single-dipole vs box-dipole relaxation")
for J in (0.0, 0.4, -0.4):
    GK, eGK, t1, tM, Ms = res[J]
    print(f"   J = {J:+.1f}: single tau1 = {t1:.3f}   box dipole tau_M = {tM:.3f}   ratio {tM / t1:.2f}   (G_K = {GK:.2f})")

# ---------- Step 5: convergence ----------
print("\\nStep 5: how precisely do runs of different length give G_K? (J = +0.4, pieces of the 200 000-step run)")
Ms = res[0.4][4]
for L_ in (1000, 2000, 5000):
    pcs = [Ms[i:i + L_] for i in range(0, len(Ms) - L_ + 1, L_)]; g = [(p**2).sum(1).mean() / N for p in pcs]
    print(f"   runs of {L_ * every:7d} steps ({L_ * fdt:5.0f} time units = {L_ * fdt / res[0.4][3]:5.0f} tau_M): G_K spread {np.std(g, ddof=1):.3f} ({100 * np.std(g, ddof=1) / np.mean(g):.0f} %)  -> eps spread {y3 * np.std(g, ddof=1):.1f}   ({len(pcs)} pieces)")
`,
  },

  exercises: [
    { q: 'A model of methanol (μ = 1.7 D, n = 14.9 nm⁻³, 298 K) gives ⟨δM²⟩/(N μ²) = 2.6 in a tin-foil Ewald run. What ε does it predict?', a: 'μ = 1.7 × 3.336×10⁻³⁰ = 5.671×10⁻³⁰ C m; μ² = 3.216×10⁻⁵⁹. n μ²/(3ε₀k_BT) = 1.49×10²⁸ × 3.216×10⁻⁵⁹ / (3 × 8.854×10⁻¹² × 4.114×10⁻²¹) = 4.792×10⁻³¹ / 1.0928×10⁻³¹ = 4.39. ε = 1 + 4.39 × 2.6 = 12.4 (experiment: about 33, so this model would underestimate the orientational correlation).' },
    { q: 'NMR gives τ₂ = 1.0 ps and dielectric relaxation of the same liquid gives 8 ps. Can you conclude that τ₁/τ₂ = 8?', a: 'No. Dielectric relaxation measures the box dipole M, whose relaxation time can be much longer than the single-molecule τ₁ when molecules are correlated (step 4: 2.7 times longer for G_K ≈ 3). The single-molecule τ₁ has to come from a simulation (or from theory relating τ_M to τ₁) before taking the ratio.' },
    { q: 'With rotational diffusion D_r = 0.1 ps⁻¹, what are τ₁ and τ₂? If a measurement finds τ₂ = 4.5 ps for this liquid, what does it suggest?', a: 'τ₁ = 1/(2D_r) = 5 ps, τ₂ = 1/(6D_r) = 1.67 ps. A measured τ₂ of 4.5 ps (ratio τ₁/τ₂ ≈ 1.1) points to large-angle jumps rather than small-step rotational diffusion.' },
  ],
};
