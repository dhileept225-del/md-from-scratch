import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-free-energy',
  lab: true,
  summary: 'Most research questions are about free energies: binding, solvation, which conformation wins, how high a barrier is. A free energy is not an average of the energy, so it needs special methods. We compute a potential of mean force from g(r), do free-energy perturbation on six numbers, charge an ion with thermodynamic integration, and stitch umbrella windows together by hand before letting WHAM do it, always on systems whose exact answer we know.',
  minutes: 40,
  sim: 'lab9',
  learn: [
    'Turn g(r) into a potential of mean force, and see why it differs from the bare potential',
    'Do free-energy perturbation and thermodynamic integration by hand, and measure the overlap problem with numbers',
    'Unbias umbrella windows by hand, combine them with WHAM, and see how many windows a barrier needs',
  ],
  simHint: 'FEP: pick how much stiffer the spring becomes and the number of samples; compare forward and backward estimates with the exact line, and watch the backward one fail when the two states hardly overlap. Umbrella: choose the number of windows and the spring; the coloured histograms must overlap for WHAM (amber) to recover the exact barrier (white).',

  setup: `
<p>Why a special method? The free energy F = −k_BT ln Z counts how many states are accessible, not just how low they are. An ordinary simulation gives averages like ⟨U⟩; the free energy <i>difference</i> between two states can be written as an average too, but of an exponential, which brings its own problems.</p>
<div class="calc">step 1     PMF from g(r)  (the 9-atom LJ box, k_BT = 0.72)
steps 2–3  a spring that gets stiffer (k_BT = 1)
           exact:  ΔF = (k_BT/2) ln(k₁/k₀)
step 4     charging an ion, solvent modelled by a fluctuating
           potential  (k_BT = 2.494 kJ/mol)
steps 5–6  the double well of Lab 7.4, U = 5(s² − 1)²
           exact barrier 5 k_BT</div>`,

  steps: [
    {
      title: 'Potential of mean force from g(r)',
      why: 'g(r) says how much more often two atoms are found at distance r than in an ideal gas. Turning this probability into an energy, w(r) = −k_BT ln g(r), gives the free energy of the pair at that distance, averaged over everything else: the potential of mean force (PMF). For ion pairs it tells whether contact pairs or solvent-separated pairs are preferred.',
      work: `<div class="calc">dilute gas (no other atoms): g = e^(−U/k_BT)
r = 1.12:  g = e^(0.9998/0.72) = 4.009
           w = −0.72 ln 4.009 = <b>−0.9998</b> = U(r)

dense liquid (our 9-atom box, ρ = 0.69):
r       g        w = −0.72 ln g     bare U
1.025   1.453    −0.2690           −0.4750
1.125   2.687    <b>−0.7117</b>           −0.9998
1.475   0.929    +0.0531           −0.3507
1.675   0.612    <b>+0.3536</b>           −0.1729</div>`,
      see: 'In a gas the PMF is just the pair potential. In the liquid it is not: the well is shallower (−0.71 instead of −1.00), and at 1.675, where the bare potential is still attractive, the PMF is <i>repulsive</i> (+0.35), because the pair at that distance has to push away the neighbours that would like to sit between them. The PMF contains the solvent. Merit: it costs nothing beyond g(r). Limitation: it only works where g(r) is well sampled; at high PMF values (rare distances) g is tiny and noisy, which is why barriers need the biased methods of steps 5–6.',
    },
    {
      title: 'Free-energy perturbation by hand',
      why: 'Zwanzig\'s formula gives the free-energy change from state 0 to state 1 using samples of state 0 only: ΔF = −k_BT ln ⟨e^(−ΔU/k_BT)⟩₀, where ΔU = U₁ − U₀ is evaluated on each sample. Test it on a spring whose stiffness doubles, where the exact answer is known.',
      work: `<div class="calc">U₀ = ½ x²  ->  U₁ = ½ (2) x²,  k_BT = 1
exact ΔF = ½ ln 2 = <b>0.3466</b>

six samples of state 0:
x       ΔU = ½x²   e^(−ΔU)
 0.30   0.0450     0.9560
−1.20   0.7200     0.4868
 0.80   0.3200     0.7261
−0.10   0.0050     0.9950
 1.90   1.8050     0.1645
−0.60   0.1800     0.8353
average e^(−ΔU) = 0.6939
ΔF = −ln 0.6939 = <b>0.3654</b>      100 000 samples: 0.3488

energy change: ⟨U₁⟩₁ − ⟨U₀⟩₀ = ½ − ½ = 0
so T ΔS = ΔU − ΔF = <b>−0.3466</b></div>`,
      see: 'Six samples already give 0.365 against the exact 0.347. The free energy rises although the average energy does not change at all (½k_BT in both states, by equipartition): the stiffer spring confines the particle to a narrower range, and the whole ΔF is entropy. This is why free energies cannot be computed from energy averages alone. Note the sample at x = 1.90 contributes only 0.16 to the average: FEP is dominated by the samples where U₁ is low.',
    },
    {
      title: 'The overlap problem',
      why: 'FEP only works if the samples of state 0 include the configurations that matter for state 1. Make the spring 100 times stiffer and run FEP in both directions: forward (sample the soft spring, evaluate the stiff one) and backward (sample the stiff spring, evaluate the soft one).',
      work: `<div class="calc">exact ΔF = ½ ln 100 = <b>2.3026</b>     (20 repeats each)
samples    forward (sample k = 1)   backward (sample k = 100)
100        2.358 ± 0.287            0.862 ± 0.269
1 000      2.298 ± 0.094            1.105 ± 0.245
100 000    2.301 ± 0.008            <b>1.345 ± 0.127</b>

backward in 4 stages (k = 1, 3.16, 10, 31.6, 100),
1000 samples each:  <b>2.180</b>

k = 1 has 32 % of its weight at |x| > 1;
samples of k = 100 found there: 0.0000 %</div>`,
      see: 'Forward works: the soft spring visits the narrow region the stiff one needs. Backward fails badly, and more samples barely help (1.35 after 100 000 samples, still 40 % low), because the stiff spring never visits |x| > 1, where the soft spring has a third of its weight: the missing part is never seen, so the estimate is biased, not just noisy. The cures: go in the direction from the broad to the narrow distribution, split a large change into stages so that neighbours overlap (2.18 with 4 stages), and use both directions together (the Bennett acceptance ratio, BAR/MBAR, the standard in GROMACS and alchemical tools).',
    },
    {
      title: 'Charging an ion: FEP fails, TI works',
      why: 'Solvation free energies are computed by slowly switching on an ion\'s charge. A simple model of the solvent: the electrostatic potential φ at the ion site fluctuates like a Gaussian with spread s; the energy of charge q is qφ. Linear response gives the exact answer, so we can test one-step FEP and thermodynamic integration (TI), which integrates the average slope ⟨∂U/∂λ⟩ along a switching parameter λ.',
      work: `<div class="calc">s = 50 kJ/mol/e,  k_BT = 2.494 kJ/mol,  q: 0 -> +1
exact ΔF = −q² s² / (2 k_BT) = −2500 / 4.988 = <b>−501.2 kJ/mol</b>

one-step FEP from the uncharged state, 10⁶ samples:
−k_BT ln ⟨e^(−φ/k_BT)⟩ = <b>−199.4</b>
needed: φ near −s²/k_BT = −1002 = −20 standard deviations
most negative sample: −232 (−4.6 s)

TI: U = λ q φ,  ∂U/∂λ = φ,  ΔF = ∫₀¹ ⟨φ⟩_λ dλ
λ = 0.0   ⟨φ⟩ = −1.2       (1000 samples)
λ = 0.5   ⟨φ⟩ = −498.8
λ = 1.0   ⟨φ⟩ = −1002.1
trapezoid: 0.5 × (−1.2 − 498.8)/2 + 0.5 × (−498.8 − 1002.1)/2
         = <b>−500.2 kJ/mol</b>
FEP in 10 stages (1000 samples each): <b>−502.1</b></div>`,
      see: 'One-step FEP is off by 300 kJ/mol even with a million samples: the charged ion needs solvent configurations 20 standard deviations away from what the neutral solute ever sees. TI with just three λ points gets −500.2 against −501.2, because it only needs ordinary averages at each λ, where the samples are typical. Staged FEP works for the same reason. In this model ⟨φ⟩ is a straight line in λ, so three points suffice; for real ions (and especially for switching off LJ, where soft-core potentials are needed) the curve bends and needs 10–20 λ windows. A value of −500 kJ/mol is the right order of magnitude for the hydration of a small ion such as Na⁺ (about −400 kJ/mol), and a real calculation also needs the charged-box correction of Lab 5.5.',
    },
    {
      title: 'Umbrella sampling by hand',
      why: 'A barrier top is visited too rarely to measure its probability. Umbrella sampling adds a spring w(s) = ½k(s − s₀)² that holds the system near s₀, measures the histogram there, and removes the known bias afterwards: F(s) = −k_BT ln P(s) − w(s) + constant.',
      work: `<div class="calc">window at s₀ = 0 (the barrier top), k = 30, 10⁶ steps
s       count     P        −ln P    − w(s)   = F + C₀   exact U
−0.15   119134    0.1191   2.128    −0.338   1.790      4.778
−0.05   139700    0.1397   1.968    −0.038   1.931      4.975
 0.05   146083    0.1461   1.924    −0.038   1.886      4.975
 0.15   130550    0.1305   2.036    −0.338   1.698      4.778
 0.25    97314    0.0973   2.330    −0.938   1.392      4.395
 0.35    71569    0.0716   2.637    −1.837   <b>0.800</b>      3.850

window at s₀ = 0.5: at the shared bin s = 0.35 it gives 3.297
shift it by 0.800 − 3.297 = −2.498 to match window 0
then F(0.75) − F(0.05) = <b>−3.943</b>   (exact −4.018)</div>`,
      see: 'Without the spring, the barrier top would be visited a fraction e^(−5) of the time; with it, 10⁵ samples per bin. After removing the bias, the window reproduces the shape of the barrier (the F + C₀ column falls from 1.93 to 0.80 as the exact U falls from 4.98 to 3.85), but only up to an unknown constant C₀. Neighbouring windows are connected through bins they both sample. By hand we matched them at one bin and got −3.94 against −4.02. Limitation: the matching relies on the overlap; a bin seen rarely by one of the windows makes the shift noisy.',
    },
    {
      title: 'WHAM: all windows at once',
      why: 'Matching windows one bin at a time wastes data and accumulates errors. The weighted histogram analysis method (WHAM) finds the constants of all windows together, so that every bin uses all the windows that visited it, weighted by how much they saw. How many windows, and how long, does a 5 k_BT barrier need?',
      work: `<div class="tbl"><table>${head(['windows (s₀ from −1.2 to 1.2)', 'steps each', 'total steps', 'barrier (exact 5.0)'])}
${row(['5', '100 000', '500 000', '5.32 ± 0.25'])}
${row(['9', '20 000', '180 000', '5.18 ± 0.37'])}
${row(['9', '100 000', '900 000', '5.14 ± 0.29'])}
${row(['17', '50 000', '850 000', '<b>5.07 ± 0.14</b>'])}</table></div>
<div class="calc">mean ± spread over 4 independent repeats
one plain run of 900 000 steps: <b>3</b> barrier crossings</div>`,
      see: 'Umbrella sampling with WHAM recovers the 5 k_BT barrier, which a plain run of the same length crosses only three times (far too few to measure anything). With 5 windows spaced 0.6 apart the overlap is poor and the result is high by 0.3; with 17 windows the error drops to ±0.14 for the same total cost as 9 long windows: overlap matters more than length. Merits: every window is an ordinary, independent simulation (easy to run in parallel), and the method is the same as metadynamics\' result in Lab 7.4 (5.00), from a different route. Limitations: the collective variable must be chosen in advance, as in metadynamics, and slow motions orthogonal to it make windows converge slowly. In practice: PLUMED for the biases, <code>gmx wham</code> or WHAM by Grossfield for the analysis.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)'])}
${row(['PMF from g(r)', 'free: −0.71 well, +0.35 solvent barrier', 'only where g(r) is well sampled'])}
${row(['FEP (Zwanzig)', '6 samples: 0.365 vs 0.347', 'needs overlap: backward 1.35 vs 2.30'])}
${row(['staged FEP / BAR', '4 stages: 2.18; 10 stages: −502 vs −501', 'more windows to run'])}
${row(['TI', '3 λ points: −500.2 vs −501.2', 'curved ⟨∂U/∂λ⟩ needs 10–20 λ'])}
${row(['umbrella + WHAM', 'barrier 5.07 ± 0.14 (plain run: 3 crossings)', 'needs overlapping windows and a good CV'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the PMF of the dilute gas and of the 9-atom liquid, FEP by hand and the overlap test, the ion-charging FEP and TI, the umbrella window by hand and the WHAM runs. Run <code>pip install numpy</code>, then <code>python lab_free_energy.py</code> (about 20 s).',
    src: `# Lab 7.5: free energies by hand. PMF from g(r), free-energy perturbation, thermodynamic integration, umbrella sampling + WHAM.
import numpy as np, time
rng = np.random.default_rng(0)

# ---------- Step 1: PMF from g(r) ----------
T = 0.72
U = lambda r: 4 * (r**-12 - r**-6)
print("Step 1: potential of mean force w(r) = -kT ln g(r), kT = 0.72 (reduced LJ)")
for r in (1.0, 1.12, 1.5):
    g = np.exp(-U(r) / T); print(f"   dilute gas, r = {r}: g = exp(-U/kT) = exp(-({U(r):.4f})/0.72) = {g:.3f}  ->  w = -0.72 ln {g:.3f} = {-T * np.log(g):.4f} = U")
L, N, RC, dt = 3.6, 9, 1.8, 0.005
LAT = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(3) for i in range(3)])
def forces(x):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < RC * RC
    ir6 = np.where(m, r2**-3.0, 0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0)
    return -(f[:, :, None] * d).sum(1)
x = LAT.copy(); v = rng.normal(size=(N, 2)) * np.sqrt(T); f = forces(x); c = np.exp(-dt); sn = np.sqrt((1 - c * c) * T)
hist = np.zeros(36); nfr = 0
for s in range(202000):
    v += 0.5 * dt * f; x += dt * v; x %= L; f = forces(x); v += 0.5 * dt * f; v = c * v + sn * rng.normal(size=v.shape)
    if s > 2000 and s % 50 == 0:
        d = x[None] - x[:, None]; d -= L * np.round(d / L); rr = np.sqrt((d**2).sum(-1))[np.triu_indices(N, 1)]
        hist += np.histogram(rr, bins=36, range=(0, 1.8))[0]; nfr += 1
edges = np.linspace(0, 1.8, 37); mid = 0.5 * (edges[1:] + edges[:-1]); rho = N / L**2
g = hist / nfr / (rho * np.pi * (edges[1:]**2 - edges[:-1]**2) * N / 2)
print("   dense liquid (the 9-atom box, rho = 0.69):")
for rr in (1.025, 1.125, 1.475, 1.675):
    k = np.argmin(np.abs(mid - rr)); print(f"   r = {mid[k]:.3f}: g = {g[k]:.3f}  ->  w = -0.72 ln g = {-T * np.log(g[k]):+.4f}   (bare U = {U(mid[k]):+.4f})")

# ---------- Step 2: FEP by hand: stiffening a spring ----------
kT = 1.0
xs = np.array([0.3, -1.2, 0.8, -0.1, 1.9, -0.6])
print("\\nStep 2: FEP, U0 = 1/2 x^2  ->  U1 = 1/2 (2) x^2  (kT = 1); exact dF = (kT/2) ln(k1/k0) = 0.5 ln 2 = %.4f" % (0.5 * np.log(2)))
dU = 0.5 * xs**2; e = np.exp(-dU / kT)
for xi, du, ei in zip(xs, dU, e): print(f"   x = {xi:5.2f}: dU = 1/2 x^2 = {du:.4f}   exp(-dU) = {ei:.4f}")
print(f"   average exp(-dU) = {e.mean():.4f}  ->  dF = -ln({e.mean():.4f}) = {-np.log(e.mean()):.4f}  (6 samples)")
big = rng.normal(size=100000); print(f"   100 000 samples: dF = {-np.log(np.mean(np.exp(-0.5 * big**2))):.4f}")
print(f"   energy change <U1>_1 - <U0>_0 = kT/2 - kT/2 = 0  ->  T dS = dU - dF = {-0.5 * np.log(2):.4f}: the whole change is entropy")

# ---------- Step 3: overlap ----------
print("\\nStep 3: overlap. Make the spring 100 times stiffer: exact dF = 0.5 ln 100 = %.4f" % (0.5 * np.log(100)))
def fep(k0, k1, n, seed):
    r = np.random.default_rng(seed); xx = r.normal(size=n) / np.sqrt(k0); return -np.log(np.mean(np.exp(-0.5 * (k1 - k0) * xx**2)))
for n in (100, 1000, 100000):
    fw = [fep(1, 100, n, s) for s in range(20)]; bw = [-fep(100, 1, n, s) for s in range(20)]
    print(f"   {n:6d} samples: forward (sample k=1) {np.mean(fw):.3f} +- {np.std(fw):.3f}    backward (sample k=100) {np.mean(bw):.3f} +- {np.std(bw):.3f}")
ks = [1, 3.16, 10, 31.6, 100]; tot = 0
for a_, b_ in zip(ks[:-1], ks[1:]): tot += np.mean([-fep(b_, a_, 1000, s) for s in range(20)])
print(f"   backward in 4 stages k = 1, 3.16, 10, 31.6, 100 (1000 samples each): {tot:.3f}")
print(f"   fraction of k=100 samples with |x| > 1 (where k=1 has {np.mean(np.abs(big) > 1) * 100:.0f} % of its weight): {np.mean(np.abs(big / 10) > 1) * 100:.4f} %")

# ---------- Step 4: charging an ion: FEP vs TI ----------
kTk = 2.494; sphi = 50.0                    # kJ/mol, std of the solvent potential at the ion site in kJ/mol/e
exact = -sphi**2 / (2 * kTk)
print(f"\\nStep 4: charging a site 0 -> +1 e. Solvent potential phi at the site: Gaussian, mean 0, std {sphi} kJ/mol/e (kT = {kTk})")
print(f"   exact (linear response): dF = -q^2 s^2 / (2 kT) = -{sphi}^2 / (2 x {kTk}) = {exact:.1f} kJ/mol")
ph = rng.normal(size=1000000) * sphi
print(f"   one-step FEP from the uncharged state, 10^6 samples: -kT ln <exp(-phi/kT)> = {-kTk * np.log(np.mean(np.exp(-ph / kTk))):.1f} kJ/mol")
print(f"   the important samples have phi near -s^2/kT = {-sphi**2 / kTk:.0f} = {-sphi / kTk:.0f} standard deviations; most negative sampled: {ph.min():.0f} ({ph.min() / sphi:.1f} s)")
print("   TI: dF = integral_0^1 <dU/dlambda>_lambda dlambda,  U = lambda q phi  ->  <dU/dlambda> = <phi>_lambda")
lams = [0.0, 0.5, 1.0]; vals = []
for lam in lams:
    # sampling at coupling lambda: phi is shifted to mean -lambda s^2/kT (linear response); draw 1000 samples
    smp = rng.normal(size=1000) * sphi - lam * sphi**2 / kTk; vals.append(smp.mean())
    print(f"   lambda = {lam}: <phi> from 1000 samples = {smp.mean():8.1f} kJ/mol/e  (exact {-lam * sphi**2 / kTk:.1f})")
ti = 0.5 * (vals[0] + vals[1]) * 0.5 + 0.5 * (vals[1] + vals[2]) * 0.5
print(f"   trapezoid: 0.5 x ({vals[0]:.1f} + {vals[1]:.1f})/2 + 0.5 x ({vals[1]:.1f} + {vals[2]:.1f})/2 = {ti:.1f} kJ/mol  (exact {exact:.1f})")
print(f"   FEP in 10 stages (1000 samples each):", end=" ")
tot = 0
for lam in np.linspace(0, 0.9, 10):
    smp = rng.normal(size=1000) * sphi - lam * sphi**2 / kTk; tot += -kTk * np.log(np.mean(np.exp(-0.1 * smp / kTk)))
print(f"{tot:.1f} kJ/mol")

# ---------- Steps 5-6: umbrella sampling and WHAM on the double well ----------
Uw = lambda s: 5.0 * (s**2 - 1)**2; dUw = lambda s: 20.0 * s * (s**2 - 1)
def window(s0, kw, nsteps, seed, dt=0.001, gamma=5.0):
    r = np.random.default_rng(seed); s = s0; v = 0.0; c = np.exp(-gamma * dt); sn = np.sqrt(1 - c * c); out = np.empty(nsteps); noise = r.normal(size=nsteps)
    F = lambda s: -dUw(s) - kw * (s - s0)
    for k in range(nsteps):
        v += 0.5 * dt * F(s); s += dt * v; v += 0.5 * dt * F(s); v = c * v + sn * noise[k]; out[k] = s
    return out[2000:]
kw = 30.0; bins = np.linspace(-1.6, 1.6, 33); bm = 0.5 * (bins[1:] + bins[:-1])
w0 = window(0.0, kw, 1002000, 3); h0 = np.histogram(w0, bins)[0]
print(f"\\nStep 5: umbrella window (10^6 steps) at s0 = 0 with a spring w(s) = 1/2 x {kw} x (s - s0)^2 holds the system on the barrier top")
print("   bin centre   count    P = count/total   -ln P     - w(s)    = F + const")
sel = [i for i in range(len(bm)) if -0.2 < bm[i] < 0.4]
for i in sel:
    P = h0[i] / h0.sum(); print(f"   {bm[i]:6.2f}   {h0[i]:7d}    {P:.4f}           {-np.log(P):6.3f}    {-0.5 * kw * bm[i]**2:6.3f}    {-np.log(P) - 0.5 * kw * bm[i]**2:6.3f}   (exact U = {Uw(bm[i]):.3f})")
w1 = window(0.5, kw, 1002000, 3); h1 = np.histogram(w1, bins)[0]
Fa = lambda h, s0: -np.log(np.where(h > 0, h / h.sum(), np.nan)) - 0.5 * kw * (bm - s0)**2
F0, F1 = Fa(h0, 0.0), Fa(h1, 0.5); j = np.argmin(np.abs(bm - 0.35))
print(f"   second window at s0 = 0.5: at the shared bin s = {bm[j]:.2f} window 0 says {F0[j]:.3f}, window 0.5 says {F1[j]:.3f}  ->  shift window 0.5 by {F0[j] - F1[j]:+.3f}")
k2 = np.argmin(np.abs(bm - 0.75))
print(f"   then F(0.75) - F(0.05) = {F1[k2] + F0[j] - F1[j] - F0[np.argmin(np.abs(bm - 0.05))]:.3f}  (exact {Uw(bm[k2]) - Uw(bm[np.argmin(np.abs(bm - 0.05))]):.3f})")

def wham(centres, nsteps, seed0=10, tol=1e-7):
    H = []; Nw = []
    for i, s0 in enumerate(centres):
        h = np.histogram(window(s0, kw, nsteps + 2000, seed0 + i), bins)[0]; H.append(h); Nw.append(h.sum())
    H = np.array(H, float); Nw = np.array(Nw, float); bias = np.array([0.5 * kw * (bm - s0)**2 for s0 in centres]); fK = np.zeros(len(centres))
    for it in range(20000):
        den = (Nw[:, None] * np.exp(-(bias - fK[:, None]))).sum(0); P = H.sum(0) / den
        fn = -np.log((P[None] * np.exp(-bias)).sum(1)); fn -= fn[0]
        if np.abs(fn - fK).max() < tol: break
        fK = fn
    Fw = -np.log(np.where(P > 0, P, np.nan)); left = Fw[(bm > -1.3) & (bm < -0.7)]; top = Fw[np.argmin(np.abs(bm))]
    return top - np.nanmin(left), it, H
t0 = time.time()
print(f"\\nStep 6: WHAM over all windows (spring {kw}), barrier F(0) - F(left minimum), exact 5.0; mean +- spread over 4 repeats")
for name, cen, ns in (("5 windows, 100 000 steps each", np.linspace(-1.2, 1.2, 5), 100000), ("9 windows, 20 000 steps each", np.linspace(-1.2, 1.2, 9), 20000),
                      ("9 windows, 100 000 steps each", np.linspace(-1.2, 1.2, 9), 100000), ("17 windows, 50 000 steps each", np.linspace(-1.2, 1.2, 17), 50000)):
    res = [wham(cen, ns, seed0=100 * rep + 10)[0] for rep in range(4)]
    print(f"   {name:30s}: barrier {np.mean(res):.3f} +- {np.std(res, ddof=1):.3f}   (total {len(cen) * ns:,} steps)")
print(f"   ({time.time() - t0:.0f} s)")
r0 = np.random.default_rng(5); s = -1.0; v = 0.0; cross = 0; side = -1; c = np.exp(-0.005); sn = np.sqrt(1 - c * c)
for k in range(9 * 100000):
    v += -0.0005 * dUw(s); s += 0.001 * v; v += -0.0005 * dUw(s); v = c * v + sn * r0.normal()
    if s * side < 0 and abs(s) > 0.8: cross += 1; side = -side
print(f"   for comparison, one plain run of the same total length (900 000 steps): {cross} barrier crossings")
`,
  },

  exercises: [
    { q: 'An ion pair has g = 3.5 at contact (0.28 nm) and g = 0.6 at the barrier between contact and solvent-separated pairs (0.40 nm), at 300 K. What is the barrier for the contact pair to separate, in kJ/mol?', a: 'w = −k_BT ln g with k_BT = 2.494 kJ/mol. w(contact) = −2.494 ln 3.5 = −3.12; w(barrier) = −2.494 ln 0.6 = +1.27. Barrier = 1.27 − (−3.12) = 4.39 kJ/mol (1.8 k_BT).' },
    { q: 'Three FEP samples of ΔU (k_BT = 1): 0.5, 1.0, 3.0. What is ΔF? Which sample dominates?', a: 'e^(−0.5) = 0.6065, e^(−1) = 0.3679, e^(−3) = 0.0498; average 0.3414; ΔF = −ln 0.3414 = 1.075. The lowest ΔU (0.5) dominates: the exponential average weighs the samples where state 1 is happiest.' },
    { q: 'TI on a charging curve gives ⟨∂U/∂λ⟩ = −10, −180, −420 kJ/mol at λ = 0, 0.5, 1. Estimate ΔF with the trapezoid rule. Why might Simpson\'s rule differ?', a: 'Trapezoid: 0.5 × (−10 − 180)/2 + 0.5 × (−180 − 420)/2 = −47.5 − 150 = −197.5 kJ/mol. Simpson: (1/6)(−10 + 4 × (−180) − 420) = −191.7. They differ because the curve bends (−170 then −240 per half-step): the trapezoid assumes straight lines, so a curved integrand needs more λ points.' },
  ],
};
