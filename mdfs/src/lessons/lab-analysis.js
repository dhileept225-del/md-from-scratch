import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-analysis',
  lab: true,
  summary: 'A trajectory is only a long list of positions; every number we report comes from analysing it. We compute the radial distribution function, the mean squared displacement, the velocity autocorrelation, the van Hove function, the conductivity, the Onsager coefficients and the Haven ratio from trajectories small enough to do by hand. The numbers show what each recipe really averages over, why some quantities converge fast and others very slowly, and how to turn simulation units into S/m.',
  minutes: 40,
  sim: 'lab6',
  learn: [
    'Compute g(r), the MSD with all time origins and the VACF from a few frames, and check D two ways',
    'Read the van Hove function and the non-Gaussian parameter, and see what a few fast jumps do to them',
    'Get the conductivity, Onsager coefficients, transference number and Haven ratio (σ_NE/σ) for 2 + 2 ions, and see why collective quantities need much longer runs',
  ],
  simHint: 'MSD mode: a random walker; compare the MSD from one time origin with the MSD averaged over all origins, and watch the long-lag end get noisy. Haven mode: 50 cations and 50 anions in 1D; set the fraction bound in neutral pairs, run, and compare the self MSD with the collective (charge) MSD. Press New run several times: D_self barely changes, the Haven ratio jumps around.',

  setup: `
<p>Two systems. <b>Part A</b> (steps 1–4) uses the 9-atom 2D Lennard-Jones box of Labs 2 and 3 (reduced units, T = 0.72, a gentle Bussi thermostat so the dynamics are realistic), plus a few made-up numbers small enough to do by hand. <b>Part B</b> (steps 5–7) uses 2 cations and 2 anions moving in 1D (nm, ps), and a Brownian model of 100 ions.</p>
<div class="calc">quantity           what it averages
g(r)               pairs at distance r, over frames
MSD(t)             displacement², over atoms and origins
VACF(t)            v(0)·v(t), over atoms and origins
G_s(r, t)          histogram of displacements
σ (conductivity)   (Σ q Δr)², the whole system as one
Haven ratio        σ_NE / σ  (self vs collective)</div>
<p>All positions used for displacements must be <b>unwrapped</b>: the atom is followed across the periodic walls instead of being put back in the box (step 2 shows why).</p>`,

  steps: [
    {
      title: 'Radial distribution function: counting pairs',
      why: 'g(r) tells how the density around an atom differs from the average: g = 1 means "as in an ideal gas", a peak means a shell of neighbours. It is the bridge to X-ray and neutron scattering, and its first peak gives the coordination number. The recipe is counting, plus one division.',
      work: `<div class="calc">one frame: the lattice of Lab 2, ρ = 9/3.6² = 0.6944
36 pairs (minimum image): 18 at r = 1.200, 18 at 1.697

bin [1.10, 1.30]:
neighbours per atom = 2 × 18 / 9 = 4.000
ideal gas expects   ρ × 2πr × dr = 0.6944 × 2π × 1.2 × 0.2
                                 = 1.0472
g = 4.000 / 1.0472 = <b>3.82</b>

bin [1.60, 1.80]:  4.000 / 1.4809 = <b>2.70</b>
everything else:   0</div>
<div class="calc">moving liquid, 4000 frames, bins of 0.05:
first peak g = <b>2.71</b> at r = 1.125;  g = 0 below 0.85
one single frame, bins in the peak region:
  0.0  1.0  1.0  3.8  2.7  0.9  1.7  1.6
coordination number (neighbours closer than 1.45)
  = <b>4.36</b></div>`,
      see: 'The factor 2 counts each pair for both atoms; the ideal-gas count 2πr·dr·ρ grows with r, which is why g → 1 at large distance even though there are more and more neighbours. A single frame is useless (36 pairs over 36 bins: mostly 0 or 3), so g(r) is always averaged over many frames, and it converges quickly because every frame contributes N(N−1)/2 pairs. Limitations: only r < L/2 is meaningful (here 1.8), and the bin width trades resolution against noise. In 3D the shell volume is 4πr²dr.',
    },
    {
      title: 'Mean squared displacement: all time origins',
      why: 'The diffusion coefficient D follows from how far atoms wander: in d dimensions MSD = 2dDt at long times (Einstein). A trajectory of n frames gives many displacements for each time lag, if we use every frame as a possible starting point (time origin).',
      work: `<div class="calc">x (nm) at t = 0, 1, 2, 3, 4, 5 ps:
0.00  0.30  0.10  0.60  0.80  1.20

lag 1:  0.3, −0.2, 0.5, 0.2, 0.4
        squares 0.09, 0.04, 0.25, 0.04, 0.16
        MSD = <b>0.1160</b>  (5 samples; origin t=0 only: 0.0900)
lag 2:  0.1, 0.3, 0.7, 0.6
        MSD = <b>0.2375</b>  (4 samples; t=0 only: 0.0100)
lag 3:  MSD = <b>0.6067</b>  (3 samples; t=0 only: 0.3600)
lag 4:  MSD = 0.7250  (2 samples)
lag 5:  MSD = 1.4400  (1 sample)

slope over lags 1–3 = 0.2453 nm²/ps
1D: D = slope / 2 = 0.1227 nm²/ps = 1.23×10⁻⁷ m²/s

wrapped in a box L = 1 nm: 0.00 0.30 0.10 0.60 0.80 <b>0.20</b>
lag-5 MSD would be 0.04 instead of 1.44</div>
<div class="calc">9-atom LJ (2D):  MSD at t = 1, 2, 5, 10:
0.227, 0.462, 1.183, 2.370
D = slope(t = 2–10) / 4 = <b>0.0598</b></div>`,
      see: 'All origins turn one trajectory into many samples, so lags 1–3 already give a steady slope, while the single-origin MSD jumps around (0.09, 0.01, 0.36). The price is the long-lag end: lag 5 has a single sample, so the last part of any MSD curve is noise and must not be fitted. Also fit only the straight (diffusive) part, after the short-time ballistic region. With wrapped coordinates the MSD falls from 1.44 to 0.04 the moment the atom crosses a wall: always unwrap (LAMMPS dumps <code>xu yu zu</code>, GROMACS <code>-pbc nojump</code>).',
    },
    {
      title: 'Velocity autocorrelation: the same D a second way',
      why: 'Displacement is the time integral of velocity, so D can also come from how long an atom remembers its velocity: D = (1/d) ∫ ⟨v(0)·v(t)⟩ dt (Green–Kubo). Getting the same D from MSD and VACF is a strong check of the analysis.',
      work: `<div class="calc">v (nm/ps) at t = 0..5 ps:  1.0 0.7 0.5 0.4 0.1 0.2
lag 0:  (1 + 0.49 + 0.25 + 0.16 + 0.01 + 0.04)/6 = 0.3250
lag 1:  (0.70 + 0.35 + 0.20 + 0.04 + 0.02)/5    = 0.2620
lag 2:  (0.50 + 0.28 + 0.05 + 0.08)/4           = 0.2275
lag 3:  (0.40 + 0.07 + 0.10)/3                  = 0.1900
lag 4:  (0.10 + 0.14)/2                         = 0.1200
1D trapezoid, lags 0–4, dt = 1 ps:
D = ½(0.3250) + 0.2620 + 0.2275 + 0.1900 + ½(0.1200)
  = <b>0.902 nm²/ps</b></div>
<div class="calc">9-atom LJ:  C(0) = ⟨v·v⟩ = N_f T / N = 1.278
C(t) turns negative at t = 0.25, minimum −0.11 C(0) at 0.35
D = ½ ∫ C dt  up to t = 1: 0.0591   2: 0.0587   5: 0.0605
D from the MSD (step 2):           0.0598</div>`,
      see: 'For the liquid, the two routes agree within 1–2 % (0.059–0.061 vs 0.0598): the analysis is consistent. The VACF also shows physics the MSD hides: it goes negative after 0.25 time units, the atom bouncing back from its cage of neighbours. Limitations: the integral must be stopped somewhere, and the noisy tail makes the answer drift (0.0587 → 0.0605); and velocities must be saved very often (here every 10 steps). The toy series shows the other danger: its VACF has not decayed by lag 4, so its "D" is not yet meaningful.',
    },
    {
      title: 'Van Hove function: how displacements are distributed',
      why: 'The MSD is just one number per time lag, the average of r². The self part of the van Hove function, G_s(r, t), is the whole histogram of displacements. If every atom diffuses the same way it is a Gaussian; if some atoms jump while others are stuck, it grows tails. The non-Gaussian parameter measures that with two averages.',
      work: `<div class="calc">10 displacements after time t (1D, nm):
−0.5 −0.2 −0.1 0.0 0.1 0.1 0.2 0.3 0.9 −1.1

⟨x²⟩ = 0.2470       ⟨x⁴⟩ = 0.21943
Gaussian would have ⟨x⁴⟩ = 3⟨x²⟩² = 0.18303
α₂ = ⟨x⁴⟩ / (3⟨x²⟩²) − 1 = <b>0.20</b>
(in 3D: α₂ = 3⟨r⁴⟩ / (5⟨r²⟩²) − 1;  in 2D: ⟨r⁴⟩/(2⟨r²⟩²) − 1)

without the two jumps (0.9, −1.1):
⟨x²⟩ = 0.0563,   α₂ = <b>−0.02</b></div>
<div class="tbl"><table>${head(['9-atom LJ, t', '0.25', '1', '5', '20'])}
${row(['α₂', '0.038', '0.033', '0.042', '0.030'])}</table></div>`,
      see: 'Two "hoppers" out of ten carry 80 % of ⟨x²⟩ and push α₂ from about 0 to 0.20. The MSD alone cannot tell this situation (a few fast ions, many stuck) from uniform diffusion; G_s and α₂ can. Our simple LJ fluid stays nearly Gaussian (α₂ ≈ 0.03–0.04). In ionic liquids and glassy electrolytes, α₂ peaks at intermediate times, which is the signature of cages and hopping. The distinct part G_d(r, t), where a <i>different</i> atom is found at time t, is the time-dependent g(r) and shows how neighbours replace each other. Limitation: the tails are rare events, so G_s and α₂ need many more samples than the MSD.',
    },
    {
      title: 'Conductivity, Onsager coefficients and the Haven ratio: 2 + 2 ions',
      why: 'Conductivity is a collective property: it depends on how the total charge moves, ΣqΔx, not on how far each ion goes. If a cation and an anion travel together, they carry no current, even though both diffuse. The Nernst–Einstein estimate σ_NE uses only the self terms; the Haven ratio σ_NE/σ measures how much correlations change the result.',
      work: `<div class="calc">x (nm) at t = 0..4 ps; c1 and a1 move together
c1 (+1)  0.0  0.2  0.5  0.6  0.9
a1 (−1)  0.4  0.6  0.9  1.0  1.3
c2 (+1)  2.0  1.8  2.1  2.3  2.2
a2 (−1)  3.0  3.3  3.1  3.0  3.4

lag-1 displacements (4 origins), mean square:
c1  0.2  0.3  0.1  0.3    0.0575
a1  0.2  0.3  0.1  0.3    0.0575
c2 −0.2  0.3  0.2 −0.1    0.0450
a2  0.3 −0.2 −0.1  0.4    0.0750
self sum (-> σ_NE)                    <b>0.2350</b>

charge displacement ΔM = Σ q Δx:
−0.5  0.5  0.3  −0.5      ⟨ΔM²⟩ = <b>0.2100</b>  (-> σ)

Haven ratio = σ_NE / σ = 0.2350 / 0.2100 = <b>1.12</b>

Onsager-type sums (cations P, anions M):
L₊₊ ~ ⟨P²⟩ = 0.1225   L₋₋ ~ ⟨M²⟩ = 0.1875
L₊₋ ~ ⟨P M⟩ = 0.0500
L₊₊ + L₋₋ − 2L₊₋ = 0.2100  (= ⟨ΔM²⟩)
t₊ = (L₊₊ − L₊₋)/(L₊₊ + L₋₋ − 2L₊₋) = <b>0.345</b>
t₊ from self diffusion only: 0.1025/0.2350 = 0.436</div>`,
      see: 'The pair c1–a1 adds 0.115 to the self sum and exactly 0 to ΔM: a neutral pair diffuses but carries no current, so σ is lower than σ_NE and the Haven ratio is above 1. The Onsager decomposition shows where it comes from: L₊₋ > 0 means cations and anions tend to move the same way. It also changes the transference number, from 0.436 (self diffusion only) to 0.345 (with correlations), which matters for battery electrolytes. These are lag-1 numbers for illustration; in practice σ and the L\'s come from the long-time <i>slopes</i>, exactly as D comes from the MSD slope. The identity L₊₊ + L₋₋ − 2L₊₋ = ⟨ΔM²⟩ is a good check of any implementation.',
    },
    {
      title: 'Units: from nm²/ps to S/m',
      why: 'Simulation output comes in nm, ps and e; experiments report S/m or mS/cm. A unit slip here gives errors of factors of 10⁶, so it is worth one careful calculation.',
      work: `<div class="calc">Nernst–Einstein (3D):  σ_NE = e²/(k_B T) × n × (D₊ + D₋)
n = 2.5 ion pairs/nm³ = 2.5×10²⁷ /m³
D₊ = 2.0×10⁻¹¹, D₋ = 1.5×10⁻¹¹ m²/s, T = 300 K
e²/(k_B T) = 2.5670×10⁻³⁸ / 4.1419×10⁻²¹ = 6.1975×10⁻¹⁸
σ_NE = 6.1975×10⁻¹⁸ × 2.5×10²⁷ × 3.5×10⁻¹¹
     = <b>0.542 S/m = 5.42 mS/cm</b>
Haven ratio 1.5:  σ = 0.542/1.5 = <b>0.362 S/m = 3.62 mS/cm</b>

from a simulation:  σ = e² / (6 V k_B T) × slope of ⟨|Σ q Δr|²⟩
V = 125 nm³: this σ needs a slope of 4.375×10⁻⁸ m²/s
           = <b>0.0437 e² nm²/ps</b>
conversion: 1 m²/s = 10¹⁸ nm² / 10¹² ps = 10⁶ nm²/ps
            D₊ = 2.0×10⁻¹¹ m²/s = 2.0×10⁻⁵ nm²/ps</div>`,
      see: 'Every factor is visible: the charge squared, the thermal energy, the ion density and the sum of diffusion coefficients. Typical ionic-liquid numbers (D ~ 10⁻¹¹ m²/s, a few ion pairs per nm³) give a few mS/cm, as measured. Remember the conversion 1 m²/s = 10⁶ nm²/ps, and that σ (from ⟨|ΣqΔr|²⟩) has a 6 (3D) in the denominator just like D from the MSD. Here n counts ion pairs, so each species has density n; for a salt like Li₂X the charges and densities must be counted species by species.',
    },
    {
      title: 'Statistics: why collective quantities converge so slowly',
      why: 'D is an average over every ion; σ comes from one quantity per time origin, the total charge displacement of the whole box. With N ions, the self average has about N times more independent samples. How much noisier is σ in practice?',
      work: `<div class="calc">Brownian model: 50 cations + 50 anions in 1D, D = 1 each
34 of the 100 ions are in neutral pairs (move together)
expected:  D_self = 1
           D_σ (per ion) = 1 − 0.34 = 0.66
           Haven = 1 / 0.66 = <b>1.515</b>

40 independent runs of 2000 steps:
D_self = 1.003 ± 0.024   (2.4 %)
D_σ    = 0.665 ± 0.171   (25.8 %)
Haven  = 1.61 ± 0.43     (range 0.91 to 3.05)
noise ratio collective / self = <b>10.8</b>   (√100 = 10)</div>`,
      see: 'The same runs give D to 2.4 % but σ only to 26 %, about √N = 10 times worse, and a single run can report a Haven ratio anywhere from 0.91 to 3.05 when the true value is 1.52. To make σ as precise as D you need about N times more sampling (longer runs or several independent runs; a bigger box does not help σ, because ΔM is still one number per time origin). That is why conductivities and Onsager coefficients from short simulations often scatter, and why papers should report error bars from independent runs or block averaging. The self quantities (D, t₊ from self diffusion) are the cheap ones; the collective ones carry the correlations but cost far more.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['quantity', 'merit (number)', 'limitation (number)'])}
${row(['g(r)', 'converges fast: N(N−1)/2 pairs per frame; CN = 4.36', 'only r < L/2; one frame useless'])}
${row(['MSD, all origins', 'D = 0.0598 from a short run', 'long lags have few samples; unwrap (1.44 vs 0.04)'])}
${row(['VACF / Green–Kubo', 'D = 0.059–0.061, shows the cage', 'where to stop the integral; frequent output'])}
${row(['G_s, α₂', 'sees hoppers MSD cannot (α₂ 0.20 vs −0.02)', 'tails need many samples'])}
${row(['σ, Onsager L\'s', 'include correlations (t₊ 0.345 vs 0.436)', '~√N noisier: 26 % vs 2.4 %'])}
${row(['Haven ratio σ_NE/σ', 'measures ion pairing (1.12, 1.52)', 'single runs range 0.91–3.05'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the RDF by hand and from the run, the MSD and VACF by hand and for the LJ liquid, the van Hove moments, the 2 + 2 ion conductivity and Onsager sums, the unit conversion and the Brownian statistics test. Run <code>pip install numpy</code>, then <code>python lab_analysis.py</code> (about 15 s).',
    src: `# Lab 6.7: analysis by hand. Part A uses the 9-atom 2D LJ box of Labs 2-3 (reduced units).
# Part B uses a made-up trajectory of 2 cations + 2 anions (nm, ps) and a simple Brownian model.
import numpy as np, time
rng = np.random.default_rng(0)
L, N, RC, dt = 3.6, 9, 1.8, 0.005
LAT = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(3) for i in range(3)])

# ---------- Step 1: RDF by hand (lattice frame) ----------
rho = N / L**2; d = LAT[None] - LAT[:, None]; d -= L * np.round(d / L); r = np.sqrt((d**2).sum(-1))[np.triu_indices(N, 1)]
print(f"Step 1: RDF of one frame (the lattice). rho = 9 / 3.6^2 = {rho:.4f}. All 36 pair distances (minimum image):")
vals, cnt = np.unique(np.round(r, 3), return_counts=True)
for v_, c_ in zip(vals, cnt): print(f"   r = {v_:.3f}: {c_} pairs")
dr = 0.2
for rc_, npair in ((1.2, 18), (1.697, 18)):
    per = 2 * npair / N; ideal = rho * 2 * np.pi * rc_ * dr
    print(f"   bin [{rc_ - 0.1:.2f}, {rc_ + 0.1:.2f}]: neighbours per atom = 2 x {npair} / 9 = {per:.3f};  ideal gas expects rho 2 pi r dr = {rho:.4f} x 2 pi x {rc_} x {dr} = {ideal:.4f};  g = {per / ideal:.2f}")

def forces(x):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < RC * RC
    ir6 = np.where(m, r2**-3.0, 0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0)
    return -(f[:, :, None] * d).sum(1)
# Bussi thermostat at T = 0.72 (tau = 0.1 to equilibrate, then a gentle tau = 1 that barely disturbs the dynamics)
x = LAT.copy(); v = rng.normal(size=(N, 2)); v -= v.mean(0); NF = 16; T0 = 0.72; f = forces(x)
def bussi(v, tau):
    K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; c = np.exp(-dt / tau); R1 = rng.normal(); S = (rng.normal(size=NF - 1)**2).sum()
    a2 = c + (1 - c) * Kt / (NF * K) * (R1**2 + S) + 2 * R1 * np.sqrt(c * (1 - c) * Kt / (NF * K))
    return v * np.sqrt(a2) * np.sign(R1 + np.sqrt(c * NF * K / ((1 - c) * Kt)))
for s in range(20000):
    v += 0.5 * dt * f; x += dt * v; x %= L; f = forces(x); v += 0.5 * dt * f; v = bussi(v, 0.1)
t0 = time.time(); nsteps = 400000; every = 10
xs = np.zeros((nsteps // every, N, 2)); vs = np.zeros((nsteps // every, N, 2)); u = x.copy(); hist = np.zeros(36); nfr = 0
for s in range(nsteps):
    v += 0.5 * dt * f; x += dt * v; u += dt * v; x %= L; f = forces(x); v += 0.5 * dt * f; v = bussi(v, 1.0)
    if s % every == 0:
        xs[s // every] = u; vs[s // every] = v
        if s % 100 == 0:
            d = x[None] - x[:, None]; d -= L * np.round(d / L); rr = np.sqrt((d**2).sum(-1))[np.triu_indices(N, 1)]
            hist += np.histogram(rr, bins=36, range=(0, 1.8))[0]; nfr += 1
T = (vs**2).sum((1, 2)).mean() / NF
edges = np.linspace(0, 1.8, 37); mid = 0.5 * (edges[1:] + edges[:-1]); ideal = rho * np.pi * (edges[1:]**2 - edges[:-1]**2) * N / 2
g = hist / nfr / ideal
print(f"   run of 2000 time units ({time.time() - t0:.0f} s, T = {T:.3f}), RDF from {nfr} frames, bins of 0.05:")
k = np.argmax(g); print(f"   first peak g = {g[k]:.2f} at r = {mid[k]:.3f};  g at r < 0.85: {g[mid < 0.85].max():.2f}")
one = np.histogram(np.sqrt(((lambda d: (d - L * np.round(d / L)))(xs[-1][None] % L - xs[-1][:, None] % L)**2).sum(-1))[np.triu_indices(N, 1)], bins=36, range=(0, 1.8))[0] / ideal
print(f"   a single frame: g values in the first-peak region: {' '.join(f'{q:.1f}' for q in one[18:26])}  (mostly 0 or large: 36 pairs in 36 bins)")
cn = (hist / nfr * 2 / N)[mid < 1.45].sum()
print(f"   coordination number (neighbours closer than 1.45, the first minimum) = {cn:.2f}")

# ---------- Step 2: MSD by hand ----------
xt = np.array([0.00, 0.30, 0.10, 0.60, 0.80, 1.20])
print("\\nStep 2: MSD of one particle, x (nm) at t = 0..5 ps:", xt)
for lag in range(1, 6):
    dd = xt[lag:] - xt[:-lag]; print(f"   lag {lag} ps: displacements {np.round(dd, 2)}  squares {np.round(dd**2, 4)}  -> MSD = {np.mean(dd**2):.4f}  ({len(dd)} samples; single origin t=0: {dd[0]**2:.4f})")
lags = np.arange(1, 4); msd = np.array([np.mean((xt[l:] - xt[:-l])**2) for l in lags]); sl = np.polyfit(lags, msd, 1)[0]
print(f"   slope of MSD for lags 1-3 = {sl:.4f} nm^2/ps;  1D: MSD = 2 D t  ->  D = {sl / 2:.4f} nm^2/ps = {sl / 2 * 1e-6:.2e} m^2/s")
xw = xt % 1.0
print(f"   wrapped into a box L = 1 nm: {np.round(xw, 2)};  lag-5 MSD would be {(xw[5] - xw[0])**2:.4f} instead of {(xt[5] - xt[0])**2:.4f}")
def msd_all(X, maxlag):
    return np.array([((X[l:] - X[:-l])**2).sum(-1).mean() for l in range(1, maxlag + 1)])
tl = np.arange(1, 1001) * every * dt; M = msd_all(xs, 1000)
sel = (tl > 2) & (tl < 10); Dm = np.polyfit(tl[sel], M[sel], 1)[0] / 4
print(f"   9-atom LJ: MSD at t = 1, 2, 5, 10: {M[19]:.3f}, {M[39]:.3f}, {M[99]:.3f}, {M[199]:.3f};  D from slope (t = 2-10), 2D: MSD = 4 D t -> D = {Dm:.4f}")

# ---------- Step 3: VACF ----------
vt = np.array([1.0, 0.7, 0.5, 0.4, 0.1, 0.2])
print("\\nStep 3: VACF of one velocity component, v (nm/ps) at t = 0..5 ps:", vt)
C = [np.mean(vt[l:] * vt[:len(vt) - l]) for l in range(6)]
for l in range(6): print(f"   lag {l}: products {np.round(vt[l:] * vt[:len(vt) - l], 3)}  -> C = {C[l]:.4f}")
I = 1.0 * (0.5 * C[0] + sum(C[1:4]) + 0.5 * C[4])
print(f"   Green-Kubo, 1D: D = integral of C dt (trapezoid, lags 0-4, dt = 1 ps) = {I:.4f} nm^2/ps")
nv = len(vs); Vf = np.fft.rfft(vs, n=2 * nv, axis=0); ac = np.fft.irfft(Vf * np.conj(Vf), axis=0)[:nv].sum((1, 2)) / (np.arange(nv, 0, -1) * N)
tv = np.arange(nv) * every * dt; dtv = every * dt
for tm in (1.0, 2.0, 5.0):
    n_ = int(tm / dtv); print(f"   9-atom LJ: D = (1/2) integral of <v(0).v(t)> up to t = {tm}: {0.5 * (ac[:n_ + 1].sum() - 0.5 * (ac[0] + ac[n_])) * dtv:.4f}")
zc = np.argmax(ac < 0); print(f"   VACF: C(0) = <v.v> = N_f T / N = {ac[0]:.3f};  first becomes negative at t = {tv[zc]:.3f};  minimum {ac[:200].min() / ac[0]:.3f} x C(0) at t = {tv[np.argmin(ac[:200])]:.3f}")

# ---------- Step 4: van Hove self part and the non-Gaussian parameter ----------
dx = np.array([-0.5, -0.2, -0.1, 0.0, 0.1, 0.1, 0.2, 0.3, 0.9, -1.1])
m2, m4 = np.mean(dx**2), np.mean(dx**4); a2 = m4 / (3 * m2**2) - 1
print("\\nStep 4: 10 displacements after time t (1D, nm):", dx)
print(f"   <x^2> = {m2:.4f}, <x^4> = {m4:.5f};  Gaussian would have <x^4> = 3<x^2>^2 = {3 * m2**2:.5f};  alpha2 = <x^4>/(3<x^2>^2) - 1 = {a2:.2f}")
core = dx[np.abs(dx) < 0.6]; print(f"   without the two jumps (0.9, -1.1): <x^2> = {np.mean(core**2):.4f}, alpha2 = {np.mean(core**4) / (3 * np.mean(core**2)**2) - 1:.2f}")
for tm in (0.25, 1.0, 5.0, 20.0):
    l = int(tm / dtv); D2 = (xs[l:] - xs[:-l]).reshape(-1, 2); r2 = (D2**2).sum(-1)
    a2s = np.mean(r2**2) / (2 * np.mean(r2)**2) - 1                      # 2D: alpha2 = <r^4>/(2<r^2>^2) - 1
    print(f"   9-atom LJ, t = {tm:5.2f}: alpha2 = {a2s:.3f}")

# ---------- Step 5: conductivity and Haven ratio, 2 cations + 2 anions ----------
traj = {'c1': [0.0, 0.2, 0.5, 0.6, 0.9], 'a1': [0.4, 0.6, 0.9, 1.0, 1.3], 'c2': [2.0, 1.8, 2.1, 2.3, 2.2], 'a2': [3.0, 3.3, 3.1, 3.0, 3.4]}
q = {'c1': 1, 'a1': -1, 'c2': 1, 'a2': -1}; X = {k: np.array(vv) for k, vv in traj.items()}
print("\\nStep 5: x (nm) of 4 ions at t = 0..4 ps; c1 and a1 travel together as a pair")
for k_ in X: print(f"   {k_} (q = {q[k_]:+d}): {X[k_]}")
dX = {k_: X[k_][1:] - X[k_][:-1] for k_ in X}
for k_ in X: print(f"   lag-1 displacements {k_}: {np.round(dX[k_], 2)}  -> mean square {np.mean(dX[k_]**2):.4f}")
selfsum = sum(np.mean(dX[k_]**2) for k_ in X)
dM = sum(q[k_] * dX[k_] for k_ in X)
print(f"   charge displacement dM = sum q dx per origin: {np.round(dM, 2)}  -> <dM^2> = {np.mean(dM**2):.4f}")
print(f"   sum of self terms = {selfsum:.4f};  Haven ratio = sigma_NE / sigma = {selfsum:.4f} / {np.mean(dM**2):.4f} = {selfsum / np.mean(dM**2):.2f}")
P = dX['c1'] + dX['c2']; Mn = dX['a1'] + dX['a2']
Lpp, Lmm, Lpm = np.mean(P**2), np.mean(Mn**2), np.mean(P * Mn)
print(f"   Onsager-type sums (lag 1): L++ ~ <(dx_c1 + dx_c2)^2> = {Lpp:.4f}, L-- ~ {Lmm:.4f}, L+- ~ <(sum cations)(sum anions)> = {Lpm:.4f}")
print(f"   check: L++ + L-- - 2 L+- = {Lpp + Lmm - 2 * Lpm:.4f} (= <dM^2>);  transference t+ = (L++ - L+-)/(L++ + L-- - 2L+-) = {(Lpp - Lpm) / (Lpp + Lmm - 2 * Lpm):.3f}")
print(f"   self part of L++: {np.mean(dX['c1']**2) + np.mean(dX['c2']**2):.4f}, distinct (c1-c2 cross) part: {2 * np.mean(dX['c1'] * dX['c2']):.4f}")
tNE = (np.mean(dX['c1']**2) + np.mean(dX['c2']**2)) / selfsum; print(f"   transference from self diffusion only: t+ = {tNE:.3f}")

# ---------- Step 6: units ----------
e, kB = 1.602177e-19, 1.380649e-23; Tk = 300.0; n = 2.5e27; Dp, Dm_ = 2.0e-11, 1.5e-11
sNE = e**2 / (kB * Tk) * n * (Dp + Dm_)
print(f"\\nStep 6: Nernst-Einstein, n = 2.5 ion pairs/nm^3 = {n:.1e} /m^3, D+ = {Dp:.1e}, D- = {Dm_:.1e} m^2/s, T = {Tk} K")
print(f"   e^2/(kB T) = {e**2:.4e} / {kB * Tk:.4e} = {e**2 / (kB * Tk):.4e};  sigma_NE = {e**2 / (kB * Tk):.4e} x {n:.1e} x {Dp + Dm_:.1e} = {sNE:.3f} S/m = {10 * sNE:.2f} mS/cm")
H = 1.5; print(f"   with Haven ratio {H}: sigma = sigma_NE / H = {sNE / H:.3f} S/m = {10 * sNE / H:.2f} mS/cm")
V = 125e-27; slope = sNE / H * 6 * V * kB * Tk / e**2
print(f"   from a simulation: V = 125 nm^3, sigma = e^2 / (6 V kB T) x slope of <|sum q dr|^2>;  the slope that gives this sigma is {slope:.3e} m^2/s = {slope * 1e6:.4f} e^2 nm^2/ps")
print(f"   1 m^2/s = 10^18 nm^2 / 10^12 ps = 10^6 nm^2/ps;  e.g. D+ = {Dp:.1e} m^2/s = {Dp * 1e6:.1e} nm^2/ps")

# ---------- Step 7: statistics, self vs collective (Brownian model) ----------
print("\\nStep 7: Brownian model, 50 cations + 50 anions in 1D, D = 1 for every ion; one third of the ions are in bound neutral pairs (move together)")
def run(seed, npair=50, nb=17, steps=2000, h=0.01, lag=100):
    r = np.random.default_rng(seed); free = r.normal(size=(steps, 2 * npair)) * np.sqrt(2 * h)
    free[:, npair:npair + nb] = free[:, :nb]                     # anions 0..nb-1 move with cations 0..nb-1
    X_ = np.cumsum(free, 0); qv = np.r_[np.ones(npair), -np.ones(npair)]
    d_ = X_[lag:] - X_[:-lag]; selfm = (d_**2).mean(); Mm = ((d_ * qv).sum(1)**2).mean()
    return selfm / (2 * lag * h), Mm / (2 * lag * h) / (2 * npair)
res = np.array([run(s) for s in range(40)])
Dself, Dsig = res[:, 0], res[:, 1]
print(f"   expected: D_self = 1, D_sigma (per ion) = 1 - 2 x 17/100 = {1 - 34 / 100:.2f}, Haven = {1 / (1 - 0.34):.3f}")
print(f"   40 runs of 2000 steps: D_self = {Dself.mean():.3f} +- {Dself.std():.3f} ({100 * Dself.std() / Dself.mean():.1f} %);  D_sigma = {Dsig.mean():.3f} +- {Dsig.std():.3f} ({100 * Dsig.std() / Dsig.mean():.1f} %)")
Hs = Dself / Dsig; print(f"   Haven ratio from single runs: {Hs.mean():.2f} +- {Hs.std():.2f} (range {Hs.min():.2f} to {Hs.max():.2f})")
print(f"   noise ratio collective / self = {(Dsig.std() / Dsig.mean()) / (Dself.std() / Dself.mean()):.1f}  (sqrt of 100 ions = 10)")
`,
  },

  exercises: [
    { q: 'In a 3D box of 500 atoms with L = 4 nm, a bin from 0.30 to 0.32 nm contains, on average per frame, 1200 pairs. What is g(r) there?', a: 'ρ = 500/64 = 7.8125 nm⁻³. Neighbours per atom = 2 × 1200 / 500 = 4.8. Ideal shell = ρ × 4πr² × dr = 7.8125 × 4π × 0.31² × 0.02 = 0.1887. g = 4.8 / 0.1887 = 25.4. That is suspiciously high for a liquid; check whether the count was per frame or summed over frames.' },
    { q: 'A 3D MSD rises from 0.30 nm² at 100 ps to 1.50 nm² at 500 ps. What is D in m²/s?', a: 'Slope = 1.20 nm² / 400 ps = 0.003 nm²/ps. D = slope/6 = 5.0×10⁻⁴ nm²/ps = 5.0×10⁻¹⁰ m²/s (divide by 10⁶).' },
    { q: 'In the 2 + 2 example, suppose c2 and a2 also became a pair (moving together). What would ⟨ΔM²⟩, σ and the Haven ratio be?', a: 'Every ΔM would be 0, so σ = 0 and the Haven ratio would be infinite: a liquid made only of neutral pairs does not conduct, however fast its ions diffuse. Real electrolytes are in between: pairs and clusters form and break, which is why the Haven ratio is usually above 1, often 1–3.' },
  ],
};
