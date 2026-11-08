import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-finite-size',
  lab: true,
  summary: 'A simulation box holds a few hundred to a few hundred thousand molecules, not 10²³, and it repeats itself periodically. Some results do not care; others depend on the box size in a predictable way and must be corrected. We compute the 1/N offset in g(r) by hand, what it does to Kirkwood–Buff integrals, the Yeh–Hummer correction for diffusion, and then measure D in four box sizes to check the 1/L law.',
  minutes: 35,
  sim: 'lab10',
  learn: [
    'See why g(r) tends to 1 − 1/N instead of 1, and correct Kirkwood–Buff integrals for it',
    'Apply the Yeh–Hummer correction by hand, and measure the 1/L dependence of D in four box sizes',
    'Know which quantities need a size correction (D, σ_NE, charged-box energies, long integrals) and which do not',
  ],
  simHint: 'g(r) tail: choose N and watch random (ideal-gas) frames build a g(r) that settles at 1 − 1/N, not 1; switch on the Ganguly–van der Vegt correction and the running Kirkwood–Buff integral goes back to 0. Diffusion: the dots are the four box sizes of step 4, the line is Yeh–Hummer; move η and L to see how large the correction is for water or an ionic liquid.',

  setup: `
<p>Two kinds of finite-size effects appear:</p>
<div class="calc">counting effects   a fixed number N of particles:
                   averages differ from the infinite system
                   by terms of order 1/N  (g(r) tail, integrals)
periodic effects   every particle interacts with its images:
                   long-range flows and fields feel the box
                   size L  (diffusion ~1/L, charged boxes ~1/L)</div>
<p>Systems: an ideal gas (exact answers), the 2D LJ liquid of earlier labs in boxes of 9, 25 and 100 atoms, and a 3D LJ liquid (ρ = 0.8, T = 1.0) with 32 to 500 atoms. For real liquids we use water (η = 0.85 mPa s) and an ionic liquid (η = 50 mPa s).</p>`,

  steps: [
    {
      title: 'The g(r) tail is 1 − 1/N, not 1',
      why: 'g(r) compares the number of pairs at distance r with an ideal gas of the same density. With a fixed number of particles, each one sees only N − 1 others, while the usual normalization uses the density N/V. That leaves a small offset at large r, even for an ideal gas, where the exact g is 1.',
      work: `<div class="calc">pairs counted:       N(N − 1)/2
normalization uses:  N²/2
g(tail) = (N − 1)/N = 1 − 1/N
N = 10:  <b>0.900</b>        N = 100:  <b>0.990</b>

random positions (ideal gas), g in 5 shells:
N = 10:   0.900 0.899 0.900 0.900 0.900
N = 100:  0.992 0.990 0.990 0.990 0.989</div>
<div class="tbl"><table>${head(['2D LJ liquid, ρ = 0.69', 'L', 'first peak g', 'g at first minimum'])}
${row(['N = 9', '3.6', '2.62', '0.985'])}
${row(['N = 25', '6.0', '2.92', '0.788'])}
${row(['N = 100', '12.0', '2.95', '0.784'])}</table></div>`,
      see: 'The ideal gas confirms the formula exactly: the tail sits at 1 − 1/N everywhere. In a dense liquid the offset is smaller, about S(0)/N, where S(0) = ρk_BTκ_T is small for a stiff liquid. The first peak and minimum are what usually matter, and they hardly change from 25 to 100 atoms (2.92 vs 2.95). The 9-atom box is too small: its half-width 1.8 cuts right through the first minimum, so its structure is distorted (peak 2.62, minimum 0.985). Rule of thumb: the box must be several times larger than the structure you want to see.',
    },
    {
      title: 'Integrals of g(r): Kirkwood–Buff and the Ganguly–van der Vegt correction',
      why: 'Kirkwood–Buff integrals G = ∫(g − 1) dV over all space connect g(r) to compressibilities and partial molar volumes, and running coordination numbers are integrals too. A tiny offset in g(r), multiplied by a large volume, becomes a large error.',
      work: `<div class="calc">ideal gas, N = 100, L = 10 (2D);  true G = 0
g − 1 = −1/N = −0.01 everywhere
G(R) = −0.01 × πR²:
R = 1:  −0.0314     R = 2:  −0.1257     R = 4:  <b>−0.5027</b>

Ganguly–van der Vegt correction at R = 4:
V_R/V = πR²/L² = 0.5027
excess neighbours  ΔN = ρ G = −0.5027
factor = N(1 − V_R/V) / (N(1 − V_R/V) − ΔN − 1)
       = 49.735 / (49.735 + 0.503 − 1) = 1.01010
corrected g = 0.99 × 1.01010 = <b>1.0000</b>

3D: 500 molecules, L = 2.5 nm, integrate to R = 1.2 nm:
sphere 7.238 nm³ × 1/500 = 0.0145 nm³ = <b>8.7 cm³/mol</b></div>`,
      see: 'The running integral drifts away linearly with the area (or volume) instead of levelling off, and the drift is pure artefact. The Ganguly–van der Vegt factor accounts for the particles that are "missing" because they are already counted within R; for the ideal gas it restores g = 1 exactly. In 3D with 500 molecules, the error in a Kirkwood–Buff integral is about 9 cm³/mol, comparable to the integrals themselves (often tens of cm³/mol). Limitations: the correction assumes the box is large compared with the range of correlations, and integrals still need R well beyond the last peak. Coordination numbers of the first shell are not affected noticeably.',
    },
    {
      title: 'Diffusion: the Yeh–Hummer correction by hand',
      why: 'A moving particle drags the surrounding liquid along. In a periodic box, the flow it creates also pushes on its own images, and the periodic flow field slows it down. Hydrodynamics gives the size of the effect: D in the box is lower than the true D by ξk_BT/(6πηL), where η is the viscosity and ξ = 2.837297 for a cubic box.',
      work: `<div class="calc">D_∞ = D_PBC + ξ k_B T / (6π η L)

water, 300 K, η = 0.85 mPa s, L = 3 nm:
2.837297 × 4.1419×10⁻²¹ / (6π × 0.00085 × 3×10⁻⁹)
= <b>2.445×10⁻¹⁰ m²/s</b>   vs D = 2.3×10⁻⁹:  +11 %

ionic liquid, 300 K, η = 50 mPa s, L = 5 nm:
2.837297 × 4.1419×10⁻²¹ / (6π × 0.05 × 5×10⁻⁹)
= <b>2.494×10⁻¹² m²/s</b>   vs D = 1.0×10⁻¹¹: <b>+25 %</b>
same liquid, L = 10 nm: +12 %,  L = 20 nm: +6 %</div>`,
      see: 'The correction does not depend on the particle: every species gets the same absolute shift. It falls only as 1/L, so doubling the box (8 times the atoms) only halves it: for the ionic liquid a 20 nm box, about a million atoms, still underestimates D by 6 %. It matters most where it is least expected: the absolute correction is smaller in a viscous liquid, but D is smaller too, so the relative error is larger (25 % vs 11 %). Limitation: it needs the viscosity, ideally computed with the same force field (next lab), and it assumes a cubic box and a simple liquid.',
    },
    {
      title: 'Measuring the 1/L law: four box sizes',
      why: 'Is the correction real, and is it the right size? Run the same 3D LJ liquid in boxes of 32, 108, 256 and 500 atoms, measure D in each with the MSD (several independent runs for the noisy small boxes), and plot D against 1/L.',
      work: `<div class="tbl"><table>${head(['N', 'L', '1/L', 'D (MSD)', 'runs'])}
${row(['32', '3.420', '0.2924', '0.0440 ± 0.0014', '12'])}
${row(['108', '5.130', '0.1949', '0.0565 ± 0.0013', '4'])}
${row(['256', '6.840', '0.1462', '0.0597 ± 0.0021', '2'])}
${row(['500', '8.550', '0.1170', '0.0626', '1'])}</table></div>
<div class="calc">straight line:  D = 0.0755 − 0.1054 / L
D_∞ = <b>0.0755</b>
slope = −ξT/(6πη)  ->  η = 2.837 × 1.0 / (6π × 0.1054) = <b>1.43</b>

D_PBC below D_∞:  32 atoms −42 %, 108: −25 %, 256: −21 %, 500: −17 %
Yeh–Hummer applied to each box:
0.0749, 0.0771, 0.0751, 0.0750   (all ≈ 0.0755)</div>`,
      see: 'D grows steadily with box size and falls on a straight line in 1/L, as hydrodynamics predicts. Even the 500-atom box underestimates D by 17 %. Two ways to fix it: run several box sizes and extrapolate to 1/L = 0 (as here), or run one box and add the Yeh–Hummer term with a known viscosity: corrected, all four boxes agree on 0.075–0.077. The slope itself gives the viscosity (1.43 in LJ units for this truncated potential), a useful cross-check of a Green–Kubo viscosity. Limitation: the small-box points need many runs, because D from 32 atoms is noisy.',
    },
    {
      title: 'Consequences for conductivity and the Haven ratio',
      why: 'Transport numbers built from self-diffusion inherit the Yeh–Hummer error. The Nernst–Einstein conductivity σ_NE uses D₊ and D₋ directly, and the Haven ratio σ_NE/σ compares it with the collective conductivity.',
      work: `<div class="calc">ionic liquid of Lab 6.7: n = 2.5 /nm³,
D₊ = 2.0×10⁻¹¹, D₋ = 1.5×10⁻¹¹ m²/s, 5 nm box, η = 50 mPa s
correction per ion 2.49×10⁻¹² m²/s
D₊ -> 2.25×10⁻¹¹,  D₋ -> 1.75×10⁻¹¹

σ_NE:  0.542 -> 0.620 S/m   (<b>+14 %</b>)
if σ (collective) = 0.362 S/m needs no correction:
Haven ratio  1.50 (box)  ->  <b>1.71</b> (corrected)</div>`,
      see: 'A 14 % change in σ_NE moves the Haven ratio from 1.50 to 1.71, larger than many differences reported between force fields. The collective conductivity is expected to need a much smaller correction, because the hydrodynamic flow carries cations and anions together and largely cancels in the charge current; but this depends on the system (masses, sizes, how the box momentum is handled), so check it by running two box sizes. Either way: state the box size with any D, σ_NE or Haven ratio, and say whether D was corrected.',
    },
    {
      title: 'Choosing the box: a checklist with numbers',
      why: 'Box size is a trade-off between cost (N grows as L³) and the size-dependent errors above. Putting the main ones side by side helps choose L before a project starts.',
      work: `<div class="tbl"><table>${head(['L (nm)', 'largest cutoff L/2', 'smallest k = 2π/L', 'IL: correction to D_PBC', 'one-ion charged-box energy'])}
${row(['3', '1.5 nm', '2.09 /nm', '+42 %', '−66 kJ/mol'])}
${row(['5', '2.5 nm', '1.26 /nm', '+25 %', '−39 kJ/mol'])}
${row(['10', '5.0 nm', '0.63 /nm', '+12 %', '−20 kJ/mol'])}</table></div>`,
      see: 'Some limits are hard: the cutoff (and any structure you want in g(r)) must stay below L/2, and no wave longer than L exists, so S(k) below 2π/L and slow collective modes are missing. Others are corrections you can calculate: D (Yeh–Hummer), charged boxes (Lab 5.5), long integrals of g(r) (Ganguly–van der Vegt). And some need checking: conductivity, viscosity and anything involving clusters or interfaces. A practical habit: run the key observable in two box sizes; if it changes by more than its error bar, either correct it or go bigger.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['quantity', 'size dependence (number)', 'what to do'])}
${row(['g(r) first peak', '2.92 (25 atoms) vs 2.95 (100)', 'nothing, if L/2 > the structure'])}
${row(['g(r) tail', '1 − 1/N (0.990 for 100)', 'Ganguly–van der Vegt for integrals'])}
${row(['Kirkwood–Buff integrals', '~9 cm³/mol error for 500 molecules', 'correct, and integrate far enough'])}
${row(['self-diffusion D', '−17 % even with 500 atoms; ∝ 1/L', 'Yeh–Hummer or extrapolate in 1/L'])}
${row(['σ_NE, Haven ratio', '+14 % σ_NE, Haven 1.50 -> 1.71', 'correct D; check σ with two boxes'])}
${row(['charged box', '−39 kJ/mol per ion at 5 nm', 'analytic correction (Lab 5.5)'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the ideal-gas g(r) tail, the 2D LJ boxes, the Kirkwood–Buff example and its correction, the Yeh–Hummer numbers, the four 3D box sizes (about 40 s) and the transport and checklist tables. Run <code>pip install numpy</code>, then <code>python lab_finite_size.py</code> (about 75 s).',
    src: `# Lab: finite-size effects by hand. How the box size changes g(r), its integrals, diffusion and transport numbers.
import numpy as np, time
rng = np.random.default_rng(0)

# ---------- Step 1: the g(r) tail is 1 - 1/N ----------
print("Step 1: ideal gas (no forces) in a 2D periodic box, density 1. Pairs per area: N(N-1)/2 / A, ideal-gas normalisation N^2/2 / A")
for N in (10, 100):
    print(f"   N = {N:3d}: g(tail) = N(N-1)/N^2 = 1 - 1/N = {1 - 1 / N:.3f}")
def g_ideal(N, frames, bins):
    L = np.sqrt(N); h = np.zeros(len(bins) - 1)
    for _ in range(frames):
        x = rng.random((N, 2)) * L; d = x[None] - x[:, None]; d -= L * np.round(d / L); r = np.sqrt((d**2).sum(-1))[np.triu_indices(N, 1)]
        h += np.histogram(r, bins)[0]
    ideal = N / L**2 * np.pi * (bins[1:]**2 - bins[:-1]**2) * N / 2
    return h / frames / ideal
for N, fr in ((10, 200000), (100, 2000)):
    b = np.linspace(0.5, np.sqrt(N) / 2, 6); g = g_ideal(N, fr, b); print(f"   N = {N:3d}, random positions, {fr} frames: g in 5 shells = {' '.join(f'{q:.3f}' for q in g)}")

# LJ liquid in 2D, same density, three box sizes: does the first peak move?
def lj_box(n, nsteps, seed, T=0.72, rc=1.8, dt=0.005):
    N = n * n; L = 1.2 * n; x = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(n) for i in range(n)]); r = np.random.default_rng(seed)
    v = r.normal(size=(N, 2)) * np.sqrt(T); c = np.exp(-dt); sn = np.sqrt((1 - c * c) * T)
    def F(x):
        d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
        ir6 = np.where(m, r2**-3.0, 0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0); return -(f[:, :, None] * d).sum(1), r2
    f, _ = F(x); bins = np.linspace(0, 1.8, 37); h = np.zeros(36); nf = 0
    for s in range(nsteps):
        v += 0.5 * dt * f; x += dt * v; x %= L; f, r2 = F(x); v += 0.5 * dt * f; v = c * v + sn * r.normal(size=v.shape)
        if s > 2000 and s % 50 == 0: h += np.histogram(np.sqrt(r2[np.triu_indices(N, 1)]), bins)[0]; nf += 1
    ideal = N / L**2 * np.pi * (bins[1:]**2 - bins[:-1]**2) * N / 2; g = h / nf / ideal; mid = 0.5 * (bins[1:] + bins[:-1])
    k = np.argmax(g); return g[k], mid[k], g[(mid > 1.4) & (mid < 1.5)].mean()
t0 = time.time()
print("   2D LJ liquid, rho = 0.69, T = 0.72:")
for n, ns in ((3, 202000), (5, 102000), (10, 42000)):
    pk, rp, mn = lj_box(n, ns, 1); print(f"   N = {n * n:3d} (L = {1.2 * n:4.1f}): first peak g = {pk:.2f} at r = {rp:.3f};  g at the first minimum (1.4-1.5) = {mn:.3f}")
print(f"   ({time.time() - t0:.0f} s)")

# ---------- Step 2: integrals of g(r): Kirkwood-Buff and the Ganguly-van der Vegt correction ----------
N, L = 100, 10.0; A = L * L
print(f"\\nStep 2: running Kirkwood-Buff integral G(R) = integral_0^R (g - 1) 2 pi r dr (2D), ideal gas, N = {N}, L = {L} (true G = 0)")
for R in (1.0, 2.0, 4.0):
    G = -(1 / N) * np.pi * R**2; print(f"   R = {R}: (g - 1) = -1/N = -0.01 everywhere -> G(R) = -0.01 x pi x {R}^2 = {G:.4f}")
R = 4.0; x = np.pi * R**2 / A; dN = (N / A) * (-1 / N) * np.pi * R**2; num = N * (1 - x); den = N * (1 - x) - dN - 1
print(f"   Ganguly-van der Vegt at R = {R}: V_R/V = pi R^2 / L^2 = {x:.4f};  excess neighbours dN = rho x G = {dN:.4f}")
print(f"   factor N(1 - V_R/V) / (N(1 - V_R/V) - dN - 1) = {num:.3f} / ({num:.3f} + {-dN:.3f} - 1) = {num / den:.5f};  corrected g = 0.99 x {num / den:.5f} = {0.99 * num / den:.4f}")
Ng = 500; Lg = 2.5; VR = 4 / 3 * np.pi * 1.2**3
print(f"   3D example: N = {Ng} molecules, L = {Lg} nm, integrate to R = 1.2 nm: sphere {VR:.3f} nm^3;  offset 1/N = {1 / Ng:.3f}  ->  error in G = {VR / Ng:.4f} nm^3 = {VR / Ng * 602.2:.1f} cm^3/mol")

# ---------- Step 3: Yeh-Hummer by hand ----------
kB = 1.380649e-23; xi = 2.837297
def yh(T, eta, L): return xi * kB * T / (6 * np.pi * eta * L)
print(f"\\nStep 3: Yeh-Hummer: D_inf = D_PBC + xi kB T / (6 pi eta L),  xi = {xi} (cubic box)")
for name, T, eta, L, D in (("water, 300 K, eta = 0.85 mPa s, L = 3 nm", 300, 0.85e-3, 3e-9, 2.3e-9), ("ionic liquid, 300 K, eta = 50 mPa s, L = 5 nm", 300, 0.05, 5e-9, 1.0e-11)):
    c = yh(T, eta, L); print(f"   {name}: correction = {xi} x {kB * T:.4e} / (6 pi x {eta} x {L:.0e}) = {c:.3e} m^2/s;  vs D_PBC = {D:.1e}: +{100 * c / D:.0f} %")
print(f"   ionic liquid in a 10 nm box: +{100 * yh(300, 0.05, 1e-8) / 1e-11:.0f} %;  a 20 nm box: +{100 * yh(300, 0.05, 2e-8) / 1e-11:.0f} %")

# ---------- Step 4: Yeh-Hummer measured: 3D LJ liquid in four box sizes ----------
RC, SK, rho, T0, dt = 1.6, 0.3, 0.8, 1.0, 0.004
def fcc(nc, a):
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]); g = np.array([[i, j, k] for i in range(nc) for j in range(nc) for k in range(nc)])
    return ((g[:, None, :] + b[None]) * a).reshape(-1, 3)
def pairs(x, L):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); return np.nonzero(np.triu(r2 < (RC + SK)**2, 1))
def forces(x, L, i, j, N):
    d = x[j] - x[i]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); m = r2 < RC * RC
    ir6 = np.where(m, 1 / r2**3, 0); ff = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0)[:, None] * d; f = np.zeros((N, 3))
    for k in range(3): f[:, k] = np.bincount(j, ff[:, k], N) - np.bincount(i, ff[:, k], N)
    return f
def diffusion(nc, nsteps, seed, tau=1.0):
    N = 4 * nc**3; L = (N / rho)**(1 / 3); x = fcc(nc, L / nc); r = np.random.default_rng(seed); NF = 3 * N - 3
    v = r.normal(size=(N, 3)); v -= v.mean(0); v *= np.sqrt(NF * T0 / (v**2).sum()); i, j = pairs(x, L); xl = x.copy(); f = forces(x, L, i, j, N)
    u = x.copy(); traj = []; Ts = []; c = np.exp(-dt / tau)
    for s in range(-3000, nsteps):
        v += 0.5 * dt * f; x += dt * v; u += dt * v
        if np.abs((x - xl + L / 2) % L - L / 2).max() > SK / 2: x %= L; i, j = pairs(x, L); xl = x.copy()
        f = forces(x, L, i, j, N); v += 0.5 * dt * f
        K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = r.normal(); S2 = r.chisquare(NF - 1)          # Bussi thermostat, gentle (tau = 1)
        v *= np.sqrt(max(c + (1 - c) * Kt / (NF * K) * (R1**2 + S2) + 2 * R1 * np.sqrt(c * (1 - c) * Kt / (NF * K)), 0))
        if s >= 0 and s % 10 == 0: traj.append(u - u.mean(0)); Ts.append((v**2).sum() / NF)
    traj = np.array(traj); lags = np.arange(50, 151, 25); msd = [((traj[l:] - traj[:-l])**2).sum(-1).mean() for l in lags]
    return L, np.polyfit(lags * 10 * dt, msd, 1)[0] / 6, np.mean(Ts)
t0 = time.time()
print("\\nStep 4: 3D LJ liquid, rho = 0.8, T = 1.0 (cutoff 1.6), four box sizes; D from the MSD, several independent runs")
res = []
for nc, nrun, ns in ((2, 12, 10000), (3, 4, 10000), (4, 2, 10000), (5, 1, 12000)):
    out = [diffusion(nc, ns, 10 * nc + k) for k in range(nrun)]; L = out[0][0]; Ds = np.array([o[1] for o in out]); Tm = np.mean([o[2] for o in out])
    err = Ds.std(ddof=1) / np.sqrt(nrun) if nrun > 1 else np.nan; res.append((L, Ds.mean(), err))
    print(f"   N = {4 * nc**3:3d}, L = {L:.3f}, 1/L = {1 / L:.4f}: D = {Ds.mean():.4f}" + (f" +- {err:.4f}" if nrun > 1 else "") + f"   ({nrun} runs, T = {Tm:.3f})")
Ls = np.array([q[0] for q in res]); Dm = np.array([q[1] for q in res]); p = np.polyfit(1 / Ls, Dm, 1)
eta = xi * T0 / (6 * np.pi * -p[0])
print(f"   straight line D = {p[1]:.4f} + ({p[0]:.4f}) / L  ->  D_inf = {p[1]:.4f};  slope = -xi T/(6 pi eta)  ->  eta = {eta:.2f} (reduced)")
for L, D, _ in res: print(f"   N = {round(rho * L**3):3d}: D_PBC is {100 * (1 - D / p[1]):.0f} % below D_inf;  Yeh-Hummer with this eta adds {xi * T0 / (6 * np.pi * eta * L):.4f} -> {D + xi * T0 / (6 * np.pi * eta * L):.4f}")
print(f"   ({time.time() - t0:.0f} s)")

# ---------- Step 5: consequences for conductivity and the Haven ratio ----------
e = 1.602177e-19; T = 300; n = 2.5e27; Dp, Dm_ = 2.0e-11, 1.5e-11; eta_il = 0.05; L = 5e-9
c = yh(T, eta_il, L); s_pbc = e**2 / (kB * T) * n * (Dp + Dm_); s_inf = e**2 / (kB * T) * n * (Dp + Dm_ + 2 * c)
print(f"\\nStep 5: ionic liquid from Lab 6.7 (n = 2.5 /nm^3, D+ = 2.0e-11, D- = 1.5e-11 m^2/s) in a 5 nm box, eta = 50 mPa s")
print(f"   correction per ion {c:.2e} m^2/s:  D+ -> {Dp + c:.2e}, D- -> {Dm_ + c:.2e}")
print(f"   sigma_NE: {s_pbc:.3f} -> {s_inf:.3f} S/m (+{100 * (s_inf / s_pbc - 1):.0f} %)")
sig = s_pbc / 1.5
print(f"   if sigma (collective) = {sig:.3f} S/m needs no correction: Haven ratio {s_pbc / sig:.2f} (box) -> {s_inf / sig:.2f} (corrected)")

# ---------- Step 6: choosing the box ----------
print("\\nStep 6: checklist for a box of side L")
for Lb in (3.0, 5.0, 10.0):
    print(f"   L = {Lb:4.1f} nm: largest cutoff L/2 = {Lb / 2:.1f} nm;  smallest wavevector 2 pi/L = {2 * np.pi / Lb:.2f} /nm;  "
          f"IL diffusion correction {100 * yh(300, 0.05, Lb * 1e-9) / 1e-11:.0f} %;  charged-box energy of one ion {-xi * 138.935 / (2 * Lb):.0f} kJ/mol")
`,
  },

  exercises: [
    { q: 'A water box of 4 nm gives D = 2.10×10⁻⁹ m²/s at 298 K. With η = 0.89 mPa s, what is the Yeh–Hummer corrected value?', a: 'Correction = 2.837297 × 1.380649×10⁻²³ × 298 / (6π × 0.00089 × 4×10⁻⁹) = 1.1673×10⁻²⁰ / 6.710×10⁻¹¹ = 1.74×10⁻¹⁰ m²/s. D_∞ = 2.27×10⁻⁹ m²/s (+8 %).' },
    { q: 'You run 2000 molecules and integrate g(r) to R = 2 nm in a box of 5 nm. Roughly how large is the 1/N error in the Kirkwood–Buff integral?', a: 'Sphere volume 4/3 π 2³ = 33.5 nm³; offset 1/2000 = 0.0005; error ≈ 0.0168 nm³ ≈ 10 cm³/mol. Bigger boxes reduce the offset, but integrating further makes the volume grow, so the correction is still needed.' },
    { q: 'Two runs of the same ionic liquid give D₊ = 1.62×10⁻¹¹ (L = 4 nm) and 1.86×10⁻¹¹ m²/s (L = 8 nm). Extrapolate to an infinite box and estimate the viscosity at 300 K.', a: 'D = D_∞ − a/L: 1.62 = D_∞ − a/4, 1.86 = D_∞ − a/8 (×10⁻¹¹, L in nm). Subtracting: a/8 = 0.24, a = 1.92, D_∞ = 2.10×10⁻¹¹ m²/s. a = ξk_BT/(6πη) = 1.92×10⁻¹¹ m²/s × 10⁻⁹ m = 1.92×10⁻²⁰ m³/s, so η = 2.837297 × 4.1419×10⁻²¹ / (6π × 1.92×10⁻²⁰) = 0.0324 Pa s = 32 mPa s.' },
  ],
};
