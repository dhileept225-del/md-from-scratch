import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-viscosity',
  lab: true,
  summary: 'Viscosity says how hard it is to shear a liquid, and it controls diffusion, conductivity and flow. In equilibrium MD there is no flow, yet the viscosity is hidden in the small, random shear stresses of the box. We compute the shear stress by hand, do the Green–Kubo and Einstein routes on ten numbers, then measure the viscosity of the LJ liquid of Lab 6.9, check it against the Yeh–Hummer slope and Stokes–Einstein, and see why it converges so slowly.',
  minutes: 35,
  sim: 'lab11',
  learn: [
    'Compute an off-diagonal pressure (shear stress) from velocities and forces by hand',
    'Get the viscosity from the stress autocorrelation (Green–Kubo) and from the Einstein form, and see where to stop the integral',
    'Check the result against the Yeh–Hummer slope and Stokes–Einstein, and estimate how long a viscosity run must be',
  ],
  simHint: 'The cyan curve is the running Green–Kubo integral of the LJ liquid; it rises to a plateau and then wanders off as noise accumulates. Pick a run length and the number of stress components: the grey curves are independent pieces, the spread of their plateaus is the error bar. Averaging the 5 shear components narrows the spread about √5 times.',

  setup: `
<p>A shear flow (layers sliding past each other) is resisted by a shear stress, the off-diagonal pressure P_xy. At rest its average is zero, but it fluctuates. The time these fluctuations take to die away is the viscosity:</p>
<div class="calc">P_xy = (1/V) [ Σ m vₓ v_y  +  Σ_pairs xᵢⱼ F_y,ᵢⱼ ]

Green–Kubo:  η = V/(k_BT) ∫₀^∞ ⟨P_xy(0) P_xy(t)⟩ dt
Einstein:    η = V/(2k_BT) × slope of ⟨[∫P_xy dt]²⟩ in time</div>
<p>System: the 3D LJ liquid of Lab 6.9 (256 atoms, ρ = 0.8, T = 1.0, cutoff 1.6, V = 320), 100 000 steps with a gentle Bussi thermostat. In an isotropic liquid five combinations of the pressure tensor carry the same information: P_xy, P_xz, P_yz, (P_xx − P_yy)/2 and (P_yy − P_zz)/2.</p>`,

  steps: [
    {
      title: 'The shear stress by hand',
      why: 'P_xy measures how much x-momentum crosses a plane perpendicular to y, per unit area and time. Atoms carry momentum across by moving (kinetic part) and push across by forces (virial part). Both are needed, and both come from quantities the simulation already has.',
      work: `<div class="calc">two atoms, V = 27, mass 1
r₁ = (1, 1, 1),   v₁ = (0.5, −0.3, 0)
r₂ = (2, 1.5, 1), v₂ = (−0.2, 0.4, 0.1)

kinetic:  Σ m vₓ v_y = 0.5 × (−0.3) + (−0.2) × 0.4 = <b>−0.2300</b>

d = r₂ − r₁ = (1, 0.5, 0),  r = 1.1180
F(r) = 24(2/r¹³ − 1/r⁷) = 0.2638 (repulsive)
force on atom 2 = 0.2638 × d/r = (0.2359, 0.1180, 0)
virial:  dₓ F_y = 1.0 × 0.1180 = <b>0.1180</b>

P_xy = (−0.2300 + 0.1180) / 27 = <b>−0.00415</b></div>`,
      see: 'One configuration gives a small number of either sign. In the LJ liquid P_xy averages to zero (−0.02 ± 0.27 per step, the mean is pure noise) and fluctuates by ±0.27, larger than any shear stress a flow would create in a simulation. The viscosity is not in the average but in how long a fluctuation lasts. Note that the diagonal pressure of Lab 3.8 is the same formula with x and x instead of x and y.',
    },
    {
      title: 'Green–Kubo by hand',
      why: 'The stress autocorrelation C(t) = ⟨P_xy(0)P_xy(t)⟩ starts at the size of the fluctuations and decays as the liquid relaxes. The longer it takes, the more viscous the liquid. Its integral, times V/k_BT, is η: the same recipe as D from the VACF in Lab 6.7.',
      work: `<div class="calc">P_xy every 0.1 time units (V = 100, k_BT = 1):
0.30 0.24 0.20 0.10 0.02 −0.06 −0.10 −0.08 −0.02 0.04

lag 0:  C = 0.02200   (mean of 10 products)
lag 1:  C = 0.01729   (9)
lag 2:  C = 0.01045   (8)
lag 3:  C = 0.00120   (7)
lag 4:  C = −0.00653  (6)

trapezoid, lags 0–4:
0.1 × (0.01100 + 0.01729 + 0.01045 + 0.00120 − 0.00327)
  = 0.003667
η = V/k_BT × 0.003667 = 100 × 0.003667 = <b>0.367</b>

Einstein: A(t) = ∫P_xy dt = 0, 0.030, 0.054, 0.074, 0.084, ...
⟨ΔA²⟩ at lags 1, 2, 3:  0.000220, 0.000733, 0.001332
slope (lags 2–3) = 0.00600  ->  η = V/(2k_BT) × 0.00600 = <b>0.300</b></div>`,
      see: 'Both routes are a few multiplications. The factor V is the important difference from diffusion: the viscosity is a property of the whole box, one number per time origin, like the conductivity in Lab 6.7, so it is a collective quantity and as noisy as one. With only ten numbers the two routes disagree (0.367 vs 0.300); the question of where to stop the integral is already visible, since C turns negative at lag 4.',
    },
    {
      title: 'The LJ liquid: where to stop the integral',
      why: 'Now the real data: 100 000 steps, all five stress components averaged. The running integral η(t) should rise and level off once C(t) has decayed. In practice the noise in C(t) keeps adding to the integral, so the plateau has to be chosen.',
      work: `<div class="calc">C(0) = 0.07218;  C(t)/C(0) = 0.307 (t = 0.1), 0.101 (0.2), 0.030 (0.5)</div>
<div class="tbl"><table>${head(['upper limit t', '0.05', '0.1', '0.2', '0.5', '1.0', '1.5', '2.0', '4.0', '8.0'])}
${row(['running η(t)', '0.97', '1.52', '1.95', '2.27', '2.40', '<b>2.43</b>', '2.46', '2.86', '3.66'])}</table></div>
<div class="calc">5 blocks of 20 000 steps, integral to t = 1.5:
2.06  2.75  2.22  2.57  2.51   ->  η = <b>2.42 ± 0.12</b>
Einstein (slope of ⟨ΔA²⟩ for t = 1–3):  η = <b>2.44</b></div>`,
      see: 'Most of the viscosity comes from the first 0.2 time units, where the stress correlation decays; by t = 1–2 the integral sits at 2.40–2.46, and both routes agree on about 2.42. Beyond that the integral drifts (2.86 at t = 4, 3.66 at t = 8), not because the liquid has a long memory but because the noisy tail of C(t) is integrated too: the random walk of the noise wins. Practical rules: plot η(t), choose the plateau where C(t) has decayed to the noise, and take the error bar from independent blocks or runs (here ±0.12). In viscous liquids such as ionic liquids the plateau comes much later, and fitting C(t) or using many independent runs becomes necessary.',
    },
    {
      title: 'Cross-checks: Yeh–Hummer and Stokes–Einstein',
      why: 'A viscosity should be consistent with diffusion. Two checks: the box-size dependence of D in Lab 6.9 must follow the Yeh–Hummer slope ξk_BT/(6πη), and Stokes–Einstein, D = k_BT/(cπηR), should give a sensible particle radius R.',
      work: `<div class="calc">predicted Yeh–Hummer slope: ξ T / (6π η) = 2.837 / (6π × 2.42)
                          = <b>0.0621</b>
D from Lab 6.9 + 0.0621/L:
N = 32:   0.0440 + 0.0182 = 0.0622   (box too small)
N = 108:  0.0567 + 0.0121 = <b>0.0688</b>
N = 256:  0.0596 + 0.0091 = <b>0.0687</b>
N = 500:  0.0615 + 0.0073 = <b>0.0688</b>
(fit of the three larger boxes alone gave η = 2.43)

Stokes–Einstein with D_∞ = 0.0688:
stick (c = 6):  R = 1/(6π × 2.42 × 0.0688) = 0.318 σ  (diameter 0.64σ)
slip  (c = 4):  R = 1/(4π × 2.42 × 0.0688) = 0.478 σ  (diameter <b>0.96σ</b>)</div>`,
      see: 'Two independent routes agree: the viscosity from stress fluctuations (2.42) predicts exactly how D changes with box size (the three larger boxes collapse onto 0.0688), and the slope of D vs 1/L gives back 2.43. Stokes–Einstein, a formula for a sphere in a continuum, works surprisingly well for atoms if the "slip" boundary condition is used: it gives a diameter of 0.96σ, almost exactly the LJ size. With "stick" the particle comes out too small. Limitation: Stokes–Einstein is a continuum law; it fails for small solutes in large solvents, in supercooled liquids, and for ions with strong local structure, which is exactly where deviations from it are interesting.',
    },
    {
      title: 'Convergence: why viscosity needs long runs',
      why: 'Viscosity, like conductivity, is a property of the whole box, so it does not get more samples from more atoms. How long must the run be for a given precision, and how much do the five equivalent stress components help?',
      work: `<div class="tbl"><table>${head(['run length (time units)', 'spread of η: 5 components', 'spread: 1 component', 'pieces'])}
${row(['20', '0.72', '1.19', '20'])}
${row(['40', '0.38', '0.74', '10'])}
${row(['80', '0.28', '0.68', '5'])}
${row(['200', '0.07', '0.38', '2 (unreliable)'])}</table></div>
<div class="calc">target 5 % (±0.12): about 400 time units with 5 components
(100 000 steps here; spread falls roughly as 1/√length)
1 component instead of 5: spread about √5 ≈ 2.2 times larger</div>`,
      see: 'A run of 20 time units, plenty for D of this liquid, gives η only to ±30 %; it takes about 400 time units (100 000 steps) to reach ±5 %, and the error falls only as one over the square root of the run length. Using all five shear components is free and is worth five times more simulation. Larger boxes do not help much: the fluctuations shrink as 1/V while the V in front of the integral grows, so the relative noise is about the same. For viscous liquids, where the stress correlation lasts 100–1000 times longer, people use many independent runs, the "time-decomposition" method, or nonequilibrium shear (NEMD) instead.',
    },
    {
      title: 'Units: from LJ to mPa s',
      why: 'Like every reduced quantity, the viscosity becomes physical once we pick a substance. For argon the LJ unit of viscosity is √(mε)/σ².',
      work: `<div class="calc">argon: m = 6.6335×10⁻²⁶ kg, ε = 1.6540×10⁻²¹ J, σ = 3.405×10⁻¹⁰ m
unit = √(m ε) / σ² = √(1.0972×10⁻⁴⁶) / 1.1594×10⁻¹⁹
     = <b>9.035×10⁻⁵ Pa s</b>

η* = 2.42  ->  η = 2.42 × 9.035×10⁻⁵ = <b>0.219 mPa s</b>
(at T = 119.8 K and 0.8/σ³ = 1.344 g/cm³)

real units (GROMACS/LAMMPS): pressure in bar or atm,
V in nm³ or Å³, time in ps or fs: check every factor,
1 bar = 10⁵ Pa,  1 mPa s = 1 cP</div>`,
      see: 'A liquid of LJ argon at these conditions comes out at about 0.22 mPa s, a quarter of water\'s viscosity, a sensible value for a simple liquid of small atoms (with our shortened cutoff it is a model, not a precise prediction for argon). In real-unit codes the Green–Kubo prefactor V/k_BT and the pressure units are the usual sources of factor-of-10⁵ mistakes; checking a code\'s viscosity on the LJ liquid against a published value is a good first test. LAMMPS: <code>fix ave/correlate</code> on the pressure tensor; GROMACS: <code>gmx energy -vis</code>.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method / check', 'merit (number)', 'limitation (number)'])}
${row(['Green–Kubo (stress ACF)', 'η = 2.42 ± 0.12 from one run', 'integral drifts after t ≈ 2 (3.66 at t = 8)'])}
${row(['Einstein (∫P dt)²', '2.44, agrees', 'same noise; slope region to choose'])}
${row(['Yeh–Hummer slope', 'larger boxes give 2.43', 'tiny box (32 atoms) gives a wrong 1.50'])}
${row(['Stokes–Einstein (slip)', 'diameter 0.96σ from η and D', 'continuum law; fails for small solutes, ions'])}
${row(['5 stress components', '√5 smaller error for free', '—'])}
${row(['run length', '±5 % after ~400 time units', '±30 % after 20 (enough for D)'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the shear stress by hand, Green–Kubo and Einstein on ten numbers, the LJ liquid run (about 55 s) with the running integral, block errors and Einstein slope, the Yeh–Hummer and Stokes–Einstein checks, the convergence table and the unit conversion. Run <code>pip install numpy</code>, then <code>python lab_viscosity.py</code> (about a minute).',
    src: `# Lab: viscosity by hand. Green-Kubo and Einstein routes, the 3D LJ liquid of Lab 6.9, Stokes-Einstein, convergence, units.
import numpy as np, time

# ---------- Step 1: the pressure tensor by hand ----------
V = 27.0
r1, r2 = np.array([1.0, 1.0, 1.0]), np.array([2.0, 1.5, 1.0]); v1, v2 = np.array([0.5, -0.3, 0.0]), np.array([-0.2, 0.4, 0.1])
d = r2 - r1; r = np.linalg.norm(d); Fmag = 24 * (2 * r**-13 - r**-7); F = Fmag * d / r          # force on atom 2 (repulsive if > 0)
kin = v1[0] * v1[1] + v2[0] * v2[1]; vir = d[0] * F[1]
print(f"Step 1: two atoms in a box of volume {V}: positions {r1}, {r2}; velocities {v1}, {v2} (mass 1)")
print(f"   kinetic part  sum m vx vy = {v1[0]} x {v1[1]} + {v2[0]} x {v2[1]} = {kin:.4f}")
print(f"   separation d = {d}, r = {r:.4f}, F(r) = {Fmag:.4f}, force on atom 2 = ({F[0]:.4f}, {F[1]:.4f}, {F[2]:.4f})")
print(f"   virial part   d_x F_y = {d[0]} x {F[1]:.4f} = {vir:.4f}")
print(f"   P_xy = ({kin:.4f} + {vir:.4f}) / {V} = {(kin + vir) / V:.5f}")

# ---------- Steps 2-3: Green-Kubo and Einstein by hand ----------
Pt = np.array([0.30, 0.24, 0.20, 0.10, 0.02, -0.06, -0.10, -0.08, -0.02, 0.04])
Vh, kT, dth = 100.0, 1.0, 0.1
print(f"\\nStep 2: P_xy every {dth} time units (V = {Vh}, kT = {kT}):", Pt)
C = [np.mean(Pt[k:] * Pt[:len(Pt) - k]) for k in range(5)]
for k in range(5): print(f"   lag {k}: C = mean of {len(Pt) - k} products = {C[k]:.5f}")
I = dth * (0.5 * C[0] + sum(C[1:4]) + 0.5 * C[4])
print(f"   integral (trapezoid, lags 0-4) = {dth} x ({0.5 * C[0]:.5f} + {C[1]:.5f} + {C[2]:.5f} + {C[3]:.5f} + {0.5 * C[4]:.5f}) = {I:.6f}")
print(f"   Green-Kubo: eta = V/kT x integral = {Vh} x {I:.6f} = {Vh * I:.4f}")
A = np.concatenate([[0], np.cumsum(Pt) * dth])
print(f"\\nStep 3: Einstein form: A(t) = integral of P_xy dt = {np.round(A, 3)}")
mA = []
for lag in (1, 2, 3):
    dA = A[lag:] - A[:-lag]; mA.append(np.mean(dA**2)); print(f"   lag {lag}: <dA^2> = {mA[-1]:.6f}")
sl = (mA[2] - mA[1]) / dth
print(f"   slope between lags 2 and 3 = ({mA[2]:.6f} - {mA[1]:.6f}) / {dth} = {sl:.5f};  eta = V/(2 kT) x slope = {Vh / 2:.0f} x {sl:.5f} = {Vh / 2 * sl:.3f}")

# ---------- Step 4: the LJ liquid ----------
RC, SK, rho, T0, dt = 1.6, 0.3, 0.8, 1.0, 0.004
def fcc(nc, a):
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]); g = np.array([[i, j, k] for i in range(nc) for j in range(nc) for k in range(nc)])
    return ((g[:, None, :] + b[None]) * a).reshape(-1, 3)
def pairs(x, L):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); return np.nonzero(np.triu(r2 < (RC + SK)**2, 1))
def forces(x, L, i, j, N):
    d = x[j] - x[i]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); m = r2 < RC * RC
    ir6 = np.where(m, 1 / r2**3, 0); fs = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0); ff = fs[:, None] * d; f = np.zeros((N, 3))
    for k in range(3): f[:, k] = np.bincount(j, ff[:, k], N) - np.bincount(i, ff[:, k], N)
    return f, np.einsum('p,pa,pb->ab', fs, d, d)
def run(nc, nsteps, seed, tau=1.0):
    N = 4 * nc**3; L = (N / rho)**(1 / 3); V = L**3; x = fcc(nc, L / nc); r = np.random.default_rng(seed); NF = 3 * N - 3
    v = r.normal(size=(N, 3)); v -= v.mean(0); v *= np.sqrt(NF * T0 / (v**2).sum()); i, j = pairs(x, L); xl = x.copy(); f, W = forces(x, L, i, j, N)
    c = np.exp(-dt / tau); out = np.zeros((nsteps, 5))
    for s in range(-3000, nsteps):
        v += 0.5 * dt * f; x += dt * v
        if np.abs((x - xl + L / 2) % L - L / 2).max() > SK / 2: x %= L; i, j = pairs(x, L); xl = x.copy()
        f, W = forces(x, L, i, j, N); v += 0.5 * dt * f
        K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = r.normal(); S2 = r.chisquare(NF - 1)          # gentle Bussi thermostat
        v *= np.sqrt(max(c + (1 - c) * Kt / (NF * K) * (R1**2 + S2) + 2 * R1 * np.sqrt(c * (1 - c) * Kt / (NF * K)), 0))
        if s >= 0:
            P = (v.T @ v + W) / V; out[s] = [P[0, 1], P[0, 2], P[1, 2], 0.5 * (P[0, 0] - P[1, 1]), 0.5 * (P[1, 1] - P[2, 2])]
    return V, out
def acf(P, nmax):
    n = len(P); f = np.fft.rfft(P, 2 * n, axis=0); return np.fft.irfft(f * np.conj(f), axis=0)[:nmax] / (n - np.arange(nmax))[:, None]
t0 = time.time(); Vb, P = run(4, 100000, 1)
print(f"\\nStep 4: 3D LJ liquid of Lab 6.9 (256 atoms, rho = 0.8, T = 1.0, V = {Vb:.1f}), 100 000 steps ({time.time() - t0:.0f} s)")
print(f"   P_xy: mean {P[:, 0].mean():+.4f}, fluctuation +-{P[:, 0].std():.4f};  the 5 shear components (xy, xz, yz, (xx-yy)/2, (yy-zz)/2) all have the same statistics")
ac = acf(P, 2500).mean(1); cum = Vb / T0 * np.cumsum(ac) * dt
print(f"   stress autocorrelation C(0) = {ac[0]:.5f};  C(t)/C(0) at t = 0.1, 0.2, 0.5: {ac[25] / ac[0]:.3f}, {ac[50] / ac[0]:.3f}, {ac[125] / ac[0]:.3f}")
print("   running Green-Kubo integral eta(t) = V/kT x integral_0^t C dt':")
for tt in (0.05, 0.1, 0.2, 0.5, 1.0, 1.5, 2.0, 4.0, 8.0):
    print(f"     t = {tt:4.2f}: {cum[int(round(tt / dt)) - 1]:.3f}")
est = {}
for tt in (0.5, 1.0, 1.5, 2.0):
    k = int(round(tt / dt)); e_ = [Vb / T0 * np.cumsum(acf(P[b * 20000:(b + 1) * 20000], k).mean(1))[-1] * dt for b in range(5)]; est[tt] = e_
    print(f"   5 blocks of 20 000 steps, integral to t = {tt}: {' '.join(f'{q:.2f}' for q in e_)}  ->  eta = {np.mean(e_):.2f} +- {np.std(e_, ddof=1) / np.sqrt(5):.2f}")
Acum = np.cumsum(P, 0) * dt; lags = np.arange(250, 751, 50)
msdA = np.array([((Acum[l:] - Acum[:-l])**2).mean() for l in lags]); slope = np.polyfit(lags * dt, msdA, 1)[0]
print(f"   Einstein: slope of <dA^2> for t = 1-3 = {slope:.5f};  eta = V/(2 kT) x slope = {Vb / 2 * slope:.2f}")
eta = np.mean(est[1.5]); deta = np.std(est[1.5], ddof=1) / np.sqrt(5)

# ---------- Step 5: Yeh-Hummer and Stokes-Einstein ----------
xi = 2.837297; D_lab = {32: (3.420, 0.0440), 108: (5.130, 0.0567), 256: (6.840, 0.0596), 500: (8.550, 0.0615)}
print(f"\\nStep 5: compare with Lab 6.9 (D in four box sizes). Yeh-Hummer slope predicted from eta = {eta:.2f}: xi T / (6 pi eta) = {xi / (6 * np.pi * eta):.4f}")
for N, (L, D) in D_lab.items(): print(f"   N = {N:3d}: D_PBC + slope/L = {D:.4f} + {xi / (6 * np.pi * eta * L):.4f} = {D + xi / (6 * np.pi * eta * L):.4f}")
Dinf = np.mean([D + xi / (6 * np.pi * eta * L) for N, (L, D) in D_lab.items() if N > 32])
print(f"   average of the three larger boxes: D_inf = {Dinf:.4f}")
for name, cst in (("stick, 6 pi", 6), ("slip, 4 pi", 4)):
    R = T0 / (cst * np.pi * eta * Dinf); print(f"   Stokes-Einstein, {name}: R = kT / ({cst} pi eta D) = 1 / ({cst} pi x {eta:.2f} x {Dinf:.4f}) = {R:.3f} sigma  (diameter {2 * R:.2f})")

# ---------- Step 6: convergence ----------
print("\\nStep 6: how the error shrinks (blocks of the same run; 5 components averaged vs one component)")
for nb_len in (5000, 10000, 20000, 50000):
    nb = len(P) // nb_len; k = int(1.5 / dt)
    e5 = [Vb * np.cumsum(acf(P[b * nb_len:(b + 1) * nb_len], k).mean(1))[-1] * dt for b in range(nb)]
    e1 = [Vb * np.cumsum(acf(P[b * nb_len:(b + 1) * nb_len, :1], k).mean(1))[-1] * dt for b in range(nb)]
    print(f"   runs of {nb_len:6d} steps ({nb_len * dt:5.0f} time units): spread of eta, 5 components {np.std(e5, ddof=1):.2f}, one component {np.std(e1, ddof=1):.2f}   ({nb} pieces)")

# ---------- Step 7: units ----------
m, eps, sig = 39.948 * 1.66054e-27, 119.8 * 1.380649e-23, 3.405e-10
unit = np.sqrt(m * eps) / sig**2
print(f"\\nStep 7: LJ viscosity unit sqrt(m eps)/sigma^2 for argon = sqrt({m:.4e} x {eps:.4e}) / ({sig:.3e})^2 = {unit:.3e} Pa s")
print(f"   eta* = {eta:.2f}  ->  eta = {eta:.2f} x {unit:.3e} = {eta * unit * 1e3:.3f} mPa s  (at T = 119.8 K, rho = 0.8/sigma^3 = {0.8 / (sig**3) / 6.022e23 * 39.948 / 1e6:.3f} g/cm^3)")
`,
  },

  exercises: [
    { q: 'Water at 298 K has η = 0.89 mPa s and D = 2.3×10⁻⁹ m²/s. What radius does Stokes–Einstein (stick) give? Is it reasonable?', a: 'R = k_BT/(6πηD) = 4.114×10⁻²¹ / (6π × 0.00089 × 2.3×10⁻⁹) = 4.114×10⁻²¹ / 3.858×10⁻¹¹ = 1.07×10⁻¹⁰ m = 0.107 nm. A water molecule has a radius of about 0.14 nm, so stick gives a slightly small value (slip, c = 4, gives 0.16 nm): Stokes–Einstein works roughly even for the solvent itself.' },
    { q: 'A Green–Kubo integral for an ionic liquid keeps rising: 25, 38, 44, 47, 51 mPa s at upper limits of 0.1, 0.5, 1, 2, 4 ns, with block errors of ±5 mPa s at 1 ns and ±12 at 4 ns. What do you report?', a: 'Beyond about 1–2 ns the changes (44 → 47 → 51) are within the growing noise, so a plateau around 45–47 mPa s with an error of at least ±5–8 mPa s is the honest result. Better: more independent runs, or fit C(t) with a stretched exponential and integrate the fit, and state the method.' },
    { q: 'Why does the viscosity not converge faster in a box with 8 times more atoms, while D does?', a: 'D is averaged over all atoms, so 8 times more atoms gives about 8 times more samples. P_xy is one number for the whole box: its fluctuations shrink as 1/V, but the factor V in η = V/k_BT ∫C dt cancels that, so the relative noise per unit time is about the same. Only longer (or more) runs help.' },
  ],
};
