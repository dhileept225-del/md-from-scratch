import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-barostats',
  lab: true,
  summary: 'A barostat changes the box size so that the pressure matches a target. Like a thermostat, it has two jobs: reach the target, and then let the box breathe by exactly the right amount. We compute the pressure of our 9-atom box by hand, work out how much its area should fluctuate, and test Berendsen, C-rescale and the Parrinello–Rahman piston against that number.',
  minutes: 35,
  sim: 'lab3b',
  learn: [
    'Compute the pressure of a configuration by hand: kinetic part plus virial part',
    'Get the correct size of the box fluctuations from the equation of state, using only addition and an exponential',
    'Apply Berendsen, C-rescale and a piston by hand, and use long runs to see which ones sample the right volumes',
  ],
  simHint: 'Pick a barostat and press Step to see one box update at a time. Press Run: the bars show how often each box area occurs, the amber line is the correct distribution from step 2. Berendsen gives a narrow peak shifted to the left; C-rescale and the piston follow the line. Start from A = 14 to compare equilibration: the light piston overshoots, the heavy one swings slowly.',

  setup: `
<p>Same <b>9 atoms</b> in a 2D periodic box. The temperature is held at <b>T₀ = 0.72</b> by a Langevin thermostat on the atoms (Lab 3.7 showed it gives the right fluctuations for any system). Now the box side L is allowed to change, so that the pressure is <b>P₀ = 3.0</b>.</p>
<div class="calc">atoms          N = 9, 2D, periodic, square box L × L
area           A = L²   (the 2D "volume")
temperature    T₀ = 0.72, Langevin γ = 1 on the atoms
target         P₀ = 3.0
time step      dt = 0.005
cutoff         r_c = 1.5 (truncated and shifted)
               the box may shrink to L = 2 r_c = 3.0
               before the minimum image breaks (Lab 2)</div>
<p>Why P₀ = 3.0 and not 0? In two dimensions, Lennard-Jones atoms at T = 0.72 are above the critical temperature: at zero pressure they would fly apart like a gas and the box would grow without limit. A pressure of 3 keeps them a dense fluid.</p>`,

  steps: [
    {
      title: 'Pressure of one configuration, by hand',
      why: 'A barostat can only react to a pressure it can measure. In a simulation there are no walls to push on, so the pressure is computed from two things we already have: how fast the atoms move (kinetic part) and how hard they push or pull on each other (virial part).',
      work: `<div class="calc">2D:  P = ( Σ m v²  +  Σ_pairs r·F ) / (2A)
average of Σ m v² = N T₀ (Langevin: N_f = 2N)

lattice, L = 3.6, A = 12.96
kinetic part  N T₀ / A = 9 × 0.72 / 12.96 = <b>0.500</b>

pairs inside r_c = 1.5: 18 at r = 1.2
(the diagonal 1.697 is beyond the cutoff)
F(1.2) = 24(2/1.2¹³ − 1/1.2⁷) = −2.2117  (attraction)
virial part  ½ × 18 × 1.2 × (−2.2117) / 12.96
           = <b>−1.843</b>

P = 0.500 − 1.843 = <b>−1.343</b></div>
<p>Now let the atoms move in the same box (NVT, 290 time units):</p>
<div class="calc">average P                     = <b>1.074</b>
instantaneous P fluctuates by ± <b>1.189</b>
averages over 1 time unit (200 steps) still ± 0.708</div>`,
      see: 'The frozen lattice says −1.343 (the box wants to shrink); the moving atoms say +1.074 (it wants to grow). A single configuration is not "the pressure". Worse, the instantaneous pressure jumps by more than its own average every step. A barostat that reacted fully to each step would shake the box violently, so every barostat responds slowly, through a coupling time τ_P or a piston mass W. It also means that a pressure from a short run is unreliable: always quote averages over long runs.',
    },
    {
      title: 'How much should the box breathe?',
      why: 'At constant pressure the box area is not fixed; it fluctuates. As for temperature in Lab 3.7, we need the correct amount before we can judge any barostat. It follows from the equation of state P(A): how the average pressure depends on the box area.',
      work: `<p>Run NVT at a series of areas (each 300 time units):</p>
<div class="tbl"><table>${head(['A', 'L', 'average P'])}
${row(['9.6', '3.098', '10.512'])}
${row(['10.0', '3.162', '6.946'])}
${row(['10.4', '3.225', '4.513'])}
${row(['10.8', '3.286', '3.046'])}
${row(['11.0', '3.317', '2.673'])}
${row(['11.2', '3.347', '2.080'])}
${row(['11.6', '3.406', '1.715'])}
${row(['12.0', '3.464', '1.316'])}
${row(['12.4', '3.521', '1.193'])}</table></div>
<div class="calc">rule: growing the box by dA does work (P − P₀) dA.
probability changes by  exp( (P − P₀) dA / T₀ )

A 11.0 -> 11.2:  mean P = (2.673 + 2.080)/2 = 2.377
(2.377 − 3.0) × 0.2 / 0.72 = −0.173
probability × e^(−0.173) = × <b>0.841</b>

chain these factors over the whole table:
correct  ⟨A⟩ = <b>11.006</b>,   σ_A = <b>0.564</b>

quick estimate from the slope near P = P₀ (A = 10.8):
dP/dA = −2.55
compressibility κ = −1/(A dP/dA) = 0.0363
σ_A = √(T₀ A κ) = √(0.72 × 10.8 × 0.0363) = 0.532</div>`,
      see: 'Where the box pressure is above P₀ the box tends to grow, where it is below it tends to shrink, and the exponential tells us by how much: this gives the full distribution of areas, the amber curve in the widget. The box should fluctuate by <b>σ_A = 0.564</b> (5 % of its area) around <b>⟨A⟩ = 11.006</b>. Note that the curve is lopsided: P(A) is steep at small A and flat at large A, so the average area is larger than the area where P = P₀ (about 10.83). This is the test for every barostat below, and the size of the fluctuation is how compressibilities are computed from simulations.',
    },
    {
      title: 'Berendsen barostat',
      why: 'The simplest idea, the same as the Berendsen thermostat: scale the box (and all coordinates) a little each step, in the direction that brings P towards P₀. The user supplies a coupling time τ_P and a guess of the compressibility β.',
      work: `<div class="calc">μ² = 1 − β (dt/τ_P)(P₀ − P)       (area factor)
β = 0.03,  τ_P = 0.5  ->  β dt/τ_P = 0.0003

from the lattice, P = −1.343:
μ² = 1 − 0.0003 × (3.0 + 1.343) = <b>0.998697</b>
μ  = 0.999348
L: 3.6 -> 3.59765     A: 12.96 -> <b>12.9431</b>
every x and y is multiplied by μ too</div>
<div class="tbl"><table>${head(['long run (700 time units)', '⟨A⟩', 'σ_A'])}
${row(['correct', '11.006', '0.564'])}
${row(['Berendsen τ_P = 0.5', '<b>10.832</b>', '<b>0.247</b>'])}
${row(['Berendsen τ_P = 5', '<b>10.846</b>', '<b>0.110</b>'])}</table></div>`,
      see: 'Merit: smooth and safe, it never overshoots, so it is a good way to bring the density to the right value (step 6). Limitation 1: the box fluctuates by 0.247 (τ_P = 0.5) or 0.110 (τ_P = 5) instead of 0.564, so any compressibility computed from it would be 5 to 25 times too small (σ² enters). Limitation 2: the average area is wrong too, 10.83 instead of 11.01, a 1.6 % error in density. Berendsen holds the box near the point where P = P₀ and never explores the lopsided tail of step 2. Weaker coupling makes the fluctuation <i>smaller</i>, not better. Not for production.',
    },
    {
      title: 'C-rescale: Berendsen plus the missing randomness',
      why: 'Exactly the fix of the Bussi thermostat, applied to the box (Bernetti and Bussi, 2020). Keep Berendsen\'s smooth pull, add a random kick of the right size, and a small correction T₀/A. Working with ln A makes the update a simple sum.',
      work: `<div class="calc">d(ln A) = −(β/τ_P)(P₀ − P − T₀/A) dt
          + √(2 T₀ β dt / (A τ_P)) × R

drift  = −(0.03/0.5)(3.0 + 1.343 − 0.0556) × 0.005
       = <b>−0.001286</b>      (= Berendsen)
noise  = √(2 × 0.72 × 0.03 × 0.005 / (12.96 × 0.5))
       = <b>0.005774</b>

R = +1:  d(ln A) = +0.004487  ->  A = 13.0183
R = −1:  d(ln A) = −0.007060  ->  A = 12.8688
R =  0:  d(ln A) = −0.001286  ->  A = 12.9433
coordinates × e^(d(ln A)/2), velocities × e^(−d(ln A)/2)</div>
<div class="tbl"><table>${head(['long run', '⟨A⟩', 'σ_A'])}
${row(['correct', '11.006', '0.564'])}
${row(['Berendsen τ_P = 0.5', '10.832', '0.247'])}
${row(['C-rescale τ_P = 0.5', '<b>11.040</b>', '<b>0.587</b>'])}</table></div>`,
      see: 'Merit: the noise is 4.5 times larger than the drift in a single step, yet on average the drift wins, so equilibration is as smooth as Berendsen\'s, and the fluctuation (0.587) and average (11.04) now agree with the correct values within the noise of a 700-time-unit run. It is first order (no piston to oscillate) and works for any τ_P. Limitation: like Bussi, it is stochastic, so trajectories are reproducible only with the same seed, and the compressibility β is still a user guess (it sets the speed, not the result). GROMACS: <code>pcoupl = C-rescale</code> (available since version 2021).',
    },
    {
      title: 'The piston: Parrinello–Rahman (Andersen in 1 dimension)',
      why: 'Instead of nudging the box, give it a mass W and let it move like a piston under Newton\'s law: the force on the piston is the pressure difference. This is deterministic, like Nosé–Hoover, and with a correctly thermostatted system it samples the right volumes. The new question is how heavy the piston should be.',
      work: `<div class="calc">W d²A/dt² = P − P₀

from the lattice (P = −1.343), W = 0.1:
box acceleration = (−1.343 − 3.0)/0.1 = <b>−43.43</b>
box velocity after 1 step   −0.217
             after 10 steps −2.17  (it keeps speeding up)

near equilibrium the piston oscillates like a spring:
period = 2π √(W κ A),  κ A = 0.0363 × 10.8 = 0.392</div>
<div class="tbl"><table>${head(['W', 'period: hand', 'period: measured', '⟨A⟩', 'σ_A'])}
${row(['correct', '—', '—', '11.006', '0.564'])}
${row(['0.01', '0.39', '0.29', '11.052', '0.583'])}
${row(['0.1', '1.24', '1.02', '11.054', '0.574'])}
${row(['2', '5.57', '4.76', '11.050', '0.559'])}</table></div>`,
      see: 'Merit: all three masses give the correct fluctuation (0.56–0.58) and average (11.05), and the dynamics are smooth and deterministic. The hand formula gets the oscillation period within 20 % (the rest is because P(A) is not a straight line). The mass only sets <i>how fast</i> the box breathes: W = 0.01 oscillates every 58 steps, which is too close to the atomic motions and disturbs them; W = 2 needs about 5 time units per breath, so a run must be many periods long to sample the volume. A common starting point: LAMMPS suggests a damping time of about 1000 time steps. The real limitation appears when you start far from equilibrium (next step). The full Parrinello–Rahman method lets the box <i>shape</i> change too (solids, phase transitions); MTK (LAMMPS <code>fix npt</code>) is the version with Nosé–Hoover chains and exact equations.',
    },
    {
      title: 'Starting far away: equilibration',
      why: 'Production barostats are judged by the fluctuations; equilibration barostats by how safely they reach the target. We start from a box that is far too big, A = 14 (P ≈ 1.1, well below 3.0), and watch the area.',
      work: `<div class="tbl"><table>${head(['time', '0.5', '1', '2', '5', '10', '20'])}
${row(['Berendsen τ_P 0.5', '12.62', '11.53', '10.98', '10.89', '10.92', '10.68'])}
${row(['C-rescale τ_P 0.5', '11.17', '11.16', '10.54', '10.55', '11.94', '11.49'])}
${row(['piston W 0.1', '<b>9.99</b>', '12.25', '11.53', '10.82', '10.88', '11.25'])}
${row(['piston W 2', '13.76', '13.12', '11.14', '<b>12.76</b>', '11.09', '11.23'])}</table></div>
<div class="calc">correct range: 11.0 ± 0.56  ->  about 10.4 to 11.6
A = 9.99:  P ≈ 7 from the table in step 2
           more than twice the target</div>`,
      see: 'Berendsen glides smoothly into the right range within 1–2 time units and stays there (with the wrong, too-small fluctuations). C-rescale is just as fast and then fluctuates correctly. The light piston (W = 0.1) gathers speed and overshoots to A = 9.99 within half a time unit, where the pressure is more than twice the target; in a real system with a big initial mismatch this can crash the run or melt a solid. The heavy piston (W = 2) swings back to 12.76 at t = 5 and takes many periods to settle. Practical rule: equilibrate with Berendsen or C-rescale, then switch to a piston for production, or simply use C-rescale for both.',
    },
    {
      title: 'What the cutoff leaves out of the pressure',
      why: 'Our forces stop at r_c = 1.5, but the real attraction continues beyond it. Every pair beyond the cutoff would pull the box inward a little. Because a barostat sets the box size from the pressure, a missing piece of pressure becomes a wrong density.',
      work: `<div class="calc">2D tail pressure (atoms beyond r_c spread evenly):
P_tail = −(π ρ² / 2) ∫ r² U′(r) dr,  from r_c to ∞
∫ = −24 (0.2 r_c⁻¹⁰ − 0.25 r_c⁻⁴) = 1.1019

ρ = 9 / 11.01 = 0.818
P_tail = −(π × 0.818² / 2) × 1.1019 = <b>−1.157</b>

effect on the box:  dA = κ A P_tail
  = 0.0363 × 11.01 × (−1.157) = <b>−0.46</b>  (−4.2 %)

3D LJ liquid, ρ = 0.8, r_c = 2.5:
P_tail = (16/3) π ρ² [ (2/3) r_c⁻⁹ − r_c⁻³ ] = <b>−0.684</b></div>`,
      see: 'With our short cutoff the missing attraction is −1.157, a third of the target pressure, and the box comes out 4.2 % too large. Even the standard 3D choice r_c = 2.5 misses −0.68 in a dense liquid, which is comparable to the pressure itself. Merit of the tail correction: it is one formula, added once. Limitation: it assumes the atoms beyond r_c are spread evenly, which fails at interfaces (liquid–vapour, membranes); there use a longer cutoff or a long-range method (LAMMPS <code>pair_modify tail yes</code>, GROMACS <code>DispCorr = EnerPres</code>).',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)', 'use it for'])}
${row(['Berendsen', 'smooth: 14 -> 11 in 2 time units', 'σ_A 0.247 vs 0.564; ⟨A⟩ 1.6 % off', 'equilibration only'])}
${row(['C-rescale', 'smooth and σ_A 0.587, ⟨A⟩ 11.04', 'stochastic; needs a seed', 'equilibration and production'])}
${row(['piston (PR / MTK)', 'σ_A 0.56–0.58, deterministic', 'overshoots from far away (A = 9.99)', 'production, after equilibration'])}
${row(['piston mass W', 'period 2π√(W κ A)', 'too light: shakes atoms; too heavy: slow', 'start near 1000 steps (LAMMPS)'])}
${row(['tail correction', 'restores −1.157 of pressure', 'wrong at interfaces', 'homogeneous liquids'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the pressure by hand, the equation of state and the correct area distribution, the first update of each barostat, the long runs (100 000+ steps each), the equilibration test and the tail correction. Run <code>pip install numpy</code>, then <code>python lab_barostats.py</code> (about a minute).',
    src: `# Lab 3.8: barostats by hand. 9 LJ atoms in a 2D periodic box, reduced units; temperature held at T0 = 0.72 by a
# Langevin thermostat (gamma = 1) on the atoms, target pressure P0 = 3.0. Cutoff rc = 1.5 (truncated and shifted),
# so the box may shrink down to L = 2 rc = 3.0 and the minimum image still works.
import numpy as np, time
N, RC, T0, P0, dt = 9, 1.5, 0.72, 3.0, 0.005
BETA = 0.03                                     # compressibility we tell the Berendsen and C-rescale barostats

def forces(x, L):                               # forces and the pair virial sum r.F (positive = repulsion)
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < RC * RC
    ir6 = np.where(m, 1 / np.where(m, r2, 1)**3, 0); ff = np.where(m, 24 * (2 * ir6**2 - ir6) / np.where(m, r2, 1), 0)
    return -(ff[:, :, None] * d).sum(1), 0.5 * (ff * np.where(m, r2, 0)).sum()

def lattice(A):
    a = np.sqrt(A) / 3; return np.array([[(i + .5) * a, (j + .5) * a] for j in range(3) for i in range(3)])

pressure = lambda v, w, A: ((v**2).sum() + w) / (2 * A)          # 2D: P = (sum m v^2 + sum r.F) / (2 A)
c = np.exp(-1.0 * dt); s = np.sqrt((1 - c * c) * T0)                # Langevin, gamma = 1

# ---------- Step 1: pressure of one configuration by hand ----------
F = lambda r: 24 * (2 * r**-13 - r**-7)
A = 3.6**2; kin = N * T0 / A; vir = 0.5 * 18 * 1.2 * F(1.2) / A
print(f"Step 1: lattice, L = 3.6, A = {A:.2f}. Kinetic part N T0 / A = {N} x {T0} / {A:.2f} = {kin:.3f}")
print(f"   18 pairs at r = 1.2 (1.697 is beyond rc = 1.5): F(1.2) = {F(1.2):.4f};  virial part (1/2)(18 x 1.2 x {F(1.2):.4f}) / {A:.2f} = {vir:.3f}")
print(f"   P = {kin:.3f} + ({vir:.3f}) = {kin + vir:.3f}")

def nvt(A, n=60000, seed=0):
    r = np.random.default_rng(seed); L = np.sqrt(A); x = lattice(A); v = r.normal(size=(N, 2)) * np.sqrt(T0); f, w = forces(x, L); P = []
    for k in range(n):
        v += 0.5 * dt * f; x += dt * v; x %= L; f, w = forces(x, L); v += 0.5 * dt * f; v = c * v + s * r.normal(size=v.shape)
        if k >= 2000: P.append(pressure(v, w, A))
    P = np.array(P); return P.mean(), P.std(), P

m, sd, Ptr = nvt(12.96)
print(f"   moving atoms at the same box (NVT, 290 time units): average P = {m:.3f}, instantaneous P fluctuates by +-{sd:.3f}")
print(f"   one step tells you almost nothing: averages over 1 time unit (200 steps) still scatter by +-{Ptr[:len(Ptr) // 200 * 200].reshape(-1, 200).mean(1).std():.3f}")

# ---------- Step 2: how much should the box fluctuate? Equation of state P(A) from NVT runs ----------
t0 = time.time()
grid = np.round(np.arange(9.2, 13.01, 0.2), 1); Pg = np.array([nvt(a)[0] for a in grid])
print(f"\\nStep 2: equation of state from NVT runs ({time.time() - t0:.0f} s)")
for a, p in zip(grid, Pg): print(f"   A = {a:5.1f}  (L = {np.sqrt(a):.3f})   P = {p:7.3f}")
i = int(np.argmin(abs(grid - 11.0)))
print(f"   step 11.0 -> 11.2: mean P = ({Pg[i]:.3f} + {Pg[i + 1]:.3f})/2 = {(Pg[i] + Pg[i + 1]) / 2:.3f};  (P - P0) dA / T0 = ({(Pg[i] + Pg[i + 1]) / 2:.3f} - {P0}) x 0.2 / {T0} = {((Pg[i] + Pg[i + 1]) / 2 - P0) * 0.2 / T0:.3f}  ->  probability x {np.exp(((Pg[i] + Pg[i + 1]) / 2 - P0) * 0.2 / T0):.3f}")
lnp = np.concatenate([[0], np.cumsum(0.5 * (Pg[1:] + Pg[:-1] - 2 * P0) * 0.2)]) / T0
p = np.exp(lnp - lnp.max()); p /= p.sum(); Aex = (p * grid).sum(); sex = np.sqrt((p * (grid - Aex)**2).sum())
print(f"   correct NPT distribution: <A> = {Aex:.3f}, sigma_A = {sex:.3f}")
j = int(np.argmin(abs(Pg - P0))); slope = (Pg[j + 1] - Pg[j - 1]) / 0.4; kap = -1 / (grid[j] * slope)
print(f"   simple estimate: P = P0 near A = {grid[j]}, slope dP/dA = {slope:.2f}, kappa_T = -1/(A dP/dA) = {kap:.4f}, sigma_A = sqrt(T0 A kappa) = {np.sqrt(T0 * grid[j] * kap):.3f}")
print("   p(A) on the grid: " + " ".join(f"{q:.4f}" for q in p))

# ---------- Steps 3-5: hand updates of each barostat, from the lattice at L = 3.6 (P = -1.343) ----------
P = kin + vir; A = 12.96; tau = 0.5
mu2 = 1 - BETA * dt / tau * (P0 - P)
print(f"\\nStep 3: Berendsen, beta = {BETA}, tau_P = {tau}: mu^2 = 1 - beta (dt/tau_P)(P0 - P) = 1 - {BETA} x {dt / tau} x ({P0} - ({P:.3f})) = {mu2:.6f}")
print(f"   mu = {np.sqrt(mu2):.6f}: L = 3.6 -> {3.6 * np.sqrt(mu2):.5f},  A = {A} -> {A * mu2:.4f}")
drift = -(BETA / tau) * (P0 - P - T0 / A) * dt; noise = np.sqrt(2 * T0 * BETA * dt / (A * tau))
print(f"\\nStep 4: C-rescale: d(ln A) = -(beta/tau_P)(P0 - P - T0/A) dt + sqrt(2 T0 beta dt / (A tau_P)) R")
print(f"   drift = -({BETA}/{tau}) x ({P0} - ({P:.3f}) - {T0 / A:.4f}) x {dt} = {drift:.6f};  noise size = {noise:.6f}")
for R in (1.0, -1.0, 0.0): print(f"   R = {R:+.1f}: d(ln A) = {drift + noise * R:+.6f}  ->  A = {A * np.exp(drift + noise * R):.4f}")
W = 0.1
print(f"\\nStep 5: piston (Andersen / isotropic Parrinello-Rahman), mass W = {W}: W d2A/dt2 = P - P0")
print(f"   box acceleration = ({P:.3f} - {P0}) / {W} = {(P - P0) / W:.2f};  after one step the box velocity is {(P - P0) / W * dt:.4f}, after 10 steps (same P) {(P - P0) / W * dt * 10:.3f}")
for Wv in (0.01, 0.1, 2.0): print(f"   W = {Wv:5.2f}: oscillation period 2 pi sqrt(W kappa A) = 2 pi sqrt({Wv} x {kap:.4f} x {grid[j]}) = {2 * np.pi * np.sqrt(Wv * kap * grid[j]):.2f}")

# ---------- long runs with each barostat ----------
def run(method, n=150000, tau=0.5, W=0.1, A0=11.0, seed=1, record=None):
    r = np.random.default_rng(seed); A = A0; L = np.sqrt(A); x = lattice(A); v = r.normal(size=(N, 2)) * np.sqrt(T0)
    f, w = forces(x, L); P = pressure(v, w, A); Pi = 0.0; As, Ps, rec = [], [], []
    for k in range(n):
        if method == 'piston': Pi += 0.5 * dt * (P - P0); v *= np.exp(-0.25 * dt * Pi / W / A)
        v += 0.5 * dt * f; x += dt * v
        if method == 'piston': An = A + dt * Pi / W; x *= np.sqrt(An / A); A = An; L = np.sqrt(A)
        x %= L; f, w = forces(x, L); v += 0.5 * dt * f
        if method == 'piston': v *= np.exp(-0.25 * dt * Pi / W / A)
        v = c * v + s * r.normal(size=v.shape); P = pressure(v, w, A)
        if method == 'piston': Pi += 0.5 * dt * (P - P0)
        if method == 'ber': mu2 = 1 - BETA * dt / tau * (P0 - P); x *= np.sqrt(mu2); A *= mu2; L = np.sqrt(A)
        if method == 'cres':
            de = -(BETA / tau) * (P0 - P - T0 / A) * dt + np.sqrt(2 * T0 * BETA * dt / (A * tau)) * r.normal()
            x *= np.exp(de / 2); v *= np.exp(-de / 2); A *= np.exp(de); L = np.sqrt(A)
        if method != 'piston': f, w = forces(x, L); P = pressure(v, w, A)
        if record and (k + 1) in record: rec.append(A)
        if k >= 10000: As.append(A)
    return (np.array(As), rec)

print(f"\\nLong runs (700 time units each), starting at A = 11. Correct: <A> = {Aex:.3f}, sigma_A = {sex:.3f}")
print("   barostat               <A>      sigma_A   volume period")
t0 = time.time()
for name, m, kw in [("Berendsen tau 0.5", 'ber', {}), ("Berendsen tau 5", 'ber', dict(tau=5.0)), ("C-rescale tau 0.5", 'cres', {}),
                    ("piston W 0.01", 'piston', dict(W=0.01)), ("piston W 0.1", 'piston', dict(W=0.1)), ("piston W 2", 'piston', dict(W=2.0))]:
    As, _ = run(m, **kw); per = "-"
    if m == 'piston':                            # period from the box velocity: its autocorrelation is lowest after half a period
        u = np.diff(As); ac = np.array([np.dot(u[:len(u) - k], u[k:]) for k in range(1, 2000)]); per = f"{2 * (np.argmin(ac) + 1) * dt:.2f}"
    print(f"   {name:20s} {As.mean():7.3f}   {As.std():7.3f}     {per}")
print(f"   ({time.time() - t0:.0f} s)")

# ---------- Step 6: equilibration from a box that is far too big ----------
marks = [int(tt / dt) for tt in (0.5, 1, 2, 5, 10, 20)]
print(f"\\nStep 6: start at A = 14.0 (P about 1.1, far below P0). Box area A at t = 0.5, 1, 2, 5, 10, 20:")
for name, m, kw in [("Berendsen tau 0.5", 'ber', {}), ("C-rescale tau 0.5", 'cres', {}), ("piston W 0.1", 'piston', dict(W=0.1)), ("piston W 2", 'piston', dict(W=2.0))]:
    _, rec = run(m, n=marks[-1], A0=14.0, record=set(marks), **kw)
    print(f"   {name:20s} " + "  ".join(f"{a:6.2f}" for a in rec))

# ---------- Step 7: what the cutoff leaves out ----------
rho = N / Aex; I2 = -24 * (0.2 * RC**-10 - 0.25 * RC**-4)          # integral of r^2 U'(r) from rc to infinity
Pt = -(np.pi * rho**2 / 2) * I2
print(f"\\nStep 7: tail pressure in 2D: P_tail = -(pi rho^2 / 2) x integral(r^2 U' dr, rc..inf);  integral = {I2:.4f}")
print(f"   rho = {N}/{Aex:.2f} = {rho:.3f}:  P_tail = {Pt:.3f}  ->  box change dA = kappa A P_tail = {kap:.4f} x {Aex:.2f} x {Pt:.3f} = {kap * Aex * Pt:.2f} ({100 * kap * Pt:.1f} %)")
r3 = 2.5; Pt3 = 16 / 3 * np.pi * 0.8**2 * (2 / 3 * r3**-9 - r3**-3)
print(f"   3D LJ liquid, rho = 0.8, rc = 2.5: P_tail = (16/3) pi rho^2 [(2/3) rc^-9 - rc^-3] = {Pt3:.3f}")
`,
  },

  exercises: [
    { q: 'For the lattice at L = 3.6 with the longer cutoff r_c = 1.8 (Lab 3.7), the 18 diagonal pairs at r = 1.697 also count, with F(1.697) = −0.5425. What is the pressure now?', a: 'Extra virial ½ × 18 × 1.697 × (−0.5425) / 12.96 = −0.639. P = 0.500 − 1.843 − 0.639 = −1.982, the value found in Lab 3.7. A longer cutoff includes more attraction and lowers the pressure.' },
    { q: 'Using the table in step 2, by what factor does the probability change going from A = 10.8 to 11.0?', a: 'Mean P = (3.046 + 2.673)/2 = 2.860. (2.860 − 3.0) × 0.2 / 0.72 = −0.039. Factor e^(−0.039) = 0.962: almost the same probability, because the pressure there is close to P₀. Compare 0.841 for the next step, where P is further below P₀.' },
    { q: 'A 3D water box has V = 27 nm³ at 300 K, and water\'s compressibility is 4.5×10⁻¹⁰ Pa⁻¹. How large should the volume fluctuation be? (k_B = 1.381×10⁻²³ J/K)', a: 'σ_V² = k_B T V κ = 1.381×10⁻²³ × 300 × 27×10⁻²⁷ × 4.5×10⁻¹⁰ = 5.03×10⁻⁵⁶ m⁶, so σ_V = 2.24×10⁻²⁸ m³ = 0.224 nm³, about 0.8 % of the volume. A Berendsen barostat would give less than this, so a compressibility computed from it would come out too small.' },
  ],
};
