import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-thermostats',
  lab: true,
  summary: 'A thermostat must do two jobs: bring the system to the target temperature, and then let the temperature fluctuate exactly as it would in a real system touching a heat bath. We take the 9-atom box of Lab 2, start it too hot, and apply every common thermostat by hand. The numbers show which ones do both jobs, which ones only do the first, and what each one costs.',
  minutes: 35,
  sim: 'lab3',
  learn: [
    'Compute the natural size of temperature fluctuations, and why a correct thermostat must keep them',
    'Apply rescaling, Berendsen, Nosé–Hoover, Bussi and Langevin by hand for a few steps',
    'Use long runs to measure what each thermostat gets right or wrong: fluctuations, energy sharing, diffusion',
  ],
  simHint: 'Pick a thermostat and press Step to see the thermostat act step by step from the hot start. The numbers differ from the hand tables because the widget also moves the atoms, so the forces change T too. Press Run and wait: the bars should follow the amber line if the thermostat is correct. Berendsen with τ = 0.1 gives a narrow histogram from the start; keep it running to time ≈ 400 and watch the sliding-mode share climb towards 1 while the atoms stop vibrating (the flying ice cube). The readouts average over the last 100 time units. Switch to "1 oscillator": Nosé–Hoover leaves holes in phase space, the NH chain and Langevin fill it.',

  setup: `
<p>We reuse the <b>9 atoms</b> of Lab 2 on the 3×3 lattice in the periodic box <b>L = 3.6</b>, with the same Gaussian velocities (drift removed). This time the velocities are scaled so that the system starts <b>too hot</b>, at T = 2.0, and the thermostat has to bring it down to <b>T₀ = 0.72</b>.</p>
<div class="calc">atoms          N = 9, 2D, periodic, L = 3.6
start          lattice, T = 2.0 (too hot)
target         T₀ = 0.72
time step      dt = 0.005
coupling time  τ = 0.1  (20 time steps)
degrees of     N_f = 2N − 2 = 16
freedom        (momentum is conserved, Lab 2)

temperature    T = Σ v² / N_f   (mass 1, k_B = 1)</div>
<p>One idea runs through the whole lab. Temperature is an <i>average</i> property; the kinetic energy of a finite system at a fixed temperature is <i>not</i> constant, it fluctuates. A thermostat that removes these fluctuations gives the right average with the wrong physics.</p>`,

  steps: [
    {
      title: 'How much should the temperature fluctuate?',
      why: 'Before we can judge a thermostat, we need the right answer. A system in contact with a heat bath (the canonical, NVT ensemble) exchanges energy with it all the time, so its kinetic energy goes up and down. How much it goes up and down is fixed by statistics, and we can compute it by hand.',
      work: `<div class="calc">hot start:  Σ v² = 32.000,  T = 32 / 16 = <b>2.000</b>

each velocity component is Gaussian; for a Gaussian
   ⟨v⁴⟩ = 3 ⟨v²⟩²   (check from 10⁶ samples: 3.002)
so one component has
   variance of v² = ⟨v⁴⟩ − ⟨v²⟩² = 2 ⟨v²⟩²

K is a sum of N_f independent v²/2 terms:
   mean of K      = N_f ⟨v²⟩/2
   std. dev. of K = √N_f × √2 ⟨v²⟩/2

relative fluctuation = std / mean = √(2 / N_f)

9 atoms, 2D     √(2/16)   = <b>0.354</b>  (35 %)
1000 atoms, 3D  √(2/2997) = <b>0.026</b>  (2.6 %)</div>`,
      see: 'For our 9 atoms at T₀ = 0.72, the instantaneous temperature should wander by about ±0.25 (35 %) around 0.72. That is not noise to be removed: it is the correct physics. In a large system it is smaller (2.6 % for 1000 atoms) but still there, and quantities like the heat capacity are computed from it. This number, <b>0.354</b>, is the test every thermostat below must pass.',
    },
    {
      title: 'Simple velocity rescaling: right temperature, no fluctuations',
      why: 'The most obvious idea: if T is 2.0 and we want 0.72, multiply every velocity by the same factor λ so that T becomes exactly 0.72. Since T is proportional to v², we need λ² = T₀/T.',
      work: `<div class="calc">λ = √(T₀ / T) = √(0.72 / 2.000) = <b>0.6000</b>

new Σ v² = 0.6² × 32.000 = 11.52
new T    = 11.52 / 16    = <b>0.720</b>  exactly

repeat every step  ->  T = 0.720 at every step
KE fluctuation (long run) = <b>0.000</b>
correct value             = 0.354</div>`,
      see: 'Merit: it reaches the target in one step, which is useful to remove a huge excess of heat at the very start. Limitation: applied every step, the temperature never fluctuates (0.000 instead of 0.354), so the system is <b>not</b> in the canonical ensemble. The sudden velocity jumps also distort the dynamics. Use it only for a first rough cooling, never for collecting data.',
    },
    {
      title: 'Berendsen: a gentle version of rescaling',
      why: 'Instead of jumping to T₀ in one step, Berendsen moves T only a small fraction dt/τ of the way each step. τ is a coupling time chosen by the user. This removes the violent jumps, but we must check whether the fluctuations come back.',
      work: `<div class="calc">rule      T_new = T + (dt/τ)(T₀ − T),  dt/τ = 0.05
factor    λ = √(1 + (dt/τ)(T₀/T − 1))

step 1  λ = √(1 + 0.05(0.72/2.0000 − 1)) = 0.98387
        T: 2.0000 -> <b>1.9360</b>
step 2  λ = √(1 + 0.05(0.72/1.9360 − 1)) = 0.98417
        T: 1.9360 -> <b>1.8752</b>
step 3  λ = 0.98448     T: 1.8752 -> <b>1.8174</b>
step 4  λ = 0.98479     T: 1.8174 -> <b>1.7626</b>

T − T₀:  1.2800, 1.2160, 1.1552, 1.0974, 1.0426
each step × (1 − dt/τ) = × 0.95
after τ/dt = 20 steps:  0.95²⁰ = 0.358 ≈ 1/e</div>
<p>(These hand steps show only the thermostat part of each step. In a real run the forces also move energy between kinetic and potential energy between thermostat steps; the long runs below include everything.) Now let it run, and watch two numbers in windows of 100 time units (20 000 steps): the kinetic-energy fluctuation, and the share of the kinetic energy in the <i>sliding modes</i>.</p>
<div class="calc">sliding modes: whole rows moving together in x,
whole columns moving together in y  ->  6 modes,
2 of them fixed by zero total momentum -> 4 free
fair share of KE = 4 / 16 = <b>0.25</b></div>
<div class="tbl"><table>${head(['time window', 'Berendsen τ = 0.1', 'Berendsen τ = 1.0', 'Bussi τ = 0.1 (step 5)'])}
${row(['correct', '0.25 / 0.354', '0.25 / 0.354', '0.25 / 0.354'])}
${row(['0–100', '0.253 / 0.150', '0.256 / 0.194', '0.251 / 0.342'])}
${row(['100–200', '0.257 / 0.148', '0.251 / 0.171', '0.261 / 0.351'])}
${row(['200–300', '<b>0.426</b> / 0.119', '0.252 / 0.170', '0.242 / 0.341'])}
${row(['300–400', '<b>0.921</b> / 0.022', '0.245 / 0.166', '0.244 / 0.354'])}
${row(['400–500', '<b>0.972</b> / 0.007', '0.260 / 0.172', '0.256 / 0.358'])}
${row(['900–1000', '<b>0.993</b> / 0.002', '0.247 / 0.175', '0.271 / 0.350'])}</table></div>
<p>Each cell is <i>sliding share / KE fluctuation</i>.</p>`,
      see: 'Merit: smooth, predictable cooling: the deviation from T₀ falls to 1/e in a time τ. That makes Berendsen excellent for <b>equilibration</b>. Limitation 1: even while things look healthy, the fluctuation is 0.15 (τ = 0.1) or 0.17 (τ = 1.0) instead of 0.354, so the ensemble is wrong, only less wrong for weak coupling. Limitation 2, much worse: with τ = 0.1, after about 250 time units the kinetic energy drains into the 4 sliding modes, until <b>99 %</b> of it is there (fair share 25 %) and the fluctuation has collapsed to 0.002. The atoms inside the rows have stopped vibrating, they are frozen, while whole rows glide past each other: the thermometer still reads 0.72. Why: rescaling all velocities by a common factor every step breaks the balance of energy flow between motions of different speed (Harvey, Tan and Cheatham, 1998). The slow collective motions exchange energy with the rest very slowly, so each small bias adds up there, step after step. In big simulations this is the famous <b>flying ice cube</b>: the whole system moves as a block while it freezes inside. Weak coupling (τ = 1.0) did not show it within 1000 time units, but it still has the wrong fluctuation. This is why Berendsen must not be used for production runs.',
    },
    {
      title: 'Nosé–Hoover: a friction that has its own equation',
      why: 'We want the right average and the right fluctuations, and we do not want to force T at every step. Nosé–Hoover adds one extra variable, a friction ξ, which acts on all velocities. ξ grows when the system is too hot and shrinks (even becomes negative, i.e. heats) when it is too cold. Its "mass" Q sets how quickly it responds.',
      work: `<div class="calc">friction  dξ/dt = (Σ v² − N_f T₀) / Q
velocity  dv/dt = F − ξ v
mass      Q = N_f T₀ τ² = 16 × 0.72 × 0.1² = <b>0.1152</b>
N_f T₀ = 11.52  (the Σ v² we are aiming for)

step 1  ξ = 0 + 0.005 × (32.000 − 11.52)/0.1152
          = <b>0.8889</b>
        v × e^(−ξ dt) = × 0.99557     T = <b>1.9823</b>
step 2  ξ = 0.8889 + 0.005 × (31.717 − 11.52)/0.1152
          = <b>1.7655</b>
        v × 0.99121                   T = <b>1.9476</b>
step 3  ξ = <b>2.6180</b>   v × 0.98700    T = <b>1.8973</b>
step 4  ξ = <b>3.4356</b>   v × 0.98297    T = <b>1.8332</b></div>
<p>The friction builds up step by step, cooling gets faster, and later ξ overshoots: T drops below 0.72, ξ turns negative and heats the system back. T oscillates around T₀ instead of being forced onto it. The long runs:</p>
<div class="tbl"><table>${head(['system', 'quantity', 'correct', 'NH', 'NH chain'])}
${row(['9 atoms', 'KE fluctuation', '0.354', '0.350', '0.353'])}
${row(['9 atoms', 'sliding-mode KE', '0.25', '0.249', '0.246'])}
${row(['1 oscillator', 'time with |v| > 1.5', '0.134', '<b>0.115</b>', '0.132'])}
${row(['1 oscillator', 'KE fluctuation', '1.414', '<b>0.904</b>', '1.428'])}</table></div>`,
      see: 'Merit: for the 9 atoms, Nosé–Hoover gets the fluctuation (0.350) and the energy sharing (0.249) right, and it is deterministic, so a run can be reproduced exactly. It is mathematically proven to give the canonical ensemble <i>if</i> the motion explores all states. Limitation: a very simple system, such as one harmonic oscillator (or a stiff, nearly harmonic solid), does not explore all states. Its fluctuation is 0.904 instead of 1.414 and fast velocities are too rare: in the widget you can see unvisited holes in the (x, v) plane. The fix is a <b>chain</b>: a second friction acting on the first one (NH chain, LAMMPS default <code>tchain 3</code>), which gives 1.428 and 0.132. Also, τ must not be too small or T oscillates strongly; 100 time steps is a common choice.',
    },
    {
      title: 'Bussi (CSVR): Berendsen plus the missing randomness',
      why: 'Berendsen fails because it pushes K towards its target smoothly, without any random kicks. Bussi, Donadio and Parrinello (2007) keep the same pull and add exactly the random part a heat bath would give, sized so that the fluctuations come out right.',
      work: `<div class="calc">target K₀ = N_f T₀ / 2 = 16 × 0.72 / 2 = 5.760
now    K  = Σ v² / 2 = <b>16.000</b>

dK = (K₀ − K) dt/τ  +  2 √(K K₀ dt / (N_f τ)) × R
   =  −0.512        +  1.073 × R

R is a fresh Gaussian random number each step:
R = +1  ->  dK = +0.561  (this step heats!)
R = −1  ->  dK = −1.585
R =  0  ->  dK = −0.512  (= Berendsen)

then all velocities × √((K + dK) / K)</div>
<div class="tbl"><table>${head(['thermostat', 'KE fluctuation', 'sliding KE', 'diffusion D'])}
${row(['correct', '0.354', '0.25', '—'])}
${row(['Berendsen τ 0.1 (t > 400)', '0.002', '0.99', '—'])}
${row(['Nosé–Hoover', '0.350', '0.249', '0.0575'])}
${row(['Bussi τ 0.1', '<b>0.350</b>', '<b>0.250</b>', '<b>0.0597</b>'])}</table></div>`,
      see: 'Merit: with R = 0 it is exactly Berendsen (smooth equilibration), and the random term restores the correct fluctuation (0.350) and energy sharing (0.250). It is robust even for the oscillator-like systems where plain NH fails, and since all velocities are scaled by one common factor, the diffusion (0.0597) is close to the Nosé–Hoover value (0.0575). Limitation: the run depends on the random numbers, so it is reproducible only with the same seed; total momentum is kept, but there is no conserved "energy" in the usual sense (codes report an effective one). In GROMACS it is <code>v-rescale</code>, in LAMMPS <code>fix temp/csvr</code>.',
    },
    {
      title: 'Langevin: friction and random kicks on every atom',
      why: 'A real heat bath is made of solvent molecules that slow an atom down (friction γ) and hit it at random. Langevin dynamics adds exactly these two forces to every atom separately. Because the kicks are on each atom, it does not matter whether the motion explores all states.',
      work: `<div class="calc">friction γ = 1:   c = e^(−γ dt) = e^(−0.005) = <b>0.99501</b>
noise size        √((1 − c²) T₀) = <b>0.08464</b>

one velocity component  v = 1.20,  random R = −0.70
v_new = 0.99501 × 1.20 + 0.08464 × (−0.70)
      = 1.19401 − 0.05925 = <b>1.13477</b>

total momentum is NOT conserved (kicks are independent)
-> N_f = 2N = 18,  correct fluctuation √(2/18) = 0.333
-> sliding modes 6 of 18 = 0.333</div>
<div class="tbl"><table>${head(['thermostat', 'KE fluctuation', 'sliding KE', 'diffusion D'])}
${row(['correct (N_f = 18)', '0.333', '0.333', '—'])}
${row(['Nosé–Hoover (reference)', '0.350', '0.249', '0.0575'])}
${row(['Langevin γ = 1', '0.330', '0.330', '<b>0.0434</b>'])}
${row(['Langevin γ = 10', '0.331', '0.333', '<b>0.0148</b>'])}</table></div>`,
      see: 'Merit: the fluctuations and the energy sharing are right (0.330 and 0.330 vs 0.333), it works for any system, including the single oscillator (0.135 vs 0.134), and it is stable even far from equilibrium. Limitation: the friction slows the atoms down. Diffusion drops from 0.058 to 0.043 at γ = 1 and to 0.015 at γ = 10, four times too slow. Momentum is not conserved, so collective flow is damped too. Use Langevin for sampling structures and energies, or with a small γ; do not take diffusion coefficients or viscosities from a strongly coupled Langevin run.',
    },
    {
      title: 'Pressure and barostats: the same story for the volume',
      why: 'To simulate at constant pressure (NPT), we first need the pressure of a configuration, and then a rule that changes the box size. The pressure has a kinetic part (atoms hitting the walls) and a virial part (forces between atoms, Lesson 3.6).',
      work: `<div class="calc">2D: P = (N T + ½ Σ_pairs r F) / A,   A = L² = 12.96

lattice pairs inside r_c: 18 at r = 1.2, 18 at 1.697
F(1.2)   = −2.2117   (negative = attraction)
F(1.697) = −0.5425
Σ r F = 18 × 1.2 × (−2.2117) + 18 × 1.697 × (−0.5425)
      = <b>−64.343</b>

P = (9 × 0.72 + 0.5 × (−64.343)) / 12.96
  = (6.48 − 32.17) / 12.96 = <b>−1.982</b>

Berendsen barostat towards P₀ = 0 (β = 1, τ_P = 1):
μ = 1 − (dt/τ_P) β (P₀ − P) / 2 = <b>0.99504</b>
-> box side × 0.99504 each step: it shrinks</div>`,
      see: 'A negative pressure means the atoms pull on each other more than they push: at this density and temperature the box is too large, so the barostat shrinks it. Merit of Berendsen scaling: smooth and fast, ideal to bring the density to the right value. Limitation: exactly like the Berendsen thermostat, it suppresses the volume fluctuations, so the ensemble is not truly NPT and the compressibility comes out wrong. For production use a barostat with the correct fluctuations: C-rescale (Bernetti–Bussi, the barostat version of step 5), Parrinello–Rahman, or MTK (LAMMPS <code>fix npt</code>).',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)', 'use it for'])}
${row(['rescaling', 'T₀ in one step', 'fluctuation 0.000 vs 0.354', 'removing a huge initial excess'])}
${row(['Berendsen', 'smooth, 1/e in τ', 'fluctuation 0.15 vs 0.354; later 99 % KE in sliding modes (flying ice cube)', 'equilibration only'])}
${row(['Nosé–Hoover', 'fluctuation 0.350, deterministic', 'single oscillator 0.904 vs 1.414', 'production, as a chain'])}
${row(['NH chain', '0.353; oscillator 1.428', 'needs τ ≈ 100 dt', 'production (LAMMPS fix nvt)'])}
${row(['Bussi (CSVR)', '0.350; D = 0.0597', 'random, needs a seed', 'production (GROMACS v-rescale)'])}
${row(['Langevin', 'works for any system', 'D 0.058 -> 0.015 at γ = 10', 'sampling; small γ if dynamics matter'])}
${row(['Berendsen barostat', 'P −1.98 -> 0 smoothly', 'volume fluctuations suppressed', 'equilibrating the density'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the fluctuation estimate, the first steps of each thermostat from the hot start, the long-run tables (100 000 steps each), the flying-ice-cube windows, the single-oscillator test and the virial pressure. Run <code>pip install numpy</code>, then <code>python lab_thermostats.py</code> (about a minute; the long runs take most of it).',
    src: `# Lab 3: thermostats by hand. 9 LJ atoms in a 2D periodic box (the system of Lab 2), reduced units, target T0 = 0.72.
import numpy as np, time
L, N, RC, T0 = 3.6, 9, 1.8, 0.72
NF = 2 * N - 2                                            # 2D, total momentum fixed
LAT = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(3) for i in range(3)])
GAU = np.array([[-0.33, -0.88], [-0.66, -0.67], [0.38, -0.11], [1.48, -1.83], [-0.00, -0.89], [0.78, -2.12], [-0.34, 0.21], [-1.48, 0.99], [0.18, 1.01]])

# ---------- Step 1: temperature and its natural fluctuation ----------
v0 = GAU - GAU.mean(0); v0 *= np.sqrt(2.0 * NF / (v0**2).sum())       # start HOT: T = 2.0
print(f"Step 1: start hot. sum v^2 = {(v0**2).sum():.3f}, N_f = {NF}, T = {(v0**2).sum() / NF:.3f}")
print(f"   canonical (correct) relative fluctuation of the kinetic energy = sqrt(2/N_f) = sqrt(2/{NF}) = {np.sqrt(2 / NF):.3f}")
print(f"   for 1000 atoms in 3D: sqrt(2/2997) = {np.sqrt(2 / 2997):.3f}")
g = np.random.default_rng(0).normal(size=1_000_000)
print(f"   check for a Gaussian velocity: <v^4>/<v^2>^2 = {np.mean(g**4) / np.mean(g**2)**2:.3f} (exactly 3), so var(v^2) = 2 <v^2>^2")

# ---------- Steps 2-6: hand updates of each thermostat (the thermostat part of one step only) ----------
dt, tau = 0.005, 0.1
T = 2.0
print(f"\\nStep 2: simple rescaling: lambda = sqrt(T0/T) = sqrt({T0}/{T}) = {np.sqrt(T0 / T):.4f} -> T = {T * T0 / T:.3f} exactly, every step")
print(f"\\nStep 3: Berendsen, dt/tau = {dt / tau:.2f}: T_new = T + (dt/tau)(T0 - T)")
for n in range(4):
    lam = np.sqrt(1 + dt / tau * (T0 / T - 1)); Tn = lam**2 * T
    print(f"   step {n + 1}: lambda = sqrt(1 + 0.05 x ({T0}/{T:.4f} - 1)) = {lam:.5f};  T: {T:.4f} -> {Tn:.4f};  T - T0 = {Tn - T0:.4f}")
    T = Tn
print(f"   T - T0 shrinks by (1 - dt/tau) = 0.95 per step; after tau/dt = 20 steps: 0.95^20 = {0.95**20:.3f} (about 1/e)")
Q = NF * T0 * tau**2; T = 2.0; xi = 0.0; s2 = T * NF
print(f"\\nStep 4: Nose-Hoover, Q = N_f T0 tau^2 = {NF} x {T0} x {tau}^2 = {Q:.4f}")
for n in range(4):
    xi += dt * (s2 - NF * T0) / Q; f = np.exp(-xi * dt); s2 *= f * f
    print(f"   step {n + 1}: xi += {dt} x ({s2 / f / f:.3f} - {NF * T0:.2f}) / {Q:.4f} -> xi = {xi:.4f};  velocities x exp(-xi dt) = {f:.5f};  T = {s2 / NF:.4f}")
c = np.exp(-dt / tau); K = 0.5 * 2.0 * NF; K0 = 0.5 * T0 * NF
print(f"\\nStep 5: Bussi (stochastic Berendsen): target K0 = N_f T0 / 2 = {K0:.3f}, now K = {K:.3f}")
for R in [+1.0, -1.0, 0.0]:
    dK = (K0 - K) * dt / tau + 2 * np.sqrt(K * K0 / NF * dt / tau) * R
    print(f"   random number R = {R:+.1f}: dK = (K0 - K) dt/tau + 2 sqrt(K K0 / N_f x dt/tau) R = {(K0 - K) * dt / tau:.3f} + {2 * np.sqrt(K * K0 / NF * dt / tau) * R:+.3f} = {dK:+.3f}")
gamma = 1.0; c = np.exp(-gamma * dt)
print(f"\\nStep 6: Langevin, gamma = {gamma}: c = exp(-gamma dt) = {c:.5f};  noise size sqrt((1 - c^2) T0) = {np.sqrt((1 - c * c) * T0):.5f}")
print(f"   one velocity component v = 1.20 with random number R = -0.70: v_new = {c:.5f} x 1.20 + {np.sqrt((1 - c * c) * T0):.5f} x (-0.70) = {c * 1.2 + np.sqrt((1 - c * c) * T0) * -0.7:.5f}")

# ---------- long runs of the 9-atom system with every thermostat ----------
def forces(x):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < RC * RC
    ir6 = np.where(m, r2**-3, 0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0)
    return -(f[:, :, None] * d).sum(1)

def run(method, nsteps=100000, tau=0.1, gamma=1.0, seed=1):
    r = np.random.default_rng(seed)
    x = LAT.copy(); v = v0.copy(); f = forces(x); xi = xi2 = 0.0; Q = NF * T0 * tau**2
    Ks, slow, traj = [], [], []; u = x.copy()
    for s in range(nsteps):
        if method in ('nh', 'nhc'):
            if method == 'nhc': xi2 += 0.25 * dt * (Q * xi**2 - T0) / Q; xi *= np.exp(-0.25 * dt * xi2)
            xi += 0.5 * dt * ((v**2).sum() - NF * T0) / Q; v *= np.exp(-0.5 * dt * xi)
        v += 0.5 * dt * f; x += dt * v; u += dt * v; x %= L; f = forces(x); v += 0.5 * dt * f
        if method in ('nh', 'nhc'):
            v *= np.exp(-0.5 * dt * xi); xi += 0.5 * dt * ((v**2).sum() - NF * T0) / Q
            if method == 'nhc': xi *= np.exp(-0.25 * dt * xi2); xi2 += 0.25 * dt * (Q * xi**2 - T0) / Q
        Tn = (v**2).sum() / NF
        if method == 'rescale': v *= np.sqrt(T0 / Tn)
        if method == 'berendsen': v *= np.sqrt(1 + dt / tau * (T0 / Tn - 1))
        if method == 'langevin':
            cc = np.exp(-gamma * dt); v = cc * v + np.sqrt((1 - cc * cc) * T0) * r.normal(size=v.shape)
        if method == 'bussi':
            cc = np.exp(-dt / tau); K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0
            R1 = r.normal(); S = (r.normal(size=NF - 1)**2).sum()
            a2 = cc + (1 - cc) * Kt / (NF * K) * (R1**2 + S) + 2 * R1 * np.sqrt(cc * (1 - cc) * Kt / (NF * K))
            v *= np.sqrt(a2) * np.sign(R1 + np.sqrt(cc * NF * K / ((1 - cc) * Kt)))
        if s >= 10000:
            Ks.append(0.5 * (v**2).sum())
            rows = np.array([v[3 * k:3 * k + 3, 0].mean() for k in range(3)]); cols = np.array([v[k::3, 1].mean() for k in range(3)])
            slow.append(1.5 * ((rows**2).sum() + (cols**2).sum()))            # KE of whole rows / columns sliding
            if s % 20 == 0: traj.append(u - u.mean(0))
    K = np.array(Ks); traj = np.array(traj); lag = 50
    D = ((traj[lag:] - traj[:-lag])**2).sum(-1).mean() / (4 * lag * 20 * dt)
    nf = 2 * N if method == 'langevin' else NF                        # Langevin does not conserve momentum
    return 2 * K.mean() / nf, K.std() / K.mean(), np.mean(slow) / K.mean(), D

print(f"\\nLong runs (500 time units each). Correct: T = {T0}, KE fluctuation {np.sqrt(2 / NF):.3f} (Langevin: {np.sqrt(2 / 18):.3f}), KE in row/column sliding 0.25 (Langevin 0.33)")
print("   thermostat              <T>     KE fluct.  sliding KE   D (diffusion)")
t0 = time.time()
for name, m, kw in [("simple rescaling", 'rescale', {}), ("Berendsen tau 0.1", 'berendsen', dict(tau=0.1)), ("Berendsen tau 1.0", 'berendsen', dict(tau=1.0)),
                    ("Nose-Hoover tau 0.1", 'nh', dict(tau=0.1)), ("NH chain tau 0.1", 'nhc', dict(tau=0.1)), ("Bussi tau 0.1", 'bussi', dict(tau=0.1)),
                    ("Langevin gamma 1", 'langevin', dict(gamma=1.0)), ("Langevin gamma 10", 'langevin', dict(gamma=10.0))]:
    Tm, rs, fs, D = run(m, **kw)
    print(f"   {name:22s} {Tm:6.3f}   {rs:7.3f}    {fs:7.3f}      {D:7.4f}")
print(f"   ({time.time() - t0:.0f} s)  The Berendsen tau 0.1 line averages a healthy start with a frozen end: see the windows below.")

# ---------- the flying ice cube: Berendsen in a longer run, watched in windows of 100 time units ----------
def watch(method, tau, nsteps=200000, seed=1):
    r = np.random.default_rng(seed); x = LAT.copy(); v = v0.copy(); f = forces(x); sl, Ks, out = [], [], []
    for s in range(nsteps):
        v += 0.5 * dt * f; x += dt * v; x %= L; f = forces(x); v += 0.5 * dt * f
        if method == 'berendsen': v *= np.sqrt(1 + dt / tau * (T0 / ((v**2).sum() / NF) - 1))
        else:
            cc = np.exp(-dt / tau); K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = r.normal(); S = (r.normal(size=NF - 1)**2).sum()
            a2 = cc + (1 - cc) * Kt / (NF * K) * (R1**2 + S) + 2 * R1 * np.sqrt(cc * (1 - cc) * Kt / (NF * K))
            v *= np.sqrt(a2) * np.sign(R1 + np.sqrt(cc * NF * K / ((1 - cc) * Kt)))
        rows = np.array([v[3 * k:3 * k + 3, 0].mean() for k in range(3)]); cols = np.array([v[k::3, 1].mean() for k in range(3)])
        sl.append(1.5 * ((rows**2).sum() + (cols**2).sum())); Ks.append(0.5 * (v**2).sum())
        if (s + 1) % 20000 == 0: out.append((np.mean(sl) / np.mean(Ks), np.std(Ks) / np.mean(Ks))); sl, Ks = [], []
    return out

print("\\nFlying ice cube: sliding-mode share of KE (correct 0.25) and KE fluctuation (correct 0.354) in each window of 100 time units")
res = {name: watch(m, tau) for name, m, tau in [("Berendsen tau 0.1", 'berendsen', 0.1), ("Berendsen tau 1.0", 'berendsen', 1.0), ("Bussi tau 0.1", 'bussi', 0.1)]}
print("   window (time)   " + "   ".join(f"{n:>20s}" for n in res))
for w in range(10):
    print(f"   {100 * w:4d}-{100 * w + 100:<4d}       " + "   ".join(f"{res[n][w][0]:9.3f} / {res[n][w][1]:5.3f}" .rjust(20) for n in res))

# ---------- one harmonic oscillator: is plain Nose-Hoover ergodic? ----------
def oscillator(method, nsteps=400000, dt=0.05, seed=2):
    r = np.random.default_rng(seed); x, v, xi, xi2, Q = 1.0, 0.0, 0.0, 0.0, 1.0; vs = []
    for s in range(nsteps):
        if method in ('nh', 'nhc'):
            if method == 'nhc': xi2 += 0.25 * dt * (Q * xi**2 - 1) / Q; xi *= np.exp(-0.25 * dt * xi2)
            xi += 0.5 * dt * (v * v - 1) / Q; v *= np.exp(-0.5 * dt * xi)
        v += 0.5 * dt * (-x); x += dt * v; v += 0.5 * dt * (-x)
        if method in ('nh', 'nhc'):
            v *= np.exp(-0.5 * dt * xi); xi += 0.5 * dt * (v * v - 1) / Q
            if method == 'nhc': xi *= np.exp(-0.25 * dt * xi2); xi2 += 0.25 * dt * (Q * xi**2 - 1) / Q
        if method == 'langevin': cc = np.exp(-dt); v = cc * v + np.sqrt(1 - cc * cc) * r.normal()
        vs.append(v)
    vs = np.array(vs[nsteps // 10:]); K = 0.5 * vs**2
    return np.mean(np.abs(vs) > 1.5), K.std() / K.mean()
print("\\nOne harmonic oscillator at T = 1 (correct: fraction of time with |v| > 1.5 = 0.134, KE fluctuation sqrt(2) = 1.414)")
for name, m in [("Nose-Hoover", 'nh'), ("NH chain", 'nhc'), ("Langevin", 'langevin')]:
    fr, rs = oscillator(m); print(f"   {name:12s} |v| > 1.5: {fr:.3f}   KE fluctuation: {rs:.3f}")

# ---------- Step 7: pressure of the lattice (virial) and one barostat step ----------
U = lambda r: 4 * (r**-12 - r**-6); F = lambda r: 24 * (2 * r**-13 - r**-7)
A = L * L; W = 18 * 1.2 * F(1.2) + 18 * np.sqrt(2) * 1.2 * F(np.sqrt(2) * 1.2)
P = (N * T0 + 0.5 * W) / A
print(f"\\nStep 7: 2D virial pressure P = (N T + 1/2 sum r F) / A")
print(f"   F(1.2) = {F(1.2):.4f}, F(1.697) = {F(np.sqrt(2) * 1.2):.4f};  sum r F = 18 x 1.2 x {F(1.2):.4f} + 18 x 1.697 x {F(np.sqrt(2) * 1.2):.4f} = {W:.3f}")
print(f"   P = ({N} x {T0} + 0.5 x {W:.3f}) / {A:.2f} = {P:.3f}   (negative: the atoms pull the box inwards)")
beta, tauP = 1.0, 1.0
mu = 1 - dt / tauP * beta * (0.0 - P) / 2
print(f"   Berendsen barostat towards P0 = 0 (compressibility 1, tau_P = 1): box scale mu = 1 - (dt/tau_P) beta (P0 - P)/2 = {mu:.5f} per step (box shrinks)")
`,
  },

  exercises: [
    { q: 'A 3D system of 125 atoms with zero total momentum is held at T = 1.0. By how much (relative, and in units of T) should its instantaneous temperature fluctuate?', a: 'N_f = 3 × 125 − 3 = 372. Relative fluctuation √(2/372) = 0.073, so T should wander by about ±0.073 around 1.0. If your thermostat gives much less, it is suppressing fluctuations.' },
    { q: 'Berendsen with dt = 0.002 and τ = 0.2 starts at T = 600 K with target 300 K. What is T after one step and after 100 steps (ignore the forces)?', a: 'dt/τ = 0.01. After one step: 600 + 0.01 × (300 − 600) = 597 K. The deviation shrinks by 0.99 per step: after 100 steps 300 × 0.99¹⁰⁰ = 300 × 0.366 = 110 K, so T ≈ 410 K. 100 steps is exactly τ, so the deviation is down to about 1/e.' },
    { q: 'Langevin with dt = 0.005 and γ = 0.1 at T₀ = 0.72: compute c and the noise size, and update v = −0.50 with the random number R = +1.2.', a: 'c = e^(−0.0005) = 0.99950; noise √((1 − 0.99950²) × 0.72) = √(0.0009995 × 0.72) = 0.02683. v_new = 0.99950 × (−0.50) + 0.02683 × 1.2 = −0.49975 + 0.03219 = −0.46756. With this small γ the dynamics are barely disturbed.' },
  ],
};
