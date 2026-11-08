import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-basics',
  lab: true,
  summary: 'Two argon atoms, a calculator and three time steps. We compute the Lennard-Jones energy and force by hand, move the atoms with two different recipes (Euler and velocity Verlet), and let the numbers show why every MD code uses Verlet and why the time step must be small.',
  minutes: 25,
  sim: 'lab1',
  learn: [
    'Compute the Lennard-Jones energy and force with only powers and fractions, and check the force',
    'Do three Euler steps and three velocity-Verlet steps by hand, and compare the energies',
    'See from the numbers why Euler fails, why Verlet works, and how large the time step may be',
  ],
  simHint: 'Press "Step" to repeat the hand calculation one step at a time; the table shows every number. Both recipes run side by side so you can compare their energies. Then press "Run" and try larger time steps.',

  setup: `
<p>Two argon atoms sit on a line. We use <b>reduced units</b> (Lesson 2.1): lengths in σ (3.4 Å for argon), energies in ε (the depth of the well), masses in units of the argon mass, so <b>m = 1</b>. One time unit is about 2.2 ps.</p>
<div class="calc">atom 1:  x1 = 0.00    v1 = 0      (at rest)
atom 2:  x2 = 1.05    v2 = 0      (at rest)
distance r = x2 − x1 = 1.05
(a little closer than the atoms like)
time step dt = 0.01</div>
<p>The atoms are slightly squeezed, so they will push apart and then vibrate like two balls on a spring. Nothing is added or taken away, so the <b>total energy must stay constant</b>. That is our test of whether a recipe can be trusted.</p>`,

  steps: [
    {
      title: 'The energy between two atoms',
      why: 'MD needs a rule for how atoms push and pull. Real atoms repel strongly when squeezed and attract weakly when a little apart. The Lennard-Jones formula captures both with just two numbers (σ and ε).',
      work: `<p>In reduced units the formula is ${tex(String.raw`U(r) = 4\left(r^{-12} - r^{-6}\right)`)}. Work out r<sup>−6</sup> first, then square it to get r<sup>−12</sup>:</p>
<div class="calc">r = 1.05:  r^6 = 1.05×1.05×1.05×1.05×1.05×1.05 = 1.3401
          r^-6  = 1/1.3401 = 0.7462
          r^-12 = 0.7462 × 0.7462 = 0.5568
          U = 4 × (0.5568 − 0.7462) = 4 × (−0.1894) = <b>−0.7575</b></div>
<div class="tbl"><table>${head(['r', 'r⁻⁶', 'r⁻¹²', 'U(r)', 'meaning'])}
${row(['0.95', '1.3604', '1.8506', '+1.9610', 'squeezed: energy up'])}
${row(['1.00', '1.0000', '1.0000', '0.0000', 'r = σ'])}
${row(['1.05', '0.7462', '0.5568', '−0.7575', 'our start'])}
${row(['1.1225', '0.5000', '0.2500', '−1.0000', 'the minimum, r = 2<sup>1/6</sup>'])}
${row(['1.50', '0.0878', '0.0077', '−0.3203', 'weak attraction'])}
${row(['2.50', '0.0041', '0.0000', '−0.0163', 'almost nothing'])}</table></div>`,
      see: 'The energy is lowest at r = 1.1225 (exactly −1 ε), rises steeply below r = 1 and fades to almost zero by r = 2.5. That last fact is why we can ignore far-away pairs (cutoffs, Lesson 2.5). Limitation: one formula for simple atoms only; molecules and ions need more terms (Chapters 4 and 5).',
    },
    {
      title: 'From energy to force',
      why: 'Newton\'s law needs forces, not energies. A force always pushes "downhill": it is minus the slope of the energy. Getting the force wrong is the most common bug in MD codes, so we also check it.',
      work: `<p>The slope of r<sup>n</sup> is n·r<sup>n−1</sup>. Applying this to both terms:</p>
<div class="calc">dU/dr = 4 × (−12 r^-13 + 6 r^-7)
F(r)  = −dU/dr = 24 × (2 r^-13 − r^-7)

r = 1.05:  r^-13 = r^-12 / r = 0.5568 / 1.05 = 0.5303
           r^-7  = r^-6 / r  = 0.7462 / 1.05 = 0.7107
           F = 24 × (2×0.5303 − 0.7107)
             = 24 × 0.3499 = <b>8.399</b>
           positive: pushes the atoms apart</div>
<p><b>Check</b> with the slope measured directly between two nearby points (a "finite difference"):</p>
<div class="calc">U(1.051) − U(1.049) = −0.016799      over a distance 0.002
F ≈ −(−0.016799)/0.002 = <b>8.3996</b>
✓ agrees with 8.399</div>
<div class="tbl"><table>${head(['r', 'F (formula)', 'F (finite difference)', 'direction'])}
${row(['1.00', '24.000', '24.001', 'strong push apart'])}
${row(['1.05', '8.399', '8.400', 'push apart'])}
${row(['1.1225', '0.000', '0.000', 'no force at the minimum'])}
${row(['1.50', '−1.158', '−1.158', 'gentle pull together'])}</table></div>`,
      see: 'The force is large and repulsive when the atoms are squeezed, zero at the minimum, and weakly attractive beyond it. The finite-difference check costs one line and catches sign errors and missing factors.',
    },
    {
      title: 'Newton\'s law for the pair',
      why: 'To move atoms we need their accelerations, a = F/m. The two atoms feel equal and opposite forces (Newton\'s third law), so their total momentum stays zero.',
      work: `<div class="calc">force on atom 2 = +8.399     a2 = +8.399 / 1 = +8.399
force on atom 1 = −8.399     a1 = −8.399
total momentum  = m v1 + m v2 = 0 at the start,
                  and it must stay 0</div>`,
      see: 'In every table below, v1 = −v2 exactly. If your own code ever gives a drifting total momentum, the forces are not equal and opposite (a bug).',
    },
    {
      title: 'Recipe 1: Euler, three steps by hand',
      why: 'The simplest idea: move each atom with its current velocity, and change the velocity with its current acceleration, both for the whole step. Let us see what happens to the energy.',
      work: `<div class="calc">Euler:   x_new = x + v·dt
         v_new = v + a·dt

step 1:  x2 = 1.05 + 0×0.01 = 1.05000
         (uses the OLD velocity 0: the atom does not move!)
         v2 = 0 + 8.3991×0.01    = 0.08399
         energy: KE = ½v1² + ½v2²
                    = 0.08399² = 0.007054
                 U(1.05) = −0.757512
                 E = −0.750457      <b>up by 0.0071 in one step</b>

step 2:  x2 = 1.05 + 0.08399×0.01 = 1.05084
         x1 = −0.00084,  r = 1.05168
         v2 = 0.08399 + 8.3991×0.01 = 0.16798
         (still the force at r = 1.05)
         E = −0.743121              up again</div>
<div class="tbl"><table>${head(['n', 'x1', 'x2', 'v2', 'r', 'F', 'E'])}
${row(['0', '0.00000', '1.05000', '0.00000', '1.05000', '8.3991', '−0.757512'])}
${row(['1', '0.00000', '1.05000', '0.08399', '1.05000', '8.3991', '−0.750457'])}
${row(['2', '−0.00084', '1.05084', '0.16798', '1.05168', '8.0653', '−0.743121'])}
${row(['3', '−0.00252', '1.05252', '0.24863', '1.05504', '7.4251', '−0.735531'])}</table></div>`,
      see: 'The energy rises by about 0.007 every step. Euler uses the force from the start of the step for the whole step, but the force is getting smaller as the atoms separate, so every kick is a little too big. The error always has the same sign, so it piles up.',
    },
    {
      title: 'Recipe 2: velocity Verlet, three steps by hand',
      why: 'Verlet uses the average of the old and new forces for the velocity, like the trapezoid rule for areas. It needs the same single force calculation per step as Euler, so it costs the same.',
      work: `<div class="calc">velocity Verlet:
  1. v_half = v + ½a·dt          (half kick, old force)
  2. x_new  = x + v_half·dt      (move)
  3. compute the NEW force at x_new
  4. v_new  = v_half + ½a_new·dt (half kick, new force)

step 1:  v2_half = 0 + ½×8.3991×0.01   = 0.041996
         x2 = 1.05 + 0.041996×0.01 = 1.05042
         x1 = −0.00042,  r = 1.05084
         new force F(1.05084)            = 8.2310
         v2 = 0.041996 + ½×8.2310×0.01   = 0.083151
         energy: KE = 0.083151² = 0.006914
                 U(1.05084) = −0.764496
                 E = −0.757582
                 <b>changed by only −0.00007</b></div>
<div class="tbl"><table>${head(['n', 'x1', 'x2', 'v2', 'r', 'F', 'E'])}
${row(['0', '0.00000', '1.05000', '0.00000', '1.05000', '8.3991', '−0.757512'])}
${row(['1', '−0.00042', '1.05042', '0.08315', '1.05084', '8.2310', '−0.757582'])}
${row(['2', '−0.00166', '1.05166', '0.16304', '1.05333', '7.7472', '−0.757771'])}
${row(['3', '−0.00368', '1.05368', '0.23679', '1.05736', '7.0031', '−0.758028'])}</table></div>`,
      see: 'After three steps the Verlet energy error is 0.0005, against 0.022 for Euler: about 40 times smaller for the same work. The small errors also go up and down instead of piling up (next step).',
    },
    {
      title: 'Many steps: drift or wobble?',
      why: 'Three steps do not prove anything about a whole simulation. We let the computer repeat our hand recipe for about three vibrations (2 time units) and compare different time steps.',
      work: `<div class="tbl"><table>${head(['dt', 'Euler: E(end) − E₀', 'Verlet: largest |E − E₀|'])}
${row(['0.020', '+1.765', '0.00497'])}
${row(['0.010', '+0.849', '0.00123'])}
${row(['0.005', '+0.331', '0.00031'])}</table></div>
<div class="calc">Euler:  halving dt roughly halves the error
        -> error grows like dt   ("first order")
Verlet: 0.00497 / 0.00123 = 4.0
        0.00123 / 0.00031 = 4.0
        -> error grows like dt²  ("second order")

Euler, dt = 0.01:  E rises from −0.76 to above 0,
                   the molecule BREAKS at t = 3.35
Euler, dt = 0.005: breaks at t = 6.17
Verlet:            never breaks (checked to t = 50)</div>`,
      see: 'Euler adds energy steadily until the atoms fly apart: an unphysical result from a numerical error. Verlet keeps the energy within 0.2 % for dt = 0.01, forever. Its limitation: the energy is not exact (it wobbles), and the wobble grows as dt². This is why every MD code uses Verlet or its twin, leap-frog.',
    },
    {
      title: 'How large may the time step be?',
      why: 'A larger dt means fewer steps, so a cheaper simulation. But the time step must resolve the fastest motion in the system. We estimate that motion by hand and then test it.',
      work: `<p>Near the minimum the LJ well looks like a spring. Its stiffness is the curvature of U there, and two atoms of mass 1 vibrate like one particle with half the mass (the reduced mass μ = 1/2):</p>
<div class="calc">stiffness  k = U''(r_min) = 36 × 2^(2/3) = 57.15
frequency  ω = √(k/μ) = √(57.15/0.5) = 10.69
period     = 2π/ω = 6.2832/10.69 = <b>0.588</b>   time units</div>
<div class="tbl"><table>${head(['dt', 'steps per vibration', 'Verlet: largest |E − E₀| (20 time units)'])}
${row(['0.01', '58.8', '0.0012'])}
${row(['0.03', '19.6', '0.011'])}
${row(['0.06', '9.8', '0.047'])}
${row(['0.12', '4.9', '1.81 (molecule nearly breaks)'])}
${row(['0.20', '2.9', '29 (explodes)'])}</table></div>
<div class="calc">rule of thumb: at least ~20 steps per period
               of the fastest vibration
Verlet is unstable when ω·dt > 2,
               i.e. dt > 2/10.69 = 0.19</div>`,
      see: 'For argon, dt = 0.01 τ ≈ 22 fs is fine. A water O–H bond vibrates with a period of about 9 fs, so water needs dt ≈ 0.5–1 fs, or constraints that freeze the bond (Lesson 4.4). The time step is set by the fastest motion, not by what you want to study.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['', 'Euler', 'Velocity Verlet'])}
${row(['force calculations per step', '1', '1'])}
${row(['energy error after 3 steps (dt = 0.01)', '+0.022', '−0.0005'])}
${row(['halving dt reduces the error by', '≈ 2×', '4×'])}
${row(['long run', 'energy drifts up, molecule breaks', 'energy wobbles, stays bounded'])}
${row(['use it for MD?', 'no', 'yes (with dt ≲ period/20)'])}</table></div>
<p>The same steps with forces from many atoms, in three dimensions, are a complete MD program. Everything later in the course adds to this loop; nothing replaces it.</p>`,

  code: {
    note: 'This script prints every number on this page: the energy and force tables, the three Euler and Verlet steps, the long-run errors and the time-step test. Run <code>pip install numpy</code>, then <code>python lab_two_atoms.py</code> (a second or two).',
    src: `# Lab 1: two argon atoms by hand. Reduced units: sigma = epsilon = mass = 1. Atoms move along a line (1D).
import numpy as np

def U(r):  return 4 * (r**-12 - r**-6)              # Lennard-Jones energy
def F(r):  return 24 * (2 * r**-13 - r**-7)          # force on atom 2 (pushing it away from atom 1 when > 0)

# ---------- Step 1: the potential at a few distances ----------
print("Step 1: U(r) = 4 (r^-12 - r^-6)")
for r in [0.95, 1.0, 1.05, 2**(1 / 6), 1.5, 2.5]:
    print(f"   r = {r:6.4f}:  r^-6 = {r**-6:8.4f}  r^-12 = {r**-12:8.4f}  U = {U(r):8.4f}")

# ---------- Step 2: the force, and a check by finite differences ----------
print("\\nStep 2: F(r) = -dU/dr = 24 (2 r^-13 - r^-7)")
for r in [1.0, 1.05, 2**(1 / 6), 1.5]:
    h = 1e-3; fd = -(U(r + h) - U(r - h)) / (2 * h)
    print(f"   r = {r:6.4f}:  F = {F(r):9.4f}   finite difference = {fd:9.4f}")

# ---------- Steps 3-5: three steps of Euler and of velocity Verlet ----------
def euler(x1, x2, v1, v2, dt):
    f = F(x2 - x1)                                    # force on atom 2; atom 1 feels -f
    x1n, x2n = x1 + v1 * dt, x2 + v2 * dt             # positions use the OLD velocities
    v1n, v2n = v1 - f * dt, v2 + f * dt               # velocities use the OLD force (a = F/m, m = 1)
    return x1n, x2n, v1n, v2n

def verlet(x1, x2, v1, v2, dt):
    f = F(x2 - x1)
    v1h, v2h = v1 - 0.5 * f * dt, v2 + 0.5 * f * dt   # half kick
    x1n, x2n = x1 + v1h * dt, x2 + v2h * dt           # drift
    fn = F(x2n - x1n)                                 # NEW force
    return x1n, x2n, v1h - 0.5 * fn * dt, v2h + 0.5 * fn * dt   # second half kick

def energy(x1, x2, v1, v2): return 0.5 * v1**2 + 0.5 * v2**2 + U(x2 - x1)

start = (0.0, 1.05, 0.0, 0.0)                         # atoms at rest, 1.05 sigma apart (slightly squeezed)
for name, step in [("Euler", euler), ("Velocity Verlet", verlet)]:
    print(f"\\n{name}, dt = 0.01")
    print("   n     x1        x2        v1        v2        r        F         E")
    s = start
    for n in range(4):
        x1, x2, v1, v2 = s
        print(f"   {n}  {x1:8.5f}  {x2:8.5f}  {v1:8.5f}  {v2:8.5f}  {x2 - x1:7.5f}  {F(x2 - x1):8.4f}  {energy(*s):9.6f}")
        s = step(*s, 0.01)

# ---------- Step 6: long runs: drift (Euler) vs bounded wobble (Verlet), and the effect of dt ----------
def run(step, dt, t_end):
    s = start; E0 = energy(*s); Es = []; rmin = 9; t_break = None
    for n in range(int(round(t_end / dt))):
        s = step(*s, dt); Es.append(energy(*s)); r = s[1] - s[0]; rmin = min(rmin, r)
        if t_break is None and r > 3.0: t_break = (n + 1) * dt          # atoms have flown apart: the molecule broke
    Es = np.array(Es); return Es[-1] - E0, np.abs(Es - E0).max(), rmin, t_break
print("\\nStep 6: energy error after 2 time units (about 3 vibrations); E0 = U(1.05) = -0.7575")
print("   dt       Euler: E(end) - E0     Verlet: max |E - E0|")
for dt in [0.02, 0.01, 0.005]:
    de = run(euler, dt, 2.0)[0]; dv = run(verlet, dt, 2.0)[1]
    print(f"   {dt:5.3f}    {de:10.4f}              {dv:10.6f}")
print("   halving dt: Euler error roughly halves (first order); Verlet error drops 4x (second order)")
for dt in [0.01, 0.005]:
    tb = run(euler, dt, 50.0)[3]; tv = run(verlet, dt, 50.0)[3]
    print(f"   over 50 time units with dt = {dt}: Euler molecule breaks at t = {tb}, Verlet molecule breaks: {tv}")

# ---------- Step 7: how big may dt be? ----------
k = 36 * 2**(2 / 3)                                    # curvature U'' at the minimum r = 2^(1/6)
mu = 0.5                                               # reduced mass of two atoms of mass 1
w = np.sqrt(k / mu); period = 2 * np.pi / w
print(f"\\nStep 7: curvature at the minimum k = 36 * 2^(2/3) = {k:.2f};  omega = sqrt(k/mu) = {w:.2f};  vibration period = {period:.3f}")
for dt in [0.01, 0.03, 0.06, 0.12, 0.2]:
    de, dv, rmin, _ = run(verlet, dt, 20.0)
    print(f"   Verlet dt = {dt:4.2f} ({period / dt:5.1f} steps per period): max |E - E0| = {dv:10.4g}, closest approach {rmin:.3f}")
`,
  },

  exercises: [
    { q: 'Compute U and F at r = 1.2. Is the force attractive or repulsive?', a: 'r⁻⁶ = 1/1.2⁶ = 1/2.9860 = 0.33490, r⁻¹² = 0.11216. U = 4 × (0.11216 − 0.33490) = −0.8910. r⁻¹³ = 0.11216/1.2 = 0.09346, r⁻⁷ = 0.33490/1.2 = 0.27908, F = 24 × (2×0.09346 − 0.27908) = −2.212. Negative, so attractive: the atoms are pulled back towards the minimum at 1.1225.' },
    { q: 'Start the atoms at rest at r = 1.2 and do one Euler step with dt = 0.01. What are v2, x2 and the change in energy?', a: 'v2 = 0 + (−2.212)(0.01) = −0.02212 and v1 = +0.02212. The positions do not change in the first Euler step (old velocity = 0), so U is unchanged and the energy rises by KE = 2 × ½ × 0.02212² = 0.00049. Energy appeared from nowhere: the Euler error.' },
    { q: 'Someone uses dt = 0.06 for this molecule. How many steps per vibration is that, and would you accept it?', a: 'Period 0.588 / 0.06 = 9.8 steps per vibration. The energy error is about 0.047, which is 6 % of |E| = 0.76: too large for production. Use dt ≤ 0.03 (about 20 steps per vibration).' },
  ],
};
