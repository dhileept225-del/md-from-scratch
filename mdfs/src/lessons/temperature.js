import { tex } from '../tex.js';

export default {
  id: 'temperature',
  summary: 'In MD, temperature is not a dial you set. It is measured from how fast the atoms move. Learn how to compute it, why it fluctuates, and why it changes on its own at the start of a run.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/temperature.mp4', poster: 'videos/temperature.jpg' },
  sim: 'temperature',
  learn: [
    'The kinetic temperature and the degrees-of-freedom count 3N − 3',
    'Equipartition: each direction carries ½k<sub>B</sub>T',
    'Why instantaneous T fluctuates as 1/√N, and why it falls after a lattice start',
  ],
  simHint: 'Press Restart and watch T fall from its starting value within a few frames: kinetic energy flows into potential energy. Compare N = 36 with N = 400: the trace gets much calmer, while fluctuation × √N stays about the same. Then press Heat once. Does T rise by exactly 1.21×, and does it stay there?',

  idea: `
<p>In a real gas, temperature measures how vigorously the molecules move. MD uses the same idea: the <b>kinetic temperature</b> is defined by the average kinetic energy. Every independent way an atom can move (a <b>degree of freedom</b>) carries on average ½k<sub>B</sub>T of kinetic energy. This is <b>equipartition</b>. So we add up the kinetic energy and divide by the number of degrees of freedom.</p>
<p><b>Counting degrees of freedom.</b> N atoms in 3D have 3N velocity components. Our simulations keep the total momentum at zero (Lesson 2.3), so 3 of them are not free: N<sub>f</sub> = 3N − 3. Constraints remove more: a rigid water molecule has 6, not 9. Getting N<sub>f</sub> wrong biases T, by 1% for 100 atoms and much more for constrained molecules.</p>
<p><b>T fluctuates.</b> The kinetic energy at any instant is a sum over a finite number of atoms, so the instantaneous T jitters around its average. The relative size of the jitter shrinks as 1/√N. A 100-atom system flickers by several per cent; a million-atom system hardly at all. Only the time average is "the temperature".</p>
<p><b>T drifts at the start.</b> If you start atoms on a perfect lattice (zero potential energy stored in displacements) with velocities at T₀, energy flows from kinetic into potential as the atoms leave their sites. In NVE the total energy is fixed, so T falls, often to about half of T₀ for a solid near equilibrium. This is why we equilibrate, and why we need thermostats (Lessons 3.3–3.5) to hold T at a target.</p>
<div class="callout"><b>Key idea.</b> T = 2⟨KE⟩/(N<sub>f</sub>k<sub>B</sub>). The instantaneous value fluctuates (by about 1/√N) and is not something you set directly. Report time averages, after equilibration.</div>`,

  algorithm: `
<p>Measuring and setting the temperature:</p>
<pre class="pseudo">def temperature(v, m):
    KE = 0.5 * sum(m * v**2)           # all atoms, all directions
    Nf = 3 * N - 3 - n_constraints
    return 2 * KE / (Nf * kB)

# set an exact starting temperature
v = gaussian(0, sqrt(kB * T0 / m))     # Maxwell-Boltzmann (Lesson 2.3)
v -= mean(v)                           # zero total momentum
v *= sqrt(T0 / temperature(v))         # rescale to exactly T0</pre>
<ol class="steps">
  <li>Use the same N<sub>f</sub> everywhere. MD codes print it (GROMACS reports the "number of degrees of freedom"; LAMMPS: <code>compute_modify thermo_temp dof</code>).</li>
  <li>Rescaling all velocities by λ multiplies T by λ². This is the simplest thermostat, and the starting point for Berendsen (Lesson 3.3). In NVE, part of the added energy then flows into potential energy, so T settles below λ²T (in a solid, about halfway).</li>
  <li>Do not judge T from one frame. Average over many steps, after the initial drift has gone.</li>
  <li>After a lattice start, either start at about 2T<sub>target</sub> or, better, equilibrate with a thermostat before the production run.</li>
</ol>`,

  code: {
    note: 'A 3D Lennard-Jones NVE run from a perfect crystal: T falls from 2.0 to about 1.1, the x, y, z kinetic energies agree, and the fluctuation × √N stays constant for N = 32, 108 and 256. Run <code>pip install numpy</code>, then <code>python temperature.py</code> (about 7 s).',
    src: `import numpy as np

def fcc(n, rho):
    """n x n x n FCC lattice (4 atoms per cell) at reduced density rho."""
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return (c[:, None] + b[None]).reshape(-1, 3) * a, L

def forces(pos, L, rc=2.5):
    d = pos[:, None] - pos[None]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = np.where(r2 < rc * rc, 1 / r2**3, 0.0)
    F = (24 * (2 * inv6**2 - inv6) / r2)[..., None] * d
    return F.sum(1), 2 * (inv6**2 - inv6).sum()          # forces, potential energy

def temperature(v):
    """Kinetic temperature: KE = (Nf/2) kB T, with Nf = 3N - 3 (total momentum is fixed)."""
    N = len(v); Nf = 3 * N - 3
    return (v ** 2).sum() / Nf                            # m = kB = 1 (reduced units)

def run(n, T0, steps, dt=0.005, rho=0.8, seed=0):
    rng = np.random.default_rng(seed)
    pos, L = fcc(n, rho)
    v = rng.normal(0, np.sqrt(T0), pos.shape); v -= v.mean(0)   # Maxwell-Boltzmann, zero momentum
    v *= np.sqrt(T0 / temperature(v))                           # exactly T0 at t = 0
    F, _ = forces(pos, L); Ts, kex = [], []
    for _ in range(steps):                                      # velocity Verlet, NVE
        v += 0.5 * dt * F; pos = (pos + dt * v) % L
        F, U = forces(pos, L); v += 0.5 * dt * F
        Ts.append(temperature(v)); kex.append((v ** 2).mean(0))
    return np.array(Ts), np.array(kex)

# 1. start a perfect crystal at T = 2.0 and watch the temperature fall
Ts, kex = run(3, 2.0, 1500)
print("N = 108, started at T = 2.00 on a perfect lattice")
for s in [0, 25, 50, 100, 200, 500, 1000, 1499]:
    print(f"  step {s:5d}: T = {Ts[s]:.3f}")
print(f"  average after step 500: T = {Ts[500:].mean():.3f}   (about half the kinetic energy became potential energy)")

# 2. equipartition: every direction gets the same share
print("\\nmean v^2 in x, y, z:", np.round(kex[500:].mean(0), 3))

# 3. instantaneous T fluctuates; the fluctuation shrinks as 1/sqrt(N)
print("\\n   N   <T>    std(T)   std/<T>   std/<T> * sqrt(N)")
for n in [2, 3, 4]:
    Ts, _ = run(n, 2.0, 1500, seed=n)
    m, s = Ts[500:].mean(), Ts[500:].std(); N = 4 * n**3
    print(f"{N:4d}  {m:.3f}  {s:.4f}   {s / m:6.3f}    {s / m * np.sqrt(N):.2f}")
`,
  },

  maths: () => `
<h4>1. Kinetic temperature</h4>
${tex(String.raw`\langle K \rangle = \Big\langle \sum_i \tfrac12 m_i v_i^2 \Big\rangle = \tfrac12 N_f\,k_B T \quad\Rightarrow\quad T = \frac{2\langle K\rangle}{N_f k_B}, \qquad N_f = 3N - 3 - N_c`, true)}
<h4>2. Equipartition</h4>
${tex(String.raw`\big\langle \tfrac12 m v_x^2 \big\rangle = \big\langle \tfrac12 m v_y^2 \big\rangle = \big\langle \tfrac12 m v_z^2 \big\rangle = \tfrac12 k_B T`, true)}
<p>Each quadratic term in the energy gets ½k<sub>B</sub>T on average. It follows from the Maxwell–Boltzmann distribution (Lesson 2.3).</p>
<h4>3. Fluctuations</h4>
${tex(String.raw`\left(\frac{\sigma_T}{\langle T\rangle}\right)_{NVT} = \sqrt{\frac{2}{N_f}} \approx \sqrt{\frac{2}{3N}}`, true)}
<p>K is a sum of N<sub>f</sub> independent squared Gaussians (a χ² variable), which gives this result in the canonical (NVT) ensemble. In NVE the fluctuations are smaller, because K and U must trade off to keep E fixed, but they still scale as 1/√N, as the code tab shows (std/⟨T⟩ × √N ≈ 0.5).</p>
<h4>4. Lattice start (virial estimate)</h4>
<p>For a harmonic solid, the average potential and kinetic energies are equal. Starting with all energy kinetic, E = K₀ is shared equally once the atoms settle:</p>
${tex(String.raw`K_0 = \langle K\rangle + \langle U - U_{\text{lattice}}\rangle = 2\langle K\rangle \quad\Rightarrow\quad T \approx \tfrac12 T_0`, true)}
<h4>5. Velocity rescaling</h4>
${tex(String.raw`\mathbf v_i \to \lambda\,\mathbf v_i \;\Rightarrow\; T \to \lambda^2 T, \qquad \lambda = \sqrt{T_{\text{target}}/T}`, true)}`,

  quiz: [
    {
      q: 'You simulate 50 rigid water molecules (3 atoms each) and remove the total momentum. How many degrees of freedom are there?',
      options: ['450', '447', '297', '147'],
      answer: 2,
      why: '150 atoms give 450 components. Each rigid water has 3 constraints (two bonds and an angle): 450 − 150 = 300, then minus 3 for momentum = 297.',
    },
    {
      q: 'Going from 100 to 10,000 atoms, the relative fluctuation of the instantaneous temperature changes by about…',
      options: ['No change', '10 times smaller', '100 times smaller', '10 times larger'],
      answer: 1,
      why: 'It scales as 1/√N: √(10,000/100) = 10.',
    },
    {
      q: 'An NVE run starts on a perfect crystal with velocities at T = 300 K. A few hundred steps later T is about 155 K. What happened?',
      options: ['The integrator is broken', 'Energy leaks out through the box walls', 'Kinetic energy turned into potential energy as atoms left their lattice sites', 'The degrees of freedom were miscounted'],
      answer: 2,
      why: 'Total energy is conserved, but it is now shared between K and U. For a near-harmonic solid they end up about equal, so T roughly halves.',
    },
  ],
};
