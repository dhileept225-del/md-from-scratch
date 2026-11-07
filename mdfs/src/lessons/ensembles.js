import { tex } from '../tex.js';

export default {
  id: 'ensembles',
  summary: 'An ensemble says which quantities a simulation holds fixed. Plain MD conserves energy (NVE); a thermostat gives constant temperature (NVT); adding a barostat gives constant pressure (NPT), like most experiments.',
  minutes: 13,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/ensembles.mp4', poster: 'videos/ensembles.jpg' },
  sim: 'ensembles',
  learn: [
    'What NVE, NVT and NPT hold fixed, and what is left to fluctuate',
    'Which ensemble to use for which question, and the usual NVT → NPT workflow',
    'How fluctuations give heat capacity and compressibility',
  ],
  simHint: 'Start in NVT and let the averages settle. Switch to NVE: the energy trace goes flat, but T keeps wobbling. Switch to NPT and set P = 0.3: the box expands until the pressure matches. Press Reset averages in each ensemble and compare the ± values. Which quantity has ± 0 in each case?',

  idea: `
<p>A simulation of N atoms in a box can be connected to the outside world in different ways. Statistical mechanics calls each choice an <b>ensemble</b>, named after the three quantities held fixed.</p>
<ul>
<li><b>NVE (microcanonical).</b> An isolated box: atoms, volume and total energy are fixed. This is what Newton's equations give with no extra terms (Lessons 1.3–1.6). Temperature and pressure fluctuate.</li>
<li><b>NVT (canonical).</b> The box exchanges energy with a large <b>heat bath</b> at temperature T. A <b>thermostat</b> imitates the bath by nudging velocities. The energy now fluctuates, and T is fixed only <i>on average</i>; the instantaneous T still fluctuates (Lesson 3.1).</li>
<li><b>NPT (isothermal–isobaric).</b> A <b>barostat</b> also lets the box volume change, like a piston, so the pressure is fixed on average. The density is then an <i>output</i>. Most experiments are at fixed pressure and temperature, so this is the ensemble that compares directly with the lab.</li>
</ul>
<p>For large systems the averages agree between ensembles (at the same state point). What differs is which quantities fluctuate, and how much. Those fluctuations are physics, not noise: energy fluctuations in NVT give the heat capacity, volume fluctuations in NPT give the compressibility.</p>
<p><b>A typical workflow:</b> minimise the energy → short NVT run to reach the temperature → NPT run until the density is stable → production run (NPT, or NVT at the equilibrated density, or NVE for clean dynamics).</p>
<div class="callout"><b>Key idea.</b> Choose the ensemble by what you want to compare. Measure densities in NPT. Check the integrator in NVE. Remember that thermostats and barostats change the dynamics a little (Lessons 3.3–3.6).</div>`,

  algorithm: `
<p>One MD step, and what each ensemble adds to it:</p>
<pre class="pseudo">velocity Verlet step                    # NVE: nothing else
if NVT or NPT:
    thermostat(v, T_target)             # e.g. Andersen: some atoms get new
                                        # Maxwell-Boltzmann velocities
if NPT:
    P = (N kB T + W / 3) / V            # virial pressure, W = Σ r·F over pairs
    mu = (1 - k (P_target - P))^(1/3)   # e.g. Berendsen barostat
    positions *= mu;  L *= mu           # rescale the box and the atoms</pre>
<ol class="steps">
  <li>The simple thermostat and barostat used here (Andersen, Berendsen) are easy to understand. Lessons 3.3–3.6 compare them with better choices (Nosé–Hoover, Langevin, Parrinello–Rahman, MTTK, C-rescale).</li>
  <li>In NPT, discard the part of the run where the box is still shrinking or growing. Then average.</li>
  <li>Do not compute dynamical properties (diffusion, viscosity) with a strong thermostat that randomises velocities often: it disturbs the dynamics.</li>
  <li>In LAMMPS these are <code>fix nve</code>, <code>fix nvt</code> and <code>fix npt</code>. In GROMACS, set <code>tcoupl</code> and <code>pcoupl</code> in the .mdp file.</li>
</ol>`,

  code: {
    note: 'The same 108 Lennard-Jones atoms, equilibrated at T = 1.5, then run in NVE, NVT (Andersen thermostat) and NPT (plus a Berendsen barostat at P = 1). The table shows which quantities stay fixed and which fluctuate. Run <code>pip install numpy</code>, then <code>python ensembles.py</code> (about 6 s).',
    src: `import numpy as np

def fcc(n, rho):
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return (c[:, None] + b[None]).reshape(-1, 3) * a, L

def forces(pos, L, rc=2.5):
    """LJ forces, potential energy U and virial W = sum over pairs of r . F."""
    d = pos[:, None] - pos[None]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = np.where(r2 < rc * rc, 1 / r2**3, 0.0)
    rF = 24 * (2 * inv6**2 - inv6)                       # r . F for each pair
    return ((rF / r2)[..., None] * d).sum(1), 2 * (inv6**2 - inv6).sum(), 0.5 * rF.sum()

T0, P0, dt = 1.5, 1.0, 0.004
rng = np.random.default_rng(1)

def md(s, ensemble, steps):
    """Velocity Verlet. NVT/NPT add an Andersen thermostat; NPT adds a Berendsen barostat."""
    pos, v, L = s["pos"].copy(), s["v"].copy(), s["L"]
    N = len(pos); Nf = 3 * N - 3; F, U, W = forces(pos, L); out = []
    for _ in range(steps):
        v += 0.5 * dt * F; pos = (pos + dt * v) % L
        F, U, W = forces(pos, L); v += 0.5 * dt * F
        if ensemble in ("NVT", "NPT"):                    # random collisions with a heat bath at T0
            hit = rng.random(N) < 5.0 * dt
            v[hit] = rng.normal(0, np.sqrt(T0), (hit.sum(), 3))
        K = 0.5 * (v ** 2).sum(); T = 2 * K / Nf; V = L**3
        P = (N * T + W / 3) / V                           # virial pressure
        if ensemble == "NPT":                             # scale the box gently towards P0
            mu = (1 - 0.01 * (P0 - P)) ** (1 / 3)
            pos *= mu; L *= mu
        out.append((K + U, T, P, N / V))
    return dict(pos=pos, v=v, L=L), np.array(out)

pos, L = fcc(3, 0.8); v = rng.normal(0, np.sqrt(T0), pos.shape); v -= v.mean(0)
state, _ = md(dict(pos=pos, v=v, L=L), "NVT", 1500)       # equilibrate at T0 first
print(f"108 LJ atoms, equilibrated in NVT at T = {T0}; then 3000 steps in each ensemble\\n")
print("ensemble   <E>/N   std(E)/N    <T>   std(T)    <P>   std(P)   density  std")
for ens in ["NVE", "NVT", "NPT"]:
    _, o = md(state, ens, 3000)
    if ens == "NPT": o = o[1000:]                         # let the box finish shrinking or growing
    E, T, P, rho = o[:, 0] / 108, o[:, 1], o[:, 2], o[:, 3]
    print(f"{ens:8s}  {E.mean():6.3f}   {E.std():7.4f}   {T.mean():5.3f}  {T.std():.3f}   {P.mean():5.2f}  {P.std():5.2f}    {rho.mean():.3f}  {rho.std():.4f}")
print("\\nNVE: E is fixed, T fluctuates. NVT: T is fixed on average, E fluctuates.")
print("NPT: P is fixed on average; the density adjusts (0.80 -> about 0.6) and fluctuates.")
`,
  },

  maths: () => `
<h4>1. The three ensembles and their probabilities</h4>
${tex(String.raw`\text{NVE: } \rho \propto \delta\big(H(\mathbf r, \mathbf p) - E\big), \qquad \text{NVT: } \rho \propto e^{-H/k_BT}, \qquad \text{NPT: } \rho \propto e^{-(H + PV)/k_BT}`, true)}
<h4>2. Pressure from the virial</h4>
${tex(String.raw`P = \frac{N k_B T}{V} + \frac{1}{3V}\Big\langle \sum_{i<j} \mathbf r_{ij}\cdot\mathbf F_{ij} \Big\rangle`, true)}
<p>The first term is the ideal gas; the second comes from the forces between atoms. In 2D, use 1/2 and the area instead of 1/3 and V.</p>
<h4>3. Fluctuations give response functions</h4>
${tex(String.raw`\text{NVT: } \langle \delta E^2 \rangle = k_B T^2\, C_V, \qquad \text{NPT: } \langle \delta V^2 \rangle = k_B T\, V\, \kappa_T`, true)}
<p>From the code tab: in NVT, σ<sub>E</sub>/N = 0.217 at T = 1.5 gives C<sub>V</sub>/N = N(σ<sub>E</sub>/N)²/T² ≈ 2.3 k<sub>B</sub> per atom. (A Berendsen barostat gives volume fluctuations that are too small, so do not use it for κ<sub>T</sub>.)</p>
<h4>4. Relative size of fluctuations</h4>
${tex(String.raw`\frac{\sigma_E}{\langle E\rangle} \sim \frac{1}{\sqrt N}`, true)}
<p>This is why ensembles give the same averages for large N: the distributions become very narrow.</p>`,

  quiz: [
    {
      q: 'Which quantity is exactly constant in an ideal NVE simulation?',
      options: ['Temperature', 'Pressure', 'Total energy', 'Potential energy'],
      answer: 2,
      why: 'Newton\'s equations conserve the total energy. Kinetic and potential energy trade off, so T and P fluctuate.',
    },
    {
      q: 'You want to predict the density of a new electrolyte at 298 K and 1 bar. Which ensemble?',
      options: ['NVE', 'NVT', 'NPT', 'Any, they all give the density'],
      answer: 2,
      why: 'Only NPT lets the volume adjust to the target pressure. In NVE or NVT you choose the density yourself, so it is an input, not a result.',
    },
    {
      q: 'In an NVT run the total energy fluctuates. What does the size of these fluctuations tell you?',
      options: ['The integrator is unstable', 'The heat capacity', 'The pressure', 'Nothing useful'],
      answer: 1,
      why: '⟨δE²⟩ = k_B T² C_V. Energy fluctuations in the canonical ensemble measure the heat capacity.',
    },
  ],
};
