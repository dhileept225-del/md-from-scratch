import { tex } from '../tex.js';

export default {
  id: 'langevin',
  summary: 'Langevin dynamics gives every atom its own small heat bath: a friction force plus random kicks, balanced so the system reaches the right temperature. It is robust and always ergodic, but strong friction slows the real dynamics.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/langevin.mp4', poster: 'videos/langevin.jpg' },
  sim: 'langevin',
  learn: [
    'The Langevin equation: friction plus random force, and the fluctuation–dissipation balance',
    'How the friction γ controls temperature relaxation, and how it slows diffusion',
    'The BAOAB integrator, and when to use Langevin dynamics',
  ],
  simHint: 'Start at γ = 1 and check that the spread of T matches the canonical value. Watch the white path of the tagged atom. Then press "Add D to chart" and try γ = 0, 5, 20 and 50 (wait for D to settle each time). How does D change? At γ = 50, compare D with the free-particle value kT/γ. Why is it smaller in a liquid?',

  idea: `
<p>Paul Langevin (1908) described a pollen grain in water as feeling two forces from the water molecules: a smooth <b>friction</b> that opposes its motion, and rapid <b>random kicks</b>. Langevin dynamics applies the same idea to every atom in a simulation, as if each atom were bathed in an invisible solvent.</p>
<p>The friction alone would cool the system to zero; the kicks alone would heat it forever. For the right temperature they must be balanced: the stronger the friction γ, the stronger the kicks must be. This balance is the <b>fluctuation–dissipation theorem</b>, and it makes the atoms sample the canonical (NVT) distribution at T.</p>
<p><b>Why it is popular.</b> Each atom is thermostatted on its own (a <i>local</i> thermostat), so there is no global ringing (unlike Nosé–Hoover) and no flying ice cube. The noise makes the dynamics ergodic even for small or stiff systems. It is very stable with large time steps, and it naturally models a solvent you do not simulate explicitly (implicit solvent, coarse-grained models, Brownian dynamics).</p>
<p><b>The price.</b> The friction and noise change the dynamics. Velocities forget their past in a time 1/γ, and at strong friction atoms move by diffusion like pollen grains. A free particle has D = k<sub>B</sub>T/(mγ); in a liquid, collisions with neighbours slow it further. Use a small γ (about 0.1–1 ps⁻¹ for atomistic systems) if you care about diffusion or viscosity; a large γ is fine if you only need equilibrium averages.</p>
<div class="callout"><b>Key idea.</b> Langevin = Newton + friction + noise, with the noise set by T and γ. Any γ > 0 gives the right ensemble; γ controls how much the dynamics is disturbed. Integrate with BAOAB.</div>`,

  algorithm: `
<p>The BAOAB splitting (Leimkuhler and Matthews, 2013) is simple and gives very accurate configurational averages:</p>
<pre class="pseudo">c1 = exp(-gamma * dt)
c2 = sqrt((1 - c1**2) * kB * T / m)

every step:
    v += (dt/2) * F / m                  # B: half kick
    x += (dt/2) * v                      # A: half drift
    v  = c1 * v + c2 * gaussian()        # O: friction + noise (exact)
    x += (dt/2) * v                      # A: half drift
    F  = forces(x)
    v += (dt/2) * F / m                  # B: half kick</pre>
<ol class="steps">
  <li>The O step solves the friction-plus-noise part exactly, so it works for any γ, even very large.</li>
  <li>Use a different random number for every atom, direction and step. Total momentum is not conserved, so N<sub>f</sub> = 3N.</li>
  <li>γ = 0 gives back velocity Verlet (NVE). 1/γ plays the role of the coupling time τ of the other thermostats (the kinetic energy relaxes with time constant 1/(2γ)).</li>
  <li>LAMMPS: <code>fix langevin T T damp seed</code> (damp = 1/γ) together with <code>fix nve</code>. GROMACS: <code>integrator = sd</code> with <code>tau-t</code> = 1/γ.</li>
</ol>`,

  code: {
    note: 'Part 1: free particles under Langevin, compared with exact formulas for heating, velocity memory and the mean squared displacement. Part 2: 108 Lennard-Jones atoms with BAOAB at four frictions: correct temperature and spread, but diffusion slows as γ grows. Run <code>pip install numpy</code>, then <code>python langevin.py</code> (about 10 s).',
    src: `import numpy as np
rng = np.random.default_rng(2)

def O_step(v, gamma, kT, dt):
    """Ornstein-Uhlenbeck step: exact friction + matching random kick (fluctuation-dissipation)."""
    c1 = np.exp(-gamma * dt)
    return c1 * v + np.sqrt((1 - c1**2) * kT) * rng.normal(size=v.shape)

# ---------- Part 1: free particles (no forces) ----------
print("Part 1: 20,000 free particles, start at T = 0, Langevin bath kT = 1, gamma = 2")
gamma, kT, dt = 2.0, 1.0, 0.01
v = np.zeros((20000, 3)); x = np.zeros_like(v)
for step in range(1, 601):
    x += 0.5 * dt * v; v = O_step(v, gamma, kT, dt); x += 0.5 * dt * v      # A O A (no forces)
    if step in (25, 50, 100, 300):
        t = step * dt
        print(f"  t = {t:4.2f}: T = {(v**2).mean():.3f}   theory 1 - exp(-2 gamma t) = {1 - np.exp(-2 * gamma * t):.3f}")
    if step == 300: v0, x0 = v.copy(), x.copy()                            # equilibrated: start measuring
    if step == 350: print(f"  velocity memory after 0.5: {(v * v0).sum(1).mean() / 3:.3f}   theory exp(-gamma t) = {np.exp(-gamma * 0.5):.3f}")
t = 300 * dt; msd = ((x - x0)**2).sum(1).mean()
print(f"  MSD after t = {t}: {msd:.3f}   theory 6 kT/gamma [t - (1 - exp(-gamma t))/gamma] = {6 * kT / gamma * (t - (1 - np.exp(-gamma * t)) / gamma):.3f}")
print(f"  -> long-time diffusion coefficient D = kT/(m gamma) = {kT / gamma}")

# ---------- Part 2: Lennard-Jones liquid with the BAOAB integrator ----------
def fcc(n, rho):
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return (c[:, None] + b[None]).reshape(-1, 3) * a, L

def forces(pos, L, rc=2.5):
    d = pos[:, None] - pos[None]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = np.where(r2 < rc * rc, 1 / r2**3, 0.0)
    return ((24 * (2 * inv6**2 - inv6) / r2)[..., None] * d).sum(1)

def baoab(gamma, T0=1.5, steps=4000, dt=0.004):
    pos, L = fcc(3, 0.8); N = len(pos); Nf = 3 * N
    v = rng.normal(0, np.sqrt(T0), pos.shape); F = forces(pos, L)
    unwrapped = pos.copy(); Ts, start = [], None
    for s in range(steps):
        g = 10.0 if s < 1000 else gamma                     # equilibrate with gamma = 10, then switch
        v += 0.5 * dt * F                                   # B
        unwrapped += 0.5 * dt * v                           # A
        v = O_step(v, g, T0, dt)                            # O  (friction + noise)
        unwrapped += 0.5 * dt * v                           # A
        F = forces(unwrapped % L, L); v += 0.5 * dt * F     # B
        Ts.append((v**2).sum() / Nf)
        if s == 1000: start = unwrapped.copy()              # measure after equilibration
    D = ((unwrapped - start)**2).sum(1).mean() / (6 * (steps - 1001) * dt)
    return np.array(Ts[1000:]), D, Nf

print("\\nPart 2: 108 LJ atoms, BAOAB Langevin at T = 1.5")
print("  gamma    <T>     std(T)   canonical   diffusion D")
for g in [0.5, 2.0, 10.0, 50.0]:
    T, D, Nf = baoab(g)
    print(f"  {g:5.1f}   {T.mean():.3f}   {T.std():.3f}     {np.sqrt(2 / Nf) * 1.5:.3f}      {D:.4f}")
print("Every gamma > 0 targets the canonical T and spread (weak friction just needs longer runs to average);")
print("strong friction slows the dynamics: D falls roughly as 1/gamma at large gamma.")
`,
  },

  maths: () => `
<h4>1. The Langevin equation</h4>
${tex(String.raw`m\,\dot{\mathbf v}_i = \mathbf F_i - \gamma m\,\mathbf v_i + \mathbf R_i(t), \qquad \langle R_{i\alpha}(t)\,R_{j\beta}(t')\rangle = 2\gamma m k_B T\,\delta_{ij}\delta_{\alpha\beta}\,\delta(t - t')`, true)}
<p>The noise strength 2γmk<sub>B</sub>T is fixed by the fluctuation–dissipation theorem: it is the only choice that makes the Maxwell–Boltzmann distribution stationary.</p>
<h4>2. A free particle (F = 0)</h4>
${tex(String.raw`\langle v^2\rangle(t) = \frac{k_BT}{m}\big(1 - e^{-2\gamma t}\big) \;\;(v(0)=0), \qquad \langle v(0)\,v(t)\rangle = \frac{k_BT}{m}e^{-\gamma t}`, true)}
${tex(String.raw`\langle \Delta x^2 \rangle = \frac{2k_BT}{m\gamma}\Big[t - \frac{1 - e^{-\gamma t}}{\gamma}\Big] \;\longrightarrow\; 2Dt, \qquad D = \frac{k_BT}{m\gamma}`, true)}
<p>This is the Einstein relation. Per dimension the MSD is 2Dt, so 6Dt in 3D and 4Dt in 2D. The code tab checks all three formulas.</p>
<h4>3. The exact O step</h4>
${tex(String.raw`v(t + \Delta t) = e^{-\gamma\Delta t}\,v(t) + \sqrt{\frac{k_BT}{m}\big(1 - e^{-2\gamma\Delta t}\big)}\;\mathcal N(0,1)`, true)}
<h4>4. Overdamped limit (Brownian dynamics)</h4>
${tex(String.raw`\gamma \gg \text{ all vibrational frequencies:}\quad \dot{\mathbf x} = \frac{\mathbf F}{m\gamma} + \sqrt{2D}\,\boldsymbol\eta(t)`, true)}`,

  quiz: [
    {
      q: 'What happens if you keep the friction in Langevin dynamics but switch off the random force?',
      options: ['Nothing changes', 'The system cools towards zero temperature', 'The system heats up forever', 'Energy is exactly conserved'],
      answer: 1,
      why: 'Friction alone only removes energy. The random kicks put it back; their size is set by the fluctuation–dissipation theorem.',
    },
    {
      q: 'You want the diffusion coefficient of water from a Langevin run. What should you do?',
      options: ['Use a very large γ for stability', 'Use a small γ (or NVE / Nosé–Hoover for the production run)', 'γ does not affect diffusion', 'Use γ = 0 and no thermostat during equilibration'],
      answer: 1,
      why: 'Friction slows the motion: at strong friction D → k_BT/(mγ). A weak coupling barely disturbs the dynamics.',
    },
    {
      q: 'In the BAOAB integrator, what does the O step do?',
      options: ['Computes the forces', 'Moves the positions', 'Applies friction and random noise to the velocities, exactly', 'Rescales the box'],
      answer: 2,
      why: 'O is the Ornstein–Uhlenbeck part: v → c₁v + c₂·noise. B steps are kicks from the forces, A steps move the positions.',
    },
  ],
};
