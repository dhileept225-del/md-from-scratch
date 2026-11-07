import { tex } from '../tex.js';

export default {
  id: 'nose-hoover',
  summary: 'Nosé–Hoover adds one extra dynamical variable, a friction ξ, that brakes the atoms when they are too hot and pushes them when too cold. It is deterministic, samples the canonical ensemble, and conserves an extended energy. Chains of thermostats fix its one weakness.',
  minutes: 15,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/nose-hoover.mp4', poster: 'videos/nose-hoover.jpg' },
  sim: 'nosehoover',
  learn: [
    'The Nosé–Hoover equations: a friction ξ driven by the temperature error',
    'The thermostat mass Q (or period τ), the conserved extended energy, and ringing',
    'Why one Nosé–Hoover thermostat can fail (ergodicity) and how chains fix it',
  ],
  simHint: 'Liquid: after the warm-up, the spread of T should be close to canonical (ratio ≈ 1). Press Kick and watch T and ξ ring: try τ = 0.05 and τ = 2. Then switch to One oscillator. With plain Nosé–Hoover the phase space fills only a ring; switch to the chain and watch it fill the Gaussian.',

  idea: `
<p>Berendsen (Lesson 3.3) reaches the target temperature but gives the wrong fluctuations. Nosé (1984) and Hoover (1985) found a deterministic way to get them right. The heat bath is represented by <b>one extra dynamical variable</b>, a friction coefficient ξ, with its own equation of motion:</p>
<ul>
<li>The atoms feel an extra force −ξ·m·v. If ξ > 0 it brakes them; if ξ < 0 it pushes them.</li>
<li>ξ itself grows when the system is hotter than T₀ and shrinks when it is colder. It is like a "thermal piston" with inertia Q.</li>
</ul>
<p>Because ξ has inertia, it overshoots: T and ξ oscillate around the target with a period set by Q. Writing Q = N<sub>f</sub>k<sub>B</sub>T₀τ² makes τ that period. Typical choices are τ ≈ 0.1–1 ps (about 100 time steps or more). Too small a τ gives fast, violent oscillations; too large makes the coupling slow.</p>
<p><b>Why it works.</b> One can prove that, if the dynamics is ergodic, the atoms sample exactly the canonical distribution. There is also a <b>conserved quantity</b>, an extended energy that includes the thermostat terms. Its drift is your energy check in NVT, just like E in NVE (Lesson 1.6).</p>
<p><b>Two weaknesses.</b> (1) A big disturbance makes the thermostat <b>ring</b> for a long time; equilibrate first (Berendsen or V-rescale), then switch. (2) For small or very regular systems (a harmonic oscillator, a stiff crystal), the dynamics is <b>not ergodic</b>: trajectories stay on a narrow band and never fill the canonical distribution. The fix is a <b>Nosé–Hoover chain</b> (Martyna, Klein and Tuckerman, 1992): ξ₁ is itself thermostatted by ξ₂, and so on. A chain of 3–5 is standard.</p>
<div class="callout"><b>Key idea.</b> Nosé–Hoover chains are the standard deterministic thermostat for production NVT. Pick τ ≈ 100 time steps or more, equilibrate before switching on, and monitor the extended energy.</div>`,

  algorithm: `
<p>A symmetric step (the thermostat acts for half a step before and after velocity Verlet):</p>
<pre class="pseudo">Q = Nf * kB * T0 * tau**2

def thermostat_half(v, h):              # advance xi and scale v over time h
    xi += (h/2) * (sum(m v**2) - Nf kB T0) / Q
    v  *= exp(-xi * h)
    eta += xi * h
    xi += (h/2) * (sum(m v**2) - Nf kB T0) / Q

every step:
    thermostat_half(v, dt/2)
    velocity Verlet step
    thermostat_half(v, dt/2)
    H_ext = K + U + Q xi**2 / 2 + Nf kB T0 eta    # should stay constant</pre>
<ol class="steps">
  <li>For a chain, each ξ<sub>j</sub> is driven by the kinetic energy of the link below and damped by the link above. The code tab shows the full routine for a chain of any length.</li>
  <li>Do not start Nosé–Hoover far from the target temperature: it rings. Equilibrate with Berendsen or V-rescale first.</li>
  <li>The friction is a smooth, global force, so dynamics (diffusion, spectra) is barely disturbed when τ is large enough.</li>
  <li>LAMMPS: <code>fix nvt temp T T tau</code> (a chain of 3 by default; <code>tau</code> ≈ 100 × timestep). GROMACS: <code>tcoupl = nose-hoover</code>, <code>nsttcouple</code>, <code>nh-chain-length</code> (default 10; the leap-frog <code>md</code> integrator supports only 1).</li>
</ol>`,

  code: {
    note: 'Part 1: Nosé–Hoover on 108 Lennard-Jones atoms. Started far from the target it rings; after a Berendsen warm-up it gives canonical fluctuations and conserves the extended energy. Part 2: one harmonic oscillator, where plain NH fails and a chain of 3 works. Run <code>pip install numpy</code>, then <code>python nose_hoover.py</code> (about 8 s).',
    src: `import numpy as np

# ---------- Part 1: Nosé-Hoover on a Lennard-Jones liquid ----------
def fcc(n, rho):
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return (c[:, None] + b[None]).reshape(-1, 3) * a, L

def forces(pos, L, rc=2.5):
    d = pos[:, None] - pos[None]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = np.where(r2 < rc * rc, 1 / r2**3, 0.0)
    ush = 4 * (rc**-12 - rc**-6)
    U = 0.5 * np.where(r2 < rc * rc, 4 * (inv6**2 - inv6) - ush, 0).sum()
    return ((24 * (2 * inv6**2 - inv6) / r2)[..., None] * d).sum(1), U

def berendsen(pos, v, L, T0, steps, dt=0.004, tau=0.1):
    """Quick equilibration with Berendsen (Lesson 3.3)."""
    Nf = 3 * len(pos) - 3; F, U = forces(pos, L)
    for _ in range(steps):
        v += 0.5 * dt * F; pos = (pos + dt * v) % L; F, U = forces(pos, L); v += 0.5 * dt * F
        v *= np.sqrt(1 + dt / tau * (0.5 * Nf * T0 / (0.5 * (v**2).sum()) - 1))
    return pos, v

def nose_hoover(pos, v, L, T0, tau, steps, dt=0.004):
    N = len(pos); Nf = 3 * N - 3
    Q = Nf * T0 * tau**2                                      # thermostat "mass" from its period tau
    xi = eta = 0.0                                            # friction and its time integral
    F, U = forces(pos, L); out = []
    for _ in range(steps):
        K = 0.5 * (v**2).sum(); xi += 0.5 * dt * (2 * K - Nf * T0) / Q; eta += 0.5 * dt * xi
        v = v * np.exp(-0.5 * dt * xi)                        # friction (negative xi = heating)
        v += 0.5 * dt * F; pos = (pos + dt * v) % L
        F, U = forces(pos, L); v += 0.5 * dt * F
        v *= np.exp(-0.5 * dt * xi)
        K = 0.5 * (v**2).sum(); eta += 0.5 * dt * xi; xi += 0.5 * dt * (2 * K - Nf * T0) / Q
        H_ext = K + U + 0.5 * Q * xi**2 + Nf * T0 * eta      # conserved "extended energy"
        out.append((2 * K / Nf, H_ext))
    return np.array(out)

rng = np.random.default_rng(0); T0 = 1.5
pos0, L = fcc(3, 0.8); Nf = 3 * len(pos0) - 3
v0 = rng.normal(0, 1.0, pos0.shape); v0 -= v0.mean(0)        # T = 1 on a perfect lattice: far from T0

print(f"Part 1: Nose-Hoover, 108 LJ atoms, target T = {T0}; canonical std(T) = {np.sqrt(2 / Nf) * T0:.3f}")
o = nose_hoover(pos0, v0.copy(), L, T0, 0.5, 3000)
print(f"  a) straight from the lattice:   <T> = {o[1000:, 0].mean():.3f}   std(T) = {o[1000:, 0].std():.3f}   <- still ringing")
pos, v = berendsen(pos0, v0.copy(), L, T0, 1000)
o = nose_hoover(pos, v, L, T0, 0.2, 3000)
print(f"  b) Berendsen first, then NH:   <T> = {o[:, 0].mean():.3f}   std(T) = {o[:, 0].std():.3f}   <- canonical")
print(f"     extended energy drift: {o[-1, 1] - o[0, 1]:+.4f}   (kinetic energy ~ {0.5 * Nf * T0:.0f})")

# ---------- Part 2: a single harmonic oscillator, NH vs NH chain ----------
def oscillator(M, steps=200000, dt=0.02, Q=1.0, kT=1.0):
    """1D harmonic oscillator (m = omega = 1) with a Nose-Hoover chain of length M (M = 1: plain NH)."""
    x, p = 1.0, 0.0; xi = np.zeros(M); xs = np.empty(steps)
    def G(j):                                                  # force on thermostat j
        return ((p * p if j == 0 else Q * xi[j - 1] ** 2) - kT) / Q
    def chain_half(h):                                         # integrate the chain for time h
        nonlocal p
        for j in reversed(range(M)):
            if j < M - 1: xi[j] *= np.exp(-xi[j + 1] * h / 4)
            xi[j] += G(j) * h / 2
            if j < M - 1: xi[j] *= np.exp(-xi[j + 1] * h / 4)
        p *= np.exp(-xi[0] * h)
        for j in range(M):
            if j < M - 1: xi[j] *= np.exp(-xi[j + 1] * h / 4)
            xi[j] += G(j) * h / 2
            if j < M - 1: xi[j] *= np.exp(-xi[j + 1] * h / 4)
    for s in range(steps):
        chain_half(dt / 2); p -= 0.5 * dt * x; x += dt * p; p -= 0.5 * dt * x; chain_half(dt / 2)
        xs[s] = x
    return xs

print("\\nPart 2: one harmonic oscillator at kT = 1 (canonical: <x^2> = 1, kurtosis ratio = 1)")
for M, name in [(1, "Nose-Hoover"), (3, "NH chain, M = 3")]:
    xs = oscillator(M)[10000:]
    x2 = (xs**2).mean(); kurt = (xs**4).mean() / (3 * x2**2)
    print(f"  {name:16s}  <x^2> = {x2:.3f}   <x^4>/(3<x^2>^2) = {kurt:.3f}")
print("A single NH thermostat cannot make one oscillator canonical (not ergodic); a chain fixes it.")
`,
  },

  maths: () => `
<h4>1. Nosé–Hoover equations</h4>
${tex(String.raw`\dot{\mathbf r}_i = \frac{\mathbf p_i}{m_i}, \qquad \dot{\mathbf p}_i = \mathbf F_i - \xi\,\mathbf p_i, \qquad \dot\xi = \frac{1}{Q}\Big(\sum_i \frac{p_i^2}{m_i} - N_f k_B T_0\Big)`, true)}
<h4>2. Conserved extended energy</h4>
${tex(String.raw`H' = K + U + \tfrac12 Q\,\xi^2 + N_f k_B T_0\,\eta, \qquad \dot\eta = \xi`, true)}
<p>Differentiate: the friction removes 2Kξ from K, the ξ² term returns (2K − N<sub>f</sub>k<sub>B</sub>T₀)ξ, and the η term supplies N<sub>f</sub>k<sub>B</sub>T₀ξ. They cancel, so dH′/dt = 0.</p>
<h4>3. Thermostat mass and period</h4>
${tex(String.raw`Q = N_f k_B T_0\,\tau^2 \qquad\Rightarrow\qquad \text{oscillation period of } T \sim \tau`, true)}
<p>For a small deviation, linearising gives an oscillator for δT with frequency ω ≈ √2/τ (for a system with no potential energy).</p>
<h4>4. Why it is canonical (if ergodic)</h4>
<p>The equations conserve H′ and have a phase-space compressibility of −N<sub>f</sub>ξ. Together these imply that the stationary distribution of the physical variables is</p>
${tex(String.raw`\rho(\mathbf r, \mathbf p) \propto e^{-(K + U)/k_B T_0}`, true)}
<h4>5. Nosé–Hoover chain</h4>
${tex(String.raw`\dot\xi_1 = \frac{2K - N_f k_B T_0}{Q_1} - \xi_1\xi_2, \qquad \dot\xi_j = \frac{Q_{j-1}\xi_{j-1}^2 - k_B T_0}{Q_j} - \xi_j\xi_{j+1}, \qquad \dot\xi_M = \frac{Q_{M-1}\xi_{M-1}^2 - k_B T_0}{Q_M}`, true)}
<p>Usually Q₁ = N<sub>f</sub>k<sub>B</sub>T₀τ² and Q<sub>j</sub> = k<sub>B</sub>T₀τ² for j > 1.</p>`,

  quiz: [
    {
      q: 'In Nosé–Hoover, what happens to the friction ξ when the system is hotter than the target?',
      options: ['ξ is set to zero', 'ξ increases, so the atoms are slowed down', 'ξ decreases, so the atoms speed up', 'The time step is reduced'],
      answer: 1,
      why: 'dξ/dt ∝ (2K − N_f k_B T₀), which is positive when T > T₀. A positive ξ acts as friction −ξp.',
    },
    {
      q: 'Your NVT run uses Nosé–Hoover. What should you monitor instead of the total energy E to check the integration?',
      options: ['The pressure', 'The extended energy H′ = K + U + Qξ²/2 + N_f k_B T₀ η', 'The friction ξ', 'Nothing, NVT has no conserved quantity'],
      answer: 1,
      why: 'E is not conserved in NVT, but H′ is. A drifting H′ means the time step is too large or something is wrong.',
    },
    {
      q: 'Why are Nosé–Hoover chains used instead of a single Nosé–Hoover thermostat?',
      options: ['They are faster to compute', 'A single thermostat can be non-ergodic for small or stiff systems; chains fill phase space properly', 'Chains also control the pressure', 'They avoid the need for a time step'],
      answer: 1,
      why: 'For a harmonic oscillator, a single NH thermostat traps the trajectory on a narrow band. Thermostatting the thermostat makes the motion chaotic enough to sample the canonical distribution.',
    },
  ],
};
