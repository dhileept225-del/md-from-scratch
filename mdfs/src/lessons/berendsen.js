import { tex } from '../tex.js';

export default {
  id: 'berendsen',
  summary: 'The simplest way to control temperature is to rescale all velocities. Berendsen does it gently, with a coupling time τ. It reaches the target smoothly, but it squeezes the natural fluctuations, so it is a tool for equilibration, not for production.',
  minutes: 13,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/berendsen.mp4', poster: 'videos/berendsen.jpg' },
  sim: 'berendsen',
  learn: [
    'Velocity rescaling, and why rescaling exactly every step is too harsh',
    'The Berendsen thermostat and its coupling time τ',
    'Why Berendsen gives the wrong fluctuations, and how stochastic rescaling (V-rescale) fixes it',
  ],
  simHint: 'With Berendsen at τ = 0.5, move the target from 1.0 to 1.5 and read the 63% time. Try τ = 0.05 and τ = 5. Then compare the spread of T with the canonical value for Rescale, Berendsen and Bussi. Which one gets it right? (None, the NVE case, is also too small, but for a different reason: there the total energy is fixed.)',

  idea: `
<p>In NVT we want the temperature to be T₀ on average. The temperature is set by the kinetic energy (Lesson 3.1), and multiplying every velocity by λ multiplies T by λ². So the simplest thermostat is: measure T, then scale all velocities.</p>
<p><b>Exact rescaling</b> picks λ = √(T₀/T) at every step. T is then exactly T₀ at all times. That sounds ideal but is wrong: a real system in contact with a heat bath has a fluctuating kinetic energy (spread √(2/N<sub>f</sub>)·T). Exact rescaling removes the fluctuations entirely and disturbs the dynamics.</p>
<p><b>Berendsen (1984)</b> scales more gently, so that T relaxes towards T₀ exponentially with a time constant τ:</p>
<p style="text-align:center">dT/dt = (T₀ − T)/τ.</p>
<p>A small τ means strong coupling (fast, close to exact rescaling); a large τ means weak coupling (close to NVE). Typical values: τ = 0.1–1 ps. It is smooth, robust, and great for bringing a system to temperature.</p>
<p><b>The catch.</b> Berendsen does not sample the canonical ensemble. The kinetic energy fluctuations come out too small, and averages that depend on them (heat capacity, some free energies) are wrong. In long runs it can also slowly pump energy from fast motions into slow ones, such as the overall translation or rotation of a molecule: the "<b>flying ice cube</b>".</p>
<p><b>The fix: stochastic velocity rescaling</b> (Bussi, Donadio and Parrinello, 2007). Like Berendsen, but with a carefully chosen random term added to the kinetic energy, so it produces exactly the canonical distribution. It keeps Berendsen's simplicity. GROMACS calls it <code>V-rescale</code>; LAMMPS has <code>fix temp/csvr</code>.</p>
<div class="callout"><b>Key idea.</b> Berendsen relaxes T to the target with time constant τ. Use it to equilibrate. For production runs, use V-rescale (Bussi) or Nosé–Hoover (Lesson 3.4), which give correct fluctuations.</div>`,

  algorithm: `
<p>Every step, after the velocity Verlet update:</p>
<pre class="pseudo">K  = 0.5 * sum(m * v**2)
K0 = 0.5 * Nf * kB * T0                 # target kinetic energy

# exact rescaling (too harsh)
lam = sqrt(K0 / K)

# Berendsen
lam = sqrt(1 + dt / tau * (K0 / K - 1))

# Bussi (V-rescale): Berendsen-like, plus a random term
c   = exp(-dt / tau)
R1  = gaussian();  S = chi_square(Nf - 1)
lam = sqrt(c + (1 - c) * K0 / (Nf * K) * (R1**2 + S)
           + 2 * R1 * sqrt(c * (1 - c) * K0 / (Nf * K)))

v *= lam</pre>
<ol class="steps">
  <li>τ must be much larger than Δt (τ/Δt ≳ 50). Otherwise Berendsen is just exact rescaling.</li>
  <li>The measured relaxation is a bit slower than τ, because some of the added energy flows into potential energy.</li>
  <li>Thermostat groups separately when they exchange heat slowly (e.g. solute and solvent), and remove centre-of-mass motion to avoid a flying ice cube.</li>
  <li>GROMACS: <code>tcoupl = berendsen</code> or <code>v-rescale</code>, with <code>tau_t</code> and <code>ref_t</code>. LAMMPS: <code>fix temp/berendsen</code> or <code>fix temp/csvr</code>, plus <code>fix nve</code>.</li>
</ol>`,

  code: {
    note: '108 Lennard-Jones atoms are equilibrated at T = 1.0, then the target jumps to 1.5. It compares exact rescaling, Berendsen (τ = 0.1 and 1.0) and Bussi: how fast each reaches the target, and how big the fluctuations are. Run <code>pip install numpy</code>, then <code>python berendsen.py</code> (about 6 s).',
    src: `import numpy as np

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

rng = np.random.default_rng(3)
dt = 0.004

def lam_rescale(K, K0, Nf, tau):          # exact rescaling to the target every step
    return np.sqrt(K0 / K)

def lam_berendsen(K, K0, Nf, tau):        # Berendsen: dT/dt = (T0 - T) / tau
    return np.sqrt(1 + dt / tau * (K0 / K - 1))

def lam_bussi(K, K0, Nf, tau):            # Bussi-Donadio-Parrinello stochastic velocity rescaling
    c = np.exp(-dt / tau); R1 = rng.normal(); S = rng.chisquare(Nf - 1)
    a2 = c + (1 - c) * K0 / (Nf * K) * (R1**2 + S) + 2 * R1 * np.sqrt(c * (1 - c) * K0 / (Nf * K))
    sign = np.sign(R1 + np.sqrt(c * Nf * K / ((1 - c) * K0)))
    return sign * np.sqrt(a2)

def md(pos, v, L, thermostat, T0, tau, steps):
    N = len(pos); Nf = 3 * N - 3; K0 = 0.5 * Nf * T0
    F = forces(pos, L); Ts = []
    for _ in range(steps):
        v += 0.5 * dt * F; pos = (pos + dt * v) % L
        F = forces(pos, L); v += 0.5 * dt * F
        K = 0.5 * (v ** 2).sum()
        v *= thermostat(K, K0, Nf, tau)          # the thermostat: scale all velocities by lambda
        Ts.append((v ** 2).sum() / Nf)
    return pos, v, np.array(Ts)

pos, L = fcc(3, 0.8); N = len(pos); Nf = 3 * N - 3
v = rng.normal(0, 1, pos.shape); v -= v.mean(0)
pos, v, _ = md(pos, v, L, lam_bussi, 1.0, 0.1, 1000)       # equilibrate at T = 1.0

print("Equilibrated at T = 1.0, then the target is raised to T0 = 1.5\\n")
print(f"canonical spread of instantaneous T: sqrt(2/Nf) * T0 = {np.sqrt(2 / Nf) * 1.5:.3f}\\n")
print("thermostat          tau    time to close 63% of gap   <T>     std(T)")
for name, th, tau in [("rescale every step", lam_rescale, 1.0), ("Berendsen", lam_berendsen, 0.1),
                      ("Berendsen", lam_berendsen, 1.0), ("Bussi (V-rescale)", lam_bussi, 0.1)]:
    _, _, T = md(pos.copy(), v.copy(), L, th, 1.5, tau, 2500)
    smooth = np.convolve(T, np.ones(25) / 25, mode="same")
    t63 = np.argmax(smooth > 1.0 + 0.63 * 0.5) * dt
    print(f"{name:18s}  {tau:4.1f}   {t63:20.2f}       {T[1000:].mean():.3f}   {T[1000:].std():.3f}")
print("\\nRescaling and Berendsen give too-small fluctuations; Bussi matches the canonical value.")
`,
  },

  maths: () => `
<h4>1. Berendsen scaling factor</h4>
<p>We want dT/dt = (T₀ − T)/τ. Over one step, T → λ²T, so</p>
${tex(String.raw`\lambda^2 T - T = \frac{\Delta t}{\tau}(T_0 - T) \quad\Rightarrow\quad \lambda = \sqrt{1 + \frac{\Delta t}{\tau}\left(\frac{T_0}{T} - 1\right)}`, true)}
<p>The solution relaxes exponentially: T(t) = T₀ + (T(0) − T₀) e<sup>−t/τ</sup>. With τ = Δt this becomes exact rescaling.</p>
<h4>2. The canonical fluctuation it should reproduce</h4>
${tex(String.raw`\frac{\sigma_K}{\langle K\rangle} = \frac{\sigma_T}{\langle T\rangle} = \sqrt{\frac{2}{N_f}}`, true)}
<p>For 108 atoms (N<sub>f</sub> = 321) at T = 1.5 this is 0.118. The code tab gives 0.054 with Berendsen (τ = 0.1) and 0.114 with Bussi.</p>
<h4>3. Stochastic velocity rescaling (Bussi)</h4>
${tex(String.raw`dK = \frac{K_0 - K}{\tau}\,dt + 2\sqrt{\frac{K K_0}{N_f \tau}}\,dW`, true)}
<p>The first term is Berendsen. The second, a random (Wiener) term, restores exactly the canonical distribution of K. For large τ it approaches NVE; for small τ, strong coupling.</p>
<h4>4. Choosing τ</h4>
${tex(String.raw`\Delta t \ll \tau, \qquad \tau \approx 0.1\text{–}1\ \text{ps for molecular systems}`, true)}`,

  quiz: [
    {
      q: 'Berendsen with a very small τ (close to the time step) behaves like…',
      options: ['NVE', 'Exact velocity rescaling every step', 'Langevin dynamics', 'A barostat'],
      answer: 1,
      why: 'λ² = 1 + (Δt/τ)(T₀/T − 1). With τ = Δt this gives λ² = T₀/T, which is exact rescaling.',
    },
    {
      q: 'Why should you not use Berendsen to calculate a heat capacity from energy fluctuations?',
      options: ['It is too slow', 'It does not sample the canonical ensemble: the fluctuations are too small', 'It changes the number of atoms', 'It only works for liquids'],
      answer: 1,
      why: 'C_V comes from ⟨δE²⟩ in NVT. Berendsen suppresses kinetic-energy fluctuations, so the result would be wrong. Use V-rescale or Nosé–Hoover.',
    },
    {
      q: 'What does the Bussi (V-rescale) thermostat add to Berendsen?',
      options: ['A piston for the volume', 'Friction on every atom', 'A random term that gives the correct canonical fluctuations', 'A second time step'],
      answer: 2,
      why: 'It rescales all velocities like Berendsen, but adds a stochastic term to the kinetic energy chosen so the distribution is exactly canonical.',
    },
  ],
};
