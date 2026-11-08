import { tex } from '../tex.js';

export default {
  id: 'metadynamics',
  summary: 'Many important events (a reaction, a conformational change, an ion leaving its solvation shell) are separated by free-energy barriers that ordinary MD almost never crosses. Metadynamics adds a history-dependent bias along a few collective variables that fills the free-energy wells, so the system escapes; the accumulated bias then gives the free-energy surface.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/metadynamics.mp4', poster: 'videos/metadynamics.jpg' },
  sim: 'metad',
  learn: [
    'Why rare events need enhanced sampling, and what a collective variable (CV) is',
    'How metadynamics builds a bias from Gaussian hills and why the bias estimates the free energy',
    'Well-tempered metadynamics, convergence, and the importance of choosing a good CV',
  ],
  simHint: 'Start with Plain MD: the walker stays in the left basin. Switch to Metadynamics or Well-tempered with CV = x: the bias (red) fills the basins, the walker hops, and the amber estimate approaches the exact free energy (cyan). Plain metadynamics keeps fluctuating; well-tempered settles. Now choose CV = y: crossings become rarer and the estimate converges much more slowly, because y is only partly related to the real barrier. Try bigger hills (faster but rougher) and a higher barrier.',

  idea: `
<p>MD samples states with probability ∝ e<sup>−F/kT</sup>. A barrier of ΔF slows transitions by roughly e<sup>ΔF/kT</sup>. For a 10 kT barrier, that is about 20 000 times slower than a barrier-free motion, and 20 kT barriers (about 50 kJ/mol at room temperature) are effectively never crossed in a few-nanosecond run (code tab, Part 1). <b>Enhanced sampling</b> methods change the simulation so that rare events happen, and then correct for the change.</p>
<p><b>Collective variables.</b> We first choose one or two <b>collective variables</b> s(<b>r</b>): low-dimensional coordinates that distinguish the states and describe the slow motion. Examples: a distance (Li⁺–anion), a coordination number (how many oxygens surround Li⁺), a dihedral, a path coordinate between two structures.</p>
<p><b>Metadynamics</b> (Laio and Parrinello 2002) then adds a bias V(s) that grows in time: every few hundred steps a small Gaussian "hill" is placed at the current value of s. Visited regions get filled, so the system is pushed to explore new values, crosses barriers and eventually diffuses freely along s. When the wells are filled the bias mirrors the free energy, F(s) ≈ −V(s) + constant (Part 2).</p>
<p><b>Well-tempered metadynamics</b> (Barducci, Bussi and Parrinello 2008) makes the hills smaller where the bias is already large: w = w₀ e<sup>−V(s)/kΔT</sup>. The bias converges smoothly instead of fluctuating, and F(s) = −(T + ΔT)/ΔT · V(s). The <b>bias factor</b> γ = (T + ΔT)/T sets how strongly the barriers are flattened (Part 3).</p>
<p><b>The CV is everything.</b> If the chosen CV misses the real slow motion, the bias pushes the system in an irrelevant direction and the barrier is never crossed, or only with hysteresis (Part 4 and the widget with CV = y). Signs of a poor CV: few back-and-forth transitions, an estimate that keeps drifting, and different answers from different starting points.</p>
<p>Related methods: umbrella sampling (fixed biases in windows + WHAM/MBAR), adaptive biasing force, replica exchange (temperature or Hamiltonian), and OPES (on-the-fly probability enhanced sampling), the newer method from the Parrinello group.</p>
<div class="callout"><b>Key idea.</b> Add Gaussian hills along a good collective variable until the wells are filled; the negative bias is the free energy. Use well-tempered metadynamics, check convergence over time, and question your CV.</div>`,

  algorithm: `
<pre class="pseudo">choose CV s(r) and its range; hill height w0 (~0.5-1 kT), width σ (~1/3 of the CV fluctuation), stride (~500 steps)
V(s) = 0
loop over MD steps:
    F_total = F_physical - ∂V/∂s · ∇s(r)          # bias force through the chain rule
    integrate
    every stride steps:
        s_now = s(r)
        w = w0 · exp(-V(s_now) / (kB ΔT))           # well-tempered (ΔT = (γ - 1) T); plain: w = w0
        V(s) += w · exp(-(s - s_now)² / (2σ²))
F(s) = -(γ / (γ - 1)) · V(s)   (well-tempered)   or   -V(s) (plain, once the wells are filled)
check: F(s) vs time (converged?), many transitions, reweighting for other observables</pre>
<ol class="steps">
  <li>In practice use <b>PLUMED</b>, which plugs into LAMMPS, GROMACS, CP2K, OpenMM, AMBER, NAMD and others. A two-line input defines a CV (<code>d: DISTANCE ATOMS=1,2</code>) and the bias (<code>METAD ARG=d SIGMA=0.05 HEIGHT=1.0 PACE=500 BIASFACTOR=10</code>).</li>
  <li>Use one or two CVs; the cost of filling grows quickly with dimension.</li>
  <li>Compute other averages by reweighting each frame with e<sup>V/kT</sup> (with the time-dependent correction for well-tempered runs), or run unbiased simulations in each basin.</li>
  <li>Multiple walkers sharing one bias speed up convergence on parallel machines.</li>
</ol>`,

  code: {
    note: 'Brownian dynamics in a double well with an 8 kT barrier. Part 1: plain dynamics never crosses (as the Kramers rate predicts). Part 2: metadynamics crosses over 100 times and −V gives the free energy. Part 3: well-tempered metadynamics converges as the run gets longer. Part 4: on a 2D surface, biasing the wrong coordinate gives no crossings at all. Run <code>pip install numpy</code>, then <code>python metadynamics.py</code> (about 6 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(9)
# Overdamped Langevin (Brownian) dynamics in units kT = 1, D = 1:  x += -U'(x) dt + sqrt(2 dt) * noise.
H = 8.0                                        # barrier height of the double well, in kT
U  = lambda x: H * (x**2 - 1) ** 2
dU = lambda x: 4 * H * x * (x**2 - 1)
GRID = np.linspace(-2.0, 2.0, 401); DG = GRID[1] - GRID[0]

def run(nsteps, dt=1e-3, meta=None, x=-1.0):
    """meta = dict(w=height, s=width, every=steps, dT=None (plain) or bias temperature for well-tempered).
    The bias is kept on a grid (value and slope), as PLUMED does with GRID_ keywords."""
    V = np.zeros_like(GRID); dV = np.zeros_like(GRID); traj = np.empty(nsteps); hills = 0
    for s in range(nsteps):
        k = min(max(int((x - GRID[0]) / DG), 0), len(GRID) - 2); f = min(max((x - GRID[k]) / DG, 0.0), 1.0)
        bias_force = -(dV[k] * (1 - f) + dV[k + 1] * f) if GRID[0] <= x <= GRID[-1] else 0.0   # no bias outside the grid
        x += (-dU(x) + bias_force) * dt + np.sqrt(2 * dt) * rng.normal()
        traj[s] = x
        if meta and s % meta["every"] == 0:                                # deposit a Gaussian hill at the current x
            Vx = V[k] * (1 - f) + V[k + 1] * f
            w = meta["w"] * (np.exp(-Vx / meta["dT"]) if meta.get("dT") else 1.0)   # well-tempered: smaller hills where bias is high
            g = w * np.exp(-(GRID - x) ** 2 / (2 * meta["s"] ** 2))
            V += g; dV += -g * (GRID - x) / meta["s"] ** 2; hills += 1
    return traj, V, hills

def crossings(traj):                           # count transitions between the wells (|x| > 0.8 on both sides)
    side = np.where(traj > 0.8, 1, np.where(traj < -0.8, -1, 0)); side = side[side != 0]
    return int(np.sum(side[1:] != side[:-1]))

def fe_error(F):
    """RMS deviation from the exact free energy (here = U, up to a constant) on -1.4 < x < 1.4, and the barrier F(0) - F(min)."""
    m = (GRID > -1.4) & (GRID < 1.4); ref = U(GRID)
    a = F[m] - F[m].min(); b = ref[m] - ref[m].min()
    return np.sqrt(np.mean((a - b) ** 2)), F[np.argmin(np.abs(GRID))] - F[m].min()

N = 300000
t0 = time.time(); plain, _, _ = run(N)
print(f"Part 1: plain dynamics in a double well with an {H:.0f} kT barrier ({N} steps, {time.time() - t0:.0f} s)")
print(f"   transitions between the wells: {crossings(plain)}   (Kramers estimate for this run: ~{N * 1e-3 * 4 * H / (np.pi * np.sqrt(2)) * np.exp(-H):.1f})")

print("\\nPart 2: metadynamics on x (hills: height 0.5 kT, width 0.1, every 250 steps)")
t0 = time.time(); traj, V, nh = run(N, meta=dict(w=0.5, s=0.1, every=250))
err, barrier = fe_error(-V)
print(f"   transitions: {crossings(traj)}   hills: {nh}   ({time.time() - t0:.0f} s)")
print(f"   free energy F = -V: RMS error {err:.2f} kT, barrier {barrier:.2f} kT (exact {H:.2f})")

print("\\nPart 3: well-tempered metadynamics (bias temperature dT = 7 kT, so bias factor gamma = 8)")
for nsteps in [100000, 300000, 600000]:
    traj, V, nh = run(nsteps, meta=dict(w=0.5, s=0.1, every=250, dT=7.0))
    F = -(1 + 7.0) / 7.0 * V                                              # F = -(T + dT)/dT * V
    err, barrier = fe_error(F)
    print(f"   {nsteps:6d} steps: transitions {crossings(traj):3d}, RMS error {err:.2f} kT, barrier {barrier:.2f} kT")
print("   hills shrink as the bias grows, so the estimate converges instead of fluctuating")

print("\\nPart 4: choosing the collective variable. 2D surface: the barrier is along y; x is an easy, irrelevant direction")
def run2d(nsteps, cv, dt=1e-3, w=0.5, s=0.1, every=250):
    p = np.array([0.0, -1.0]); V = np.zeros_like(GRID); dV = np.zeros_like(GRID); ys = np.empty(nsteps)
    for st in range(nsteps):
        q = p[cv]; k = min(max(int((q - GRID[0]) / DG), 0), len(GRID) - 2); f = min(max((q - GRID[k]) / DG, 0.0), 1.0)
        force = np.array([-2 * p[0], -dU(p[1])])                         # U = x^2 + H (y^2 - 1)^2
        if GRID[0] <= q <= GRID[-1]: force[cv] -= dV[k] * (1 - f) + dV[k + 1] * f
        p += force * dt + np.sqrt(2 * dt) * rng.normal(size=2); ys[st] = p[1]
        if st % every == 0 and GRID[0] <= q <= GRID[-1]:
            g = w * np.exp(-(GRID - q) ** 2 / (2 * s * s)); V += g; dV += -g * (GRID - q) / (s * s)
    return ys
for name, cv in [("bias on x (wrong CV)", 0), ("bias on y (right CV)", 1)]:
    print(f"   {name}: transitions over the barrier in y = {crossings(run2d(200000, cv))}")
`,
  },

  maths: () => `
<h4>1. Free energy along a collective variable</h4>
${tex(String.raw`F(s) = -k_BT\ln P(s), \qquad P(s) = \big\langle \delta\big(s - s(\mathbf r)\big)\big\rangle`, true)}
<h4>2. Metadynamics bias</h4>
${tex(String.raw`V(s, t) = \sum_{t' = \tau, 2\tau, \dots < t} w\,\exp\!\left(-\frac{\big(s - s(t')\big)^2}{2\sigma^2}\right), \qquad \lim_{t\to\infty} V(s, t) = -F(s) + C`, true)}
<h4>3. Well-tempered metadynamics</h4>
${tex(String.raw`w(t) = w_0\,e^{-V(s(t),\,t)/k_B\Delta T}, \qquad F(s) = -\frac{T + \Delta T}{\Delta T}V(s, t\to\infty) + C, \qquad \gamma = \frac{T + \Delta T}{T}`, true)}
<p>The biased system samples P(s)<sup>1/γ</sup>: barriers are lowered by a factor γ rather than removed.</p>
<h4>4. Kramers rate (why plain MD fails)</h4>
${tex(String.raw`k \approx \frac{\sqrt{U''(x_{\min})\,|U''(x_{\text{barrier}})|}}{2\pi\gamma_f}\, e^{-\Delta U/k_BT}`, true)}
<p>(overdamped limit, friction γ<sub>f</sub>). Each extra kT of barrier costs a factor e ≈ 2.7 in waiting time.</p>
<h4>5. Bias force on the atoms</h4>
${tex(String.raw`\mathbf F_i^{\text{bias}} = -\frac{\partial V}{\partial s}\,\frac{\partial s}{\partial \mathbf r_i}`, true)}`,

  quiz: [
    {
      q: 'After a long, converged plain metadynamics run, how is the free energy along the CV related to the bias V(s)?',
      options: ['F(s) = V(s)', 'F(s) ≈ −V(s) + constant', 'F(s) = 0', 'F(s) = −kT ln V(s)'],
      answer: 1,
      why: 'The hills fill the free-energy wells until the biased landscape F + V is flat, so the bias is the negative of the free energy (up to a constant). In well-tempered metadynamics there is an extra factor (T + ΔT)/ΔT.',
    },
    {
      q: 'You run metadynamics on a CV but see almost no transitions and the free-energy estimate keeps drifting. What is the most likely problem?',
      options: ['The hills are too small', 'The CV does not capture the slow motion of the real transition', 'The temperature is too low', 'Metadynamics cannot handle barriers'],
      answer: 1,
      why: 'If the slow degree of freedom is not in the CV, filling the CV does not lower the actual barrier. That gives hysteresis and drifting estimates. Rethink the CV (or add a second one).',
    },
    {
      q: 'What does the well-tempered bias factor γ control?',
      options: ['The width of the Gaussian hills', 'How strongly barriers are flattened: the biased distribution is P(s)^(1/γ)', 'The number of CVs', 'The MD time step'],
      answer: 1,
      why: 'With γ = (T + ΔT)/T the system samples P(s)^(1/γ), so barriers are reduced by a factor γ. Large γ approaches plain metadynamics; γ = 1 means no bias.',
    },
  ],
};
