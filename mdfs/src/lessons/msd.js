import { tex } from '../tex.js';

export default {
  id: 'msd',
  summary: 'The mean squared displacement (MSD) measures how far atoms wander. At long times it grows linearly, MSD = 6Dt, and its slope gives the self-diffusion coefficient D, one of the most common numbers reported from MD and directly comparable with NMR experiments.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/msd.mp4', poster: 'videos/msd.jpg' },
  sim: 'msd',
  learn: [
    'What the MSD is and how to compute it with all time origins from unwrapped coordinates',
    'Ballistic (t²) and diffusive (t) regimes, and fitting D from the linear part',
    'Error bars, the wrapped-coordinates mistake and the finite-size (Yeh–Hummer) correction',
  ],
  simHint: 'Let the liquid run until the MSD is straight and D settles near 0.03. Switch to log–log: the curve starts on the dotted ballistic line (slope 2) and bends over to slope 1. Turn on "+ wrapped" to see how wrong wrapped coordinates are. Then compare Solid (D ≈ 0, the MSD flattens) and Gas: atoms fly so far between collisions that even at 5 τ the slope is not yet 1, and the widget warns you not to trust D. Finally, try heating the liquid.',

  idea: `
<p>Follow each atom from time 0 to time t and measure how far it has moved, r(t) − r(0). The average of the square of this displacement over all atoms is the <b>mean squared displacement</b>:</p>
<p style="text-align:center">MSD(t) = ⟨ |r(t) − r(0)|² ⟩</p>
<p>Because every moment in an equilibrium run is an equally good starting point, we also average over <b>all time origins</b>. This makes the MSD much smoother than using t = 0 only.</p>
<p><b>Two regimes</b> (code tab, Part 2):</p>
<ul>
<li><b>Short times (ballistic)</b>: before an atom hits its neighbours it flies freely, r = vt, so MSD = ⟨v²⟩t² = 3(kT/m)t². On a log–log plot the slope is 2.</li>
<li><b>Long times (diffusive)</b>: after many collisions the motion is a random walk and the MSD grows linearly, <b>MSD = 6Dt</b> in 3D. On a log–log plot the slope is 1.</li>
</ul>
<p>In between, a dense liquid shows a "cage" region where atoms rattle among their neighbours. In a solid the MSD stops growing (D = 0); in glasses and polymer electrolytes the cage can last nanoseconds.</p>
<p><b>Getting D.</b> Fit a straight line to the MSD in the region where the log–log slope is 1, and divide the slope by 6 (by 4 in 2D, by 2 in 1D). For the LJ liquid at ρ* = 0.8442, T* = 0.722 this gives D* ≈ 0.03 (for argon, about 1.6 × 10⁻⁵ cm²/s).</p>
<p><b>Three classic mistakes:</b></p>
<ul>
<li><b>Wrapped coordinates.</b> An atom that leaves the box re-enters on the other side, jumping by L. Always use unwrapped coordinates (Part 3).</li>
<li><b>Fitting too early or too late.</b> Too early includes the ballistic or cage region; too late uses lag times with few origins and large noise.</li>
<li><b>Ignoring the box size.</b> Hydrodynamic interactions with periodic images slow diffusion by about k<sub>B</sub>T ξ/(6πηL). That is about 15 % for this small box and still a few percent for typical boxes (Part 4).</li>
</ul>
<p>For ionic liquids and battery electrolytes, the MSD can stay sub-diffusive for many nanoseconds, so check the log–log slope before trusting a D. Self-diffusion is what pulsed-field-gradient NMR measures.</p>
<div class="callout"><b>Key idea.</b> MSD(t) = ⟨|r(t) − r(0)|²⟩, from unwrapped coordinates and averaged over atoms and time origins. Where the log–log slope is 1, D = slope / 6. Add the finite-size correction for small boxes.</div>`,

  algorithm: `
<pre class="pseudo">save unwrapped positions every n steps (or unwrap: add L each time an atom crosses the box)
for each lag k = 1 ... K:
    MSD[k] = mean over atoms i and origins t0 of |r_i(t0 + k) - r_i(t0)|²
check the log-log slope; pick t1..t2 where it is ≈ 1
D = slope of a linear fit of MSD(t) on [t1, t2] / 6
error: split the run into blocks (or use x, y, z separately) -> D ± std / sqrt(n_blocks)
finite size: D_inf = D(L) + 2.837297 k_B T / (6 π η L)   (cubic box; η = shear viscosity)</pre>
<ol class="steps">
  <li>Remove centre-of-mass drift (or check it is zero) before computing the MSD.</li>
  <li>The lag time should not exceed about a third of the run; long lags have few independent origins.</li>
  <li>For mixtures compute one MSD per species (Li⁺, anion, solvent). Cross terms and collective motion belong to conductivity (Lesson 6.5).</li>
  <li>Fast all-origin MSD uses FFTs (cost N_frames log N_frames), as in MDAnalysis <code>EinsteinMSD(fft=True)</code>. Other tools: GROMACS <code>gmx msd</code>, LAMMPS <code>compute msd</code> (single origin!), TRAVIS.</li>
</ol>`,

  code: {
    note: 'Part 1 checks the MSD routine on random walkers with known D. Part 2 runs 256 LJ atoms for 40 τ and fits D, with block error bars and the log–log slopes. Part 3 shows the wrapped-coordinates mistake and Part 4 the Yeh–Hummer correction. Run <code>pip install numpy</code>, then <code>python msd.py</code> (about 35 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(2)

def fcc(n, rho):
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + b[None]).reshape(-1, 3)) * a, n * a

def forces(x, L, rc=2.5):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)
    r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
    ir6 = np.where(m, 1 / r2**3, 0.0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
    return (f[:, :, None] * d).sum(1)

def msd(traj, max_lag):
    """MSD(t) averaged over atoms and over all time origins. traj: (frames, N, 3), UNWRAPPED."""
    out = np.zeros(max_lag + 1)
    for lag in range(1, max_lag + 1):
        d = traj[lag:] - traj[:-lag]
        out[lag] = (d**2).sum(-1).mean()
    return out

# ---------- Part 1: random walkers: MSD = 6 D t exactly ----------
D_true, dt = 0.5, 0.01
steps = rng.normal(0, np.sqrt(2 * D_true * dt), (2000, 500, 3))      # each Cartesian step has variance 2 D dt
walk = np.cumsum(steps, 0); m = msd(walk, 200); t = np.arange(201) * dt
print(f"Part 1: 500 random walkers with D = {D_true}: slope/6 of MSD = {np.polyfit(t[20:], m[20:], 1)[0] / 6:.3f}")

# ---------- Part 2: LJ liquid near the triple point ----------
t0 = time.time(); rho, T, dt, every = 0.8442, 0.722, 0.005, 10
x, L = fcc(4, rho); N = len(x)
v = rng.normal(0, np.sqrt(T), x.shape); v -= v.mean(0); f = forces(x, L)
for step in range(1500):                                             # equilibrate with velocity rescaling
    v += 0.5 * dt * f; x = (x + dt * v) % L; f = forces(x, L); v += 0.5 * dt * f
    v *= np.sqrt(T * (3 * N - 3) / (v**2).sum())
unwrapped = x.copy(); traj, wrapped = [], []
for step in range(8000):                                             # production in NVE
    v += 0.5 * dt * f
    unwrapped += dt * v; x = (x + dt * v) % L                         # keep an unwrapped copy for the MSD!
    f = forces(x, L); v += 0.5 * dt * f
    if step % every == 0: traj.append(unwrapped.copy()); wrapped.append(x.copy())
traj, wrapped = np.array(traj), np.array(wrapped); tau = np.arange(301) * dt * every
m = msd(traj, 300)
print(f"\\nPart 2: LJ liquid, N = {N}, rho* = {rho}, T* = {T}, {len(traj)} frames over {8000 * dt:.0f} tau ({time.time() - t0:.0f} s)")
for tt in [0.05, 0.1, 0.5, 1, 2, 3]:
    k = int(round(tt / (dt * every)))
    print(f"   t = {tt:4.2f}: MSD = {m[k]:.4f}  (ballistic 3T t^2 = {3 * T * tt**2:.4f})")
fit = slice(20, 61)                                                 # lags from tau = 1 to 3
slope = np.polyfit(tau[fit], m[fit], 1)[0]; D = slope / 6
print(f"   linear fit between 1 and 3 tau: D* = slope / 6 = {D:.4f}")
Dc = []                                                              # error bar: 4 independent blocks
for blk in np.array_split(traj, 4):
    mb = msd(blk, 60); Dc.append(np.polyfit(tau[fit], mb[fit], 1)[0] / 6)
print(f"   from 4 blocks: D* = {np.mean(Dc):.4f} +- {np.std(Dc, ddof=1) / 2:.4f} (standard error)")
print(f"   local slope d ln MSD / d ln t: {np.gradient(np.log(m[1:]), np.log(tau[1:]))[0]:.2f} at the start (2 = ballistic), "
      f"{np.gradient(np.log(m[1:]), np.log(tau[1:]))[-1]:.2f} at the end (1 = diffusive)")

# ---------- Part 3: the classic mistake: wrapped coordinates ----------
mw = msd(wrapped, 300)
print(f"\\nPart 3: MSD at t = 15 tau from unwrapped coordinates: {m[300]:.3f} (= 6 D t);  from wrapped coordinates: {mw[300]:.3f}")
print(f"   wrapped atoms that cross the box jump by L = {L:.2f}, so the wrapped MSD is garbage (here fraction of atoms that crossed: "
      f"{np.mean(np.any(np.abs(np.diff(wrapped, axis=0)) > L / 2, axis=(0, 2))):.2f})")

# ---------- Part 4: finite-size correction (Yeh and Hummer 2004) ----------
eta = 3.2                                                            # shear viscosity of the LJ liquid here (literature, approx.)
corr = 2.837297 * T / (6 * np.pi * eta * L)
print(f"\\nPart 4: periodic images slow diffusion by k_B T xi / (6 pi eta L) = {corr:.4f}  (eta* ~ {eta}, L = {L:.2f})")
print(f"   D(L) = {D:.4f}  ->  D(infinite) = {D + corr:.4f}  ({corr / D * 100:.0f} % correction for this small box)")
`,
  },

  maths: () => `
<h4>1. Einstein relation</h4>
${tex(String.raw`D = \lim_{t\to\infty}\frac{1}{6t}\Big\langle |\mathbf r(t) - \mathbf r(0)|^2 \Big\rangle \qquad (3\text{D}; \text{ use } 2d\,t \text{ in } d \text{ dimensions})`, true)}
<h4>2. Short-time limit</h4>
${tex(String.raw`\mathbf r(t) \approx \mathbf r(0) + \mathbf v t \;\Rightarrow\; \text{MSD} \approx \langle v^2\rangle t^2 = \frac{3k_BT}{m}\,t^2`, true)}
<h4>3. Link to the velocity autocorrelation (Green–Kubo, next lesson)</h4>
${tex(String.raw`\text{MSD}(t) = 2\int_0^t (t - s)\,\langle \mathbf v(s)\cdot\mathbf v(0)\rangle\,ds \;\Rightarrow\; D = \frac13\int_0^\infty \langle \mathbf v(t)\cdot\mathbf v(0)\rangle\,dt`, true)}
<h4>4. All-origins estimator</h4>
${tex(String.raw`\text{MSD}(k\Delta t) = \frac{1}{N(M-k)}\sum_{i=1}^{N}\sum_{t_0=0}^{M-k-1}\big|\mathbf r_i(t_0 + k) - \mathbf r_i(t_0)\big|^2`, true)}
<h4>5. Finite-size correction (Yeh and Hummer 2004)</h4>
${tex(String.raw`D_\infty = D(L) + \frac{\xi\,k_BT}{6\pi\eta L}, \qquad \xi = 2.837297 \text{ (cubic box)}`, true)}`,

  quiz: [
    {
      q: 'On a log–log plot of the MSD of a liquid, what slope do you expect at long times?',
      options: ['2', '1', '0', '0.5'],
      answer: 1,
      why: 'At long times MSD = 6Dt, a straight line of slope 1 on log–log axes. Slope 2 is the short-time ballistic regime; slope 0 means the atoms are trapped (a solid).',
    },
    {
      q: 'In 3D the MSD has a slope of 0.18 σ²/τ in its linear region. What is D?',
      options: ['0.18 σ²/τ', '0.06 σ²/τ', '0.03 σ²/τ', '0.09 σ²/τ'],
      answer: 2,
      why: 'MSD = 6Dt in three dimensions, so D = 0.18 / 6 = 0.03 σ²/τ.',
    },
    {
      q: 'Why must the MSD be computed from unwrapped coordinates?',
      options: ['Wrapped coordinates are less precise', 'An atom that crosses the box is moved by L, which looks like a huge jump', 'Unwrapped coordinates remove the thermostat', 'It is only a convention'],
      answer: 1,
      why: 'Wrapping puts atoms back into the box, so a small real step across the boundary appears as a jump of a whole box length, ruining the MSD.',
    },
  ],
};
