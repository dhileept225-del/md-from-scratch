import { tex } from '../tex.js';

export default {
  id: 'van-hove',
  summary: 'The Van Hove function G(r, t) extends g(r) to time: how likely is it to find an atom at distance r at time t from where an atom was at time 0? Its self part shows how individual atoms move (smoothly, or by hopping), its distinct part shows how neighbour shells break up, and its Fourier transform is what neutron scattering measures.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/van-hove.mp4', poster: 'videos/van-hove.jpg' },
  sim: 'vanhove',
  learn: [
    'The self and distinct Van Hove functions and how to compute them from a trajectory',
    'The Gaussian approximation, the non-Gaussian parameter α₂ and the signature of hopping',
    'How neighbour shells relax (G_d) and the link to F_s(k, t) and neutron scattering',
  ],
  simHint: 'LJ liquid: move the lag time from 0.1 to 5 τ and watch G_s stay close to the dashed Gaussian, while the G_d peak melts from about 3 towards 1. Then switch to the Hopping model: at t ≈ 1 there are two clear peaks (stayed, hopped once) and α₂ is large. Make the waiting time longer, or the rattling larger, and see how the peaks change. "Moved less than 0.5" should follow exp(−t/wait) at short times.',

  idea: `
<p>The radial distribution function g(r) (Lesson 6.1) is a snapshot. Léon Van Hove (1954) added time: <b>G(r, t)</b> is the density of atoms at distance r and time t from the place where an atom was at time 0. It splits into two parts.</p>
<p><b>Self part G<sub>s</sub>(r, t)</b>: the same atom. 4πr²G<sub>s</sub>(r, t) dr is the probability that an atom has moved a distance between r and r + dr in time t. Its second moment is the MSD (Lesson 6.2), but its <i>shape</i> carries more information:</p>
<ul>
<li>For smooth diffusion it is a <b>Gaussian</b> whose width grows as √(6Dt). The simple LJ liquid stays close to this (code tab, Part 1).</li>
<li>The <b>non-Gaussian parameter</b> α₂ = 3⟨r⁴⟩/5⟨r²⟩² − 1 measures the deviation. It is 0 for a Gaussian; it peaks at intermediate times in liquids (about 0.13 here) and becomes large in supercooled liquids, glasses, polymers and ionic liquids, where some atoms are fast and others stuck (<b>dynamic heterogeneity</b>).</li>
<li>For <b>hopping</b> (jump diffusion), G<sub>s</sub> has separate peaks: one near r ≈ 0 for atoms still rattling in their site, and one at the jump distance a for atoms that have hopped (Part 3). For Li⁺ in electrolytes, a peak at the Li–Li or site–site distance is the classic evidence for hopping rather than vehicular transport.</li>
</ul>
<p><b>Distinct part G<sub>d</sub>(r, t)</b>: other atoms. At t = 0, G<sub>d</sub>(r, 0)/ρ = g(r). As time goes on, the first-shell peak decays (the shell breaks up) and the hole near r = 0 fills as neighbours move into the space the atom has left (Part 2). Its decay time is a measure of how long a solvation or coordination shell survives. For ion pairs (Li⁺–anion), G<sub>d</sub> tells you how long a pair stays together.</p>
<p><b>Scattering.</b> The Fourier transform of G<sub>s</sub> is the self-intermediate scattering function F<sub>s</sub>(k, t), measured by incoherent quasi-elastic neutron scattering (QENS). Its decay at the first peak of S(k) defines the structural relaxation time τ<sub>α</sub> (Part 4). The full G(r, t) transforms into the coherent F(k, t).</p>
<div class="callout"><b>Key idea.</b> G<sub>s</sub>(r, t): how far single atoms have moved (Gaussian = smooth diffusion; extra peaks = hopping; α₂ = heterogeneity). G<sub>d</sub>(r, t): how the neighbour shells relax. Their Fourier transforms connect MD to neutron scattering.</div>`,

  algorithm: `
<pre class="pseudo">save UNWRAPPED positions every n steps
self part, lag k:
    for each origin t0 and atom i:  r = |u_i(t0 + k) - u_i(t0)|;  histogram r
    4πr² G_s(r) = hist / (n_samples · dr)               # normalized: integrates to 1
    MSD = <r²>;  α2 = 3<r⁴> / (5 <r²>²) - 1
distinct part, lag k (fewer origins: costs N² per origin):
    for each origin t0, each i, each j ≠ i:
        d = minimum image of (x_j(t0 + k) - x_i(t0));  histogram |d|
    G_d(r)/ρ = hist / (n_origins · N · ρ · 4/3 π ((r+dr)³ - r³))
F_s(k, t) = < sin(k r) / (k r) >  over the same displacements (isotropic average)</pre>
<ol class="steps">
  <li>For mixtures, compute G_s per species and G_d per pair (Li–Li, Li–O, Li–anion). For distinct pairs A–B divide by N_A ρ_B.</li>
  <li>G_d needs positions consistently wrapped at t0 + t: use the wrapped position at t0 plus the unwrapped displacement, then the minimum image.</li>
  <li>Plot 4πr²G_s on a log scale of r (or log G_s vs r) to see small hopping populations.</li>
  <li>Tools: MDAnalysis (<code>van_hove</code>-style analyses), TRAVIS, LiquidLib, and your own NumPy code (as here); MDANSE computes F(k, t) for neutron comparisons.</li>
</ol>`,

  code: {
    note: 'Part 1 computes the self part and α₂ for the LJ liquid, Part 2 the distinct part, Part 3 a hopping model with two-peaked G_s, and Part 4 the self-intermediate scattering function F_s(k, t). Run <code>pip install numpy</code>, then <code>python van_hove.py</code> (about 30 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(4)

def fcc(n, rho):
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + b[None]).reshape(-1, 3)) * a, n * a

def forces(x, L, rc=2.5):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)
    r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
    ir6 = np.where(m, 1 / r2**3, 0.0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
    return (f[:, :, None] * d).sum(1)

# ---------- an LJ liquid trajectory (unwrapped positions every 0.05 tau) ----------
t0 = time.time(); rho, T, dt, every = 0.8442, 0.722, 0.005, 10
x, L = fcc(4, rho); N = len(x)
v = rng.normal(0, np.sqrt(T), x.shape); v -= v.mean(0); f = forces(x, L)
for s in range(1500):
    v += 0.5 * dt * f; x = (x + dt * v) % L; f = forces(x, L); v += 0.5 * dt * f
    v *= np.sqrt(T * (3 * N - 3) / (v**2).sum())
u = x.copy(); traj = []
for s in range(6000):
    v += 0.5 * dt * f; x = (x + dt * v) % L; u += dt * v; f = forces(x, L); v += 0.5 * dt * f
    if s % every == 0: traj.append(u.copy())
traj = np.array(traj); h = dt * every
print(f"LJ liquid: N = {N}, rho* = {rho}, T* = {T}, {len(traj)} frames over {6000 * dt:.0f} tau ({time.time() - t0:.0f} s)")

def self_vh(traj, lag, edges):
    """4 pi r^2 G_s(r, t): distribution of distances each atom moved in time t (all origins). Also <r^2>, <r^4>."""
    d = np.sqrt(((traj[lag:] - traj[:-lag])**2).sum(-1)).ravel()
    hist = np.histogram(d, edges, density=True)[0]
    return hist, (d**2).mean(), (d**4).mean()

# ---------- Part 1: self part G_s(r, t) and the non-Gaussian parameter ----------
print("\\nPart 1: self part 4 pi r^2 G_s(r, t) (where has an atom got to after time t?)")
edges = np.linspace(0, 4, 161); rc = 0.5 * (edges[1:] + edges[:-1])
for tt in [0.1, 0.5, 2.0, 5.0]:
    hist, r2, r4 = self_vh(traj, int(round(tt / h)), edges)
    a2 = 3 * r4 / (5 * r2**2) - 1                                  # non-Gaussian parameter
    print(f"   t = {tt:3.1f} tau: MSD = {r2:.3f};  peak at r = {rc[np.argmax(hist)]:.2f} (Gaussian with the same MSD: {np.sqrt(2 * r2 / 3):.2f});  alpha2 = {a2:.3f}")
print("   alpha2 = 3<r^4>/(5<r^2>^2) - 1 is 0 for a Gaussian; small here: atoms in a simple liquid move smoothly")

# ---------- Part 2: distinct part G_d(r, t) / rho: how the neighbour shell dissolves ----------
print("\\nPart 2: distinct part G_d(r, t) / rho (at t = 0 it is g(r))")
def distinct_vh(traj, lag, edges, n_origins=40):
    """Where are OTHER atoms at time t, relative to where atom i was at time 0? Normalized by the ideal gas."""
    hist = np.zeros(len(edges) - 1); orig = np.linspace(0, len(traj) - lag - 1, n_origins).astype(int)
    for o in orig:
        d = traj[o + lag][None, :, :] - traj[o][:, None, :]; d -= L * np.round(d / L)     # positions (unwrapped -> min image)
        r = np.sqrt((d**2).sum(-1)); np.fill_diagonal(r, np.inf)                         # exclude j = i (that is the self part)
        hist += np.histogram(r.ravel(), edges)[0]
    shell = 4 / 3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
    return hist / (n_origins * N * rho * shell)
edd = np.linspace(0, L / 2, 121); rd = 0.5 * (edd[1:] + edd[:-1])
for tt in [0.0, 0.5, 1.0, 3.0]:
    gd = distinct_vh(traj, int(round(tt / h)), edd)
    k1 = np.argmax(gd)
    print(f"   t = {tt:3.1f} tau: first peak G_d/rho = {gd[k1]:.2f} at r = {rd[k1]:.2f};  mean over r < 0.5: {gd[rd < 0.5].mean():.2f}")
print("   the shell peak melts away and the hole near r = 0 fills in as neighbours move into the place the atom left")

# ---------- Part 3: hopping (as for Li+ in some electrolytes): a model where G_s is far from Gaussian ----------
print("\\nPart 3: hopping model: particles rattle on sites (sd 0.1) and jump by a = 1.0 at random (mean wait 1.0)")
M, steps, hdt, a, wait = 4000, 400, 0.01, 1.0, 1.0
site = np.zeros((M, 3)); pos = []
for s in range(steps):
    jump = rng.random(M) < hdt / wait
    dirs = rng.normal(size=(M, 3)); dirs /= np.linalg.norm(dirs, axis=1, keepdims=True)
    site[jump] += a * dirs[jump]
    pos.append(site + rng.normal(0, 0.1, (M, 3)))
pos = np.array(pos); eh = np.linspace(0, 3, 121); rh = 0.5 * (eh[1:] + eh[:-1])
for tt in [0.5, 1.0, 3.0]:
    lag = int(round(tt / hdt)); d = np.sqrt(((pos[lag:] - pos[:-lag])**2).sum(-1)).ravel()
    hist = np.histogram(d, eh, density=True)[0]; r2, r4 = (d**2).mean(), (d**4).mean()
    stay, hop = hist[rh < 0.5].sum() * (rh[1] - rh[0]), hist[(rh > 0.7) & (rh < 1.3)].sum() * (rh[1] - rh[0])
    print(f"   t = {tt:3.1f}: peaks of 4 pi r^2 G_s at r = {rh[rh < 0.5][np.argmax(hist[rh < 0.5])]:.2f} (stayed, {stay:.2f} of particles) and "
          f"{rh[(rh > 0.7) & (rh < 1.3)][np.argmax(hist[(rh > 0.7) & (rh < 1.3)])]:.2f} (one hop, {hop:.2f});  alpha2 = {3 * r4 / (5 * r2**2) - 1:.2f}")
print(f"   exp(-t/wait) stayed: {np.exp(-0.5):.2f}, {np.exp(-1):.2f}, {np.exp(-3):.2f}.  Separate peaks at r ~ 0.2 (rattling) and r = a (hopped) are the fingerprint of jump diffusion")

# ---------- Part 4: self-intermediate scattering function F_s(k, t) (what neutrons measure) ----------
print("\\nPart 4: F_s(k, t) = <sin(k r)/(k r)> at k = 7.1 (first peak of S(k))")
k = 7.1; lags = np.arange(0, 121, 2); Fs = []
for lag in lags:
    d = np.sqrt(((traj[lag:] - traj[:len(traj) - lag])**2).sum(-1)).ravel() if lag else np.zeros(1)
    Fs.append(np.mean(np.sinc(k * d / np.pi)))                       # np.sinc(x) = sin(pi x)/(pi x)
Fs = np.array(Fs); ta = lags * h
tau_a = np.interp(1 / np.e, Fs[::-1], ta[::-1])
print(f"   F_s drops to 1/e at tau_alpha = {tau_a:.2f} tau (structural relaxation time);  F_s(t = 3 tau) = {Fs[ta >= 3][0]:.3f}")
`,
  },

  maths: () => `
<h4>1. Definition</h4>
${tex(String.raw`G(\mathbf r, t) = \frac{1}{N}\Big\langle \sum_{i}\sum_{j} \delta\big(\mathbf r + \mathbf r_i(0) - \mathbf r_j(t)\big)\Big\rangle = G_s(\mathbf r, t) + G_d(\mathbf r, t)`, true)}
<p>G<sub>s</sub> keeps the terms i = j, G<sub>d</sub> the terms i ≠ j. At t = 0: G<sub>s</sub> = δ(<b>r</b>) and G<sub>d</sub> = ρ g(r). At long times: G<sub>s</sub> → 0 (spread over the box) and G<sub>d</sub> → ρ.</p>
<h4>2. Gaussian approximation and α₂</h4>
${tex(String.raw`G_s(r, t) \approx \left(\frac{3}{2\pi\langle r^2(t)\rangle}\right)^{3/2}\exp\!\left(-\frac{3r^2}{2\langle r^2(t)\rangle}\right), \qquad \alpha_2(t) = \frac{3\langle r^4\rangle}{5\langle r^2\rangle^2} - 1`, true)}
<h4>3. Moments</h4>
${tex(String.raw`\int 4\pi r^2 G_s\,dr = 1, \qquad \int 4\pi r^2\, r^2\, G_s\,dr = \text{MSD}(t) \;\to\; 6Dt`, true)}
<h4>4. Intermediate scattering functions</h4>
${tex(String.raw`F_s(k, t) = \int G_s(\mathbf r, t)\, e^{-i\mathbf k\cdot\mathbf r}\, d\mathbf r = \Big\langle \frac{\sin kr}{kr}\Big\rangle \;\approx\; e^{-k^2\langle r^2\rangle/6}`, true)}
<p>The last form holds in the Gaussian approximation. Hopping gives a k-dependent relaxation rate, as in the Chudley–Elliott jump model.</p>
<h4>5. Chudley–Elliott jump diffusion (isotropic jumps of length a, mean waiting time τ₀)</h4>
${tex(String.raw`F_s(k, t) = \exp\!\left[-\frac{t}{\tau_0}\Big(1 - \frac{\sin ka}{ka}\Big)\right], \qquad D = \frac{a^2}{6\tau_0}`, true)}`,

  quiz: [
    {
      q: 'The self Van Hove function 4πr²G_s(r, t) of Li⁺ shows two peaks: one near 0.3 Å and one near 3 Å. What does this suggest?',
      options: ['Li⁺ diffuses smoothly like a simple liquid', 'Li⁺ moves by hopping between sites about 3 Å apart', 'The simulation is not equilibrated', 'Li⁺ does not move at all'],
      answer: 1,
      why: 'Separate peaks mean that ions are either still rattling near their start (small r) or have jumped by a fixed distance. That is the signature of jump (hopping) diffusion; smooth diffusion gives a single Gaussian peak.',
    },
    {
      q: 'What is G_d(r, t)/ρ at t = 0?',
      options: ['Zero everywhere', 'The radial distribution function g(r)', 'A Gaussian', 'The MSD'],
      answer: 1,
      why: 'At t = 0 the distinct part counts the other atoms around each atom at the same instant, which is ρ g(r). Over time it relaxes to ρ (uniform).',
    },
    {
      q: 'The non-Gaussian parameter α₂(t) is large at intermediate times. What does that tell you?',
      options: ['Diffusion is perfectly Gaussian', 'Dynamics are heterogeneous: some particles move much further than others in the same time', 'The temperature is too high', 'The box is too small'],
      answer: 1,
      why: 'α₂ = 0 for a Gaussian. Large α₂ means a broad, long-tailed displacement distribution: a mixture of fast and slow (or hopping and waiting) particles.',
    },
  ],
};
