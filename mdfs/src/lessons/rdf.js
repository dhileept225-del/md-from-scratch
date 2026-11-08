import { tex } from '../tex.js';

export default {
  id: 'rdf',
  summary: 'The radial distribution function g(r) tells you how likely you are to find a neighbour at distance r, compared with a random (ideal-gas) arrangement. It is the first thing to compute from any trajectory: it shows the structure, gives coordination numbers and can be compared directly with X-ray and neutron scattering.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/rdf.mp4', poster: 'videos/rdf.jpg' },
  sim: 'rdf',
  learn: [
    'What g(r) means and how to compute it from a trajectory (histogram + ideal-gas normalization)',
    'Reading g(r): excluded core, neighbour shells, gas vs liquid vs solid',
    'Coordination numbers, the potential of mean force, and the link to scattering',
  ],
  simHint: 'Watch the liquid g(r) build up from noise and approach the dashed curve from the code tab. Then press Solid, and slowly raise the temperature until the sharp peaks melt. Try Gas: the single broad peak needs many frames to become smooth. The cyan atoms in the slice are the tagged atom\'s first shell; their 3D count fluctuates around the coordination number.',

  idea: `
<p>Sit on one atom and ask: <i>how many neighbours are at distance r?</i> Count the atoms in a thin spherical shell between r and r + dr, and divide by the number you would expect if the atoms were scattered at random with the same density ρ (an ideal gas). That ratio is the <b>radial distribution function</b>:</p>
<p style="text-align:center">g(r) = (neighbours in the shell) / (ρ · 4πr² dr)</p>
<p>Average over every atom as the centre and over many frames, and g(r) becomes smooth.</p>
<p><b>Reading it</b> (code tab, Part 2; Lennard-Jones liquid near its triple point):</p>
<ul>
<li><b>g = 0 for r &lt; 0.9σ</b>: atoms cannot overlap.</li>
<li><b>A tall first peak</b> (g ≈ 3 at r ≈ 1.09σ): the first shell of neighbours, three times more likely than at random.</li>
<li><b>Smaller peaks</b> at about 2σ and 3σ, then <b>g → 1</b>: beyond a few molecular diameters a liquid has no memory of the central atom.</li>
</ul>
<p><b>Gas, liquid, solid.</b> A dilute gas has one broad peak, close to e<sup>−u(r)/kT</sup>. A liquid has a few decaying peaks. A crystal has sharp peaks at the lattice distances that never die out (Part 3).</p>
<p><b>Coordination number.</b> Integrating ρ g(r) 4πr² up to the first minimum counts the neighbours in the first shell: about 12 for the LJ liquid and exactly 12 for fcc. For battery electrolytes this is how you count how many anion oxygens or solvent molecules surround a Li⁺ ion, using the Li–O g(r).</p>
<p><b>Potential of mean force.</b> w(r) = −kT ln g(r) is the effective free energy between two atoms, averaged over everything else (Part 4). It differs from the bare pair potential because of the surrounding atoms.</p>
<p><b>Experiments.</b> X-ray and neutron scattering measure the structure factor S(k), which is a Fourier transform of g(r) − 1. That makes g(r) one of the most direct checks of a force field.</p>
<div class="callout"><b>Key idea.</b> g(r) = local density at distance r / average density. Zero inside the core, peaks for the neighbour shells, 1 far away. Integrate the first peak to get coordination numbers.</div>`,

  algorithm: `
<pre class="pseudo">choose dr (about 0.01-0.02 σ, or 0.002 nm); bins up to r_max ≤ L/2
h[:] = 0
for each saved frame:
    for each pair i < j:
        r = |minimum image of (x_i - x_j)|
        if r < r_max: h[floor(r / dr)] += 2        # the pair counts for both atoms
for each bin k (centre r_k):
    shell = 4/3 π ((r_k + dr/2)³ - (r_k - dr/2)³)
    g[k] = h[k] / (n_frames · N · ρ · shell)
n_1 = Σ_{r_k < r_min} 4π ρ r_k² g[k] dr               # coordination number</pre>
<ol class="steps">
  <li>Never go beyond half the box length: past L/2 the minimum image no longer sees every neighbour.</li>
  <li>For two species A and B, count only A–B pairs and divide by N<sub>A</sub> ρ<sub>B</sub>; the coordination number of B around A uses ρ<sub>B</sub>.</li>
  <li>Equilibrate first, and use uncorrelated frames (every 10–100 steps). The first peak converges fast; the tail and dilute species need many frames.</li>
  <li>Tools: GROMACS <code>gmx rdf</code>, LAMMPS <code>compute rdf</code> + <code>fix ave/time</code>, VMD, MDAnalysis <code>InterRDF</code>, and TRAVIS.</li>
</ol>`,

  code: {
    note: 'A complete small MD code (256 LJ atoms, velocity Verlet) plus an RDF routine. Part 1 checks the normalization on random points, Part 2 analyses the liquid, Part 3 compares gas, liquid and solid, and Part 4 computes the potential of mean force. Run <code>pip install numpy</code>, then <code>python rdf.py</code> (about 30 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(1)

def fcc(n, rho):
    """n x n x n fcc unit cells (4 atoms each) at number density rho (LJ units)."""
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + b[None]).reshape(-1, 3)) * a, n * a

def forces(x, L, rc=2.5):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)          # minimum image
    r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
    ir6 = np.where(m, 1 / r2**3, 0.0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
    return (f[:, :, None] * d).sum(1)

def rdf(frames, L, nbins=150):
    """g(r): histogram of pair distances / ideal-gas count in each shell, averaged over frames."""
    N = len(frames[0]); rmax = L / 2; edges = np.linspace(0, rmax, nbins + 1); h = np.zeros(nbins)
    iu = np.triu_indices(N, 1)
    for x in frames:
        d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)
        r = np.sqrt((d**2).sum(-1))[iu]
        h += np.histogram(r, edges)[0]
    rho = N / L**3; shell = 4 / 3 * np.pi * (edges[1:]**3 - edges[:-1]**3)
    g = 2 * h / (len(frames) * N * rho * shell)                    # factor 2: each pair counted once
    return 0.5 * (edges[1:] + edges[:-1]), g, rho

def run_md(x, L, T, nequil, nprod, dt=0.005, every=10):
    v = rng.normal(0, np.sqrt(T), x.shape); v -= v.mean(0); f = forces(x, L); frames = []
    for step in range(nequil + nprod):
        v += 0.5 * dt * f; x = (x + dt * v) % L; f = forces(x, L); v += 0.5 * dt * f
        if step < nequil: v *= np.sqrt(T / (v**2).sum() * (3 * len(x) - 3))   # crude rescaling during equilibration
        elif step % every == 0: frames.append(x.copy())
    return frames

def report(name, r, g, rho):
    i1 = np.argmax(g); j = i1 + np.argmin(g[i1:i1 + int(0.6 / (r[1] - r[0]))]); dr = r[1] - r[0]
    n1 = np.sum(4 * np.pi * rho * r[:j + 1]**2 * g[:j + 1]) * dr       # coordination number up to the first minimum
    print(f"   {name:22s} first peak r = {r[i1]:.2f}, g = {g[i1]:.2f};  first minimum r = {r[j]:.2f};  neighbours n1 = {n1:.1f}")
    return i1, j

# ---------- Part 1: an ideal gas has g(r) = 1 (checks the normalization) ----------
L = 6.0; frames = [rng.random((256, 3)) * L for _ in range(50)]
r, g, rho = rdf(frames, L, 60)
print(f"Part 1: random points (ideal gas): mean g(r) = {g[5:].mean():.3f} +- {g[5:].std():.3f}  (should be 1)")

# ---------- Part 2: LJ liquid near the triple point ----------
t0 = time.time(); x, L = fcc(4, 0.8442)
frames = run_md(x, L, 0.722, 1500, 3000)
r, g, rho = rdf(frames, L)
print(f"\\nPart 2: LJ liquid, N = 256, rho* = 0.8442, T* = 0.722 ({len(frames)} frames, {time.time() - t0:.0f} s)")
i1, j = report('liquid', r, g, rho)
gs = np.convolve(g, np.ones(5) / 5, mode="same")                    # light smoothing before peak picking
peaks = [k for k in range(3, len(g) - 3) if gs[k] == gs[k - 3:k + 4].max() and gs[k] > 1.02]
print("   peaks at r =", ", ".join(f"{r[k]:.2f}" for k in peaks[:3]), "(about one, two, three neighbour shells)")
print(f"   g(r) at r = 0.8: {g[np.searchsorted(r, 0.8)]:.3f}  (atoms cannot overlap); beyond the peaks g oscillates around 1 and decays to it")

# ---------- Part 3: same picture for a gas and a crystal ----------
print("\\nPart 3: gas, liquid and solid")
xg, Lg = fcc(4, 0.05)                                                # dilute lattice start (random starts overlap)
r_, g_, _ = rdf(run_md(xg, Lg, 2.0, 1000, 2000), Lg, 80); k = np.argmax(g_)
print(f"   {'gas (rho 0.05, T 2.0)':22s} first peak r = {r_[k]:.2f}, g = {g_[k]:.2f};  low-density theory exp(-u/kT) at the LJ minimum = {np.exp(1 / 2.0):.2f}")
report('liquid (rho 0.84, T 0.72)', r, g, rho)
xs, Ls = fcc(4, 1.1); fr = run_md(xs, Ls, 0.3, 500, 1500); report('fcc solid (rho 1.1, T 0.3)', *rdf(fr, Ls))

# ---------- Part 4: potential of mean force ----------
w = -0.722 * np.log(np.maximum(g, 1e-12))
print(f"\\nPart 4: potential of mean force w(r) = -kT ln g(r), liquid:")
print(f"   well at the first peak: w = {w[i1]:.2f} epsilon;  barrier at the first minimum: w = {w[j]:.2f} epsilon")
print(f"   (compare the bare LJ well depth, -1 epsilon at r = 1.12: the liquid around a pair changes the effective interaction)")
`,
  },

  maths: () => `
<h4>1. Definition</h4>
${tex(String.raw`g(r) = \frac{1}{\rho N}\Big\langle \sum_{i}\sum_{j\ne i} \delta(\mathbf r - \mathbf r_{ij}) \Big\rangle, \qquad \rho = \frac{N}{V}`, true)}
<p>Averaged over directions, the number of neighbours between r and r + dr is ρ g(r) 4πr² dr.</p>
<h4>2. Coordination number</h4>
${tex(String.raw`n(r) = 4\pi\rho\int_0^{r} g(r')\, r'^2\, dr', \qquad n_1 = n(r_{\min})`, true)}
<h4>3. Potential of mean force</h4>
${tex(String.raw`w(r) = -k_BT\,\ln g(r) \quad\xrightarrow{\;\rho\to 0\;}\quad u(r)`, true)}
<p>At low density g(r) → e<sup>−u(r)/kT</sup>, so the PMF becomes the pair potential (Part 3: g<sub>max</sub> = 1.61 vs e<sup>1/2</sup> = 1.65).</p>
<h4>4. Structure factor (what scattering measures)</h4>
${tex(String.raw`S(k) = 1 + 4\pi\rho\int_0^\infty \big[g(r) - 1\big]\, \frac{\sin kr}{kr}\, r^2\, dr`, true)}
<h4>5. Thermodynamics from g(r) (pairwise potentials)</h4>
${tex(String.raw`\frac{U}{N} = \tfrac32 k_BT + 2\pi\rho\int_0^\infty u(r)\,g(r)\,r^2\,dr, \qquad P = \rho k_BT - \frac{2\pi\rho^2}{3}\int_0^\infty r\,u'(r)\,g(r)\,r^2\,dr`, true)}`,

  quiz: [
    {
      q: 'Far from the central atom (several molecular diameters) in a liquid, g(r) approaches…',
      options: ['0', '1', 'the coordination number', 'infinity'],
      answer: 1,
      why: 'Far away the positions of the neighbours are uncorrelated with the central atom, so the local density equals the average density and g(r) = 1.',
    },
    {
      q: 'How do you get the number of nearest neighbours from g(r)?',
      options: ['Read the height of the first peak', 'Integrate 4πρ r² g(r) from 0 to the first minimum', 'Take the position of the first peak', 'Multiply the first peak height by ρ'],
      answer: 1,
      why: 'ρ g(r) 4πr² dr is the number of neighbours in a shell, so integrating it up to the first minimum counts the whole first shell (about 12 in the LJ liquid).',
    },
    {
      q: 'Your box is 4 nm long. Up to what distance can you safely compute g(r) with the minimum image convention?',
      options: ['8 nm', '4 nm', '2 nm', 'Any distance'],
      answer: 2,
      why: 'Beyond L/2 the minimum image no longer includes every neighbour at that distance, so g(r) would be wrong. Use r_max ≤ L/2 (or a bigger box).',
    },
  ],
};
