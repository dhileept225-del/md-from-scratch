import { tex } from '../tex.js';

export default {
  id: 'lattice-start',
  summary: 'Before the first step, every atom needs a position. Drop them at random and some will sit on top of each other, giving huge forces that blow the simulation apart. Starting from a crystal lattice avoids this, and then the system melts into the liquid you want.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/lattices.mp4', poster: 'videos/lattices.jpg' },
  sim: 'lattice',
  learn: [
    'Why random starting positions are dangerous',
    'How to build simple cubic, BCC and FCC lattices',
    'How the box size follows from N and the density',
  ],
  simHint: 'Start with FCC and check the energy. Switch to Random: red atoms overlap and the energy explodes. Tick "Reject overlaps" and compare. Then raise the density until random insertion fails. Which lattice has the closest pair furthest apart at the same density?',

  idea: `
<p>Choosing random coordinates sounds easy, but at liquid density some atoms will land almost on top of each other. The Lennard-Jones repulsion grows as r<sup>−12</sup>, so one close pair gives an enormous force. On the first step those atoms are shot apart, the energy blows up, and the run crashes.</p>
<p>The standard fix is to start from a <b>crystal lattice</b>. Atoms on a regular grid are evenly spaced, so no two are too close.</p>
<ul>
  <li><b>Simple cubic (SC):</b> one atom per cubic cell. Easy, but loosely packed.</li>
  <li><b>Body-centred cubic (BCC):</b> 2 atoms per cell.</li>
  <li><b>Face-centred cubic (FCC):</b> 4 atoms per cell (8 corners × ⅛ + 6 faces × ½). This is the crystal LJ atoms actually form, and the closest packing. It is the usual choice, so N = 4n³: 32, 108, 256, 500, 864…</li>
</ul>
<p>The <b>density</b> fixes the box: for N atoms at number density ρ, the box length is L = (N/ρ)<sup>1/3</sup>. Then we give the atoms velocities (Lesson 2.3) and let them run. At a liquid temperature the crystal <b>melts</b> within a few picoseconds, and from then on it has no memory of its start.</p>
<p>For molecules such as water or ions in a solvent, place each molecule on a grid point with a <b>random orientation</b>, or use a packing tool such as Packmol. Then run an <b>energy minimisation</b> to remove any remaining close contacts before MD.</p>
<div class="callout"><b>Key idea.</b> The start does not have to look like the final liquid. It only has to be safe: no overlaps, and a sensible density. Equilibration does the rest.</div>`,

  algorithm: `
<p>Building an FCC start for N = 4n³ atoms at reduced density ρ:</p>
<ol class="steps">
  <li>Box length L = (N/ρ)<sup>1/3</sup>, and cell size a = L/n.</li>
  <li>Loop over cells (i, j, k) from 0 to n−1.</li>
  <li>In each cell, place the 4 basis atoms at (0,0,0), (½,½,0), (½,0,½) and (0,½,½), in units of a.</li>
  <li>Check: the closest pair distance should be a/√2, well above σ at liquid densities.</li>
</ol>
<pre class="pseudo">basis = [(0,0,0), (.5,.5,0), (.5,0,.5), (0,.5,.5)]
for i, j, k in cells:
    for b in basis:
        pos.append(((i, j, k) + b) * a)</pre>
<p><b>Random with rejection</b> (for mixtures and molecules): propose a random position and accept it only if it is at least d<sub>min</sub> from every atom already placed. This gets very slow, or fails, near liquid density, which is why packing tools and minimisation exist.</p>`,

  code: {
    note: 'Builds a 256-atom FCC lattice at liquid-argon density, then compares its closest pair and LJ energy with pure random placement and random placement with rejection. Run <code>pip install numpy</code>, then <code>python lattice.py</code>.',
    src: `import numpy as np

def fcc_lattice(n_cells, rho):
    """n_cells^3 fcc unit cells (4 atoms each) at reduced density rho."""
    N = 4 * n_cells**3
    L = (N / rho) ** (1 / 3)                     # box length from the density
    a = L / n_cells                              # unit-cell size
    basis = np.array([[0, 0, 0], [0.5, 0.5, 0], [0.5, 0, 0.5], [0, 0.5, 0.5]])
    cells = np.array([[i, j, k] for i in range(n_cells)
                      for j in range(n_cells) for k in range(n_cells)])
    pos = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 3) * a
    return pos, L

def check(pos, L, rc=2.5):
    """LJ energy per atom and the closest pair distance (periodic box)."""
    d = pos[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    inv6 = 1 / r2[r2 < rc**2] ** 3
    return 2 * (inv6**2 - inv6).sum() / len(pos), np.sqrt(r2.min())

def random_no_overlap(N, L, rng, dmin=0.8):
    """Random positions, rejecting any atom closer than dmin to another."""
    pos = []
    while len(pos) < N:
        p = rng.random(3) * L
        if pos:
            d = np.array(pos) - p
            d -= L * np.round(d / L)
            if (d ** 2).sum(1).min() < dmin**2:
                continue
        pos.append(p)
    return np.array(pos)

rng = np.random.default_rng(0)
pos, L = fcc_lattice(4, 0.8442)                 # 256 atoms, liquid-argon density
print(f"N = {len(pos)}, box L = {L:.3f} sigma")
for name, p in [("fcc lattice", pos),
                ("pure random", rng.random(pos.shape) * L),
                ("random, d_min=0.8", random_no_overlap(len(pos), L, rng))]:
    E, dmin = check(p, L)
    print(f"{name:18s} closest pair {dmin:6.3f} sigma   energy/atom {E:12.3e}")`,
  },

  maths: () => `
<h4>1. Box size from density</h4>
${tex(String.raw`\rho = \frac{N}{L^3} \quad\Longrightarrow\quad L = \left(\frac{N}{\rho}\right)^{1/3}, \qquad a = \frac{L}{n}`, true)}
<h4>2. Atoms per cell and nearest-neighbour distance</h4>
${tex(String.raw`\begin{array}{lccc} & \text{SC} & \text{BCC} & \text{FCC} \\ \text{atoms per cell} & 1 & 2 & 4 \\ \text{nearest neighbour} & a & \tfrac{\sqrt3}{2}a & \tfrac{a}{\sqrt2} \\ \text{packing fraction} & 0.52 & 0.68 & 0.74 \end{array}`, true)}
<p>For FCC at ρ* = 0.8442 with N = 256: L = 6.72σ, a = 1.68σ, nearest neighbour 1.19σ, close to the LJ minimum 2<sup>1/6</sup>σ ≈ 1.12σ.</p>
<h4>3. Why one overlap is fatal</h4>
<p>At r = 0.2σ the LJ energy is</p>
${tex(String.raw`V(0.2\sigma) = 4\varepsilon\left[5^{12} - 5^{6}\right] \approx 10^{9}\,\varepsilon`, true)}
<p>That one pair alone carries billions of times the thermal energy (about 1 ε). The first step converts it into kinetic energy, and the system explodes.</p>
<h4>4. Random insertion has a limit</h4>
<p>Random sequential insertion of spheres of diameter d jams near a packing fraction of about 0.38. With d = 0.8σ, it becomes hopeless above ρ* ≈ 1.4. In practice it is already very slow around liquid density.</p>`,

  quiz: [
    {
      q: 'Why not just place atoms at random positions in the box?',
      options: ['Random numbers are too slow', 'Some atoms land almost on top of each other, causing huge forces that blow up the run', 'Atoms must always be in a crystal', 'Random positions break periodic boundaries'],
      answer: 1,
      why: 'At liquid density, random placement almost always creates a few very close pairs. The r⁻¹² repulsion makes their energy astronomically large.',
    },
    {
      q: 'How many atoms does an FCC lattice of 5 × 5 × 5 unit cells contain?',
      options: ['125', '250', '500', '1000'],
      answer: 2,
      why: 'FCC has 4 atoms per cell: 4 × 5³ = 500. That is why LJ simulations often use 32, 108, 256, 500 or 864 atoms.',
    },
    {
      q: 'You want N = 500 atoms at ρ* = 0.8. What box length do you need?',
      options: ['L ≈ 5.0σ', 'L ≈ 8.5σ', 'L ≈ 62.5σ', 'L ≈ 400σ'],
      answer: 1,
      why: 'L = (N/ρ)^(1/3) = (500/0.8)^(1/3) = 625^(1/3) ≈ 8.55σ.',
    },
  ],
};
