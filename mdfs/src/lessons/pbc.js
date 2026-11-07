import { tex } from '../tex.js';

export default {
  id: 'pbc',
  summary: 'A simulation box holds only thousands of atoms, so a large fraction would sit at the walls. Periodic boundary conditions remove the walls: the box is surrounded by copies of itself, and an atom leaving one side comes back in through the other.',
  minutes: 13,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/pbc.mp4', poster: 'videos/pbc.jpg' },
  sim: 'pbc',
  learn: [
    'Why small simulation boxes need periodic boundary conditions',
    'How to wrap positions and use the minimum-image distance',
    'Two rules: keep the cutoff below L/2, and keep unwrapped positions for diffusion',
  ],
  simHint: 'Pick an atom near an edge. With "Minimum image" on, its neighbours include atoms on the other side of the box. Turn it off and count again. Push the cutoff past 0.5 L and watch the warning. The unwrapped path can leave the box and end on one of the copies: that is where the atom really is. Then switch the path to wrapped: the trail breaks into jumps, and "distance travelled" becomes wrong.',

  idea: `
<p>Even a million atoms is a very small piece of matter. In a cube of 1000 atoms, about half of them touch a wall or surface. Walls change everything: density, structure and dynamics would all be wrong for a bulk liquid.</p>
<p><b>Periodic boundary conditions (PBC)</b> remove the walls. Imagine the box copied in every direction to fill all space. Every copy moves exactly like the original, so we only simulate one box, but its atoms feel neighbours on all sides.</p>
<ul>
  <li><b>Wrapping:</b> an atom that leaves through the right face re-enters through the left face. The number of atoms stays constant.</li>
  <li><b>Minimum image:</b> when computing forces, each atom interacts with the <b>nearest copy</b> of every other atom. Two atoms at opposite edges are really neighbours.</li>
  <li><b>Cutoff rule:</b> the cutoff must be at most <b>L/2</b>. Otherwise an atom could interact with two copies of the same neighbour, or with its own image.</li>
  <li><b>Unwrapped coordinates:</b> for diffusion (Lesson 6.2) you need the true distance travelled. Store image counters, or the "unwrapped" positions, alongside the wrapped ones.</li>
</ul>
<p>PBC is not perfect. The system is still artificially periodic, so very long-range effects such as large-scale fluctuations, slow hydrodynamic flows and some electrostatics feel the box size. Always check important results with a bigger box.</p>
<div class="callout"><b>Key idea.</b> PBC turns a small box into an infinite, wall-free material. Wrap positions, use the nearest copy for distances, keep r<sub>c</sub> ≤ L/2, and save unwrapped positions.</div>`,

  algorithm: `
<p>Two one-line functions do almost all the work in a cubic box of side L:</p>
<pre class="pseudo">def wrap(x):                 # keep positions in [0, L)
    return x - L * floor(x / L)

def min_image(d):            # shortest vector to the nearest copy
    return d - L * round(d / L)</pre>
<ol class="steps">
  <li><b>After each move</b>, wrap the positions. Optionally count crossings in an integer "image" array, so the unwrapped position is wrap(x) + image·L.</li>
  <li><b>In the force loop</b>, compute d = r<sub>i</sub> − r<sub>j</sub> and apply min_image before taking the length.</li>
  <li><b>Check</b> that r<sub>c</sub> ≤ L/2 before running. LAMMPS and GROMACS stop with an error if it is not.</li>
  <li><b>For analysis</b>, use unwrapped coordinates for diffusion, and minimum-image distances for structure (RDF).</li>
</ol>
<p class="muted">Non-cubic (triclinic) boxes need a slightly more general version, using the box matrix. The idea is identical.</p>`,

  code: {
    note: 'The wrap and minimum-image functions, a check that two atoms on opposite edges are really 1.0 apart, and a demonstration that wrapped coordinates give a wrong mean squared displacement. Run <code>pip install numpy</code>, then <code>python pbc.py</code>.',
    src: `import numpy as np

L = 10.0                                   # box length

def wrap(x):
    """Put positions back inside the box [0, L)."""
    return x - L * np.floor(x / L)

def min_image(d):
    """Shortest vector between two atoms, using the nearest periodic copy."""
    return d - L * np.round(d / L)

a = np.array([0.5, 5.0, 5.0])
b = np.array([9.5, 5.0, 5.0])
print(f"naive distance {np.linalg.norm(a - b):.2f}   minimum-image distance {np.linalg.norm(min_image(a - b)):.2f}")

# random walkers: wrapped positions hide how far atoms really travelled
rng = np.random.default_rng(0)
N, steps, s = 500, 10000, 0.1
x = rng.random((N, 3)) * L                 # unwrapped positions
x0 = x.copy()
for _ in range(steps):
    x += rng.normal(0, s, (N, 3))
xw = wrap(x)                               # what a wrapped trajectory file stores
print(f"true MSD {((x - x0) ** 2).sum(1).mean():6.1f}   "
      f"MSD from wrapped coords {((xw - x0) ** 2).sum(1).mean():6.1f}   "
      f"expected {3 * steps * s**2:6.1f}")

# image counters recover the unwrapped positions from the wrapped ones
image = np.floor(x / L).astype(int)
print("unwrapped = wrapped + image*L ?", np.allclose(xw + image * L, x))`,
  },

  maths: () => `
<h4>1. Periodic copies</h4>
<p>The infinite system is the box plus all its translations by integer multiples of L:</p>
${tex(String.raw`\mathbf r_i^{(\mathbf n)} = \mathbf r_i + \mathbf n L, \qquad \mathbf n \in \mathbb Z^3`, true)}
<h4>2. Minimum image</h4>
${tex(String.raw`\mathbf d_{ij} = \mathbf r_i - \mathbf r_j - L\,\operatorname{round}\!\left(\frac{\mathbf r_i - \mathbf r_j}{L}\right), \qquad |d_{ij,\alpha}| \le \frac L2`, true)}
<p>Each component of the minimum-image vector lies between −L/2 and L/2. So a cutoff r<sub>c</sub> ≤ L/2 guarantees each pair is counted exactly once.</p>
<h4>3. Surface fraction without PBC</h4>
<p>For a cube of N atoms, about N<sup>1/3</sup> per side, the fraction in the outer layer is roughly</p>
${tex(String.raw`f_{\text{surface}} \approx 1 - \left(1 - \frac{2}{N^{1/3}}\right)^3 \qquad (N = 1000:\ \approx 49\%,\quad N = 10^6:\ \approx 6\%)`, true)}
<h4>4. Unwrapped coordinates</h4>
${tex(String.raw`\mathbf r_i^{\text{unwrapped}} = \mathbf r_i^{\text{wrapped}} + \mathbf n_i L`, true)}
<p>The mean squared displacement must use unwrapped positions. Computed from wrapped ones, it saturates near L²/4 per dimension (L²/2 in total for a uniform spread), which is the plateau the code tab shows.</p>`,

  quiz: [
    {
      q: 'Atom A is at x = 0.3 and atom B at x = 9.8 in a box of length L = 10. What is their minimum-image distance along x?',
      options: ['9.5', '0.5', '10.1', '4.75'],
      answer: 1,
      why: 'd = 0.3 − 9.8 = −9.5; −9.5 − 10·round(−0.95) = −9.5 + 10 = 0.5. Through the periodic boundary they are close neighbours.',
    },
    {
      q: 'Why must the cutoff be at most half the box length?',
      options: ['To save memory', 'Otherwise an atom can interact with two copies of the same neighbour, or with its own image', 'Because forces are zero beyond L/2', 'It does not matter'],
      answer: 1,
      why: 'Beyond L/2, the cutoff sphere reaches past the nearest copy, so some pairs would be counted twice and atoms could see themselves.',
    },
    {
      q: 'You compute diffusion from a trajectory that stores wrapped positions. What goes wrong?',
      options: ['Nothing', 'The mean squared displacement saturates, because jumps across the box are hidden, so D comes out too small', 'D comes out too large', 'The simulation crashes'],
      answer: 1,
      why: 'Wrapped positions never move further than the box size, so the mean squared displacement levels off instead of growing linearly. Use unwrapped coordinates or image counters.',
    },
  ],
};
