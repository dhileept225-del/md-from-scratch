import { tex } from '../tex.js';

export default {
  id: 'neighbour-lists',
  summary: 'Even with a cutoff, checking every pair still costs N². Neighbour lists remember who is close, and cell lists find them in linear time. Together they make large simulations possible.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/neighbour-lists.mp4', poster: 'videos/neighbour-lists.jpg' },
  sim: 'neighbours',
  learn: [
    'Why a cutoff alone still leaves N² distance tests',
    'Verlet lists: a skin, and the skin/2 rule for rebuilding',
    'Cell lists: finding neighbours in O(N) time, and combining the two',
  ],
  simHint: 'Start with All pairs and raise N: the tests (and time) grow as N². Switch to Cell list: only the 9 shaded cells are searched. Then pick Verlet list and move the skin slider. A bigger skin means fewer rebuilds but a longer list. Set the skin to 0: the list is rebuilt every step and saves nothing. Finally tick “Ignore the skin/2 rule” with a small skin and press Check: why are pairs now missed?',

  idea: `
<p>A cutoff (Lesson 2.5) tells us which pairs to <i>compute</i>, but to know which pairs are inside r<sub>c</sub> we still measure all N(N−1)/2 distances. At N = 100,000 that is 5 × 10<sup>9</sup> tests per step, almost all of them for atoms far apart.</p>
<p><b>Verlet list.</b> Atoms move only about 0.01σ per step. So once in a while, list every pair closer than r<sub>c</sub> + <b>skin</b>, and for the next steps test only those pairs. The list is still correct as long as no atom has moved more than skin/2 since it was built. (Two atoms moving towards each other by skin/2 each close the gap of one skin.) When any atom passes skin/2, rebuild.</p>
<p><b>Cell list.</b> Divide the box into cells at least r<sub>c</sub> wide. A neighbour of an atom must then be in its own cell or one of the adjacent cells: 9 cells in 2D, 27 in 3D. Each cell holds a fixed number of atoms at fixed density, so the search costs O(N) instead of O(N²).</p>
<p><b>Both together.</b> Real codes build the Verlet list <i>with</i> a cell list (cells of size r<sub>c</sub> + skin), then reuse the list for 10–20 steps. LAMMPS: <code>neighbor 0.3 bin</code> and <code>neigh_modify every 1 delay 0 check yes</code>. GROMACS: <code>cutoff-scheme = Verlet</code>, with <code>nstlist</code> and an automatic buffer (<code>rlist</code>).</p>
<div class="callout"><b>Key idea.</b> Build a list of pairs within r<sub>c</sub> + skin using cells (O(N)), reuse it until some atom moves skin/2, then rebuild. A skin of about 0.3σ (0.1–0.2 nm in force fields) is typical.</div>`,

  algorithm: `
<p>Inside the MD loop:</p>
<pre class="pseudo">every step:
    if list is empty or max |pos - pos_at_build| > skin / 2:
        put atoms into cells of size ≥ rc + skin
        list = []
        for each cell, for each of its 27 neighbour cells:
            for each pair i < j in these two cells:
                if |min_image(pos[i] - pos[j])| < rc + skin:
                    list.append((i, j))
        pos_at_build = pos
    for (i, j) in list:                  # forces only from the list
        if r2(i, j) < rc2: add LJ force</pre>
<ol class="steps">
  <li>Cells need at least 3 per side, so the box must be at least 3(r<sub>c</sub> + skin) wide. For small boxes, fall back to all pairs.</li>
  <li>Measure displacements with the minimum image, since atoms wrap through the walls.</li>
  <li>Too small a skin: rebuilds every few steps. Too large: a long list full of pairs outside r<sub>c</sub>. Test a few values and pick the fastest.</li>
  <li>Always check: forces with the list must match an all-pairs calculation. A skin that is too small, or skipping the check, silently loses pairs.</li>
</ol>`,

  code: {
    note: 'Cell list against brute force (same pairs, far fewer tests), how the work scales with N, and a Verlet list that rebuilds itself with the skin/2 rule. Run <code>pip install numpy</code>, then <code>python neighbours.py</code> (about 2 s).',
    src: `
import numpy as np
from collections import defaultdict

rng = np.random.default_rng(0)
rho, rc, skin = 0.8, 2.5, 0.3

def brute_force(pos, L, rcut):
    """Check every pair: N(N-1)/2 distance tests."""
    d = pos[:, None] - pos[None]
    d -= L * np.round(d / L)
    i, j = np.where(np.triu((d ** 2).sum(-1) < rcut**2, 1))
    return set(zip(i.tolist(), j.tolist()))

def cell_list(pos, L, rcut):
    """Sort atoms into cells of size >= rcut; test only atoms in the 27 neighbouring cells."""
    nc = int(L // rcut)                                   # cells per side (needs nc >= 3)
    cell_of = np.floor(pos / (L / nc)).astype(int) % nc
    cells = defaultdict(list)
    for a, c in enumerate(map(tuple, cell_of)):
        cells[c].append(a)
    cells = {c: np.array(m) for c, m in cells.items()}
    pairs, tests = [], 0
    for (x, y, z), A in cells.items():
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    B = cells.get(((x + dx) % nc, (y + dy) % nc, (z + dz) % nc))
                    if B is None:
                        continue
                    d = pos[A][:, None] - pos[B][None]
                    d -= L * np.round(d / L)
                    keep = A[:, None] < B[None]                # count each pair once
                    tests += keep.sum()
                    hit = keep & ((d ** 2).sum(-1) < rcut**2)
                    ia, ib = np.nonzero(hit)
                    pairs += list(zip(A[ia].tolist(), B[ib].tolist()))
    return set(pairs), int(tests)

# 1. same answer, fewer tests
N = 1000; L = (N / rho) ** (1 / 3); pos = rng.random((N, 3)) * L
pb = brute_force(pos, L, rc); pc, tc = cell_list(pos, L, rc)
print(f"N={N}: brute force and cell list find the same {len(pb)} pairs: {pb == pc}")

# 2. how the work grows with N
print("\\n     N   brute-force tests   cell-list tests   ratio")
for N in [1000, 4000, 16000]:
    L = (N / rho) ** (1 / 3); pos = rng.random((N, 3)) * L
    _, tc = cell_list(pos, L, rc); tb = N * (N - 1) // 2
    print(f"{N:6d}   {tb:17,d}   {tc:15,d}   {tb / tc:5.1f}x")

# 3. Verlet list: build with rc + skin, reuse until any atom has moved skin/2
N = 4000; L = (N / rho) ** (1 / 3); pos = rng.random((N, 3)) * L
vlist, _ = cell_list(pos, L, rc + skin); ref, rebuilds, steps = pos.copy(), 1, 200
for step in range(steps):
    pos = (pos + rng.normal(0, 0.005, pos.shape)) % L    # small moves, like one MD step
    disp = pos - ref; disp -= L * np.round(disp / L)
    if np.sqrt((disp ** 2).sum(1)).max() > skin / 2:
        vlist, _ = cell_list(pos, L, rc + skin); ref = pos.copy(); rebuilds += 1
print(f"\\nVerlet list (N={N}): {len(vlist):,} pairs per step, rebuilt {rebuilds} times in {steps} steps")
`,
  },

  maths: () => `
<h4>1. Cost per step</h4>
${tex(String.raw`\text{all pairs: } \tfrac12 N(N-1) \sim N^2, \qquad \text{cells: } \tfrac12 N \cdot 27\,\rho\,r_{\text{cell}}^3 \sim N`, true)}
<p>With r<sub>cell</sub> = r<sub>c</sub>, each atom tests 27ρr<sub>c</sub>³ candidates against 4/3πρr<sub>c</sub>³ ≈ 4.2ρr<sub>c</sub>³ true neighbours: about 16% of tests are hits.</p>
<h4>2. Verlet list size</h4>
${tex(String.raw`N_{\text{list}} = \tfrac12 N\,\rho\,\tfrac43\pi (r_c + r_s)^3, \qquad \frac{N_{\text{list}}}{N_{\text{pairs}}} = \left(1 + \frac{r_s}{r_c}\right)^3`, true)}
<p>For r<sub>c</sub> = 2.5σ and skin r<sub>s</sub> = 0.3σ the list is 1.12³ ≈ 1.4 times longer than the true pair count.</p>
<h4>3. Safe reuse</h4>
${tex(String.raw`|\Delta \mathbf r_i| + |\Delta \mathbf r_j| \le r_s \quad\Leftarrow\quad \max_i |\Delta \mathbf r_i| \le \frac{r_s}{2}`, true)}
<h4>4. Rebuild interval</h4>
${tex(String.raw`n_{\text{rebuild}} \approx \frac{r_s/2}{v_{\max}\,\Delta t}`, true)}
<p>With r<sub>s</sub> = 0.3σ, Δt = 0.005τ and v<sub>max</sub> ≈ 3σ/τ (the fastest atoms), about 10 steps. Atoms do not move in straight lines, so in practice it is often 10–20.</p>`,

  quiz: [
    {
      q: 'Why may a Verlet list be reused only until an atom has moved half the skin?',
      options: ['Because the cell size is half the skin', 'Two atoms can each move skin/2 towards each other, closing the whole skin', 'Because the force halves at the skin', 'It is just a convention with no reason'],
      answer: 1,
      why: 'A pair just outside rc + skin can enter rc only if the two atoms together move more than one skin. If no atom moves more than skin/2, no pair can slip in unseen.',
    },
    {
      q: 'Cells are r_c wide in 3D. How many cells must be searched for each atom?',
      options: ['1', '8', '9', '27'],
      answer: 3,
      why: 'Its own cell plus all neighbours: 3 × 3 × 3 = 27 (9 in 2D). Any atom within rc must be in one of them.',
    },
    {
      q: 'You make the skin very large. What happens?',
      options: ['Faster in every way', 'Fewer rebuilds, but each step tests many pairs outside rc, so it gets slower', 'Pairs get missed', 'The energy is no longer conserved'],
      answer: 1,
      why: 'The list grows as (rc + skin)³. Rebuilds become rare, but every step wastes time on far-away pairs. The best skin balances the two (about 0.3σ).',
    },
  ],
};
