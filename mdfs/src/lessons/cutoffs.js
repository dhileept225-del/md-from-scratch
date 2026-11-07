import { tex } from '../tex.js';

export default {
  id: 'cutoffs',
  summary: 'Computing every pair of atoms costs N², which quickly becomes impossible. The Lennard-Jones interaction fades so fast that pairs beyond a cutoff can be ignored, if you handle the edge carefully and add back what you left out.',
  minutes: 13,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/cutoffs.mp4', poster: 'videos/cutoffs.jpg' },
  sim: 'cutoff',
  learn: [
    'Why a cutoff makes MD affordable, and how its cost grows as rc³',
    'Truncated, shifted, force-shifted and switched potentials',
    'Tail corrections, and why charges cannot simply be cut off',
  ],
  simHint: 'Start with Truncate and see the jump at rc in both V and F. Switch to Shift: the energy jump vanishes but the force still jumps. Force-shift and Switch make both smooth. Then pick Truncate, tick the tail correction and compare rc = 2.5σ with rc = 1.5σ. Why does the correction work so much better at 2.5σ?',

  idea: `
<p>N atoms have N(N−1)/2 pairs. For 100,000 atoms that is 5 × 10<sup>9</sup> pair calculations every step: far too many. Luckily the Lennard-Jones attraction falls off as r<sup>−6</sup>. At 2.5σ it is only −0.016ε, about 1.6% of the well depth. So we introduce a <b>cutoff</b> r<sub>c</sub> and simply ignore pairs beyond it.</p>
<p>Each atom then only interacts with neighbours inside a sphere of radius r<sub>c</sub>. That is about ρ·(4/3)πr<sub>c</sub>³ atoms, roughly 55 at 2.5σ in a liquid. The cost grows as <b>r<sub>c</sub>³</b>: doubling the cutoff makes the run eight times slower.</p>
<p>How you cut matters:</p>
<ul>
  <li><b>Truncate:</b> V = 0 beyond r<sub>c</sub>. The energy jumps when a pair crosses r<sub>c</sub>, which spoils energy conservation (Lesson 1.6).</li>
  <li><b>Shift:</b> subtract V(r<sub>c</sub>) so the energy is continuous. The force still jumps at r<sub>c</sub>.</li>
  <li><b>Force-shift:</b> also subtract the force at r<sub>c</sub>, so both energy and force go smoothly to zero. It changes the potential slightly everywhere.</li>
  <li><b>Switch:</b> leave V untouched up to r<sub>on</sub>, then smoothly bend it to zero by r<sub>c</sub>. Common in force fields (CHARMM, GROMACS).</li>
</ul>
<p>Pairs beyond r<sub>c</sub> still add a small attraction. The <b>tail correction</b> adds it back analytically, assuming the liquid is uniform beyond r<sub>c</sub>. It matters for energy and especially for <b>pressure</b>.</p>
<div class="callout"><b>Key idea.</b> Cut the Lennard-Jones interaction at about 2.5σ (or 1.0–1.2 nm in force fields), make the edge smooth, and add tail corrections. Coulomb forces decay only as 1/r, so they cannot be cut. They need Ewald summation (Chapter 5).</div>`,

  algorithm: `
<p>The force loop with a cutoff, using a shifted potential:</p>
<pre class="pseudo">rc2 = rc * rc
V_rc = 4 * (rc**-12 - rc**-6)
for each pair i < j:
    d = min_image(pos[i] - pos[j])
    r2 = dot(d, d)
    if r2 >= rc2: continue              # the cutoff
    inv6 = 1 / r2**3
    U += 4 * (inv6**2 - inv6) - V_rc    # shifted energy
    f = 24 * (2 * inv6**2 - inv6) / r2 * d
    F[i] += f;  F[j] -= f</pre>
<ol class="steps">
  <li>Compare <b>r²</b> with <b>r<sub>c</sub>²</b> and avoid the square root, which is slow.</li>
  <li>Keep r<sub>c</sub> ≤ L/2 (Lesson 2.4).</li>
  <li>Report energies with the tail correction, or state clearly which cutoff and scheme you used. Results are not comparable otherwise.</li>
  <li>The loop above still checks all N² pairs. Lesson 2.6 uses <b>neighbour lists</b> so each atom only checks atoms that are nearby.</li>
</ol>`,

  code: {
    note: 'A 500-atom LJ configuration: the energy per atom for several cutoffs, truncated, with tail correction, and shifted, compared with a near-complete reference. Run <code>pip install numpy</code>, then <code>python cutoffs.py</code>.',
    src: `import numpy as np

def fcc(n, rho):
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((c[:, None] + b[None]).reshape(-1, 3) * a), L

rng = np.random.default_rng(1)
pos, L = fcc(5, 0.8442)                       # 500 atoms, L = 8.39 sigma
pos += rng.normal(0, 0.08, pos.shape)         # shake it off the perfect lattice
N, rho = len(pos), len(pos) / L**3

d = pos[:, None] - pos[None]
d -= L * np.round(d / L)                      # minimum image
r = np.linalg.norm(d, axis=-1)[np.triu_indices(N, 1)]   # every pair once

lj = lambda r: 4 * (r**-12 - r**-6)
def tail(rc):
    """Energy per atom from all pairs beyond rc, assuming g(r) = 1 there."""
    return 8 / 3 * np.pi * rho * (rc**-9 / 3 - rc**-3)

ref = lj(r[r < 4.15]).sum() / N + tail(4.15)  # best estimate (rc close to L/2)
print(f"reference energy per atom {ref:.4f}\\n")
print("  rc   neighbours  truncated  +tail corr.  shifted")
for rc in [1.5, 2.0, 2.5, 3.0, 3.5]:
    m = r < rc
    trunc = lj(r[m]).sum() / N
    shifted = (lj(r[m]) - lj(rc)).sum() / N
    print(f"{rc:4.1f}  {2 * m.sum() / N:9.1f}  {trunc:9.4f}  {trunc + tail(rc):10.4f}  {shifted:8.4f}")`,
  },

  maths: () => `
<h4>1. Cost of a cutoff</h4>
${tex(String.raw`N_{\text{neigh}} \approx \rho\,\frac{4}{3}\pi r_c^3 \quad(\rho^* = 0.84,\ r_c = 2.5\sigma:\ \approx 55), \qquad \text{cost} \propto N\,r_c^3`, true)}
<h4>2. Shifted and force-shifted potentials</h4>
${tex(String.raw`V_{\text{sh}}(r) = V(r) - V(r_c), \qquad V_{\text{fs}}(r) = V(r) - V(r_c) - (r - r_c)\,V'(r_c) \qquad (r < r_c)`, true)}
<h4>3. Switching function (CHARMM form)</h4>
${tex(String.raw`S(r) = \frac{(r_c^2 - r^2)^2\,(r_c^2 + 2r^2 - 3r_{\text{on}}^2)}{(r_c^2 - r_{\text{on}}^2)^3}, \qquad r_{\text{on}} \le r \le r_c`, true)}
<p>S goes from 1 to 0 with zero slope at both ends, so V·S and its force stay smooth.</p>
<h4>4. Tail corrections (assume g(r) = 1 beyond r<sub>c</sub>)</h4>
${tex(String.raw`\frac{U_{\text{tail}}}{N} = \frac{8}{3}\pi\rho\varepsilon\sigma^3\left[\frac13\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3\right]`, true)}
${tex(String.raw`P_{\text{tail}} = \frac{16}{3}\pi\rho^2\varepsilon\sigma^3\left[\frac23\left(\frac{\sigma}{r_c}\right)^9 - \left(\frac{\sigma}{r_c}\right)^3\right]`, true)}
<p>At ρ* = 0.84 and r<sub>c</sub> = 2.5σ, U<sub>tail</sub>/N ≈ −0.45ε, the same size as the missing energy in the code tab. The correction assumes the liquid is structureless (g(r) = 1) beyond r<sub>c</sub>. That holds at 2.5σ but not at 1.5σ, where the first neighbour shells are still strongly structured, so a short cutoff cannot be rescued by the tail correction.</p>
<h4>5. Why Coulomb cannot be cut</h4>
<p>The number of neighbours in a shell grows as r², so the energy from shells at distance r scales as r²·V(r). For LJ this is r<sup>−4</sup> and converges fast. For charges it is r²·r<sup>−1</sup> = r, which grows. Truncating gives meaningless results, which is why Ewald summation exists.</p>`,

  quiz: [
    {
      q: 'Doubling the cutoff from 2.5σ to 5σ makes the force calculation roughly how much more expensive?',
      options: ['2 times', '4 times', '8 times', 'No change'],
      answer: 2,
      why: 'Neighbours per atom grow as rc³: (5/2.5)³ = 8.',
    },
    {
      q: 'What does shifting the potential (subtracting V(rc)) fix, and what does it leave broken?',
      options: ['Fixes the force jump, leaves the energy jump', 'Fixes the energy jump, leaves a force jump at rc', 'Fixes both', 'Fixes neither'],
      answer: 1,
      why: 'Shifting makes V continuous at rc, but the force (the slope) still drops suddenly to zero there. Force-shifting or switching smooths the force too.',
    },
    {
      q: 'Why can\'t we use a simple cutoff for charges (Coulomb)?',
      options: ['Charges are too small', '1/r decays so slowly that distant shells contribute more and more; truncation gives wrong results', 'Coulomb forces are always repulsive', 'It is too cheap to need one'],
      answer: 1,
      why: 'The number of atoms in a shell grows as r², so r²·(1/r) grows with distance. The sum does not converge. Ewald summation (Chapter 5) handles this properly.',
    },
  ],
};
