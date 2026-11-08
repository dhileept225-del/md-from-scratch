import { tex } from '../tex.js';

export default {
  id: 'ewald',
  summary: 'Ewald summation computes the full periodic Coulomb energy exactly. It splits 1/r into a short-ranged part, summed directly with a cutoff, and a smooth part, summed in Fourier space with a few wave vectors. A self-energy term completes it.',
  minutes: 18,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/ewald.mp4', poster: 'videos/ewald.jpg' },
  sim: 'ewald',
  learn: [
    'The Gaussian screening idea behind the erfc/erf split',
    'The real-space, reciprocal-space and self terms, and their forces',
    'How to choose α, the cutoff and k<sub>max</sub>; why the total does not depend on α',
  ],
  simHint: 'Move α and watch the real and reciprocal parts trade off while the total stays (almost) the same. Find the bottom of the U-shaped error curve. Then shrink the cutoff r_c: the curve moves right (you need a larger α). Lower k_max: the right side of the U rises. This is exactly how MD codes choose their Ewald parameters.',

  idea: `
<p>Lesson 5.1 showed that the Coulomb sum over a periodic system converges only conditionally. Ewald (1921) found a way to rewrite it as two sums that both converge fast and absolutely.</p>
<p><b>The trick.</b> Surround every point charge with a Gaussian cloud of the opposite charge. The point charge plus its cloud is <b>screened</b>: from outside it looks neutral, so its potential, erfc(αr)/r, dies off quickly. Sum these screened interactions directly in real space, with a cutoff and neighbour lists, as for Lennard-Jones.</p>
<p>Then add back Gaussian clouds of the original sign. Their potential, erf(αr)/r, is <b>smooth</b> everywhere: no singularity at r = 0. A smooth periodic function needs only a few Fourier components, so this part is summed in <b>reciprocal space</b> over wave vectors <b>k</b> = 2π<b>m</b>/L, using the structure factor S(<b>k</b>) = Σ q<sub>j</sub> e<sup>i<b>k</b>·<b>r</b><sub>j</sub></sup>.</p>
<p>Finally, each point charge has been made to interact with its own Gaussian in the Fourier part. That artificial <b>self-energy</b>, −(α/√π) Σ q<sub>i</sub>², is subtracted. (With the usual "tin-foil" boundary condition there is no extra surface term.)</p>
<p><b>The parameter α</b> decides how the work is divided. A large α makes the real-space part short-ranged (small cutoff) but needs more k vectors; a small α does the opposite. The <i>total</i> does not depend on α: a powerful check. In practice r<sub>c</sub> is set (about 1 nm), α is chosen so that erfc(αr<sub>c</sub>) ≈ 10⁻⁵, and k<sub>max</sub> is raised until the reciprocal error matches. The cost scales as N<sup>3/2</sup> at best; PME (Lesson 5.3) brings it to N log N.</p>
<div class="callout"><b>Key idea.</b> E = E<sub>real</sub> (screened, short-ranged) + E<sub>recip</sub> (smooth, Fourier) + E<sub>self</sub>. The split parameter α moves work between the two sums but does not change the answer.</div>`,

  algorithm: `
<p>Energy and forces for a neutral cubic box of side L (units with 1/(4πε₀) = 1):</p>
<pre class="pseudo"># real space (with cutoff rc and neighbour lists)
for pairs i < j with r = |r_ij (minimum image)| < rc:
    E_real += q_i q_j erfc(α r) / r
    F_ij    = q_i q_j [erfc(α r)/r + 2α/√π exp(-α² r²)] r_ij / r²

# reciprocal space
for k = 2π m / L, m ≠ 0, |m_x|,|m_y|,|m_z| ≤ k_max:
    A(k) = exp(-k² / 4α²) / k²
    S(k) = Σ_j q_j exp(i k·r_j)
    E_recip += (2π/V) A(k) |S(k)|²
    F_i     += (4π q_i/V) A(k) k Im[S*(k) exp(i k·r_i)]

E_self = -(α/√π) Σ_i q_i²
E = E_real + E_recip + E_self</pre>
<ol class="steps">
  <li>Use k and −k symmetry to halve the reciprocal work, and update S(k) once per step for all atoms.</li>
  <li>Excluded pairs (bonded 1–2, 1–3) are inside the reciprocal sum anyway, so their erf(αr)/r part must be subtracted in real space.</li>
  <li>Check: the energy must be independent of α, and the forces must match a numerical derivative (both in the code tab).</li>
  <li>GROMACS sets α from <code>rcoulomb</code> and <code>ewald-rtol</code> (default 10⁻⁵); LAMMPS from <code>kspace_style ewald 1e-5</code>. Plain Ewald is fine for a few hundred atoms; use PME/PPPM beyond that.</li>
</ol>`,

  code: {
    note: 'A complete Ewald implementation (energy and forces). Part 1 reproduces the NaCl Madelung constant to six digits with four different α values. Part 2 shows convergence with k<sub>max</sub> for 40 random ions and checks the forces numerically. Run <code>pip install numpy</code>, then <code>python ewald.py</code> (about 1 s).',
    src: `import numpy as np
from math import erfc, sqrt, pi
from itertools import product
erfc_v = np.vectorize(erfc)

def ewald(pos, q, L, alpha, kmax, nimg=1, rc=None, forces=False):
    """Ewald energy (units with 1/(4 pi eps0) = 1) for a neutral cubic box of side L.
    Returns real, reciprocal and self parts (and forces if asked)."""
    N = len(q); V = L**3; F = np.zeros_like(pos)
    # real space: erfc(alpha r)/r over image cells -nimg..nimg (within rc if given)
    Er = 0.0
    for n in product(range(-nimg, nimg + 1), repeat=3):
        shift = np.array(n) * L
        d = pos[:, None, :] - pos[None, :, :] + shift
        r = np.linalg.norm(d, axis=2)
        mask = r > 1e-12
        if rc is not None: mask &= r < rc
        rr = np.where(mask, r, 1.0)
        e = np.where(mask, erfc_v(alpha * rr) / rr, 0.0)
        Er += 0.5 * np.sum(q[:, None] * q[None, :] * e)
        if forces:
            g = np.where(mask, (erfc_v(alpha * rr) / rr + 2 * alpha / sqrt(pi) * np.exp(-(alpha * rr)**2)) / rr**2, 0.0)
            F += np.sum((q[:, None] * q[None, :] * g)[..., None] * d, axis=1)
    # reciprocal space: (2 pi / V) sum_k exp(-k^2/4a^2)/k^2 |S(k)|^2
    Ek = 0.0
    for m in product(range(-kmax, kmax + 1), repeat=3):
        if m == (0, 0, 0): continue
        k = 2 * pi * np.array(m) / L; k2 = k @ k
        A = np.exp(-k2 / (4 * alpha**2)) / k2
        phase = np.exp(1j * pos @ k); S = np.sum(q * phase)
        Ek += 2 * pi / V * A * abs(S)**2
        if forces:
            F += (4 * pi / V) * A * (q * np.imag(np.conj(S) * phase))[:, None] * k
    Es = -alpha / sqrt(pi) * np.sum(q**2)
    return (Er, Ek, Es, F) if forces else (Er, Ek, Es)

# ---------- Part 1: NaCl unit cell (8 ions, cube of side 2 nearest-neighbour distances) ----------
pos = np.array([p for p in product([0, 1], repeat=3)], float); q = np.array([(-1.0) ** sum(p) for p in pos])
print("Part 1: rock salt, 8-ion cell. Madelung constant M = -2 E / N  (exact 1.747565)")
print("   alpha    real part    recip part    self part      M")
for alpha in [1.0, 1.5, 2.0, 3.0]:
    Er, Ek, Es = ewald(pos, q, 2.0, alpha, kmax=6, nimg=3)
    print(f"   {alpha:4.1f}   {Er:10.5f}   {Ek:10.5f}   {Es:10.5f}   {-2 * (Er + Ek + Es) / 8:.6f}")
print("   The split changes with alpha; the total does not.")

# ---------- Part 2: a random ionic configuration: forces and convergence ----------
rng = np.random.default_rng(0); N, L = 40, 4.0
pos = rng.uniform(0, L, (N, 3)); q = np.array([1.0, -1.0] * (N // 2))
alpha, rc = 2.2, 2.0                                    # erfc(alpha rc) ~ 5e-10, see below
print(f"\\nPart 2: {N} random ions in a box of side {L}; real-space cutoff rc = {rc}")
print("   kmax    E_total        (alpha = 2.2)")
ref = None
for kmax in [2, 4, 6, 8, 10]:
    E = sum(ewald(pos, q, L, alpha, kmax, nimg=1, rc=rc))
    print(f"   {kmax:4d}   {E:12.6f}"); ref = E
print("   real-space accuracy depends on erfc(alpha rc):", ", ".join(f"alpha={a}: {erfc(a * rc):.0e}" for a in (1.5, 2.2, 3.0)))

Er, Ek, Es, F = ewald(pos, q, L, alpha, 10, nimg=1, rc=rc, forces=True)
h = 1e-5; Fn = np.zeros(3)
for d in range(3):
    p1, p2 = pos.copy(), pos.copy(); p1[0, d] += h; p2[0, d] -= h
    Fn[d] = -(sum(ewald(p1, q, L, alpha, 10, 1, rc)) - sum(ewald(p2, q, L, alpha, 10, 1, rc))) / (2 * h)
print(f"\\n   force on ion 0, analytic: {np.round(F[0], 5)}")
print(f"   force on ion 0, numerical: {np.round(Fn, 5)}")
print(f"   sum of all forces: {np.abs(F.sum(0)).max():.1e}")
`,
  },

  maths: () => `
<h4>1. The split</h4>
${tex(String.raw`\frac{1}{r} = \frac{\operatorname{erfc}(\alpha r)}{r} + \frac{\operatorname{erf}(\alpha r)}{r}`, true)}
<p>A Gaussian charge cloud ρ(r) = q(α/√π)³ e<sup>−α²r²</sup> produces the potential q erf(αr)/r.</p>
<h4>2. The Ewald energy</h4>
${tex(String.raw`E = \frac12\sum_{\mathbf n}{}'\sum_{i,j} q_iq_j\frac{\operatorname{erfc}(\alpha|\mathbf r_{ij} + \mathbf nL|)}{|\mathbf r_{ij} + \mathbf nL|} + \frac{2\pi}{V}\sum_{\mathbf k\ne 0}\frac{e^{-k^2/4\alpha^2}}{k^2}\,|S(\mathbf k)|^2 - \frac{\alpha}{\sqrt\pi}\sum_i q_i^2`, true)}
${tex(String.raw`S(\mathbf k) = \sum_j q_j e^{i\mathbf k\cdot\mathbf r_j}, \qquad \mathbf k = \frac{2\pi}{L}(m_x, m_y, m_z)`, true)}
<h4>3. Errors and cost</h4>
${tex(String.raw`\delta E_{\text{real}} \sim \operatorname{erfc}(\alpha r_c), \qquad \delta E_{\text{recip}} \sim e^{-k_{\max}^2/4\alpha^2}, \qquad \text{cost} \propto N^{3/2} \text{ at the optimal } \alpha`, true)}
<p>Keeping both errors equal while N grows gives α ∝ N<sup>1/6</sup>, and the total work grows as N<sup>3/2</sup>.</p>`,

  quiz: [
    {
      q: 'Why is the reciprocal-space part of Ewald summation efficient?',
      options: ['It ignores distant charges', 'The erf(αr)/r potential is smooth, so only a few Fourier components are needed', 'It uses a cutoff in k-space of zero', 'It only includes nearest neighbours'],
      answer: 1,
      why: 'The Gaussian clouds remove the 1/r singularity. A smooth periodic function has rapidly decaying Fourier coefficients, exp(−k²/4α²).',
    },
    {
      q: 'You double α but keep the cutoff and k_max fixed. What happens?',
      options: ['The total energy doubles', 'The real-space part becomes shorter-ranged (more accurate), but more k vectors are needed for the same reciprocal accuracy', 'Nothing changes at all', 'The self term disappears'],
      answer: 1,
      why: 'α shifts work between the two sums. The exact total is independent of α, but with fixed k_max the reciprocal error grows.',
    },
    {
      q: 'What is the self-energy term in Ewald summation?',
      options: ['The energy of the box with itself', 'A correction removing each point charge\'s interaction with its own Gaussian cloud', 'The Lennard-Jones energy', 'An error from the cutoff'],
      answer: 1,
      why: 'The reciprocal sum includes i = j terms: each charge interacting with its own smooth Gaussian. E_self = −(α/√π) Σ q_i² removes them.',
    },
  ],
};
