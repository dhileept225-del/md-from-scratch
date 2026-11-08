import { tex } from '../tex.js';

export default {
  id: 'coulomb-pbc',
  summary: 'Coulomb interactions decay so slowly that a plain cutoff gives wrong answers, and in a periodic box every ion interacts with infinitely many images. The rock-salt Madelung constant shows the problem, and why neutral groups and damping help.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/coulomb-pbc.mp4', poster: 'videos/coulomb-pbc.jpg' },
  sim: 'coulomb',
  learn: [
    'Why the 1/r Coulomb interaction cannot simply be cut off like Lennard-Jones',
    'Conditional convergence: the result depends on the order of summation (Madelung constant)',
    'Neutral groups, reaction field and damped-shifted (Wolf) sums as approximations, leading to Ewald',
  ],
  simHint: 'With Sphere cutoff, slide R from 1.5 to 12 and watch the estimate jump wildly as whole shells of + or − ions enter. Check the net charge. Switch to Neutral cubes: it converges to 1.7476 within a few shells. Then try the Wolf sum with α = 0 and α = 0.3. What does the damping do?',

  idea: `
<p>For Lennard-Jones we used a cutoff of about 2.5σ and a small tail correction (Lesson 2.5). That works because the number of neighbours at distance r grows as r², while the LJ energy falls as r⁻⁶: the far region contributes almost nothing. Coulomb energy falls only as 1/r. The far shells hold more and more charges, and their total contribution does not die away at all. The sum is finite only because positive and negative charges cancel.</p>
<p><b>Conditional convergence.</b> The classic example is the Madelung constant of rock salt, M = 1.7476: the electrostatic energy of an ion in an infinite NaCl crystal. Add up the ions inside a growing sphere and the answer swings wildly (2.1, −4.2, 8.0, …) and never settles, because each new shell is not neutral. Add them in neutral blocks (expanding cubes with boundary ions shared, Evjen's method) and it converges within a few shells. For such sums, the answer depends on <i>how</i> you add the terms.</p>
<p><b>In a simulation</b> with periodic boundaries (Lesson 2.4), each charge interacts with every other charge and with all their periodic images. A plain cutoff of 1–1.2 nm makes the included region randomly charged and produces known artefacts: wrong ion–ion structure, artificial ordering of water near the cutoff, and wrong dielectric behaviour.</p>
<p><b>Quick fixes</b> keep a cutoff but soften its effect:</p>
<ul>
<li><b>Charge groups</b>: include or exclude whole neutral groups (old GROMACS, GROMOS).</li>
<li><b>Reaction field</b>: treat everything beyond the cutoff as a continuum with dielectric constant ε<sub>RF</sub>.</li>
<li><b>Damped, shifted sums</b> (Wolf, DSF): damp 1/r with erfc(αr), shift it to zero at the cutoff, and neutralise the cutoff sphere. Cheap and often accurate for dense systems.</li>
</ul>
<p>The proper solution, used in almost all modern simulations, is <b>Ewald summation</b> and its fast version PME (Lessons 5.2–5.3). It splits 1/r into a short-ranged part summed directly and a smooth part summed exactly in Fourier space.</p>
<div class="callout"><b>Key idea.</b> Never use a plain cutoff for charged or polar systems. The Coulomb sum over a periodic system is only conditionally convergent; use Ewald/PME, or at least a damped-shifted or reaction-field method.</div>`,

  algorithm: `
<p>The Madelung constant three ways (reference ion +1 at the origin, nearest-neighbour distance 1):</p>
<pre class="pseudo"># 1. plain sphere: does not converge
M = -Σ_{r_j < R} q_j / r_j

# 2. neutral cubes (Evjen): converges fast
M = -Σ_{|x|,|y|,|z| ≤ L} w_j q_j / r_j
    w = 1 inside, 1/2 on faces, 1/4 on edges, 1/8 at corners

# 3. Wolf: damp, shift and neutralise
e_i = ½ Σ_{r_j < Rc} q_i q_j [erfc(α r_j)/r_j - erfc(α Rc)/Rc]
      - [erfc(α Rc)/(2 Rc) + α/√π] q_i²</pre>
<ol class="steps">
  <li>The shift term makes each pair interaction zero at Rc. Together with the self term, the cutoff sphere behaves as if it were neutral.</li>
  <li>α trades accuracy at short range against convergence at long range; α ≈ 0.2–0.3 Å⁻¹ with Rc ≈ 10–12 Å is common for ionic liquids and molten salts.</li>
  <li>In GROMACS, <code>coulombtype = Reaction-Field</code> or <code>PME</code>; plain <code>Cut-off</code> is a reaction field with ε<sub>RF</sub> = 1. LAMMPS: <code>pair_style coul/wolf</code>, <code>coul/dsf</code>, or <code>coul/long</code> with <code>kspace_style pppm</code>.</li>
</ol>`,

  code: {
    note: 'Computes the rock-salt Madelung constant with a plain spherical cutoff (never converges), with neutral Evjen cubes (converges to 1.74757), and with the Wolf damped-shifted sum. Run <code>pip install numpy</code>, then <code>python coulomb_pbc.py</code> (under 1 s).',
    src: `import numpy as np
from math import erfc, sqrt, pi
M_EXACT = 1.747565                     # Madelung constant of rock salt (NaCl), nearest-neighbour distance 1

# rock-salt lattice around a reference ion (+1) at the origin: charge (-1)^(i+j+k)
n = 14
g = np.arange(-n, n + 1)
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
pos = np.stack([X.ravel(), Y.ravel(), Z.ravel()], 1).astype(float)
q = (-1.0) ** (np.abs(pos).sum(1))
keep = (pos != 0).any(1); pos, q = pos[keep], q[keep]
r = np.linalg.norm(pos, axis=1)

# ---------- 1. plain spherical cutoff: does not converge ----------
print("Part 1: Madelung constant from a plain spherical cutoff  M(R) = -sum_{r<R} q_j / r_j")
for R in [2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0, 12.0]:
    m = r < R
    print(f"   R = {R:5.1f}:  M = {-np.sum(q[m] / r[m]):7.3f}   net charge inside = {q[m].sum() + 1:+.0f}")
print(f"   exact: {M_EXACT}")

# ---------- 2. neutral cubes (Evjen): converges ----------
print("\\nPart 2: neutral expanding cubes with Evjen weights (1/2 face, 1/4 edge, 1/8 corner)")
for L in [1, 2, 3, 4, 6, 8]:
    a = np.abs(pos); inside = (a <= L).all(1)
    w = np.ones(len(pos)) * 0.5 ** (a == L).sum(1)               # share boundary ions with neighbouring cubes
    print(f"   cube half-width {L}:  M = {-np.sum((w * q / r)[inside]):.6f}")

# ---------- 3. Wolf summation: damp, shift, neutralise ----------
def wolf(Rc, alpha):
    m = r < Rc
    erfcv = np.array([erfc(alpha * x) for x in r[m]])
    pair = np.sum(q[m] * (erfcv / r[m] - erfc(alpha * Rc) / Rc))             # q_ref = +1
    self_term = erfc(alpha * Rc) / (2 * Rc) + alpha / sqrt(pi)
    return 0.5 * pair - self_term                                            # energy per ion
print("\\nPart 3: Wolf (damped, shifted, charge-neutralised) sum; exact energy per ion = -M/2 = -0.873783")
for Rc in [2.0, 3.0, 4.0, 6.0, 8.0]:
    print(f"   Rc = {Rc:4.1f}:  alpha = 0  -> {wolf(Rc, 0.0):8.4f}     alpha = 0.3 -> {wolf(Rc, 0.3):9.6f}")

# ---------- 4. why Coulomb is long-ranged ----------
print("\\nPart 4: energy carried by pairs beyond a distance R (per unit density, 3D)")
print("   LJ:      integral of 4 pi r^2 * r^-6 dr from R to infinity = 4 pi / (3 R^3)  -> finite, falls fast")
print("   Coulomb: integral of 4 pi r^2 * r^-1 dr from R to infinity = infinite       -> only cancellation between")
print("            + and - charges makes the sum finite, and the answer depends on how you add the terms up")
`,
  },

  maths: () => `
<h4>1. Coulomb energy of a periodic system</h4>
${tex(String.raw`U = \frac{1}{2}\sum_{\mathbf n}{}' \sum_{i,j} \frac{q_i q_j}{4\pi\varepsilon_0\,|\mathbf r_{ij} + \mathbf n L|}`, true)}
<p>The prime means i = j is skipped for <b>n</b> = 0. This sum is <i>conditionally convergent</i>: its value depends on the order of summation (and, for a sphere of images, on the dielectric surrounding it).</p>
<h4>2. Why the tail matters</h4>
${tex(String.raw`\int_R^\infty 4\pi r^2 \rho\, r^{-6}\, dr = \frac{4\pi\rho}{3R^3} \quad\text{(finite)}, \qquad \int_R^\infty 4\pi r^2 \rho\, r^{-1}\, dr = \infty`, true)}
<h4>3. Madelung constant</h4>
${tex(String.raw`U_{\text{ion}} = -\frac{M}{2}\,\frac{q^2}{4\pi\varepsilon_0 a}, \qquad M_{\text{NaCl}} = 1.747565`, true)}
<h4>4. Wolf (damped, shifted) pair energy</h4>
${tex(String.raw`E = \frac12\sum_i\sum_{j\ne i}^{r_{ij}<R_c} q_iq_j\left[\frac{\operatorname{erfc}(\alpha r_{ij})}{r_{ij}} - \frac{\operatorname{erfc}(\alpha R_c)}{R_c}\right] - \left[\frac{\operatorname{erfc}(\alpha R_c)}{2R_c} + \frac{\alpha}{\sqrt\pi}\right]\sum_i q_i^2`, true)}
<p>(in units where 1/(4πε₀) = 1). With α = 0.3 and R<sub>c</sub> = 8 the code tab gives −0.87364 per ion against the exact −0.87378.</p>`,

  quiz: [
    {
      q: 'Why can Lennard-Jones be cut off at 2.5σ but Coulomb cannot be cut off safely?',
      options: ['Coulomb forces are always repulsive', 'The number of neighbours grows as r² while Coulomb falls only as 1/r, so distant shells never become negligible', 'Coulomb is cheaper to compute', 'LJ is long-ranged'],
      answer: 1,
      why: 'For LJ, r² × r⁻⁶ dies quickly. For Coulomb, r² × r⁻¹ grows; only cancellation between charges keeps the sum finite.',
    },
    {
      q: 'Summing the NaCl Madelung constant over growing spheres gives values that jump around and never settle. Why?',
      options: ['Rounding errors', 'Each sphere contains a net charge, which changes from shell to shell', 'The lattice is not periodic', 'The ions are too close'],
      answer: 1,
      why: 'Spheres cut through shells of like charges. A neutral summation order (Evjen cubes, Ewald) converges.',
    },
    {
      q: 'What does the Wolf method do to make a cutoff sum of charges behave well?',
      options: ['Ignores all charges beyond 0.5 nm', 'Damps 1/r with erfc(αr), shifts it to zero at the cutoff, and neutralises the cutoff sphere', 'Adds a Fourier-space sum', 'Scales all charges by 0.5'],
      answer: 1,
      why: 'Damping makes the remaining interaction short-ranged, and the shift plus self term removes the effect of the net charge inside the cutoff sphere.',
    },
  ],
};
