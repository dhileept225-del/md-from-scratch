import { tex } from '../tex.js';

export default {
  id: 'lennard-jones',
  summary: 'Two neutral atoms still feel each other: a gentle pull when they are near, a hard push when they overlap. One simple curve captures both.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/lj-potential.mp4', poster: 'videos/lj-potential.jpg' },
  sim: 'pair',
  learn: [
    'Why neutral atoms attract at a distance and repel up close',
    'What the two parameters σ and ε mean',
    'How to compute the Lennard-Jones energy and force in code',
  ],
  simHint: 'Drag the red atom and watch the arrows flip from pull to push. Then set r = 1.05 and press Release: the atoms vibrate in the well. Now try r = 0.95. Why do they fly apart?',

  idea: `
<p>Argon atoms have no charge, yet argon still turns into a liquid when it is cold enough. Something must hold the atoms together.</p>
<p>When two atoms are <b>far apart</b>, they hardly notice each other. A bit <b>closer</b>, the electron clouds of each atom shift slightly and create tiny, fluctuating dipoles that attract each other. This is the <b>van der Waals</b> (dispersion) pull. If they come <b>too close</b>, the electron clouds overlap and a very strong <b>repulsion</b> pushes them apart.</p>
<p>The <b>Lennard-Jones potential</b> puts all of this into one curve V(r), the energy of the pair at distance r. It has just two numbers:</p>
<ul>
  <li><b>σ</b> (sigma), the size: the distance where V = 0. Atoms cannot get much closer than this.</li>
  <li><b>ε</b> (epsilon), the stickiness: the depth of the energy well. A bigger ε means atoms hold on to each other more strongly.</li>
</ul>
<p>The bottom of the well, at r ≈ 1.12σ, is where the push and pull balance and the force is zero. Atoms in a liquid or solid jiggle around this distance.</p>
<div class="callout"><b>Key idea.</b> The force is the slope of the energy curve. Downhill to the left of the minimum means push, uphill to the right means pull. An MD code only needs this one rule, applied to every pair of atoms.</div>`,

  algorithm: `
<p>To get the force on every atom, add up the Lennard-Jones force from every other atom:</p>
<ol class="steps">
  <li>Set all forces to zero.</li>
  <li>Loop over every <b>pair</b> i &lt; j, so each pair is counted once.</li>
  <li>Find the distance vector <b>d</b> = r<sub>i</sub> − r<sub>j</sub> and its length r.</li>
  <li>If r is larger than the <b>cutoff</b> (usually 2.5σ), skip the pair. The force there is almost zero.</li>
  <li>Compute the force size F(r) and add F·<b>d</b>/r to atom i.</li>
  <li>Subtract the same amount from atom j. Newton's third law: equal and opposite.</li>
</ol>
<pre class="pseudo">forces = zeros(N, 3)
for i in range(N):
    for j in range(i + 1, N):
        d = pos[i] - pos[j]
        r = length(d)
        if r < r_cut:
            f = F(r) * d / r      # vector along the line between the atoms
            forces[i] += f
            forces[j] -= f</pre>
<p>This double loop costs about N²/2 pair checks. That is fine for a few hundred atoms. Lesson 2.6 shows how <b>neighbour lists</b> make it much faster.</p>`,

  code: {
    note: 'Two parts: plot V(r) and F(r), then let two atoms vibrate in the well using the same velocity Verlet loop as Lesson 1.1. Run <code>pip install numpy matplotlib</code>, then <code>python lj.py</code>. The energy E should stay constant while r oscillates.',
    src: `import numpy as np
import matplotlib.pyplot as plt

eps, sigma = 1.0, 1.0

def lj(r):
    """Lennard-Jones energy V and force F for two atoms a distance r apart."""
    sr6 = (sigma / r) ** 6
    V = 4 * eps * (sr6**2 - sr6)
    F = 24 * eps * (2 * sr6**2 - sr6) / r    # F = -dV/dr  (positive = push apart)
    return V, F

# --- 1. plot the curves ---
r = np.linspace(0.9, 3.0, 400)
V, F = lj(r)
r_min = 2 ** (1/6) * sigma
print(f"minimum at r = {r_min:.4f} sigma, V = {lj(r_min)[0]:.4f} eps")

plt.plot(r, V, label="V(r)")
plt.plot(r, F, "--", label="F(r)")
plt.axhline(0, color="gray", lw=0.8)
plt.ylim(-3, 3); plt.xlabel("r / sigma"); plt.legend()
plt.show()

# --- 2. let two atoms vibrate (velocity Verlet on the distance r) ---
mu, dt = 0.5, 0.002        # reduced mass of two atoms with mass 1
r, v = 1.05, 0.0           # start a little inside the minimum, at rest
f = lj(r)[1]
for step in range(1500):
    v += 0.5 * dt * f / mu
    r += dt * v
    V, f = lj(r)
    v += 0.5 * dt * f / mu
    if step % 250 == 0:
        E = 0.5 * mu * v**2 + V
        print(f"t = {step*dt:4.2f}   r = {r:.3f}   E = {E:.5f}")`,
  },

  maths: () => `
<h4>1. The potential</h4>
${tex(String.raw`V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]`, true)}
<p>The <b>r<sup>−6</sup></b> term is the attraction. It comes from London dispersion between fluctuating dipoles. The <b>r<sup>−12</sup></b> term is the repulsion. Its exact form is a convenient choice (it is the square of r<sup>−6</sup>, which is cheap to compute), not a law of nature.</p>
<h4>2. The force</h4>
${tex(String.raw`F(r) = -\frac{dV}{dr} = \frac{24\varepsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]`, true)}
<p>F &gt; 0 pushes the atoms apart; F &lt; 0 pulls them together. As a vector on atom i:</p>
${tex(String.raw`\mathbf{F}_{ij} = F(r_{ij})\,\frac{\mathbf{r}_i - \mathbf{r}_j}{r_{ij}}, \qquad \mathbf{F}_{ji} = -\mathbf{F}_{ij}`, true)}
<h4>3. Key points on the curve</h4>
<p>V = 0 when σ/r = 1, so at <b>r = σ</b>. The minimum is where F = 0:</p>
${tex(String.raw`2\left(\frac{\sigma}{r}\right)^{12} = \left(\frac{\sigma}{r}\right)^{6} \;\Rightarrow\; r_{\min} = 2^{1/6}\sigma \approx 1.122\,\sigma, \qquad V(r_{\min}) = -\varepsilon`, true)}
<h4>4. Many atoms</h4>
<p>The total energy is a sum over all pairs:</p>
${tex(String.raw`U = \sum_{i<j} V(r_{ij})`, true)}
<p class="muted">Typical values for argon: σ ≈ 3.4 Å and ε/k<sub>B</sub> ≈ 120 K. In reduced units, σ = ε = 1, so every LJ system looks the same, and results can be scaled to any noble gas.</p>`,

  quiz: [
    {
      q: 'At what distance is the Lennard-Jones force zero?',
      options: ['r = 0', 'r = σ', 'r = 2^(1/6) σ ≈ 1.12 σ', 'r = 2.5 σ'],
      answer: 2,
      why: 'The force is zero at the bottom of the well, r = 2^(1/6) σ. At r = σ the energy is zero, but the curve is steep there, so the force is strongly repulsive.',
    },
    {
      q: 'What does ε control?',
      options: ['The size of the atoms', 'The depth of the energy well: how strongly atoms stick', 'The mass of the atoms', 'The time step'],
      answer: 1,
      why: 'ε is the well depth. σ sets the size. A bigger ε means atoms hold together more strongly, which gives higher melting and boiling points.',
    },
    {
      q: 'Which term makes atoms repel when they are very close?',
      options: ['The (σ/r)¹² term', 'The (σ/r)⁶ term', 'The 4ε in front', 'The cutoff'],
      answer: 0,
      why: 'At small r, (σ/r)¹² grows much faster than (σ/r)⁶, so the positive repulsive term wins and the energy shoots up.',
    },
  ],
};
