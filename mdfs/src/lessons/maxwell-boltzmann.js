import { tex } from '../tex.js';

export default {
  id: 'maxwell-boltzmann',
  summary: 'Temperature is atomic motion. To start a simulation at a chosen temperature, we give each atom a random velocity from the Maxwell–Boltzmann distribution, remove any net drift, and rescale to hit the target exactly.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/maxwell-boltzmann.mp4', poster: 'videos/maxwell-boltzmann.jpg' },
  sim: 'velocities',
  learn: [
    'How temperature and atomic velocities are linked',
    'The Maxwell–Boltzmann distribution, and why heavy atoms move slower',
    'The three-step recipe for initial velocities',
  ],
  simHint: 'It starts with every atom at the same speed: right temperature, wrong shape. Tick "Collisions on" and watch it relax to the Maxwell–Boltzmann curve. Then try the recipe: press 1, 2, 3 and watch the temperature and drift readouts. Switch to heavy atoms: what happens to the temperature?',

  idea: `
<p>In MD, <b>temperature is the average kinetic energy</b> of the atoms. Hotter means faster. But at a given temperature the atoms do <i>not</i> all move at the same speed. Some are slow, some fast, and the speeds follow the <b>Maxwell–Boltzmann distribution</b>.</p>
<ul>
  <li>Each velocity <b>component</b> (v<sub>x</sub>, v<sub>y</sub>, v<sub>z</sub>) is a <b>Gaussian</b> centred on zero, with variance k<sub>B</sub>T/m.</li>
  <li>The <b>speed</b> |v| then follows the skewed Maxwell–Boltzmann curve: no atoms at zero speed, a peak, and a long tail of fast atoms.</li>
  <li><b>Heavier atoms move slower</b> at the same temperature: typical speed ∝ √(T/m).</li>
</ul>
<p>The three-step recipe used by MD codes:</p>
<ol class="steps">
  <li><b>Draw</b> each component from a Gaussian with variance k<sub>B</sub>T/m.</li>
  <li><b>Remove the drift</b>: subtract the centre-of-mass velocity so the whole box does not fly off.</li>
  <li><b>Rescale</b> all velocities by √(T<sub>target</sub>/T<sub>now</sub>) so the temperature is exactly right.</li>
</ol>
<p>Getting the exact shape right is less critical than it sounds. Collisions redistribute energy within a few picoseconds, so any reasonable start relaxes to Maxwell–Boltzmann. The Gaussian draw just gets you there immediately.</p>
<div class="callout"><b>Key idea.</b> Temperature fixes the <i>width</i> of the velocity distribution, not the speed of each atom. Draw from a Gaussian, remove the drift, rescale.</div>`,

  algorithm: `
<p>The recipe with masses, for N atoms:</p>
<pre class="pseudo">for each atom i, each direction d:
    v[i, d] = sqrt(kB * T / m[i]) * gaussian()
v_cm = sum(m[i] * v[i]) / sum(m)        # mass-weighted
v -= v_cm                               # no net momentum
T_now = sum(m[i] * |v[i]|^2) / ((3N - 3) * kB)
v *= sqrt(T_target / T_now)</pre>
<ol class="steps">
  <li>Use the <b>mass-weighted</b> centre-of-mass velocity when atoms have different masses.</li>
  <li>Count <b>3N − 3</b> degrees of freedom after removing the centre-of-mass motion. Subtract more if bonds are constrained (Lesson 4.4).</li>
  <li>Use a <b>new random seed</b> for each independent run. Different seeds give different trajectories, which is how you estimate error bars.</li>
</ol>
<p>In LAMMPS this is <code>velocity all create 300.0 4928459 dist gaussian mom yes</code>. In GROMACS it is <code>gen_vel = yes</code>, <code>gen_temp = 300</code>.</p>`,

  code: {
    note: 'The three-step recipe, followed by a check against the Maxwell–Boltzmann formula: a text histogram of sampled speeds next to theory. Run <code>pip install numpy</code>, then <code>python velocities.py</code>.',
    src: `import numpy as np

rng = np.random.default_rng(42)
N, T, m = 1000, 1.5, 1.0                 # reduced units (kB = 1)

# 1. draw each velocity component from a Gaussian with variance kT/m
v = rng.normal(0.0, np.sqrt(T / m), size=(N, 3))
print("drift before:", v.mean(axis=0).round(4))

# 2. remove the centre-of-mass velocity (no drifting box)
v -= v.mean(axis=0)

# 3. rescale so the temperature is exactly T  (3N - 3 degrees of freedom)
T_now = m * (v ** 2).sum() / (3 * N - 3)
v *= np.sqrt(T / T_now)
print(f"T before scaling {T_now:.4f}, after {m * (v ** 2).sum() / (3 * N - 3):.4f}")

# compare the speeds with the Maxwell-Boltzmann distribution
speed = np.linalg.norm(v, axis=1)
hist, edges = np.histogram(speed, bins=12, range=(0, 4 * np.sqrt(T / m)), density=True)
c = 0.5 * (edges[1:] + edges[:-1])
mb = 4 * np.pi * (m / (2 * np.pi * T)) ** 1.5 * c**2 * np.exp(-m * c**2 / (2 * T))
print("\\n speed  sampled  theory")
for ci, h, f in zip(c, hist, mb):
    print(f"{ci:5.2f}  {h:7.3f}  {f:6.3f}  " + "#" * int(h * 60))
print(f"\\nmean speed {speed.mean():.3f}   theory sqrt(8kT/pi m) = {np.sqrt(8 * T / (np.pi * m)):.3f}")`,
  },

  maths: () => `
<h4>1. One component: a Gaussian</h4>
${tex(String.raw`p(v_x) = \sqrt{\frac{m}{2\pi k_B T}}\;\exp\!\left(-\frac{m v_x^2}{2 k_B T}\right), \qquad \langle v_x^2\rangle = \frac{k_B T}{m}`, true)}
<h4>2. The speed: Maxwell–Boltzmann</h4>
<p>Combining three independent Gaussians and counting the shell of radius v gives</p>
${tex(String.raw`f(v) = 4\pi\left(\frac{m}{2\pi k_B T}\right)^{3/2} v^2\,\exp\!\left(-\frac{m v^2}{2 k_B T}\right)`, true)}
${tex(String.raw`v_{\text{peak}} = \sqrt{\frac{2k_BT}{m}},\qquad \langle v\rangle = \sqrt{\frac{8k_BT}{\pi m}},\qquad v_{\text{rms}} = \sqrt{\frac{3k_BT}{m}}`, true)}
<h4>3. Temperature from velocities (equipartition)</h4>
${tex(String.raw`\tfrac12 m\langle v_x^2\rangle = \tfrac12 k_B T \quad\Rightarrow\quad T = \frac{1}{N_{\text{dof}}\,k_B}\sum_i m_i |\mathbf v_i|^2, \qquad N_{\text{dof}} = 3N - 3`, true)}
<h4>4. Rescaling</h4>
${tex(String.raw`\mathbf v_i \;\to\; \lambda\,\mathbf v_i, \qquad \lambda = \sqrt{T_{\text{target}}/T_{\text{now}}}`, true)}
<p class="muted">Example: argon at 300 K has v<sub>rms</sub> = √(3k<sub>B</sub>T/m) ≈ 433 m/s, roughly the speed of a passenger jet. Xenon at the same temperature: about 239 m/s.</p>`,

  quiz: [
    {
      q: 'From which distribution do we draw each velocity component (vx, vy, vz)?',
      options: ['A uniform distribution between −1 and 1', 'A Gaussian centred on zero with variance kBT/m', 'The Maxwell–Boltzmann speed curve', 'All components equal to √(kBT/m)'],
      answer: 1,
      why: 'Each component is Gaussian with mean zero and variance kBT/m. Combining three of them gives the Maxwell–Boltzmann distribution of speeds.',
    },
    {
      q: 'Why remove the centre-of-mass velocity after drawing velocities?',
      options: ['To make the run faster', 'Otherwise the whole box drifts, and that motion wrongly counts as temperature', 'Because velocities must be positive', 'To change the density'],
      answer: 1,
      why: 'Random draws leave a small net momentum. The box would slowly fly off, and that collective motion would be counted as thermal energy.',
    },
    {
      q: 'Xenon is about 3.3 times heavier than argon. At the same temperature, how does xenon\'s rms speed compare?',
      options: ['3.3 times faster', 'The same', 'About 0.55 times (1/√3.3)', 'About 0.30 times (1/3.3)'],
      answer: 2,
      why: 'v_rms = √(3kBT/m) scales as 1/√m: 1/√3.3 ≈ 0.55.',
    },
  ],
};
