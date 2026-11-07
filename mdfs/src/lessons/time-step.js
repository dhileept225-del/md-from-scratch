import { tex } from '../tex.js';

export default {
  id: 'time-step',
  summary: 'The time step Δt is the most important number in an MD run. Too big and the simulation goes wrong or explodes; too small and it takes forever. The fastest vibration in your system decides it.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/time-step.mp4', poster: 'videos/time-step.jpg' },
  sim: 'timestep',
  learn: [
    'Why the fastest vibration sets the time step',
    'The rule of thumb: 10–20 steps per vibration',
    'How to check your choice with energy conservation',
  ],
  simHint: 'Start with C–H at Δt = 1 fs. Slowly raise Δt and watch the computed steps (cyan) drift away from the exact curve (amber). Where does it explode? Then switch to C–C: the same Δt is suddenly very safe. Why?',

  idea: `
<p>Every MD step assumes the force stays <b>constant</b> during Δt. That is only true if the atoms barely move in one step. So Δt must be short compared with the <b>fastest motion</b> in the system.</p>
<p>The fastest motions are bond vibrations involving <b>hydrogen</b>, the lightest atom. A C–H bond vibrates back and forth about once every <b>11 fs</b>, and an O–H bond about every 9 fs. To follow a vibration well you need about <b>10–20 steps per period</b>, which gives the familiar <b>Δt ≈ 1 fs</b> for all-atom MD.</p>
<ul>
  <li><b>Too big:</b> the steps overshoot. The energy drifts, and past a limit (about T/π) the simulation explodes.</li>
  <li><b>Too small:</b> everything is accurate, but you need more steps. 1 ns at 1 fs is already a million steps.</li>
</ul>
<p>Common choices in practice:</p>
<ul>
  <li><b>1 fs:</b> all-atom, flexible bonds.</li>
  <li><b>2 fs:</b> bonds to hydrogen held fixed with constraints (SHAKE/LINCS, Lesson 4.4).</li>
  <li><b>4 fs:</b> hydrogen mass repartitioning, which makes H heavier.</li>
  <li><b>10–40 fs:</b> coarse-grained models with no light atoms at all.</li>
</ul>
<div class="callout"><b>Key idea.</b> The lightest, stiffest bond sets the time step for the whole system. Remove or slow down that motion, and you can take bigger steps.</div>`,

  algorithm: `
<p>How to choose and check Δt for a new system:</p>
<ol class="steps">
  <li><b>Find the fastest motion.</b> Usually the stiffest bond to a light atom: O–H, N–H, C–H.</li>
  <li><b>Estimate its period</b> from its vibrational wavenumber: T = 1/(c ν̃). For 3000 cm⁻¹, T ≈ 11 fs.</li>
  <li><b>Pick Δt ≈ T/20 to T/10.</b> Never go near T/π, where the update becomes unstable.</li>
  <li><b>Run a short NVE test</b> (no thermostat) and plot the total energy. It should fluctuate but not drift.</li>
  <li><b>If it drifts</b>, halve Δt or constrain the bonds to hydrogen, then test again.</li>
</ol>
<pre class="pseudo">for dt in [0.5, 1.0, 2.0]:          # fs
    run_nve(dt, n_steps=10_000)
    plot(total_energy)            # flat = good, drifting = too big</pre>
<p>The code tab does exactly this test for a single C–H bond.</p>`,

  code: {
    note: 'A C–H bond as a spring in real units (fs), integrated with velocity Verlet at several time steps. Run <code>pip install numpy</code>, then <code>python timestep.py</code>. Watch the error grow with Δt, and the explosion past T/π ≈ 3.5 fs.',
    src: `import numpy as np

# A C-H bond as a spring, in real units (fs)
wavenumber = 3000.0                 # cm^-1, typical C-H stretch
c = 2.998e-5                        # speed of light in cm/fs
T = 1 / (wavenumber * c)            # vibration period in fs
omega = 2 * np.pi / T
print(f"C-H vibration period = {T:.1f} fs")

def run(dt, t_total=200.0):
    """Integrate the bond for t_total fs; return the worst relative energy error."""
    x, v = 0.1, 0.0                 # stretched by 0.1 A, at rest
    E0 = 0.5 * (omega * x) ** 2     # energy per unit mass
    worst = 0.0
    a = -omega**2 * x               # a = F/m for a spring
    for _ in range(int(t_total / dt)):
        v += 0.5 * a * dt           # velocity Verlet (the update from Lesson 1.1)
        x += v * dt
        a = -omega**2 * x
        v += 0.5 * a * dt
        E = 0.5 * v**2 + 0.5 * (omega * x) ** 2
        worst = max(worst, abs(E - E0) / E0)
        if worst > 10:
            return "blows up!"
    return f"max energy error {100 * worst:6.1f} %"

for dt in [0.25, 0.5, 1.0, 2.0, 3.0, 4.0]:
    print(f"dt = {dt:4.2f} fs  ({T / dt:5.1f} steps per vibration):  {run(dt)}")`,
  },

  maths: () => `
<h4>1. Period of a bond vibration</h4>
<p>From the spectroscopic wavenumber ν̃ (in cm⁻¹) and the speed of light c:</p>
${tex(String.raw`T = \frac{1}{c\,\tilde\nu}, \qquad \omega = \frac{2\pi}{T} = \sqrt{\frac{k}{\mu}}, \qquad \mu = \frac{m_1 m_2}{m_1 + m_2}`, true)}
<p>For C–H, ν̃ ≈ 3000 cm⁻¹ gives T ≈ 11.1 fs. The small reduced mass μ (hydrogen is light) is what makes it so fast.</p>
<h4>2. Stability limit</h4>
<p>For a harmonic vibration, velocity Verlet (and the update from Lesson 1.3) stays stable only if</p>
${tex(String.raw`\omega\,\Delta t < 2 \quad\Longleftrightarrow\quad \Delta t < \frac{T}{\pi}`, true)}
<p>For C–H this is about 3.5 fs. Beyond it, errors grow by a constant factor every step: the explosion you saw.</p>
<h4>3. Accuracy</h4>
<p>Velocity Verlet's energy error for one vibration scales as the square of ωΔt:</p>
${tex(String.raw`\frac{\delta E}{E} \sim \frac{(\omega\,\Delta t)^2}{8}`, true)}
<p>Halving Δt cuts the error about four times. With 20 steps per period, ωΔt ≈ 0.31 and the error is about 1%.</p>
<h4>4. Cost</h4>
${tex(String.raw`N_{\text{steps}} = \frac{t_{\text{sim}}}{\Delta t}: \qquad 100\ \text{ns} \,/\, 2\ \text{fs} = 5\times10^{7}\ \text{steps}`, true)}`,

  quiz: [
    {
      q: 'A C–H bond vibrates with a period of about 11 fs. Which time step is a sensible choice?',
      options: ['10 fs', '5 fs', '1 fs', '0.001 fs'],
      answer: 2,
      why: '1 fs gives about 11 steps per vibration, which is enough. 5 and 10 fs are past the stability limit (T/π ≈ 3.5 fs). 0.001 fs works but wastes a thousand times more steps.',
    },
    {
      q: 'Why can many simulations use Δt = 2 fs?',
      options: ['Computers got faster', 'Bonds to hydrogen are held fixed, so the fastest vibration is gone', 'Hydrogen atoms are removed from the system', 'The thermostat corrects the errors'],
      answer: 1,
      why: 'Constraint algorithms such as SHAKE or LINCS freeze bond lengths to hydrogen. The fastest remaining motions are slower, so the step can double.',
    },
    {
      q: 'How many steps do you need for 100 ns with Δt = 2 fs?',
      options: ['5 × 10⁴', '5 × 10⁷', '2 × 10⁸', '5 × 10¹⁰'],
      answer: 1,
      why: '100 ns = 10⁸ fs. Divided by 2 fs per step gives 5 × 10⁷ steps.',
    },
  ],
};
