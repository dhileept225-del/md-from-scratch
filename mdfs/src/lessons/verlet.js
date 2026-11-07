import { tex } from '../tex.js';

export default {
  id: 'verlet',
  summary: 'There are many ways to turn F = ma into steps. The obvious one, Euler, slowly pumps energy into the system. Velocity Verlet costs the same but keeps energy steady and can even run backwards. That is why every MD code uses it.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/euler-vs-verlet.mp4', poster: 'videos/euler-vs-verlet.jpg' },
  sim: 'integrators',
  learn: [
    'Why the simple Euler update makes energy grow',
    'How velocity Verlet works: half kick, move, half kick',
    'Two properties that make Verlet the MD standard: no energy drift and time reversibility',
  ],
  simHint: 'Watch the three dots orbit the same spring. Red (Euler) spirals out, and its energy climbs. Let it run, then press Reverse velocities: only one dot comes back to ◎. Try a bigger Δt too.',

  idea: `
<p>In Lesson 1.3 we updated the velocity first and then moved the atom. The most obvious rule does it the other way round. It moves the atom using its <b>old velocity</b>, then updates the velocity with the <b>old force</b>. This is the <b>explicit Euler</b> method.</p>
<p>For an orbit, Euler always steps along the tangent of the circle, which points slightly outward. Every step lands a bit too far out, so the orbit <b>spirals outward</b> and the energy <b>grows</b> without limit. A simulation like that heats itself up and is useless.</p>
<p><b>Velocity Verlet</b> splits the step into three moves:</p>
<ol class="steps">
  <li><b>Half kick:</b> update the velocity with half of the current force.</li>
  <li><b>Move:</b> update the position with that half-updated velocity.</li>
  <li><b>Half kick:</b> compute the force at the new position, and use it for the other half of the velocity update.</li>
</ol>
<p>The step is <b>symmetric in time</b>. Because of that:</p>
<ul>
  <li>The energy <b>wiggles</b> slightly but does not drift, even over millions of steps.</li>
  <li>It is <b>reversible</b>: flip all velocities and the system retraces its path back to the start.</li>
  <li>It needs only <b>one force calculation per step</b>, the same as Euler, because the new force is reused in the next step.</li>
</ul>
<div class="callout"><b>Key idea.</b> A good integrator is not just "accurate for one step". It must respect the physics over millions of steps. Verlet does: same cost as Euler, far better behaviour.</div>`,

  algorithm: `
<p>Velocity Verlet in a real MD loop. Note that the force from the end of one step is reused at the start of the next:</p>
<pre class="pseudo">F = compute_forces(pos)                 # once, before the loop
for step in range(n_steps):
    vel += 0.5 * dt * F / mass          # half kick
    pos += dt * vel                     # move
    F = compute_forces(pos)             # the only force call this step
    vel += 0.5 * dt * F / mass          # half kick</pre>
<p>For comparison, explicit Euler and Euler–Cromer (Lesson 1.3):</p>
<pre class="pseudo"># explicit Euler: energy grows
F = compute_forces(pos)
pos, vel = pos + dt * vel, vel + dt * F / mass

# Euler–Cromer: stable but first-order and not reversible
F = compute_forces(pos)
vel += dt * F / mass
pos += dt * vel</pre>
<p>Velocity Verlet is equivalent to <b>leapfrog</b> (velocities at half steps), which is what GROMACS uses by default. Both are the same algorithm written differently.</p>`,

  code: {
    note: 'The three update rules on the same oscillator for 1000 steps, followed by a reversibility test: flip the velocity, run 1000 more steps, and check whether you return to x = 1. Run <code>pip install numpy</code>, then <code>python verlet.py</code>.',
    src: `import numpy as np

# A 1D oscillator (a bond): a = -x, exact period 2*pi, exact energy 0.5
def accel(x):
    return -x

def energy(x, v):
    return 0.5 * v**2 + 0.5 * x**2

def euler(x, v, dt):            # explicit Euler: old velocity, old force
    return x + v * dt, v + accel(x) * dt

def euler_cromer(x, v, dt):     # Lesson 1.3: new velocity, then move
    v = v + accel(x) * dt
    return x + v * dt, v

def verlet(x, v, dt):           # velocity Verlet: half kick, move, half kick
    v = v + 0.5 * accel(x) * dt
    x = x + v * dt
    v = v + 0.5 * accel(x) * dt # in a real code, keep this force for the next step
    return x, v

dt, n = 0.1, 1000               # 1000 steps = about 16 periods
for name, step in [("Euler", euler), ("Euler-Cromer", euler_cromer), ("Verlet", verlet)]:
    x, v = 1.0, 0.0
    for _ in range(n):
        x, v = step(x, v, dt)
    drift = (energy(x, v) - 0.5) / 0.5 * 100
    # reversibility test: flip the velocity, run n steps, see if we get home
    v = -v
    for _ in range(n):
        x, v = step(x, v, dt)
    print(f"{name:13s} energy change {drift:+10.2f} %   back at start? x = {x:+.6f} (start 1.0)")`,
  },

  maths: () => `
<h4>1. Where Verlet comes from: Taylor expansion</h4>
${tex(String.raw`\mathbf{r}(t+\Delta t) = \mathbf{r}(t) + \mathbf{v}(t)\,\Delta t + \tfrac12\,\mathbf{a}(t)\,\Delta t^2 + \mathcal{O}(\Delta t^3)`, true)}
${tex(String.raw`\mathbf{v}(t+\Delta t) = \mathbf{v}(t) + \tfrac12\left[\mathbf{a}(t) + \mathbf{a}(t+\Delta t)\right]\Delta t + \mathcal{O}(\Delta t^3)`, true)}
<p>These are exactly the three moves (half kick, move, half kick). The error after a fixed simulated time scales as <b>Δt²</b>, compared with Δt for Euler. Halving Δt makes Verlet four times more accurate.</p>
<h4>2. Why Euler's energy grows</h4>
<p>For a harmonic oscillator (a = −ω²x), one Euler step multiplies the energy by a constant factor larger than 1:</p>
${tex(String.raw`E_{n+1} = \left(1 + \omega^2\Delta t^2\right) E_n \quad\Rightarrow\quad E_n = E_0\left(1+\omega^2\Delta t^2\right)^n`, true)}
<p>With ωΔt = 0.1 and 1000 steps, the energy grows by a factor of about 2 × 10⁴. That is what the code tab shows.</p>
<h4>3. Time reversibility</h4>
<p>Replace Δt by −Δt in a Verlet step and you get exactly the inverse step. So flipping all velocities retraces the trajectory, apart from rounding errors. Euler does not have this symmetry.</p>
<h4>4. Symplectic: no long-term drift</h4>
<p>Verlet is <b>symplectic</b>: it exactly conserves a slightly modified "shadow" energy that stays within about (ωΔt)² of the true one. So the true energy wiggles but cannot drift away, however long you run.</p>
<p class="muted">The equivalent position (Störmer) Verlet form, without velocities: r(t+Δt) = 2r(t) − r(t−Δt) + a(t)Δt².</p>`,

  quiz: [
    {
      q: 'Why does the explicit Euler orbit spiral outward?',
      options: ['The force is computed wrongly', 'Each step moves along the old velocity, the tangent, which lands slightly outside the orbit, so energy grows', 'Rounding errors in the computer', 'The time step is always too small'],
      answer: 1,
      why: 'Euler uses the velocity and force from the start of the step. Stepping along the tangent always overshoots outward, so the energy grows by a factor (1 + ω²Δt²) every step.',
    },
    {
      q: 'How many force calculations does velocity Verlet need per step?',
      options: ['One', 'Two', 'Three', 'Four'],
      answer: 0,
      why: 'The force computed at the new position is used for the second half kick and reused for the first half kick of the next step. So it costs the same as Euler.',
    },
    {
      q: 'You run Verlet for 1000 steps, flip all velocities and run 1000 more. Where do the atoms end up?',
      options: ['Somewhere random', 'Back where they started (to rounding error)', 'Twice as far away', 'At rest'],
      answer: 1,
      why: 'Velocity Verlet is time-reversible: a step with −Δt exactly undoes a step with +Δt, so the trajectory retraces itself.',
    },
  ],
};
