import { tex } from '../tex.js';

export default {
  id: 'what-is-md',
  summary: 'Atoms never sit still. MD follows every one of them, step by tiny step, and turns that motion into properties you can measure.',
  minutes: 10,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v=").
  // Leave it empty to play the local MP4 instead.
  video: { youtube: '', mp4: 'videos/what-is-md.mp4', poster: 'videos/what-is-md.jpg' },
  sim: 'lj',
  learn: [
    'What a molecular dynamics simulation actually does',
    'The four-step loop at the heart of every MD code',
    'Why time steps are tiny and runs are long',
  ],
  simHint: 'Cool the target temperature to 0.3 and watch the atoms clump into liquid droplets. Heat it to 2.0 and they fly apart as a gas. Then turn the time step all the way up and see what happens.',

  idea: `
<p>Everything around you is made of atoms, and those atoms are always moving: vibrating, colliding and drifting. How they move decides whether something is a solid, a liquid or a gas, how fast ions travel in a battery, and how a protein folds.</p>
<p><b>Molecular dynamics (MD)</b> is a way to watch that motion on a computer. We put atoms in a box, work out the force on each one, and use Newton's law to move every atom forward by a very small time step. Then we do it again, millions of times.</p>
<p>The result is a <b>trajectory</b>: a movie of every atom's position over time. From it we can measure temperature, pressure, structure and diffusion, and compare them with experiments.</p>
<div class="callout"><b>Key idea.</b> MD does not guess how atoms behave. It only needs a rule for the forces between atoms. Everything else follows from Newton's law.</div>`,

  algorithm: `
<p>Every MD program, from a 40-line script to LAMMPS, runs the same loop:</p>
<ol class="steps">
  <li><b>Set up.</b> Place N atoms in a box and give them random velocities.</li>
  <li><b>Forces.</b> For every atom, add up the forces from its neighbours.</li>
  <li><b>Move.</b> Use the forces to update velocities and positions over a time step Δt.</li>
  <li><b>Measure.</b> Record energy, temperature or anything else you need.</li>
  <li><b>Repeat</b> from step 2 until you have enough data.</li>
</ol>
<pre class="pseudo">positions, velocities = set_up_box()
forces = compute_forces(positions)
for step in range(n_steps):
    velocities += 0.5 * dt * forces / mass
    positions  += dt * velocities
    forces      = compute_forces(positions)
    velocities += 0.5 * dt * forces / mass
    record(positions, velocities)</pre>
<p>That update rule is called <b>velocity Verlet</b>. You will see why it is better than the obvious approach in Lesson 1.5.</p>`,

  code: {
    note: 'A complete MD program in about 40 lines. Save it as <code>md.py</code>, run <code>pip install numpy</code>, then <code>python md.py</code>. The total energy E should stay nearly constant: that is your first check that the simulation is right.',
    src: `import numpy as np

# --- settings (reduced Lennard-Jones units: sigma = epsilon = mass = 1) ---
N, L = 64, 6.0              # number of atoms, box length
dt, steps, T = 0.005, 2000, 1.0
rc = 2.5                    # cutoff distance

# --- start: atoms on a cubic lattice, random velocities ---
n = round(N ** (1/3))
grid = np.arange(n) * (L / n)
pos = np.array(np.meshgrid(grid, grid, grid)).reshape(3, -1).T
vel = np.random.normal(0, np.sqrt(T), (N, 3))
vel -= vel.mean(axis=0)     # remove net drift

def forces(pos):
    d = pos[:, None, :] - pos[None, :, :]    # every pair vector
    d -= L * np.round(d / L)                 # minimum image (periodic box)
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)             # no self-interaction
    r2[r2 > rc**2] = np.inf                  # ignore pairs beyond the cutoff
    inv6 = 1 / r2**3
    f = (24 * (2 * inv6**2 - inv6) / r2)[..., None] * d
    inside = np.isfinite(r2)
    vc = rc**-12 - rc**-6                    # shift so V(rc) = 0
    pe = 2 * (inv6**2 - inv6 - vc)[inside].sum()
    return f.sum(axis=1), pe

# --- the MD loop (velocity Verlet) ---
f, pe = forces(pos)
for step in range(steps):
    vel += 0.5 * dt * f     # half kick
    pos += dt * vel         # drift
    pos %= L                # wrap atoms back into the box
    f, pe = forces(pos)
    vel += 0.5 * dt * f     # half kick
    if step % 250 == 0:
        ke = 0.5 * (vel ** 2).sum()
        print(f"{step:5d}  KE={ke:8.3f}  PE={pe:8.3f}  E={ke+pe:8.3f}")`,
  },

  maths: () => `
<h4>1. The force rule: Lennard-Jones</h4>
<p>Two neutral atoms a distance r apart attract weakly at long range and repel strongly when they overlap:</p>
${tex(String.raw`V(r) = 4\varepsilon\left[\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]`, true)}
<p>The force is minus the slope of the potential:</p>
${tex(String.raw`\mathbf{F}(r) = -\frac{dV}{dr}\,\hat{\mathbf{r}} = \frac{24\varepsilon}{r}\left[2\left(\frac{\sigma}{r}\right)^{12} - \left(\frac{\sigma}{r}\right)^{6}\right]\hat{\mathbf{r}}`, true)}
<h4>2. Newton's second law</h4>
${tex(String.raw`m_i\,\mathbf{a}_i = \mathbf{F}_i = \sum_{j \ne i} \mathbf{F}_{ij}`, true)}
<h4>3. Velocity Verlet update</h4>
${tex(String.raw`\mathbf{v}_i\!\left(t+\tfrac{\Delta t}{2}\right) = \mathbf{v}_i(t) + \frac{\Delta t}{2m_i}\mathbf{F}_i(t)`, true)}
${tex(String.raw`\mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) + \Delta t\,\mathbf{v}_i\!\left(t+\tfrac{\Delta t}{2}\right)`, true)}
${tex(String.raw`\mathbf{v}_i(t+\Delta t) = \mathbf{v}_i\!\left(t+\tfrac{\Delta t}{2}\right) + \frac{\Delta t}{2m_i}\mathbf{F}_i(t+\Delta t)`, true)}
<h4>4. Temperature from motion</h4>
<p>Temperature is measured from the kinetic energy (the 3 removes the box's overall drift):</p>
${tex(String.raw`T = \frac{2}{(3N-3)\,k_B}\sum_i \tfrac12 m_i v_i^2`, true)}
<p class="muted">In reduced units, σ = ε = m = k<sub>B</sub> = 1. For argon, one time unit is about 2.16 ps, so Δt = 0.005 is about 11 fs.</p>`,

  quiz: [
    {
      q: 'What does an MD program do at every step?',
      options: ['Picks random new positions for the atoms', 'Computes the force on each atom, then moves it a small step', 'Solves the Schrödinger equation for all electrons', 'Finds the lowest-energy structure'],
      answer: 1,
      why: 'MD computes forces, then uses Newton\'s law to move each atom forward by Δt. Random moves are Monte Carlo, and finding the lowest-energy structure is minimisation.',
    },
    {
      q: 'What is a typical time step for an all-atom MD simulation?',
      options: ['About 1 femtosecond', 'About 1 nanosecond', 'About 1 microsecond', 'About 1 second'],
      answer: 0,
      why: 'The fastest motions (bonds to hydrogen) vibrate in about 10 fs, so the step must be about 1 fs to follow them.',
    },
    {
      q: 'Why do MD runs need millions of steps?',
      options: ['To make the computer warm up', 'Each step covers so little time that you need millions of steps for nanoseconds of motion', 'Each step is less accurate than the one before', 'Atoms only move every thousand steps'],
      answer: 1,
      why: '10⁶ steps × 1 fs = 1 ns. Many processes, such as diffusion, need nanoseconds or more to show up.',
    },
  ],
};
