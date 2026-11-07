import { tex } from '../tex.js';

export default {
  id: 'newton-in-code',
  summary: 'A force alone does not move an atom: it changes the velocity, and the velocity changes the position. In code, that becomes two lines repeated over and over.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/newton-in-code.mp4', poster: 'videos/newton-in-code.jpg' },
  sim: 'newton',
  learn: [
    'How F = ma turns a force into motion',
    'The two update lines at the heart of every MD step',
    'How to store and update many atoms at once with arrays',
  ],
  simHint: 'Press Step a few times and watch the numbers: force, then acceleration, then velocity, then position. Then press Run. Switch to Constant push and compare a heavy atom with a light one. Finally, set the mass to 0.25 and Δt to 1.2 on the spring.',

  idea: `
<p>In Lesson 1.2 we found the force between atoms. Now we need to turn that force into <b>motion</b>.</p>
<p>Newton's second law says <b>F = ma</b>, or <b>a = F/m</b>. A force gives an atom an <b>acceleration</b>: a change in its velocity. The same push changes a light atom's velocity a lot and a heavy atom's velocity only a little.</p>
<p>A computer cannot follow time smoothly, so it moves forward in tiny jumps of length <b>Δt</b>. In each jump:</p>
<ol class="steps">
  <li>Compute the <b>force</b> on the atom.</li>
  <li>Get the <b>acceleration</b>: a = F/m.</li>
  <li>Update the <b>velocity</b>: new v = v + a·Δt.</li>
  <li>Update the <b>position</b>: new x = x + v·Δt.</li>
</ol>
<p>Repeat this many times and the dots join up into a <b>trajectory</b>. In 3D, each position, velocity and force has three numbers (x, y, z), so for N atoms we store N × 3 arrays and update them all together.</p>
<div class="callout"><b>Key idea.</b> Force changes velocity; velocity changes position. Everything in MD is this pair of updates, done for every atom at every step.</div>`,

  algorithm: `
<p>One MD step for N atoms, written with arrays of shape (N, 3):</p>
<pre class="pseudo">F = compute_forces(pos)      # (N, 3): from Lesson 1.2
a = F / mass                 # mass has shape (N, 1)
vel = vel + a * dt           # update all velocities
pos = pos + vel * dt         # move all atoms</pre>
<ol class="steps">
  <li><b>Order matters.</b> Here we update the velocity first and then use the <i>new</i> velocity to move. This version (semi-implicit Euler) is simple and surprisingly stable.</li>
  <li><b>Arrays, not loops.</b> One line updates every atom. In NumPy this is much faster than a Python loop over atoms.</li>
  <li><b>Mass as a column.</b> Storing masses as an (N, 1) column lets each atom's force row be divided by its own mass.</li>
</ol>
<p>The step size Δt is a trade-off: small steps follow the true path closely but cost more steps. Lesson 1.4 looks at how to choose it, and Lesson 1.5 replaces this update with the more accurate velocity Verlet.</p>`,

  code: {
    note: 'First one atom on a spring in 2D, then the same update applied to many atoms at once. Run <code>pip install numpy</code>, then <code>python newton.py</code>. In the last output, the two heavy atoms (mass 4) change velocity four times less.',
    src: `import numpy as np

# --- one atom on a spring (like a bond to a fixed atom), in 2D ---
k, m, dt = 1.0, 1.0, 0.05
x = np.array([1.0, 0.0])      # position
v = np.array([0.0, 0.8])      # velocity

def force(x):
    return -k * x             # Hooke's law: always pulls back toward the origin

for step in range(201):
    a = force(x) / m          # Newton: a = F / m
    v = v + a * dt            # 1. update the velocity
    x = x + v * dt            # 2. move the atom
    if step % 40 == 0:
        print(f"step {step:3d}   x = {x.round(3)}   v = {v.round(3)}")

# --- many atoms: the same two lines work on whole arrays at once ---
N = 4
rng = np.random.default_rng(0)
pos = rng.random((N, 3))                  # N rows, columns = x, y, z
vel = np.zeros((N, 3))
mass = np.array([1.0, 1.0, 4.0, 4.0])[:, None]   # shape (N, 1) so it divides each row

F = force(pos)                            # (N, 3) forces, one row per atom
vel += F / mass * dt                      # update every velocity in one line
pos += vel * dt                           # move every atom in one line
print("velocities after one step:\\n", vel.round(4))`,
  },

  maths: () => `
<h4>1. Newton's second law as two equations</h4>
${tex(String.raw`\frac{d\mathbf{v}_i}{dt} = \mathbf{a}_i = \frac{\mathbf{F}_i}{m_i}, \qquad \frac{d\mathbf{r}_i}{dt} = \mathbf{v}_i`, true)}
<p>Every atom i has its own position r<sub>i</sub>, velocity v<sub>i</sub> and mass m<sub>i</sub>. The forces F<sub>i</sub> depend on all positions, which couples the atoms together.</p>
<h4>2. From derivatives to steps</h4>
<p>Over a short time Δt, a derivative is close to "change divided by Δt". Rearranging gives the update rules:</p>
${tex(String.raw`\mathbf{v}_i(t+\Delta t) = \mathbf{v}_i(t) + \frac{\mathbf{F}_i(t)}{m_i}\,\Delta t`, true)}
${tex(String.raw`\mathbf{r}_i(t+\Delta t) = \mathbf{r}_i(t) + \mathbf{v}_i(t+\Delta t)\,\Delta t`, true)}
<p>Using the new velocity in the position update makes this the <b>semi-implicit (symplectic) Euler</b> method. Its error per unit time shrinks in proportion to Δt.</p>
<h4>3. Real units</h4>
<p>MD codes often use forces in kcal mol<sup>−1</sup> Å<sup>−1</sup>, masses in g mol<sup>−1</sup> and time in fs. Then the acceleration needs a conversion factor:</p>
${tex(String.raw`a\;[\text{Å/fs}^2] = 4.184\times10^{-4}\;\frac{F\;[\text{kcal mol}^{-1}\,\text{Å}^{-1}]}{m\;[\text{g mol}^{-1}]}`, true)}
<p class="muted">Forgetting a unit conversion like this is one of the most common bugs in a first MD code. Reduced units (Lesson 2.1) avoid it.</p>`,

  quiz: [
    {
      q: 'Atom B is 4 times heavier than atom A. Both feel the same force. How does B\'s acceleration compare?',
      options: ['4 times larger', 'The same', '4 times smaller', '16 times smaller'],
      answer: 2,
      why: 'a = F/m. With the same F and 4 times the mass, the acceleration is 4 times smaller.',
    },
    {
      q: 'In each MD step, what does the acceleration directly update?',
      options: ['The position', 'The velocity', 'The force', 'The mass'],
      answer: 1,
      why: 'Acceleration is the rate of change of velocity: new v = v + a·Δt. The position is then updated from the velocity.',
    },
    {
      q: 'What happens if the time step Δt is too large?',
      options: ['The simulation becomes more accurate', 'Nothing changes', 'The steps overshoot and drift off the true path, or even blow up', 'The atoms get heavier'],
      answer: 2,
      why: 'Each step assumes the force stays constant over Δt. With big steps that assumption fails, so the path drifts and can become unstable. Try it in the simulation with a light atom.',
    },
  ],
};
