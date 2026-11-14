import { tex } from '../tex.js';

export default {
  id: 'verlet',
  summary: 'There are many ways to turn F = ma into steps. The obvious one, Euler, slowly pumps energy into the system. Velocity Verlet costs the same but keeps energy steady and can even run backwards. Leapfrog, used by GROMACS, moves the atoms along exactly the same path; only the way it reports the energy differs.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/euler-vs-verlet.mp4', poster: 'videos/euler-vs-verlet.jpg' },
  sim: 'integrators',
  learn: [
    'Why the simple Euler update makes energy grow',
    'How velocity Verlet works: half kick, move, half kick',
    'Two properties that make Verlet the MD standard: no energy drift and time reversibility',
    'Leapfrog: the same positions as Verlet, with velocities at half steps, and what that does to the reported energy',
  ],
  simHint: 'Watch the three dots orbit the same spring. Red (Euler) spirals out, and its energy climbs. Let it run, then press Reverse velocities: only Verlet and leapfrog come back to ◎. The green leapfrog dot hides under the blue Verlet dot all the time: same path. Compare their energy readouts, switch the leapfrog kinetic-energy choice, and try a bigger Δt.',

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
<p><b>Leapfrog</b> is the same idea written with velocities at <b>half steps</b>: kick the velocity a whole step, from v(t − Δt/2) to v(t + Δt/2), then move the atom with it. Starting from the same point, it puts every atom at <b>exactly the same positions</b> as velocity Verlet. Its catch: it never has the velocity at time t, so it must <b>estimate the kinetic energy</b>, and different estimates report slightly different energies (section 5 of the maths). Neither estimate drifts.</p>
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
<p><b>Leapfrog</b> (GROMACS <code>integrator = md</code>) keeps the velocity half a step ahead of the position:</p>
<pre class="pseudo">v_half = vel - 0.5 * dt * F / mass      # v(−dt/2), once, before the loop
for step in range(n_steps):
    F = compute_forces(pos)
    v_new = v_half + dt * F / mass      # whole kick: v(t − dt/2) → v(t + dt/2)
    KE = 0.5 * mass * (v_half**2 + v_new**2) / 2   # GROMACS: average of the two half-step KEs
    pos += dt * v_new                   # move with the half-step velocity
    v_half = v_new</pre>
<p>Positions are the same as velocity Verlet's (the average of the two half-step velocities is exactly Verlet's v(t)). What is not fixed is the kinetic energy at time t: averaging the two half-step <i>velocities</i> reproduces Verlet's energy; averaging the two half-step <i>kinetic energies</i>, as GROMACS does, reports a slightly different one.</p>`,

  code: {
    note: 'Leapfrog\'s first three steps and both of its energies, then the three update rules on the same oscillator for 1000 steps, followed by a reversibility test: flip the velocity, run 1000 more steps, and check whether you return to x = 1. Run <code>pip install numpy</code>, then <code>python verlet.py</code>.',
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

# leapfrog: velocities at half steps. Start: v(-dt/2) and v(+dt/2) from x0 = 1, v0 = 0
x, vm, vp = 1.0, 0.0 - 0.5 * accel(1.0) * dt, 0.0 + 0.5 * accel(1.0) * dt
for i in range(1, n + 1):
    x = x + vp * dt                       # move with v(t + dt/2)
    vm, vp = vp, vp + accel(x) * dt       # whole kick
    e_v  = 0.5 * ((vm + vp) / 2)**2 + 0.5 * x**2     # KE of the average velocity (= Verlet)
    e_ke = 0.25 * (vm**2 + vp**2) + 0.5 * x**2       # average of the two KEs (GROMACS)
    if i <= 3:
        print(f"leapfrog step {i}: x = {x:.10g}  E(v avg) = {e_v:.8f}  E(KE avg) = {e_ke:.8f}")
print(f"leapfrog after {n} steps: E(v avg) {(e_v-0.5)/0.5*100:+.4f} %   E(KE avg) {(e_ke-0.5)/0.5*100:+.4f} %")

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
<p class="muted">The equivalent position (Störmer) Verlet form, without velocities: r(t+Δt) = 2r(t) − r(t−Δt) + a(t)Δt².</p>
<h4>5. Leapfrog by hand: same positions, two energies</h4>
<p>The code-tab oscillator (a = −x, x₀ = 1, v₀ = 0, E₀ = 0.5) with Δt = 0.1. Leapfrog starts from v(−Δt/2) = v₀ − ½Δt·a₀ = +0.05 and v(+Δt/2) = −0.05. Step 1: x₁ = 1 + 0.1 × (−0.05) = 0.995, then v(3Δt/2) = −0.05 + 0.1 × (−0.995) = −0.1495.</p>
<div class="tbl"><table><tr><th>step</th><th>x (Verlet)</th><th>x (leapfrog)</th><th>E Verlet</th><th>E leapfrog, KE of average v</th><th>E leapfrog, average of KEs</th></tr>
<tr><td>1</td><td>0.995</td><td>0.995</td><td>0.49998753</td><td>0.49998753</td><td>0.50122506</td></tr>
<tr><td>2</td><td>0.98005</td><td>0.98005</td><td>0.49995062</td><td>0.49995062</td><td>0.50115125</td></tr>
<tr><td>3</td><td>0.9552995</td><td>0.9552995</td><td>0.49989075</td><td>0.49989075</td><td>0.50103149</td></tr></table></div>
<p>The positions agree to every digit. The two leapfrog energies differ by exactly</p>
${tex(String.raw`\tfrac12 m\,\frac{v_-^2 + v_+^2}{2} - \tfrac12 m\left(\frac{v_- + v_+}{2}
ight)^2 = \frac{m}{8}\,(v_+ - v_-)^2 = \frac{m\,a^2\,\Delta t^2}{8}`, true)}
<p>Step 1: 0.1² × 0.995² / 8 = 0.0012375, the difference in the table. It is always positive and shrinks as Δt², and it follows a(t), so it goes up and down with the motion instead of piling up.</p>
<p>Over 1000 periods (62 832 steps) neither drifts:</p>
<div class="tbl"><table><tr><th>Δt</th><th>leapfrog energy</th><th>largest |E − E₀|</th><th>mean E − E₀, first 10 %</th><th>last 10 %</th></tr>
<tr><td>0.1</td><td>KE of average v (= Verlet)</td><td>1.25 × 10⁻³</td><td>−6.25 × 10⁻⁴</td><td>−6.25 × 10⁻⁴</td></tr>
<tr><td>0.1</td><td>average of KEs (GROMACS)</td><td>1.25 × 10⁻³</td><td>+4.5 × 10⁻⁷</td><td>+1.0 × 10⁻⁷</td></tr>
<tr><td>0.05</td><td>KE of average v (= Verlet)</td><td>3.1 × 10⁻⁴</td><td>−1.56 × 10⁻⁴</td><td>−1.56 × 10⁻⁴</td></tr>
<tr><td>0.05</td><td>average of KEs (GROMACS)</td><td>3.1 × 10⁻⁴</td><td>+2.3 × 10⁻⁸</td><td>+7.5 × 10⁻⁹</td></tr></table></div>
<p>For this spring the averaged-KE energy wobbles symmetrically around E₀, while Verlet's sits slightly below it; both stay put. In a molecule with stiff bonds the averaged-KE wobble is larger (a test on butane with bonds, angles and a dihedral gave about twice Verlet's), but again without drift.</p>
<p><b>Where leapfrog codes do drift.</b> With bond <b>constraints</b> (SHAKE, LINCS, SETTLE: Lesson 4.4), many codes recompute the half-step velocity from the constrained positions, v = (x_new − x_old)/Δt. In single precision this loses digits and makes the energy drift, as shown for GROMACS and AMBER by Eastwood et al., <a href="https://www.deshawresearch.com/publications/A%20common,%20avoidable%20source%20of%20error%20in%20molecular%20dynamics%20integrators.pdf" target="_blank" rel="noopener">J. Chem. Phys. 132, 2010</a>. That is a property of how constraints are added, not of leapfrog itself.</p>`,

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
    {
      q: 'GROMACS (leapfrog) and LAMMPS (velocity Verlet) run the same NVE system from the same start with the same Δt and no constraints. What differs?',
      options: ['The positions: leapfrog is less accurate', 'Nothing at all', 'The reported kinetic energy: leapfrog only has half-step velocities and must estimate it', 'Leapfrog drifts in energy, Verlet does not'],
      answer: 2,
      why: 'The positions are the same (apart from rounding). Leapfrog estimates the kinetic energy at time t from v(t − Δt/2) and v(t + Δt/2); averaging the two kinetic energies, as GROMACS does, adds m·a²Δt²/8. Neither drifts.',
    },
  ],
};
