import { tex } from '../tex.js';

export default {
  id: 'energy-check',
  summary: 'An MD code can run happily and still be wrong. The cheapest, most powerful test is energy conservation: with no thermostat, the total energy must stay flat. If it drifts or jumps, something is broken.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/energy-check.mp4', poster: 'videos/energy-check.jpg' },
  sim: 'energylab',
  learn: [
    'Why total energy is the first thing to check in any MD run',
    'How to measure drift and the fluctuation ratio',
    'The three usual culprits: time step, forces that do not match the energy, and cutoff jumps',
  ],
  simHint: 'It starts healthy. Break it in different ways, one at a time: raise Δt, untick "Shift potential" and lower the cutoff, or switch on the force bug. Press Force check each time. Which problems does the force check catch, and which does only the energy plot catch?',

  idea: `
<p>With no thermostat (the <b>NVE</b> ensemble), Newton's laws conserve the total energy E = K + U. Kinetic and potential energy trade back and forth, but their <b>sum must stay constant</b>.</p>
<p>A numerical simulation never conserves it perfectly. A healthy run shows only small <b>fluctuations</b> around a flat line. What you must never see is <b>drift</b> (a steady climb or fall) or <b>jumps</b>. These are the usual culprits:</p>
<ol class="steps">
  <li><b>Time step too big.</b> The integrator cannot follow the fastest motion (Lesson 1.4).</li>
  <li><b>Forces that do not match the energy.</b> The force must be exactly minus the slope of the potential. A typo, a missing factor or a wrong unit breaks this.</li>
  <li><b>Jumps at the cutoff.</b> If V(r) is simply cut off at r<sub>c</sub>, a pair crossing the cutoff suddenly changes the energy. Shifting the potential so that V(r<sub>c</sub>) = 0 removes the jump.</li>
</ol>
<p>Other common causes are thermostats left on, broken periodic boundaries, or single-precision arithmetic in very long runs.</p>
<div class="callout"><b>Key idea.</b> Before any production run, do a short NVE test and plot the total energy. Also check your forces against a finite difference of the energy. These two tests catch most MD bugs.</div>`,

  algorithm: `
<p><b>Test 1: energy conservation (NVE).</b></p>
<ol class="steps">
  <li>Turn the thermostat and barostat off.</li>
  <li>Run a few thousand steps and record K, U and E = K + U.</li>
  <li>Fit a straight line to E/N versus time. Its slope is the <b>drift</b>.</li>
  <li>Compute the <b>fluctuation ratio</b> σ(E)/σ(K). A healthy run gives roughly 0.01 or less.</li>
</ol>
<p><b>Test 2: forces versus energy (finite difference).</b></p>
<pre class="pseudo">h = 1e-5
x[i] += h;      U_plus  = energy(x)
x[i] -= 2*h;    U_minus = energy(x)
x[i] += h
F_numeric = -(U_plus - U_minus) / (2*h)
assert abs(F[i] - F_numeric) < 1e-6 * abs(F_numeric)</pre>
<p>The force test catches coding errors in the forces. It does <i>not</i> catch cutoff jumps or a time step that is too big. Those only show up in the energy test, so you need both.</p>`,

  code: {
    note: 'The LJ code from Lesson 1.1, with both tests: a finite-difference force check, then NVE runs with and without the cutoff shift and at larger time steps. Run <code>pip install numpy</code>, then <code>python energy_check.py</code> (it takes about a minute).',
    src: `import numpy as np

N, L, rc = 64, 6.0, 2.5                 # 64 LJ atoms (reduced units)

def energy_forces(pos, shift=True):
    d = pos[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)            # minimum image
    r2 = (d ** 2).sum(-1)
    np.fill_diagonal(r2, np.inf)
    inside = r2 < rc**2
    r2 = np.where(inside, r2, np.inf)
    inv6 = 1 / r2**3
    vc = (rc**-12 - rc**-6) if shift else 0.0
    U = 2 * (inv6**2 - inv6 - vc)[inside].sum()          # sum over pairs
    F = ((24 * (2 * inv6**2 - inv6) / r2)[..., None] * d).sum(axis=1)
    return U, F

def start(T=1.0, seed=1):
    n = round(N ** (1/3)); g = np.arange(n) * (L / n)
    pos = np.array(np.meshgrid(g, g, g)).reshape(3, -1).T
    vel = np.random.default_rng(seed).normal(0, np.sqrt(T), (N, 3))
    return pos, vel - vel.mean(axis=0)

# --- check 1: do the forces match the energy?  F = -dU/dx ---
pos, _ = start()
pos += np.random.default_rng(2).normal(0, 0.05, pos.shape)    # off the perfect lattice
U0, F = energy_forces(pos)
h, i, k = 1e-5, 0, 0                    # nudge atom 0 along x
p = pos.copy(); p[i, k] += h; Up, _ = energy_forces(p)
p[i, k] -= 2 * h;            Um, _ = energy_forces(p)
print(f"force check: analytic {F[i, k]:.6f}   numeric {-(Up - Um) / (2 * h):.6f}")

# --- check 2: is the total energy conserved in NVE? ---
def run(dt, shift, steps=4000):
    pos, vel = start()
    U, F = energy_forces(pos, shift)
    E = []
    for _ in range(steps):
        vel += 0.5 * dt * F
        pos = (pos + dt * vel) % L
        U, F = energy_forces(pos, shift)
        vel += 0.5 * dt * F
        K = 0.5 * (vel ** 2).sum()
        E.append((K, U))
    K, U = np.array(E).T
    Etot = K + U
    t = np.arange(steps) * dt
    drift = np.polyfit(t, Etot / N, 1)[0]                 # energy per atom per time unit
    ratio = Etot.std() / K.std()                         # small = healthy
    print(f"dt={dt:<6} shift={str(shift):5}  drift {drift:+.2e} per atom per time   fluctuation ratio {ratio:.3f}")

run(0.005, True)
run(0.005, False)    # unshifted cutoff: energy jumps when pairs cross rc
run(0.010, True)
run(0.020, True)     # time step too big`,
  },

  maths: () => `
<h4>1. Why energy is conserved</h4>
<p>For forces that come from a potential, F<sub>i</sub> = −∇<sub>i</sub>U, the total energy has zero time derivative:</p>
${tex(String.raw`\frac{dE}{dt} = \sum_i m_i \mathbf{v}_i\cdot\dot{\mathbf{v}}_i + \sum_i \nabla_i U\cdot\mathbf{v}_i = \sum_i \mathbf{v}_i\cdot\left(\mathbf{F}_i + \nabla_i U\right) = 0`, true)}
<p>This only works if F<sub>i</sub> really equals −∇<sub>i</sub>U. That is why a force bug shows up as energy drift.</p>
<h4>2. The cutoff jump and the shift</h4>
<p>Truncating V at r<sub>c</sub> makes U jump by V(r<sub>c</sub>) every time a pair crosses the cutoff. The shifted potential removes the jump:</p>
${tex(String.raw`V_{\text{sh}}(r) = \begin{cases} V(r) - V(r_c) & r < r_c \\ 0 & r \ge r_c \end{cases}`, true)}
<p>The force is unchanged, but it still has a small step at r<sub>c</sub>. "Shifted-force" or switching functions smooth that too.</p>
<h4>3. Two numbers to report</h4>
${tex(String.raw`\text{drift} = \frac{1}{N}\,\frac{d\langle E\rangle}{dt}, \qquad R = \frac{\sigma(E)}{\sigma(K)} \lesssim 0.01`, true)}
<p>Drift is often quoted in kJ mol⁻¹ ns⁻¹ per atom, or as k<sub>B</sub>T per ns per degree of freedom. R compares the noise in E with the natural exchange between K and U.</p>
<h4>4. Finite-difference force check</h4>
${tex(String.raw`F_{i,x} \approx -\frac{U(x_i + h) - U(x_i - h)}{2h}, \qquad \text{error} \sim \mathcal{O}(h^2)`, true)}
<p class="muted">Choose h around 10⁻⁵ (in σ units) in double precision. Much smaller and rounding error dominates.</p>`,

  quiz: [
    {
      q: 'In an NVE run, the total energy rises steadily over time. Which is NOT a likely cause?',
      options: ['The time step is too large', 'The forces do not match the potential', 'The system has too many atoms', 'A thermostat was accidentally left on'],
      answer: 2,
      why: 'System size does not break conservation. Too big a Δt, inconsistent forces, a thermostat, cutoff jumps or broken boundaries do.',
    },
    {
      q: 'Why do we shift the LJ potential so that V(rc) = 0?',
      options: ['To make the atoms bigger', 'To remove the energy jump when a pair crosses the cutoff', 'To make the force stronger', 'To speed up the force loop'],
      answer: 1,
      why: 'Without the shift, U jumps by V(rc) each time a pair crosses rc. That shows up as noise and spikes in the total energy.',
    },
    {
      q: 'Your finite-difference force check passes, but the total energy still jumps. What is the most likely cause?',
      options: ['A typo in the force formula', 'An unshifted cutoff', 'Wrong atom masses in the force', 'Nothing: that is normal'],
      answer: 1,
      why: 'A force typo would fail the force check. Cutoff jumps happen only when pairs cross rc, which a small finite-difference nudge does not see. Only the energy test catches them.',
    },
  ],
};
