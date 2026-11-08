import { tex } from '../tex.js';

export default {
  id: 'bonds-angles',
  summary: 'Atoms in a molecule are held together by bonds and angles, which force fields model as springs. Learn the potentials, how to turn them into forces on each atom, and why their fast vibrations set the time step.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/bonds-angles.mp4', poster: 'videos/bonds-angles.jpg' },
  sim: 'bonds',
  learn: [
    'Harmonic (and Morse) bond potentials, and harmonic angle potentials',
    'Deriving the forces on each atom, and checking them numerically',
    'Vibrational periods, and why flexible O–H or C–H bonds need Δt ≤ 1 fs',
  ],
  simHint: 'Drag one H away from the O and let go: watch the bond (cyan) and angle (amber) forces. Then raise the time step from 0.5 fs to 2 and 3 fs and watch the total energy change. Switch to Morse bonds and pull an H far away: the bond breaks. With harmonic bonds it never can.',

  idea: `
<p>So far our atoms only interacted through Lennard-Jones forces. In a molecule, some atoms are <b>chemically bonded</b>. A force field describes this with simple terms that depend only on the internal geometry (Lesson 4.3 puts them all together):</p>
<ul>
<li><b>Bonds</b> (two atoms): the energy rises when the distance r moves away from its equilibrium value r₀. Near r₀ a parabola, U = ½k<sub>b</sub>(r − r₀)², is accurate. Real bonds can break; the <b>Morse</b> potential captures that, but most force fields keep the harmonic form because bonds do not break at normal temperatures.</li>
<li><b>Angles</b> (three atoms, i–j–k): the energy rises when the angle θ at the central atom j moves from θ₀: U = ½k<sub>θ</sub>(θ − θ₀)².</li>
</ul>
<p><b>Forces.</b> The force is minus the gradient of the energy with respect to each atom's position. For a bond the force acts along the bond and pulls the atoms back towards r₀. For an angle the forces on the outer atoms act perpendicular to their bonds, opening or closing the angle; the central atom gets minus their sum. In every case the forces add up to zero (Newton's third law). Always test analytic forces against a numerical derivative of the energy; the code tab does this.</p>
<p><b>Fast vibrations.</b> Bonds are stiff and hydrogen is light, so an O–H or C–H bond vibrates with a period of only about 10 fs. With 10–20 steps per vibration needed (Lesson 1.4), flexible bonds to hydrogen force Δt ≈ 0.5–1 fs. This is why most simulations constrain these bonds (Lesson 4.4) and use 2 fs.</p>
<div class="callout"><b>Key idea.</b> Bonds and angles are springs on the internal coordinates. Their forces follow from the chain rule; check them numerically. The stiffest, lightest spring (X–H bonds) sets the largest stable time step.</div>`,

  algorithm: `
<p>Forces from a harmonic bond (i–j) and angle (i–j–k):</p>
<pre class="pseudo"># bond
d = r_i - r_j;  r = |d|
F_i = -k_b (r - r0) d / r;   F_j = -F_i

# angle at j
a = r_i - r_j;  b = r_k - r_j
cos_t = a·b / (|a||b|);  theta = arccos(cos_t)
g = k_theta (theta - theta0) / sin(theta)
F_i = g * ( b/(|a||b|) - cos_t a/|a|² )
F_k = g * ( a/(|a||b|) - cos_t b/|b|² )
F_j = -F_i - F_k</pre>
<ol class="steps">
  <li>Bonded atoms (1–2) and atoms sharing an angle (1–3) are usually <b>excluded</b> from the Lennard-Jones and Coulomb interactions. Their interaction is already in the bonded terms.</li>
  <li>Watch the convention! GROMACS writes ½k(r − r₀)². AMBER and LAMMPS write K(r − r₀)², so K = k/2. Units also differ: kJ mol⁻¹ nm⁻² vs kcal mol⁻¹ Å⁻².</li>
  <li>θ near 180° makes sin θ → 0 in the formula. Linear molecules need special treatment (or a cosine-based angle potential).</li>
  <li>Check: finite-difference forces must match the analytic ones to many digits, and all forces must sum to zero.</li>
</ol>`,

  code: {
    note: 'A flexible water model (flexible SPC parameters) in real units: nm, ps, amu, kJ/mol. It checks the analytic forces against numerical derivatives, predicts the O–H period, then runs the molecule in vacuum with four time steps and measures the vibration frequency. Run <code>pip install numpy</code>, then <code>python bonds_angles.py</code> (about 1 s).',
    src: `import numpy as np
# Units (GROMACS style): nm, ps, amu, kJ/mol.  1 fs = 0.001 ps.
mO, mH = 15.999, 1.008
kb, r0 = 345000.0, 0.1                    # O-H bond: kJ/mol/nm^2, nm   (flexible SPC)
kt, t0 = 383.0, np.deg2rad(109.47)        # H-O-H angle: kJ/mol/rad^2, rad

def bond(ri, rj, k, b0):
    """Harmonic bond U = 1/2 k (r - b0)^2. Returns U, force on i (force on j is minus this)."""
    d = ri - rj; r = np.linalg.norm(d)
    return 0.5 * k * (r - b0)**2, -k * (r - b0) * d / r

def angle(ri, rj, rk, k, th0):
    """Harmonic angle U = 1/2 k (theta - th0)^2 at the central atom j. Returns U, F_i, F_j, F_k."""
    a, b = ri - rj, rk - rj; la, lb = np.linalg.norm(a), np.linalg.norm(b)
    c = np.dot(a, b) / (la * lb); th = np.arccos(np.clip(c, -1, 1))
    dU = k * (th - th0); s = np.sin(th)
    Fi = dU / s * (b / (la * lb) - c * a / la**2)        # -dU/dtheta * dtheta/dr_i
    Fk = dU / s * (a / (la * lb) - c * b / lb**2)
    return 0.5 * k * (th - th0)**2, Fi, -Fi - Fk, Fk

def energy_forces(x):                     # x[0] = O, x[1], x[2] = H
    U1, F1 = bond(x[1], x[0], kb, r0); U2, F2 = bond(x[2], x[0], kb, r0)
    U3, Fa, Fo, Fb = angle(x[1], x[0], x[2], kt, t0)
    F = np.array([-F1 - F2 + Fo, F1 + Fa, F2 + Fb])
    return U1 + U2 + U3, F

# ---------- 1. check the analytic forces against a numerical derivative ----------
x = np.array([[0, 0, 0], [0.105, 0.01, 0], [-0.03, 0.098, 0.01]])
U, F = energy_forces(x); h = 1e-6; Fnum = np.zeros_like(x)
for i in range(3):
    for d in range(3):
        xp, xm = x.copy(), x.copy(); xp[i, d] += h; xm[i, d] -= h
        Fnum[i, d] = -(energy_forces(xp)[0] - energy_forces(xm)[0]) / (2 * h)
print(f"1. analytic vs numerical forces: max difference = {np.abs(F - Fnum).max():.2e} kJ/mol/nm  (forces ~ {np.abs(F).max():.0f})")
print(f"   sum of forces = {np.abs(F.sum(0)).max():.1e}  (Newton's third law)")

# ---------- 2. one O-H bond alone: period from theory ----------
mu = mO * mH / (mO + mH); period = 2 * np.pi * np.sqrt(mu / kb)
print(f"\\n2. O-H bond: reduced mass {mu:.3f} amu, period 2*pi*sqrt(mu/k) = {period * 1000:.1f} fs"
      f"  -> wavenumber {1 / (period * 1e-12 * 2.998e10):.0f} cm^-1")

# ---------- 3. a vibrating water molecule: NVE with different time steps ----------
def run(dt_fs, t_ps=2.0):
    dt = dt_fs / 1000; m = np.array([mO, mH, mH])[:, None]
    x = np.array([[0, 0, 0], [0.11, 0, 0], [0.1 * np.cos(t0), 0.1 * np.sin(t0), 0]])   # one bond stretched
    v = np.zeros_like(x); U, F = energy_forces(x); E, roh = [], []
    for _ in range(int(t_ps / dt)):
        v += 0.5 * dt * F / m; x += dt * v; U, F = energy_forces(x); v += 0.5 * dt * F / m
        E.append(U + 0.5 * (m * v**2).sum()); roh.append(np.linalg.norm(x[1] - x[0]))
    return np.array(E), np.array(roh), dt

print("\\n3. flexible water in vacuum, velocity Verlet")
for dt_fs in [0.25, 0.5, 1.0, 2.0]:
    E, roh, dt = run(dt_fs)
    print(f"   dt = {dt_fs:4.2f} fs: energy spread = {E.std() / E.mean() * 100:7.3f} %")
E, roh, dt = run(0.25)
spec = np.abs(np.fft.rfft(roh - roh.mean())); freq = np.fft.rfftfreq(len(roh), dt)   # 1/ps
peaks = freq[np.argsort(spec)[-1]] / 2.998e-2                                         # 1/ps -> cm^-1
print(f"   strongest O-H vibration in the bond-length signal: {peaks:.0f} cm^-1")
print("A C-H or O-H bond vibrates in ~10 fs, so flexible bonds need dt <= 1 fs (or constraints, Lesson 4.4).")
`,
  },

  maths: () => `
<h4>1. Bond and angle potentials</h4>
${tex(String.raw`U_{\text{bond}} = \tfrac12 k_b (r - r_0)^2, \qquad U_{\text{Morse}} = D_e\big(1 - e^{-a(r - r_0)}\big)^2, \qquad U_{\text{angle}} = \tfrac12 k_\theta(\theta - \theta_0)^2`, true)}
<p>Near r₀ the Morse potential is harmonic with k<sub>b</sub> = 2D<sub>e</sub>a².</p>
<h4>2. Bond force</h4>
${tex(String.raw`\mathbf F_i = -\frac{\partial U}{\partial \mathbf r_i} = -k_b (r - r_0)\,\frac{\mathbf r_i - \mathbf r_j}{r}, \qquad \mathbf F_j = -\mathbf F_i`, true)}
<h4>3. Angle force (chain rule)</h4>
${tex(String.raw`\cos\theta = \frac{\mathbf a\cdot\mathbf b}{|\mathbf a||\mathbf b|}, \qquad \frac{\partial\theta}{\partial \mathbf r_i} = -\frac{1}{\sin\theta}\left(\frac{\mathbf b}{|\mathbf a||\mathbf b|} - \cos\theta\,\frac{\mathbf a}{|\mathbf a|^2}\right)`, true)}
${tex(String.raw`\mathbf F_i = -k_\theta(\theta - \theta_0)\,\frac{\partial\theta}{\partial \mathbf r_i}, \qquad \mathbf F_j = -\mathbf F_i - \mathbf F_k`, true)}
<h4>4. Vibrational frequency</h4>
${tex(String.raw`\omega = \sqrt{\frac{k_b}{\mu}}, \quad \mu = \frac{m_1 m_2}{m_1 + m_2}, \qquad \tilde\nu = \frac{\omega}{2\pi c}`, true)}
<p>For O–H with k<sub>b</sub> = 345,000 kJ mol⁻¹ nm⁻² and μ = 0.948 amu: period 10.4 fs, 3200 cm⁻¹ (experiment: about 3657 cm⁻¹ in the gas phase).</p>`,

  quiz: [
    {
      q: 'A paper lists a bond constant of 450 kcal/mol/Å² in the AMBER convention, U = K(r − r₀)². What is k in the GROMACS convention, U = ½k(r − r₀)², in kcal/mol/Å²?',
      options: ['225', '450', '900', '4500'],
      answer: 2,
      why: 'K(r − r₀)² = ½k(r − r₀)² means k = 2K = 900 kcal/mol/Å². (Then convert units: × 4.184 × 100 to get kJ/mol/nm².)',
    },
    {
      q: 'In a harmonic angle term i–j–k, what is the sum of the three forces?',
      options: ['Zero', 'Along the bisector of the angle', 'Equal to the force on the central atom', 'It depends on θ₀'],
      answer: 0,
      why: 'The energy depends only on the internal geometry, so it does not change if the molecule is translated: the forces sum to zero.',
    },
    {
      q: 'Why do flexible O–H bonds need a time step of about 1 fs or less?',
      options: ['Oxygen is heavy', 'The O–H vibration period is about 10 fs, and we need 10–20 steps per period', 'Hydrogen bonds are long-ranged', 'Angle forces are large'],
      answer: 1,
      why: 'Stiff spring plus light hydrogen gives a fast vibration. Constraining X–H bonds (Lesson 4.4) removes it and allows 2 fs.',
    },
  ],
};
