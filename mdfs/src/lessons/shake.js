import { tex } from '../tex.js';

export default {
  id: 'shake',
  summary: 'The fastest motions in a molecule, bond vibrations involving hydrogen, force a tiny time step. Constraints freeze those bond lengths. SHAKE, RATTLE, SETTLE and LINCS enforce them every step, which lets simulations use 2 fs instead of 0.5 fs.',
  minutes: 15,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/shake.mp4', poster: 'videos/shake.jpg' },
  sim: 'shake',
  learn: [
    'Why constraining X–H bonds allows larger time steps',
    'How SHAKE corrects positions and RATTLE corrects velocities (Lagrange multipliers)',
    'SETTLE and LINCS, the effect on degrees of freedom, and hydrogen mass repartitioning',
  ],
  simHint: 'Start with Flexible bonds at 0.5 fs, then raise the time step to 2 and 3 fs and watch the energy. Switch to Rigid: 2 fs and even 4 fs stay stable. Now make the SHAKE tolerance loose (10⁻²): fewer sweeps per step, but what happens to the bond lengths and the energy?',

  idea: `
<p>Lesson 4.1 showed that an O–H or C–H bond vibrates with a period of about 10 fs, which limits flexible simulations to Δt ≈ 0.5–1 fs. These vibrations are also not very realistic classically: at room temperature they are almost entirely in their quantum ground state. So most simulations simply <b>freeze</b> them: the bond length is held fixed at its equilibrium value.</p>
<p><b>How.</b> A constraint is a condition, |<b>r</b><sub>i</sub> − <b>r</b><sub>j</sub>| = d<sub>ij</sub>, that must hold at every step. It is enforced by a <b>constraint force</b> along the bond, whose size is a Lagrange multiplier chosen exactly so the condition holds. In practice:</p>
<ul>
<li><b>SHAKE</b> (Ryckaert, Ciccotti and Berendsen, 1977): take the normal, unconstrained step, then move each pair of atoms back along their <i>old</i> bond direction until the length is right. Constraints share atoms, so this is repeated (sweeps) until all lengths agree within a tolerance.</li>
<li><b>RATTLE</b> (Andersen, 1983): the velocity Verlet version. It also removes the velocity components along each constrained bond, so the bonds do not try to stretch.</li>
<li><b>SETTLE</b> (Miyamoto and Kollman, 1992): an exact, non-iterative solution for rigid three-site water. Used for water in essentially all codes.</li>
<li><b>LINCS</b> (Hess et al., 1997): solves all constraints at once with a matrix expansion. Fast and parallel; the GROMACS default for bonds.</li>
</ul>
<p><b>Consequences.</b> With all bonds to hydrogen constrained (and rigid water), Δt = 2 fs is standard. Each constraint removes one degree of freedom, which must be subtracted when computing the temperature. Making hydrogens heavier by moving mass from the atoms they are bonded to (<b>hydrogen mass repartitioning</b>, HMR) slows the remaining fast motions enough for 4 fs.</p>
<div class="callout"><b>Key idea.</b> Constrain X–H bonds (and water) to remove the fastest vibrations. SHAKE/RATTLE iterate, SETTLE is exact for water, LINCS handles the rest in parallel. Count N<sub>f</sub> = 3N − N<sub>c</sub> − 3.</div>`,

  algorithm: `
<p>Velocity Verlet with RATTLE:</p>
<pre class="pseudo">v += (dt/2) * F / m
x_new = x + dt * v                         # unconstrained step
repeat until all |x_i - x_j| = d_ij within tol:      # SHAKE
    for each constraint (i, j, d):
        s = x_new_i - x_new_j;  r = x_i - x_j      # r = old bond vector
        g = (s·s - d²) / (2 (1/m_i + 1/m_j) s·r)
        x_new_i -= g r / m_i;  x_new_j += g r / m_j
v = (x_new - x) / dt;  x = x_new
F = forces(x)
v += (dt/2) * F / m
repeat until all r_ij·v_ij = 0:                       # RATTLE
    for each constraint: k = r·v_ij / ((1/m_i + 1/m_j) d²)
        v_i -= k r / m_i;  v_j += k r / m_j</pre>
<ol class="steps">
  <li>Use a tight tolerance (10⁻⁶ to 10⁻⁸ relative). A loose tolerance lets the bond lengths drift and the energy is no longer conserved.</li>
  <li>Constrain only bonds (to hydrogen, or all bonds). Constraining angles as well makes the system too rigid and is avoided, except in rigid water.</li>
  <li>Thermostats and the temperature must use the reduced N<sub>f</sub>; MD codes do this automatically from the topology.</li>
  <li>GROMACS: <code>constraints = h-bonds</code>, <code>constraint-algorithm = lincs</code> (water uses SETTLE). LAMMPS: <code>fix shake</code> or <code>fix rattle</code>. AMBER: <code>ntc = 2, ntf = 2</code> (SHAKE on bonds to H). HMR: <code>parmed</code> HMassRepartition, or the force field's HMR topology.</li>
</ol>`,

  code: {
    note: 'Part 1 shows SHAKE restoring the bond lengths after one unconstrained step, sweep by sweep. Part 2 runs an SPC water dimer in vacuum with flexible bonds and with rigid SHAKE/RATTLE water at time steps from 0.5 to 4 fs. Part 3 counts degrees of freedom. Run <code>pip install numpy</code>, then <code>python shake.py</code> (about 10 s).',
    src: `import numpy as np
# Units: nm, ps, amu, kJ/mol.  SPC water: q_O = -0.82, sigma 0.3166 nm, eps 0.650 kJ/mol, O-H 0.1 nm, HOH 109.47 deg
COUL, qO, qH, sOO, eOO = 138.935, -0.82, 0.41, 0.3166, 0.650
r0, th0 = 0.1, np.radians(109.47); dHH = 2 * r0 * np.sin(th0 / 2)
kb, kth = 345000.0, 383.0                         # flexible SPC bond and angle constants
m = np.array([15.999, 1.008, 1.008] * 2)
q = np.array([qO, qH, qH] * 2)
cons = [(0, 1, r0), (0, 2, r0), (1, 2, dHH), (3, 4, r0), (3, 5, r0), (4, 5, dHH)]   # 3 per water

def forces(x, flexible):
    F = np.zeros_like(x); U = 0.0
    for i in range(3):                                         # between the molecules: Coulomb + O-O LJ
        for j in range(3, 6):
            d = x[i] - x[j]; r = np.linalg.norm(d)
            U += COUL * q[i] * q[j] / r; f = COUL * q[i] * q[j] / r**3 * d
            if i == 0 and j == 3:
                sr6 = (sOO / r) ** 6; U += 4 * eOO * (sr6**2 - sr6); f += 24 * eOO * (2 * sr6**2 - sr6) / r**2 * d
            F[i] += f; F[j] -= f
    if flexible:                                               # inside the molecules: bonds and angles
        for o in (0, 3):
            for h in (o + 1, o + 2):
                d = x[h] - x[o]; r = np.linalg.norm(d); U += 0.5 * kb * (r - r0)**2
                f = -kb * (r - r0) * d / r; F[h] += f; F[o] -= f
            a, b = x[o + 1] - x[o], x[o + 2] - x[o]; la, lb = np.linalg.norm(a), np.linalg.norm(b)
            c = a @ b / (la * lb); th = np.arccos(c); g = kth * (th - th0) / np.sin(th)
            Fa = g * (b / (la * lb) - c * a / la**2); Fb = g * (a / (la * lb) - c * b / lb**2)
            U += 0.5 * kth * (th - th0)**2; F[o + 1] += Fa; F[o + 2] += Fb; F[o] -= Fa + Fb
    return F, U

def shake(x_new, x_old, tol=1e-10, verbose=False):
    """Iteratively move atoms along the OLD bond vectors until every constraint length is satisfied."""
    for it in range(1, 500):
        worst = 0.0
        for i, j, d in cons:
            s = x_new[i] - x_new[j]; diff = s @ s - d * d; worst = max(worst, abs(diff) / d**2)
            rij = x_old[i] - x_old[j]
            g = diff / (2 * (1 / m[i] + 1 / m[j]) * (s @ rij))
            x_new[i] -= g * rij / m[i]; x_new[j] += g * rij / m[j]
        if verbose and (it in (1, 2, 3, 5, 10, 20, 30, 40) or worst < tol): print(f"   SHAKE iteration {it:2d}: largest relative length error {worst:.1e}")
        if worst < tol: return x_new, it
    raise RuntimeError("SHAKE did not converge")

def rattle_v(x, v, tol=1e-12):
    """Remove velocity components along each constraint (second half of RATTLE)."""
    for it in range(500):
        worst = 0.0
        for i, j, d in cons:
            rij, vij = x[i] - x[j], v[i] - v[j]; k = (rij @ vij) / ((1 / m[i] + 1 / m[j]) * d * d)
            worst = max(worst, abs(rij @ vij)); v[i] -= k * rij / m[i]; v[j] += k * rij / m[j]
        if worst < tol: return v

def start():
    def water(O, ang, tilt):
        h = [np.array([np.cos(th0 / 2), s * np.sin(th0 / 2) * np.cos(tilt), s * np.sin(th0 / 2) * np.sin(tilt)]) * r0 for s in (1, -1)]
        R = np.array([[np.cos(ang), -np.sin(ang), 0], [np.sin(ang), np.cos(ang), 0], [0, 0, 1]])
        return [O] + [O + R @ hh for hh in h]
    x = np.array(water(np.zeros(3), np.pi + th0 / 2, 0) + water(np.array([0.28, 0, 0]), np.radians(-30), np.pi / 2))
    v = np.random.default_rng(5).normal(0, 1, x.shape) * np.sqrt(2.494 / m)[:, None]   # ~300 K
    return x, v

# ---------- Part 1: one SHAKE correction ----------
print("Part 1: one unconstrained step of 2 fs breaks the bond lengths; SHAKE restores them")
print("   (water's triangle of constraints couples strongly, so plain SHAKE needs many sweeps; SETTLE solves it exactly)")
x, v = start(); x, _ = shake(x.copy(), x.copy())                         # make the start exactly rigid
x_new = x + 0.002 * v
print("   O-H lengths before SHAKE:", np.round([np.linalg.norm(x_new[0] - x_new[k]) for k in (1, 2)], 5), "nm")
x_new, n = shake(x_new, x, verbose=True)
print("   O-H lengths after SHAKE: ", np.round([np.linalg.norm(x_new[0] - x_new[k]) for k in (1, 2)], 5), "nm")

# ---------- Part 2: flexible vs rigid water dimer, energy conservation ----------
def run(dt_fs, rigid, t_ps=1.0):
    dt = dt_fs / 1000; x, v = start()
    if rigid: x, _ = shake(x.copy(), x.copy()); v = rattle_v(x, v)
    F, U = forces(x, not rigid); E = []
    for _ in range(int(t_ps / dt)):
        v += 0.5 * dt * F / m[:, None]; x_new = x + dt * v
        if rigid: x_new, _ = shake(x_new, x, tol=1e-9); v = (x_new - x) / dt     # positions and velocities consistent
        x = x_new; F, U = forces(x, not rigid); v += 0.5 * dt * F / m[:, None]
        if rigid: v = rattle_v(x, v)
        E.append(U + 0.5 * (m[:, None] * v**2).sum())
        if not np.isfinite(E[-1]) or abs(E[-1] - E[0]) > 1e4: return None
    E = np.array(E); return E.std(), E[-1] - E[0]

print("\\nPart 2: SPC water dimer, 1 ps NVE. Spread and drift of the total energy (kJ/mol)")
print("   dt (fs)      flexible: spread  drift        rigid (SHAKE/RATTLE): spread  drift")
for dt in [0.5, 1.0, 2.0, 3.0, 4.0]:
    f, r = run(dt, False), run(dt, True)
    fs = "    blew up       " if f is None else f"{f[0]:8.3f} {f[1]:+8.3f}"
    print(f"   {dt:4.1f}        {fs}            {r[0]:8.4f} {r[1]:+8.4f}")

# ---------- Part 3: degrees of freedom ----------
Nw = 1000; N = 3 * Nw
print(f"\\nPart 3: {Nw} rigid waters: N_f = 3N - N_constraints - 3 = {3 * N} - {3 * Nw} - 3 = {3 * N - 3 * Nw - 3}")
`,
  },

  maths: () => `
<h4>1. Constrained equations of motion</h4>
${tex(String.raw`\sigma_k(\mathbf r) = |\mathbf r_i - \mathbf r_j|^2 - d_{ij}^2 = 0, \qquad m_i\ddot{\mathbf r}_i = \mathbf F_i - \sum_k \lambda_k \nabla_i \sigma_k`, true)}
<p>The Lagrange multipliers λ<sub>k</sub> are the constraint forces, determined by requiring σ<sub>k</sub> = 0 after the step.</p>
<h4>2. One SHAKE correction (single constraint)</h4>
${tex(String.raw`\mathbf s = \mathbf r_i' - \mathbf r_j', \quad \mathbf r_{ij} = \mathbf r_i - \mathbf r_j, \qquad g = \frac{s^2 - d^2}{2\big(\tfrac{1}{m_i} + \tfrac{1}{m_j}\big)\,\mathbf s\cdot\mathbf r_{ij}}, \qquad \mathbf r_i' \mathrel{-}= \frac{g}{m_i}\mathbf r_{ij}, \;\; \mathbf r_j' \mathrel{+}= \frac{g}{m_j}\mathbf r_{ij}`, true)}
<p>This is a Newton step for that constraint alone. Coupled constraints need repeated sweeps; the correction conserves the total momentum because it moves the atoms in opposite directions, weighted by 1/m.</p>
<h4>3. RATTLE velocity condition</h4>
${tex(String.raw`\frac{d\sigma_k}{dt} = 2\,\mathbf r_{ij}\cdot\mathbf v_{ij} = 0`, true)}
<h4>4. Degrees of freedom</h4>
${tex(String.raw`N_f = 3N - N_c - 3, \qquad T = \frac{2K}{N_f k_B}`, true)}
<p>Rigid water has 3 constraints per molecule: 1000 waters give N<sub>f</sub> = 9000 − 3000 − 3 = 5997.</p>`,

  quiz: [
    {
      q: 'What is the main reason to constrain bonds to hydrogen?',
      options: ['To make the force field more accurate', 'To remove the fastest vibrations so a larger time step (about 2 fs) can be used', 'To speed up the Coulomb calculation', 'To keep molecules from rotating'],
      answer: 1,
      why: 'X–H bonds vibrate with a period of about 10 fs. Freezing them lets the next-fastest motions set the time step.',
    },
    {
      q: 'In SHAKE, along which direction are the atoms moved to correct a bond length?',
      options: ['Along the new (unconstrained) bond vector', 'Along the old bond vector from the previous step', 'Along the force', 'Along the velocity'],
      answer: 1,
      why: 'The constraint force acts along the bond at the start of the step, so the correction is applied along the old bond vector r_ij.',
    },
    {
      q: 'A box of 500 rigid water molecules (with the total momentum removed). How many degrees of freedom?',
      options: ['4500', '4497', '2997', '1497'],
      answer: 2,
      why: '3N = 4500 for 1500 atoms, minus 3 constraints × 500 waters = 3000, minus 3 for momentum = 2997.',
    },
  ],
};
