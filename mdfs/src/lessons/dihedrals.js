import { tex } from '../tex.js';

export default {
  id: 'dihedrals',
  summary: 'A dihedral (torsion) angle describes the twist around a bond, defined by four atoms. Its energy, a short sum of cosines, controls the shapes molecules adopt: trans or gauche in butane, the backbone of proteins, the flatness of rings.',
  minutes: 15,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/dihedrals.mp4', poster: 'videos/dihedrals.jpg' },
  sim: 'dihedrals',
  learn: [
    'How the dihedral angle φ is defined and computed from four positions',
    'Torsion potentials as cosine series (and Ryckaert–Bellemans), and the butane trans/gauche example',
    'Improper dihedrals and the scaled 1–4 non-bonded interactions',
  ],
  simHint: 'Drag the red atom in the Newman projection to feel the energy profile. With the butane preset, read the trans and gauche populations at 300 K, then raise the temperature. Switch to 3-fold and Double bond and watch how the rotor (red dot) hops. Try building your own profile: what does k₁ do on its own?',

  idea: `
<p>Bonds (two atoms) and angles (three atoms) fix distances and bends. The last internal coordinate is the <b>twist</b> around a bond. Take four bonded atoms 1–2–3–4. Atoms 1, 2, 3 define one plane, and 2, 3, 4 another. The <b>dihedral angle</b> φ is the angle between these planes. Looking straight down the 2–3 bond (a <b>Newman projection</b>), φ is the angle between atom 1 in front and atom 4 behind. By convention φ = 180° is <i>trans</i> (anti), ±60° <i>gauche</i>, and 0° <i>cis</i>.</p>
<p><b>Torsion energy.</b> Rotating around a single bond is periodic, so the energy is a short Fourier series: U(φ) = Σ k<sub>n</sub>[1 + cos(nφ − δ<sub>n</sub>)], with multiplicities n = 1, 2, 3 (sometimes up to 6). An sp³–sp³ bond like ethane's has a 3-fold term; a double bond or amide a 2-fold term that keeps it planar. GROMACS also offers the equivalent <b>Ryckaert–Bellemans</b> form, a polynomial in cos ψ.</p>
<p><b>Butane.</b> The trans state is lowest. Gauche is about 3 kJ/mol higher (≈ 1 k<sub>B</sub>T at 300 K), so both are populated: about 65% trans and 35% gauche. The barrier between them is about 12 kJ/mol (≈ 5 k<sub>B</sub>T), so changes are rare events, a few per 100 ps. Torsions with higher barriers (peptide bonds, about 80 kJ/mol) essentially never flip in a normal simulation.</p>
<p><b>Impropers.</b> To keep planar groups flat (aromatic rings, sp² carbons, amides) or to stop a chiral centre inverting, force fields add an <b>improper dihedral</b>: a harmonic or cosine term on the angle defined by four atoms that are not in a chain.</p>
<p><b>1–4 interactions.</b> Atoms 1 and 4 of a dihedral are only three bonds apart. Their Lennard-Jones and Coulomb interactions are kept, but <b>scaled</b> (OPLS: ×0.5 for both; AMBER: LJ ×0.5, Coulomb ×1/1.2), because the torsion term was fitted together with them. 1–2 and 1–3 pairs are excluded entirely.</p>
<div class="callout"><b>Key idea.</b> Dihedrals control conformations. The torsion term, the scaled 1–4 non-bonded forces, and impropers together decide which shapes a molecule adopts and how often it changes between them.</div>`,

  algorithm: `
<p>Computing φ and the forces:</p>
<pre class="pseudo">b1 = r2 - r1;  b2 = r3 - r2;  b3 = r4 - r3
n1 = b1 × b2;  n2 = b2 × b3              # normals of the two planes
m1 = n1 × (b2 / |b2|)
phi = atan2(m1 · n2, n1 · n2)            # signed, -180..180, trans = 180

dU = dU/dphi                             # e.g. -Σ n k_n sin(n phi - delta_n)
F1 = -dU * |b2| / |n1|² * n1             # (sign depends on the phi convention)
F4 = +dU * |b2| / |n2|² * n2
F2, F3 from the projections of b1, b3 on b2, with F1 + F2 + F3 + F4 = 0</pre>
<ol class="steps">
  <li>Use atan2, not arccos: arccos loses the sign of φ and is inaccurate near 0° and 180°.</li>
  <li>Sign conventions differ between references. Always check analytic forces against a numerical derivative (the code tab does).</li>
  <li>A rotatable bond usually has several dihedrals (e.g. 9 H–C–C–H style terms in ethane); their sum is what matters.</li>
  <li>GROMACS: <code>[ dihedrals ]</code> function 9 (multiple cosines), 3 (Ryckaert–Bellemans), 4 or 2 (impropers); <code>fudgeLJ</code>, <code>fudgeQQ</code> set the 1–4 scaling. LAMMPS: <code>dihedral_style charmm / opls / multi/harmonic</code>, <code>special_bonds</code> for 1–4 scaling.</li>
</ol>`,

  code: {
    note: 'Builds a butane carbon chain, measures signed dihedrals, checks the analytic dihedral forces against numerical derivatives, computes trans/gauche populations from the Ryckaert–Bellemans potential, and runs 400 Langevin rotors to confirm them. Run <code>pip install numpy</code>, then <code>python dihedrals.py</code> (about 11 s).',
    src: `import numpy as np
kT = 2.494                                   # kJ/mol at 300 K

def dihedral(p0, p1, p2, p3):
    """Signed dihedral angle (radians) for atoms 0-1-2-3, IUPAC convention (trans = 180 deg)."""
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    m1 = np.cross(n1, b2 / np.linalg.norm(b2))
    return np.arctan2(np.dot(m1, n2), np.dot(n1, n2))

# Ryckaert-Bellemans torsion for butane (GROMACS units kJ/mol), psi = phi - 180 deg
C = [9.28, 12.16, -13.12, -3.06, 26.24, -31.5]
def U_rb(phi):
    c = np.cos(phi - np.pi)
    return sum(Cn * c**n for n, Cn in enumerate(C))
def dU_rb(phi):
    c, s = np.cos(phi - np.pi), np.sin(phi - np.pi)
    return sum(-n * Cn * c**(n - 1) * s for n, Cn in enumerate(C) if n > 0)

def dihedral_forces(p, dUdphi):
    """Forces F_i = -dU/dphi * dphi/dr_i. The gradients follow Blondel & Karplus (1996),
    with the sign flipped because their phi has the opposite sign to the IUPAC phi above."""
    F, G, H = p[0] - p[1], p[1] - p[2], p[3] - p[2]
    A, B = np.cross(F, G), np.cross(H, G); A2, B2, lG = A @ A, B @ B, np.linalg.norm(G)
    g0 = lG / A2 * A
    g3 = -lG / B2 * B
    g1 = -lG / A2 * A - (F @ G) / (A2 * lG) * A + (H @ G) / (B2 * lG) * B
    g2 = lG / B2 * B + (F @ G) / (A2 * lG) * A - (H @ G) / (B2 * lG) * B
    return -dUdphi * np.array([g0, g1, g2, g3])

# ---------- 1. a butane carbon skeleton ----------
def butane(phi):
    """C1-C2-C3-C4 with C-C = 0.153 nm, angle 112 deg, dihedral phi."""
    b, th = 0.153, np.deg2rad(112)
    p1, p2 = np.zeros(3), np.array([b, 0, 0])
    p0 = p1 + b * np.array([np.cos(th), np.sin(th), 0])                       # angle C1-C2-C3 = th
    d = np.array([-np.cos(th), np.sin(th) * np.cos(phi), -np.sin(th) * np.sin(phi)])   # angle C2-C3-C4 = th
    return np.array([p0, p1, p2, p2 + b * d])

for deg in [180, 60, -60, 0]:
    print(f"built with {deg:4d} deg -> measured {np.degrees(dihedral(*butane(np.radians(deg)))):7.1f} deg")

# ---------- 2. check the dihedral forces numerically ----------
p = butane(np.radians(75)) + np.random.default_rng(0).normal(0, 0.01, (4, 3))
phi = dihedral(*p); F = dihedral_forces(p, dU_rb(phi)); h = 1e-6; Fn = np.zeros_like(p)
for i in range(4):
    for d in range(3):
        q1, q2 = p.copy(), p.copy(); q1[i, d] += h; q2[i, d] -= h
        Fn[i, d] = -(U_rb(dihedral(*q1)) - U_rb(dihedral(*q2))) / (2 * h)
print(f"\\nanalytic vs numerical forces: max diff {np.abs(F - Fn).max():.1e} (forces ~ {np.abs(F).max():.0f} kJ/mol/nm), sum {np.abs(F.sum(0)).max():.1e}")

# ---------- 3. the butane torsion profile and Boltzmann populations ----------
ph = np.radians(np.arange(-180, 180, 0.5)); U = U_rb(ph)
print(f"\\nU(trans 180) = {round(U_rb(np.pi), 2) + 0.0:.2f}, U(gauche 60) = {U_rb(np.radians(60)):.2f} min at "
      f"{np.degrees(ph[(ph > 0) & (ph < np.pi / 2)][np.argmin(U[(ph > 0) & (ph < np.pi / 2)])]):.0f} deg, "
      f"barrier at 120 deg = {U_rb(np.radians(120)):.2f}, cis (0) = {U_rb(0):.2f} kJ/mol")
w = np.exp(-U / kT); trans = np.abs(ph) > np.radians(120)
print(f"Boltzmann at 300 K: trans {w[trans].sum() / w.sum() * 100:.1f} %, gauche+ and gauche- {w[~trans].sum() / w.sum() * 100:.1f} %")

# ---------- 4. let 400 rotors move with Langevin dynamics ----------
rng = np.random.default_rng(1); I, gamma, dt = 0.34, 5.0, 0.001          # amu nm^2, 1/ps, ps
x = rng.uniform(-np.pi, np.pi, 400); w_ = np.zeros(400); hops = 0; tr = []
c1 = np.exp(-gamma * dt); c2 = np.sqrt((1 - c1**2) * kT / I)
for s in range(80000):                                                      # 80 ps
    w_ += -0.5 * dt * dU_rb(x) / I; x += 0.5 * dt * w_
    w_ = c1 * w_ + c2 * rng.normal(size=400)
    x += 0.5 * dt * w_; w_ += -0.5 * dt * dU_rb(x) / I
    x = (x + np.pi) % (2 * np.pi) - np.pi
    new = np.abs(x) > np.radians(120)
    if s == 20000: state = new
    if s > 20000:                                                           # skip the first 20 ps
        hops += (new != state).sum(); state = new; tr.append(new.mean())
print(f"Langevin rotors (60 ps): trans {np.mean(tr) * 100:.1f} %, about {hops / 400 / 60:.2f} trans/gauche changes per ps per rotor")
`,
  },

  maths: () => `
<h4>1. The dihedral angle</h4>
${tex(String.raw`\mathbf n_1 = \mathbf b_1\times\mathbf b_2, \quad \mathbf n_2 = \mathbf b_2\times\mathbf b_3, \qquad \phi = \operatorname{atan2}\!\Big(\big(\mathbf n_1\times\hat{\mathbf b}_2\big)\cdot\mathbf n_2,\; \mathbf n_1\cdot\mathbf n_2\Big)`, true)}
<h4>2. Torsion potentials</h4>
${tex(String.raw`U_{\text{prop}}(\phi) = \sum_n k_n\big[1 + \cos(n\phi - \delta_n)\big], \qquad U_{\text{RB}}(\psi) = \sum_{m=0}^{5} C_m \cos^m\psi, \;\; \psi = \phi - 180^\circ`, true)}
<p>Since cos<sup>m</sup>ψ can be written with cos(nψ) for n ≤ m, the two forms are equivalent.</p>
<h4>3. Improper dihedral</h4>
${tex(String.raw`U_{\text{imp}}(\xi) = \tfrac12 k_\xi(\xi - \xi_0)^2`, true)}
<h4>4. Populations and rates</h4>
${tex(String.raw`p(\phi) = \frac{e^{-U(\phi)/k_BT}}{\int e^{-U(\phi')/k_BT}\,d\phi'}, \qquad k_{\text{hop}} \sim \nu_0\, e^{-\Delta U^\ddagger / k_BT}`, true)}
<p>For butane, ΔU ≈ 3 kJ/mol gives about 65% trans at 300 K. A barrier of 12 kJ/mol (≈ 5k<sub>B</sub>T) makes changes rare compared with vibrations.</p>`,

  quiz: [
    {
      q: 'Which atoms define the dihedral angle around the bond between atoms 2 and 3?',
      options: ['Atoms 2 and 3 only', 'Atoms 1, 2 and 3', 'Atoms 1, 2, 3 and 4', 'Any three atoms'],
      answer: 2,
      why: 'φ is the angle between the plane of atoms 1–2–3 and the plane of atoms 2–3–4.',
    },
    {
      q: 'A torsion term k[1 + cos(2φ − 180°)] has minima at…',
      options: ['60° and 180°', '0° and 180° (planar)', '90° and −90°', 'Only 180°'],
      answer: 1,
      why: '1 + cos(2φ − 180°) = 1 − cos 2φ, which is zero at φ = 0° and 180°: a 2-fold term that keeps a double bond or amide planar.',
    },
    {
      q: 'Why are 1–4 non-bonded interactions scaled rather than removed?',
      options: ['To save computer time', 'The torsion parameters were fitted together with partly scaled 1–4 interactions, so the force field needs both', 'Because Coulomb forces are not defined at short range', 'They are always removed'],
      answer: 1,
      why: 'Atoms three bonds apart still feel each other, but the torsion term already captures part of it. Each force field uses the scaling factors it was parameterised with (OPLS 0.5/0.5, AMBER 0.5/0.833).',
    },
  ],
};
