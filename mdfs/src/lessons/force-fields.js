import { tex } from '../tex.js';

export default {
  id: 'force-fields',
  summary: 'A force field is the recipe for the potential energy: a fixed set of simple functions (bonds, angles, dihedrals, Lennard-Jones, Coulomb) plus thousands of parameters for each atom type. Here is how the pieces fit together, where the numbers come from, and the rules for using them.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/force-fields.mp4', poster: 'videos/force-fields.jpg' },
  sim: 'forcefield',
  learn: [
    'The full classical energy expression: bonded plus non-bonded terms',
    'Atom types, partial charges, combining rules and exclusions',
    'How force fields are parameterised, the main families, and polarizable and machine-learned alternatives',
  ],
  simHint: 'Drag water B around water A and watch which charge pairs attract (cyan) and repel (coral). Press "Find the minimum". Then change the oxygen charge: how do the energy and distance of the hydrogen bond change? Set ε to 0: what happens at short distance without the LJ repulsion?',

  idea: `
<p>Every lesson so far needed a potential energy U(<b>r</b>₁, …, <b>r</b><sub>N</sub>): the forces are its gradient. For real molecules we cannot solve quantum mechanics at every step, so we use a <b>force field</b>: an approximate, classical U built from simple pieces.</p>
<ul>
<li><b>Bonded terms</b> (Lessons 4.1–4.2): bonds, angles, proper and improper dihedrals. They act only between atoms connected through the molecule's topology.</li>
<li><b>Non-bonded terms</b>: Lennard-Jones (repulsion plus dispersion) and Coulomb interactions between fixed <b>partial charges</b>. They act between all pairs in different molecules, and between atoms in the same molecule that are more than three bonds apart (1–4 pairs scaled, 1–2 and 1–3 excluded).</li>
</ul>
<p><b>Atom types.</b> Parameters are not per element but per <i>atom type</i>: a carbon in CH₃, in C=O and in a benzene ring each have their own LJ parameters and bonded terms. A force field is essentially a large table of these types and the parameters for every combination.</p>
<p><b>Where the numbers come from.</b> Bond and angle constants from vibrational spectra or quantum (QM) calculations. Torsions fitted to QM rotation profiles. Partial charges fitted to the QM electrostatic potential (e.g. RESP) or chosen for liquid properties. LJ parameters tuned so that the simulation reproduces densities, heats of vaporisation and solvation free energies. Every parameter is fitted <i>together</i> with all the others, the functional form, the combining rule and the water model.</p>
<p><b>Families.</b> AMBER (GAFF for small molecules), CHARMM (CGenFF), OPLS-AA and GROMOS each have their own rules (combining rule, 1–4 scaling, cutoffs). Never mix parameters from different families. Beyond fixed charges, <b>polarizable</b> force fields (Drude, AMOEBA) let charges respond to their environment, and <b>machine-learned potentials</b> (MACE, NequIP, …) learn U directly from QM data, at a higher cost.</p>
<div class="callout"><b>Key idea.</b> Force field = functional form + parameters + rules (exclusions, 1–4 scaling, combining rule, cutoff and long-range treatment). The parameters only make sense with the rest of the recipe they were fitted with.</div>`,

  algorithm: `
<p>Evaluating a classical force field each step:</p>
<pre class="pseudo">U = 0
for each bond (i, j):         U += ½ k_b (r_ij - r0)²
for each angle (i, j, k):     U += ½ k_θ (θ_ijk - θ0)²
for each dihedral (i,j,k,l):  U += Σ_n k_n [1 + cos(n φ - δ_n)]
for each improper:            U += ½ k_ξ (ξ - ξ0)²
for each non-excluded pair (i, j) within the cutoff:
    σ_ij, ε_ij from the combining rule
    s = 1–4 scaling factor (or 1)
    U += s_LJ * 4 ε_ij [(σ_ij/r)¹² - (σ_ij/r)⁶]
       + s_QQ * q_i q_j / (4π ε0 r)              # plus long-range electrostatics (Ch. 5)</pre>
<ol class="steps">
  <li><b>Combining rules:</b> Lorentz–Berthelot σ<sub>ij</sub> = (σ<sub>i</sub> + σ<sub>j</sub>)/2, ε<sub>ij</sub> = √(ε<sub>i</sub>ε<sub>j</sub>) (AMBER, CHARMM), or geometric for both (OPLS). Use the one your force field was fitted with.</li>
  <li>The <b>topology</b> lists atoms, types, charges and all bonded terms. Exclusions are generated from it automatically (GROMACS <code>nrexcl = 3</code>).</li>
  <li>Small-molecule parameters come from tools such as antechamber/GAFF, CGenFF, LigParGen (OPLS) or the ATB (GROMOS). Always check the assigned charges and torsions.</li>
  <li>Validate before use: density, heat of vaporisation, diffusion or a solvation free energy against experiment.</li>
</ol>`,

  code: {
    note: 'Part 1: a complete united-atom force field (TraPPE-UA) for butane, term by term. Part 2: how the combining rule changes the interaction between unlike sites. Part 3: the TIP3P water dimer, LJ plus Coulomb, optimised to the literature value of about −27 kJ/mol. Run <code>pip install numpy</code>, then <code>python force_field.py</code> (under 1 s).',
    src: `import numpy as np
kB = 0.0083145                                   # kJ/mol per K
COUL = 138.935                                   # kJ mol^-1 nm e^-2

# ---------- Part 1: a complete (united-atom) force field for butane: TraPPE-UA ----------
# sites: CH3 - CH2 - CH2 - CH3 ; parameters from Martin & Siepmann (1998)
eps = {"CH3": 98 * kB, "CH2": 46 * kB}           # kJ/mol
sig = {"CH3": 0.375, "CH2": 0.395}               # nm
b0, th0, kth = 0.154, np.radians(114), 62500 * kB            # fixed bond, angle (U = kth/2 (th-th0)^2)
c1, c2, c3 = 355.03 * kB, -68.19 * kB, 791.32 * kB           # torsion

def angle(a, b, c):
    u, v = a - b, c - b
    return np.arccos(np.dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v))

def dihedral(p0, p1, p2, p3):
    b1, b2, b3 = p1 - p0, p2 - p1, p3 - p2
    n1, n2 = np.cross(b1, b2), np.cross(b2, b3)
    return np.arctan2(np.dot(np.cross(n1, b2 / np.linalg.norm(b2)), n2), np.dot(n1, n2))

def butane(phi, th=th0, x0=np.zeros(3)):
    p1, p2 = np.zeros(3), np.array([b0, 0, 0])
    p0 = p1 + b0 * np.array([np.cos(th), np.sin(th), 0])                    # angle p0-p1-p2 = th
    d = np.array([-np.cos(th), np.sin(th) * np.cos(phi), -np.sin(th) * np.sin(phi)])   # angle p1-p2-p3 = th
    return np.array([p0, p1, p2, p2 + b0 * d]) + x0

types = ["CH3", "CH2", "CH2", "CH3"]
def intramolecular(x):
    """Energy terms of one molecule. TraPPE excludes all non-bonded pairs closer than 4 bonds,
    so butane has no intramolecular LJ at all."""
    E = {"bonds": 0.0}                                       # bonds are fixed in TraPPE
    E["angles"] = sum(0.5 * kth * (angle(x[i], x[i + 1], x[i + 2]) - th0) ** 2 for i in range(2))
    phi = dihedral(*x)
    E["torsion"] = c1 * (1 + np.cos(phi)) + c2 * (1 - np.cos(2 * phi)) + c3 * (1 + np.cos(3 * phi))
    return E, np.degrees(phi)

print("Part 1: TraPPE-UA butane, energy terms (kJ/mol)")
for name, phi in [("trans", np.pi), ("gauche", np.radians(65))]:
    E, deg = intramolecular(butane(phi, th=np.radians(112)))
    print(f"  {name:6s} phi = {deg:6.1f}:  " + "  ".join(f"{k} {v:6.2f}" for k, v in E.items()) + f"   total {sum(E.values()):6.2f}")

# ---------- Part 2: two molecules: intermolecular LJ and combining rules ----------
def lj(r, s, e):
    return 4 * e * ((s / r) ** 12 - (s / r) ** 6)

def inter(xa, xb, rule):
    U = 0.0
    for i, ti in enumerate(types):
        for j, tj in enumerate(types):
            e = np.sqrt(eps[ti] * eps[tj])
            s = 0.5 * (sig[ti] + sig[tj]) if rule == "LB" else np.sqrt(sig[ti] * sig[tj])
            U += lj(np.linalg.norm(xa[i] - xb[j]), s, e)
    return U

print("\\nPart 2: combining rules matter for unlike sites")
print("  CH3 (sigma 0.375 nm, eps 0.815 kJ/mol) with a small site (sigma 0.20 nm, eps 0.20 kJ/mol; illustrative)")
s1, e1, s2, e2 = sig["CH3"], eps["CH3"], 0.20, 0.20
for name, s12 in [("Lorentz-Berthelot (arithmetic)", 0.5 * (s1 + s2)), ("geometric", np.sqrt(s1 * s2))]:
    e12 = np.sqrt(e1 * e2)
    print(f"  {name:31s} sigma_ij = {s12:.4f} nm, minimum at {2**(1/6) * s12:.3f} nm, U(0.30 nm) = {lj(0.30, s12, e12):+.3f} kJ/mol")
print("  Parameters only make sense together with the combining rule (and functional form) they were fitted with.")

# ---------- Part 3: TIP3P water dimer: LJ + Coulomb ----------
qO, qH = -0.834, 0.417; sO, eO = 0.315061, 0.6364; rOH, hoh = 0.09572, np.radians(104.52)
def water(origin, axis_angle, tilt):
    """Molecule with O at origin; plane rotated by tilt about the bisector, then the bisector tilted by axis_angle (about z)."""
    h = np.array([[np.cos(hoh / 2), s * np.sin(hoh / 2) * np.cos(tilt), s * np.sin(hoh / 2) * np.sin(tilt)] for s in (1, -1)]) * rOH
    c, s_ = np.cos(axis_angle), np.sin(axis_angle)
    R = np.array([[c, -s_, 0], [s_, c, 0], [0, 0, 1]])
    return np.vstack([origin, origin + h @ R.T])

def dimer_energy(x1, x2):
    q = np.array([qO, qH, qH]); Ec = 0.0
    for i in range(3):
        for j in range(3):
            Ec += COUL * q[i] * q[j] / np.linalg.norm(x1[i] - x2[j])
    return lj(np.linalg.norm(x1[0] - x2[0]), sO, eO), Ec

best = None
for R in np.arange(0.25, 0.32, 0.0025):
    for beta in np.radians(np.arange(-80, 81, 2)):
        d = water(np.zeros(3), np.radians(180) - hoh / 2, 0.0)          # donor: one O-H along +x
        d[1] = np.array([rOH, 0, 0]); d[2] = rOH * np.array([np.cos(hoh), np.sin(hoh), 0])
        a = water(np.array([R, 0, 0]), beta, np.pi / 2)                  # acceptor: H atoms point away
        Elj, Ec = dimer_energy(d, a)
        if best is None or Elj + Ec < best[0]: best = (Elj + Ec, R, np.degrees(beta), Elj, Ec)
E, R, beta, Elj, Ec = best
print(f"\\nPart 3: TIP3P water dimer (linear hydrogen bond)")
print(f"  best O-O distance {R:.3f} nm, acceptor tilted {abs(beta):.0f} deg")
print(f"  energy {E:.1f} kJ/mol = LJ {Elj:+.1f} + Coulomb {Ec:+.1f}   (TIP3P literature: about -27 kJ/mol)")
`,
  },

  maths: () => `
<h4>1. A typical (class I) force field</h4>
${tex(String.raw`U = \sum_{\text{bonds}} \tfrac12 k_b(r - r_0)^2 + \sum_{\text{angles}} \tfrac12 k_\theta(\theta - \theta_0)^2 + \sum_{\text{dihedrals}} \sum_n k_n\big[1 + \cos(n\phi - \delta_n)\big] + \sum_{\text{impropers}} \tfrac12 k_\xi(\xi - \xi_0)^2`, true)}
${tex(String.raw`\;+\; \sum_{i<j}^{\text{non-excl.}} \left\{ 4\varepsilon_{ij}\left[\Big(\frac{\sigma_{ij}}{r_{ij}}\Big)^{12} - \Big(\frac{\sigma_{ij}}{r_{ij}}\Big)^{6}\right] + \frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}} \right\}`, true)}
<h4>2. Combining rules</h4>
${tex(String.raw`\text{Lorentz–Berthelot: } \sigma_{ij} = \frac{\sigma_i + \sigma_j}{2},\; \varepsilon_{ij} = \sqrt{\varepsilon_i\varepsilon_j} \qquad \text{geometric: } \sigma_{ij} = \sqrt{\sigma_i\sigma_j},\; \varepsilon_{ij} = \sqrt{\varepsilon_i\varepsilon_j}`, true)}
<h4>3. Coulomb constant in MD units</h4>
${tex(String.raw`\frac{1}{4\pi\varepsilon_0} = 138.935\ \text{kJ mol}^{-1}\,\text{nm}\,e^{-2} = 332.06\ \text{kcal mol}^{-1}\,\text{Å}\,e^{-2}`, true)}
<p>Two charges of ±0.4 e at 0.2 nm: 138.935 × 0.16 / 0.2 ≈ 111 kJ/mol. Electrostatics dominate hydrogen bonding, which is why charges are the most important (and most debated) parameters.</p>`,

  quiz: [
    {
      q: 'Which non-bonded pairs inside one molecule are normally excluded completely?',
      options: ['None', '1–2 and 1–3 pairs (directly bonded or sharing an angle)', 'All pairs', 'Only 1–4 pairs'],
      answer: 1,
      why: 'Their interaction is described by the bond and angle terms. 1–4 pairs are kept but scaled; pairs further apart are fully included.',
    },
    {
      q: 'You take LJ parameters for a new ion from a CHARMM paper and use them with the OPLS-AA force field. What is the main risk?',
      options: ['None, LJ parameters are universal', 'Different combining rules and fitting conditions make the parameters inconsistent', 'The simulation will run slower', 'Charges will be ignored'],
      answer: 1,
      why: 'Parameters are fitted together with a functional form, a combining rule, 1–4 scaling and a water model. Mixing families gives unpredictable errors.',
    },
    {
      q: 'In the TIP3P water dimer, which term provides most of the hydrogen-bond attraction?',
      options: ['Lennard-Jones dispersion', 'Coulomb interaction between the partial charges', 'Bond stretching', 'Angle bending'],
      answer: 1,
      why: 'Coulomb ≈ −34 kJ/mol, while LJ is +7 kJ/mol (repulsive at the H-bond distance). The total is about −27 kJ/mol.',
    },
  ],
};
