import { tex } from '../tex.js';

export default {
  id: 'water-models',
  summary: 'In most biomolecular and electrolyte simulations, the majority of atoms are water. Simple rigid models with 3, 4 or 5 interaction sites reproduce different properties well or badly. Choosing (and not mixing) the right one matters.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/water-models.mp4', poster: 'videos/water-models.jpg' },
  sim: 'water',
  learn: [
    '3-, 4- and 5-site rigid water models, and the virtual site M',
    'Why model dipoles are larger than the gas-phase value (effective polarisation)',
    'How models compare for density, diffusion, dielectric constant and melting point, and how to choose one',
  ],
  simHint: 'Switch between the models and watch the electrostatic potential: in the 4-site models the negative charge sits on M, slightly in front of the oxygen. Compare the dipole moments with the gas value (1.85 D). Then flip through the property bars: which model is closest to experiment for the dielectric constant, and which for the melting point?',

  idea: `
<p>Water is the solvent of life and of most electrolytes, and in a typical simulation box 80–95% of the atoms belong to water. A water model must be cheap (it is evaluated billions of times) yet reproduce the liquid's structure, dynamics and thermodynamics.</p>
<p><b>Rigid point-charge models</b> are the workhorses. The geometry is fixed (Lesson 4.4, SETTLE), one Lennard-Jones site sits on the oxygen, and partial charges model the electrostatics:</p>
<ul>
<li><b>3-site</b> (TIP3P, SPC, SPC/E): charges on O and on both H atoms. Simple and cheap.</li>
<li><b>4-site</b> (TIP4P, TIP4P/2005, TIP4P-Ew, OPC): the negative charge moves from O to a massless <b>virtual site M</b> on the H–O–H bisector, about 0.15 Å from O. This gives a better charge distribution (quadrupole). Forces on M are passed back to O and H.</li>
<li><b>5-site</b> (TIP5P): two negative sites at the "lone pair" positions, for a tetrahedral charge distribution.</li>
</ul>
<p><b>Effective polarisation.</b> A water molecule in the liquid is polarised by its neighbours; its dipole grows from 1.85 D (gas) to an estimated 2.6–3.0 D. Fixed-charge models cannot respond, so their charges are <b>boosted</b> to give about 2.3–2.5 D: a compromise that works in the liquid but overbinds the gas-phase dimer. Polarizable models (SWM4, AMOEBA, Drude) and machine-learned potentials treat this explicitly.</p>
<p><b>No model gets everything right.</b> TIP3P diffuses more than twice too fast and melts at 146 K, but it is the water that AMBER and CHARMM were parameterised with. SPC/E is a good all-rounder. TIP4P/2005 reproduces densities, the temperature of maximum density and the phase diagram well, but underestimates the dielectric constant. OPC was optimised for its charge distribution and gives a good dielectric constant.</p>
<div class="callout"><b>Key idea.</b> Use the water model your force field was parameterised with (TIP3P for most AMBER/CHARMM work, SPC for GROMOS, TIP4P-type for many OPLS and phase-behaviour studies), and know which properties it gets wrong.</div>`,

  algorithm: `
<p>Handling a 4-site model in an MD code:</p>
<pre class="pseudo">every step, after moving O, H1, H2:
    M = O + a * ((H1 - O) + (H2 - O))          # rebuild the virtual site
compute forces on all sites (LJ on O, Coulomb on H1, H2, M)
spread the force on M back to the real atoms:
    F_O  += (1 - 2a) * F_M
    F_H1 += a * F_M;   F_H2 += a * F_M
integrate only O, H1, H2 (M has no mass), then SETTLE</pre>
<ol class="steps">
  <li>a is fixed by the geometry: a = d<sub>OM</sub> / |(H1 − O) + (H2 − O)|. Spreading this way conserves the total force and torque (checked in the code tab).</li>
  <li>Models are fitted with particular simulation settings (cutoff, long-range electrostatics, LJ tail correction). TIP4P-Ew and OPC assume Ewald/PME.</li>
  <li>Use the same model for every water in the system, including ion and ligand parameterisation: ion parameters (e.g. Joung–Cheatham) are water-model specific.</li>
  <li>GROMACS: <code>#include "amber99sb.ff/tip3p.itp"</code> or <code>tip4p.itp</code>; build with <code>gmx solvate -cs tip4p.gro</code>. LAMMPS: <code>pair_style lj/cut/tip4p/long</code> for 4-site models.</li>
</ol>`,

  code: {
    note: 'Builds TIP3P, SPC/E, TIP4P/2005 and OPC, including the virtual site M. Computes their dipole moments and gas-phase dimer energies, checks that spreading the force on M conserves force and torque, and prints literature liquid properties. Run <code>pip install numpy</code>, then <code>python water_models.py</code> (about 1 s).',
    src: `import numpy as np
COUL = 138.935                      # kJ mol^-1 nm e^-2
DEBYE = 48.032                      # 1 e*nm in Debye
# rigid models: O-H length (nm), HOH angle (deg), q_H, O-M distance (nm, 0 = 3-site), LJ sigma (nm), eps (kJ/mol)
MODELS = {
    "TIP3P":      dict(r=0.09572, a=104.52, qH=0.4170, dM=0.0,     s=0.315061, e=0.63640),
    "SPC/E":      dict(r=0.1000,  a=109.47, qH=0.4238, dM=0.0,     s=0.316557, e=0.65019),
    "TIP4P/2005": dict(r=0.09572, a=104.52, qH=0.5564, dM=0.01546, s=0.315890, e=0.77490),
    "OPC":        dict(r=0.08724, a=103.60, qH=0.6791, dM=0.01594, s=0.316655, e=0.89036),
}

def sites(p, O=np.zeros(3), R=np.eye(3)):
    """O, H1, H2 and the charge site (O itself for 3-site models, the virtual site M otherwise)."""
    h = np.radians(p["a"]) / 2
    H1, H2 = p["r"] * np.array([np.cos(h), np.sin(h), 0]), p["r"] * np.array([np.cos(h), -np.sin(h), 0])
    M = p["dM"] * np.array([1.0, 0, 0])          # on the HOH bisector, towards the H atoms
    return [O + R @ v for v in (np.zeros(3), H1, H2, M)]

def charges(p):
    return np.array([0.0, p["qH"], p["qH"], -2 * p["qH"]])      # [O (LJ only), H, H, M or O]

print("Part 1: geometry and dipole moment (gas-phase experiment: 1.85 D; liquid: estimated 2.6-3.0 D)")
for name, p in MODELS.items():
    x, q = sites(p), charges(p)
    mu = np.linalg.norm((q[:, None] * np.array(x)).sum(0)) * DEBYE
    print(f"  {name:11s} sites = {3 if p['dM'] == 0 else 4}   q_H = {p['qH']:.4f}   dipole = {mu:.3f} D")

def rot_z(t):
    return np.array([[np.cos(t), -np.sin(t), 0], [np.sin(t), np.cos(t), 0], [0, 0, 1]])
def rot_x(t):
    return np.array([[1, 0, 0], [0, np.cos(t), -np.sin(t)], [0, np.sin(t), np.cos(t)]])

def dimer(p, R, beta):
    """Donor with its first O-H along +x, acceptor at distance R with its plane perpendicular, tilted by beta."""
    h = np.radians(p["a"]) / 2
    A = sites(p, R=rot_z(-h))                                    # donor: H1 now along +x
    B = sites(p, O=np.array([R, 0, 0]), R=rot_z(beta) @ rot_x(np.pi / 2))
    q = charges(p); E = 0.0
    for i in range(4):
        for j in range(4):
            if q[i] and q[j]: E += COUL * q[i] * q[j] / np.linalg.norm(A[i] - B[j])
    r = np.linalg.norm(A[0] - B[0]); sr6 = (p["s"] / r) ** 6
    return E + 4 * p["e"] * (sr6**2 - sr6)

print("\\nPart 2: gas-phase dimer minimum (experiment: about -21 kJ/mol at 0.298 nm)")
for name, p in MODELS.items():
    best = min((dimer(p, R, b), R, b) for R in np.arange(0.25, 0.31, 0.001) for b in np.radians(np.arange(-80, 81, 2)))
    print(f"  {name:11s} E = {best[0]:6.1f} kJ/mol at O-O = {best[1]:.3f} nm")
print("  All models over-bind the gas dimer: their charges are boosted to mimic polarisation in the liquid.")

# ---------- Part 3: forces on a virtual site ----------
p = MODELS["TIP4P/2005"]; O, H1, H2, M = sites(p)
# M = O + a (H1 - O) + a (H2 - O) with a chosen so that |M - O| = dM
a = p["dM"] / np.linalg.norm((H1 - O) + (H2 - O))
print(f"\\nPart 3: virtual site M = O + a[(H1 - O) + (H2 - O)], a = {a:.5f};  check |M - O| = {np.linalg.norm(a * ((H1 - O) + (H2 - O))):.5f} nm")
FM = np.array([10.0, -4.0, 2.0])                                 # any force acting on M
FO, FH1, FH2 = (1 - 2 * a) * FM, a * FM, a * FM                  # spread it onto the real atoms
torque = lambda r, F: np.cross(r, F)
print(f"  force on M {FM} -> O {np.round(FO, 3)}, H1 {np.round(FH1, 3)}, H2 {np.round(FH2, 3)}")
print(f"  total force kept: {np.allclose(FO + FH1 + FH2, FM)};  torque kept: "
      f"{np.allclose(torque(O, FO) + torque(H1, FH1) + torque(H2, FH2), torque(M, FM))}")

# ---------- Part 4: literature properties at 298 K, 1 bar ----------
print("\\nPart 4: liquid properties at 298 K (Vega & Abascal, PCCP 2011; Izadi et al. 2014)")
print("  model        density   diffusion (1e-5 cm2/s)   dielectric   T_melt (K)   T_max density (K)")
for row in [("TIP3P", "0.98", "~5.5-6", "94", "146", "182"), ("SPC/E", "0.994", "2.5", "71", "215", "241"),
            ("TIP4P/2005", "0.997", "2.1", "58", "252", "278"), ("OPC", "0.997", "2.3", "78", "-", "272"),
            ("experiment", "0.997", "2.3", "78.4", "273", "277")]:
    print(f"  {row[0]:11s}  {row[1]:7s}   {row[2]:22s}   {row[3]:10s}   {row[4]:10s}   {row[5]}")
`,
  },

  maths: () => `
<h4>1. The interaction between two rigid waters</h4>
${tex(String.raw`U_{ab} = 4\varepsilon_{OO}\left[\Big(\frac{\sigma_{OO}}{r_{OO}}\Big)^{12} - \Big(\frac{\sigma_{OO}}{r_{OO}}\Big)^{6}\right] + \sum_{i\in a}\sum_{j\in b}\frac{q_i q_j}{4\pi\varepsilon_0 r_{ij}}`, true)}
<p>The charge sites are H, H and O (3-site) or H, H and M (4-site).</p>
<h4>2. Virtual site and dipole moment</h4>
${tex(String.raw`\mathbf r_M = \mathbf r_O + a\big[(\mathbf r_{H_1} - \mathbf r_O) + (\mathbf r_{H_2} - \mathbf r_O)\big], \qquad \boldsymbol\mu = \sum_i q_i \mathbf r_i`, true)}
<p>1 e·nm = 48.03 D. For TIP4P/2005: q<sub>H</sub> = 0.5564 e, d<sub>OM</sub> = 0.01546 nm gives μ = 2.305 D.</p>
<h4>3. Dielectric constant from dipole fluctuations</h4>
${tex(String.raw`\varepsilon_r = 1 + \frac{\langle \mathbf M^2\rangle - \langle\mathbf M\rangle^2}{3\varepsilon_0 V k_B T}, \qquad \mathbf M = \sum_{\text{molecules}} \boldsymbol\mu`, true)}
<p>(with conducting boundary conditions, as in Ewald/PME). Larger, more correlated dipoles give a larger ε<sub>r</sub>.</p>
<h4>4. Polarisation penalty</h4>
${tex(String.raw`E_{\text{pol}} = \frac{(\mu_{\text{liq}} - \mu_{\text{gas}})^2}{2\alpha}`, true)}
<p>SPC/E adds this self-energy correction (about 5 kJ/mol) to fit the heat of vaporisation, which is what the "E" stands for.</p>`,

  quiz: [
    {
      q: 'What is the virtual site M in TIP4P-type models?',
      options: ['A fourth atom with mass', 'A massless site on the H–O–H bisector carrying the negative charge', 'A lone-pair site with LJ parameters', 'The centre of mass of the molecule'],
      answer: 1,
      why: 'M carries the oxygen\'s negative charge slightly towards the H atoms. It has no mass; its forces are distributed onto O and H.',
    },
    {
      q: 'Why are the dipole moments of fixed-charge water models (about 2.3–2.5 D) larger than the gas-phase value of 1.85 D?',
      options: ['A fitting error', 'To include, on average, the polarisation by neighbours in the liquid', 'Because the molecules are rigid', 'To make the molecule lighter'],
      answer: 1,
      why: 'In the liquid, each molecule\'s dipole is enhanced by its neighbours. Fixed charges cannot respond, so they are set to an effective, boosted value.',
    },
    {
      q: 'You run a protein simulation with the AMBER ff14SB force field. Which water model should you use by default?',
      options: ['Any model gives the same result', 'TIP3P (or a model validated with that force field, such as OPC)', 'TIP5P', 'No water: use vacuum'],
      answer: 1,
      why: 'The protein parameters were fitted together with TIP3P. Changing the water changes the balance of protein–water interactions unless the combination has been validated.',
    },
  ],
};
