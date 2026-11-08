import { tex } from '../tex.js';

export default {
  id: 'polarizable',
  summary: 'Fixed-charge force fields give every atom one charge forever. Real electron clouds shift in an electric field. Polarizable models add this response, as induced dipoles, Drude oscillators or fluctuating charges, at two to four times the cost.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/polarizable.mp4', poster: 'videos/polarizable.jpg' },
  sim: 'polarizable',
  learn: [
    'Why fixed charges miss induction, and how big the induction energy is near an ion',
    'Why polarization is many-body (not pairwise additive), and the polarization catastrophe',
    'The three common models (induced dipoles, Drude oscillators, fluctuating charges), their cost and when to use them',
  ],
  simHint: 'In "1 ion · 3 sites", drag the ion around and compare the many-body energy with the sum of pair terms. Then switch Thole off and move the slider below 1.45 Å: the dipoles run away. Switch to "3 ions · 1 site" and drag the ions to the same side, then spread them evenly: the pair sum hardly changes, but the true energy goes from large to almost zero. Try both with Drude too.',

  idea: `
<p>In a fixed-charge force field (Lessons 4.3 and 4.5) each atom carries one partial charge, chosen once. But electrons respond to their surroundings: an electric field <b>E</b> pulls the electron cloud one way and the nucleus the other, creating an <b>induced dipole</b></p>
<p style="text-align:center"><b>μ</b><sub>ind</sub> = α <b>E</b></p>
<p>where α is the <b>polarizability</b> (about 1.45 Å³ for a water molecule, 0.03 Å³ for Li⁺, 4.5 Å³ for Cl⁻). The induced dipole lowers the energy by ½αE². Next to a Li⁺ ion at 2 Å this is about −15 kcal/mol, as large as a hydrogen bond, but it falls off as r⁻⁴ and is under 0.2 kcal/mol at 6 Å (code tab, Part 1).</p>
<p><b>Polarization is many-body.</b> The induced dipole depends on the <i>total</i> field, and the energy goes as the square of it. Three ions on the same side of a molecule give −41 kcal/mol; the same three ions spread evenly around it give zero, because their fields cancel. A sum of pair terms gives −31 in both cases (Part 2). Induced dipoles also polarize each other (<b>mutual polarization</b>), which is why water's dipole grows from 1.85 D in the gas to about 2.5–3 D in the liquid.</p>
<p><b>The polarization catastrophe.</b> Two point dipoles that get too close keep reinforcing each other, and below r = (4αα′)<sup>1/6</sup> (1.43 Å for two water-like sites in line) the induced dipoles become infinite (Part 3). Real electron clouds overlap and the field is weaker, so models <b>damp</b> the short-range dipole field, most often with Thole's smeared-charge functions.</p>
<p><b>Three ways to put polarization into MD:</b></p>
<ul>
<li><b>Induced point dipoles</b> (AMOEBA, AMBER ff02/ipolq): solve μ = α(E₀ + Tμ) every step, iteratively to self-consistency (or propagate the dipoles with an extended Lagrangian).</li>
<li><b>Drude oscillators</b> (CHARMM Drude, SWM4-NDP water, OPLS/CL&amp;Pol for ionic liquids): attach a charged massive particle to each atom with a spring; α = q<sub>D</sub>²/k<sub>D</sub>. It only needs ordinary charges and springs, so any MD code can run it. The Drude particles are kept "cold" (about 1 K) with a separate thermostat, which is close to keeping them relaxed.</li>
<li><b>Fluctuating charges</b> (CHARMM FQ): charges flow within each molecule to equalize electronegativity. Cheap, but no out-of-plane response.</li>
</ul>
<p><b>Cost and when to use it.</b> Polarizable runs cost about 2–4× a fixed-charge run (more charges, SCF iterations or a smaller time step). They matter most where fields are strong or change a lot: ions and ion selectivity, interfaces (air–water, electrodes), ionic liquids and molten salts, and spectroscopy. A cheap mean-field shortcut for ionic liquids and ions is to scale the ionic charges by 0.7–0.8 (the "electronic continuum correction").</p>
<div class="callout"><b>Key idea.</b> Polarization is the response of the electron cloud to the local field: μ = αE, energy −½αE². It is many-body and short-ranged, needs damping at short range, and costs 2–4× more. Use it when the field around your molecules changes a lot.</div>`,

  algorithm: `
<p>One MD step with induced point dipoles (Drude version below):</p>
<pre class="pseudo">E0_i = field at site i from all permanent charges (with PME)
μ_i = α_i E0_i                              # initial guess (or last step's dipoles)
repeat:                                      # self-consistent field (SCF)
    E_i = E0_i + Σ_j T_ij(damped) μ_j
    μ_i_new = α_i E_i
until max |μ_new - μ| < 1e-6 D (typically 5-10 iterations)
U_pol = -½ Σ_i μ_i · E0_i
forces = permanent + dipole terms; integrate as usual

Drude (extended Lagrangian, e.g. LAMMPS fix drude, GROMACS, OpenMM):
each polarizable atom = core (q - q_D) + Drude particle (q_D, mass ~0.4 u) on a spring k_D
integrate the centre of mass of each pair at T (normal thermostat)
integrate the core-Drude relative motion at T_D ~ 1 K (second, cold thermostat)
time step ~1 fs; add a hard wall (~0.2 Å) so a Drude cannot fly off</pre>
<ol class="steps">
  <li>Use PME for the permanent and induced fields (codes extend PME to dipoles).</li>
  <li>Exclude or screen interactions between dipoles on bonded neighbours (1–2, 1–3); Thole damping for the rest.</li>
  <li>For SCF, a tight tolerance matters: a loose one leaks energy and makes NVE runs drift.</li>
  <li>Software: AMOEBA in Tinker, Tinker-HP, OpenMM, AMBER; Drude in CHARMM, NAMD, GROMACS (drude branch), OpenMM and LAMMPS (<code>fix drude</code>, <code>pair thole</code>). Parameters: CHARMM Drude FF 2019, SWM4-NDP water, CL&amp;Pol for ionic liquids.</li>
</ol>`,

  code: {
    note: 'The four pieces of physics in this lesson, with numbers: the induction energy of a Li+ ion near a water-like site, many-body non-additivity, mutual polarization with and without Thole damping, and the same response from a Drude oscillator. Run <code>pip install numpy</code>, then <code>python polarizable.py</code> (instant).',
    src: `import numpy as np
K = 332.06                     # kcal/mol * Angstrom / e^2  (Coulomb constant)
ALPHA = 1.45                   # polarizability of a water-like site, Angstrom^3

# ---------- Part 1: induction energy of an ion near a polarizable site ----------
print("Part 1: Li+ near a neutral polarizable site (alpha = 1.45 A^3)")
print("   r (A)    field E (e/A^2)   induced dipole (D)   E_pol = -1/2 alpha E^2 (kcal/mol)")
for r in [1.8, 2.0, 2.5, 3.0, 4.0, 6.0]:
    E = 1.0 / r**2; mu = ALPHA * E
    print(f"   {r:4.1f}       {E:.4f}             {mu * 4.803:5.2f}               {-0.5 * ALPHA * E**2 * K:8.2f}")
print("   (1 e*A = 4.803 D; falls off as r^-4, so short-ranged but strong at contact)")

# ---------- Part 2: many-body: induction is not pairwise additive ----------
print("\\nPart 2: three Li+ ions at 2.2 A around one polarizable site")
for name, angles in [("all on one side (60 deg apart)", [0, 60, 120]), ("symmetric (120 deg apart)", [0, 120, 240])]:
    Etot = np.zeros(2); pair = 0.0
    for a in np.radians(angles):
        R = 2.2 * np.array([np.cos(a), np.sin(a)]); Ei = -R / 2.2**3      # field at the site from this ion
        Etot += Ei; pair += -0.5 * ALPHA * (Ei @ Ei) * K
    print(f"   {name:31s}: sum of pair terms {pair:7.2f}, true many-body E_pol {-0.5 * ALPHA * (Etot @ Etot) * K:7.2f} kcal/mol")

# ---------- Part 3: mutual polarization and the polarization catastrophe ----------
def thole(r, a1, a2, a=0.39):
    u3 = (r / (a1 * a2) ** (1 / 6)) ** 3
    return 1 - np.exp(-a * u3), 1 - (1 + a * u3) * np.exp(-a * u3)        # lambda3, lambda5

def induced(r, E0, damp):
    """Two sites on the x axis in a uniform field E0 along x: solve mu = alpha (E0 + T mu) self-consistently."""
    l3, l5 = thole(r, ALPHA, ALPHA) if damp else (1.0, 1.0)
    Txx = (3 * l5 - l3) / r**3                 # field along x at one site from a unit x-dipole at the other
    # iterate (as MD codes do) -- diverges when alpha * Txx >= 1
    mu = ALPHA * E0
    for it in range(200):
        new = ALPHA * (E0 + Txx * mu)
        if abs(new) > 1e6: return np.inf
        if abs(new - mu) < 1e-12: break
        mu = new
    return mu / (ALPHA * E0)

print("\\nPart 3: two polarizable sites in a field along their axis: dipole / (alpha E0)")
print("   r (A)    no damping      Thole damping")
for r in [3.0, 2.0, 1.6, 1.45, 1.3, 1.0]:
    print(f"   {r:4.2f}     {induced(r, 0.01, False):10.3f}      {induced(r, 0.01, True):10.3f}")
print(f"   without damping the dipoles blow up below r = (2 alpha)^(1/3) = {(2 * ALPHA) ** (1 / 3):.2f} A")

# ---------- Part 4: the same physics with a Drude oscillator ----------
kD = 1000.0                                      # spring constant, kcal/mol/A^2  (U = 1/2 kD d^2)
qD = -np.sqrt(ALPHA * kD / K)                    # alpha = qD^2 K / kD
print(f"\\nPart 4: Drude oscillator with kD = {kD:.0f} kcal/mol/A^2 needs a Drude charge qD = {qD:.3f} e")
for r in [2.0, 3.0, 4.0]:
    d = 0.0                                      # relax the Drude particle in the field of Li+ (on the x axis)
    for it in range(200):
        x = -r + d                               # Li+ at origin, site at x = -r, Drude at -r + d
        F = K * qD * x / abs(x)**3               # Coulomb force on the Drude from Li+ (charge +1)
        d_new = F / kD
        if abs(d_new - d) < 1e-12: break
        d = d_new
    core = K * (-qD) / r                         # core keeps charge -qD so the site is neutral overall
    E_drude = 0.5 * kD * d**2 + K * qD / abs(-r + d) + core
    print(f"   r = {r:.1f} A: Drude moves {d * 100:6.2f} pm, E_pol(Drude) = {E_drude:7.2f}  vs  -1/2 alpha E^2 = {-0.5 * ALPHA / r**4 * K:7.2f} kcal/mol")
`,
  },

  maths: () => `
<h4>1. Induced dipole and its energy</h4>
${tex(String.raw`\boldsymbol\mu = \alpha\,\mathbf E, \qquad U_{\text{pol}} = -\tfrac12\,\alpha E^2 = -\tfrac12\,\boldsymbol\mu\cdot\mathbf E`, true)}
<p>The −½ (not −1) is because creating the dipole costs +½αE² of internal energy. For an ion of charge q at distance r, E = q/r², so U<sub>pol</sub> = −αq²/2r⁴.</p>
<h4>2. Mutual polarization</h4>
${tex(String.raw`\boldsymbol\mu_i = \alpha_i\Big(\mathbf E^0_i + \sum_{j\ne i}\mathbf T_{ij}\,\boldsymbol\mu_j\Big), \qquad \mathbf T_{ij} = \frac{3\lambda_5\,\hat{\mathbf r}\hat{\mathbf r} - \lambda_3\,\mathbf I}{r^3}`, true)}
${tex(String.raw`U_{\text{pol}} = -\tfrac12\sum_i \boldsymbol\mu_i\cdot\mathbf E^0_i`, true)}
<p>This is a linear system (α⁻¹ − T)μ = E⁰. It has a stable solution only while α⁻¹ − T is positive definite. For two equal sites in line, T = 2/r³ and μ = αE⁰/(1 − 2α/r³), which diverges at r³ = 2α: the polarization catastrophe.</p>
<h4>3. Thole damping (exponential, as in AMOEBA)</h4>
${tex(String.raw`u = \frac{r}{(\alpha_i\alpha_j)^{1/6}}, \qquad \lambda_3 = 1 - e^{-a u^3}, \qquad \lambda_5 = 1 - (1 + a u^3)\,e^{-a u^3}, \qquad a \approx 0.39`, true)}
<p>At large r both λ → 1 (no change); at short range the dipole field goes smoothly to zero.</p>
<h4>4. Drude oscillator</h4>
${tex(String.raw`k_D\,\mathbf d = q_D\,\mathbf E \;\Rightarrow\; \boldsymbol\mu = q_D\mathbf d = \frac{q_D^2}{k_D}\,\mathbf E, \qquad \alpha = \frac{q_D^2}{k_D}`, true)}
<p>With k<sub>D</sub> = 1000 kcal mol⁻¹ Å⁻² (CHARMM) and α = 1.45 Å³, q<sub>D</sub> = −2.09 e and the Drude particle moves only about 0.1–0.2 Å. In a uniform field the Drude and point-dipole models are identical; close to an ion they differ slightly because the field changes across the displacement.</p>`,

  quiz: [
    {
      q: 'Three identical cations sit at the same distance from a polarizable molecule, spread evenly around it (120° apart). What is the induction energy of the molecule?',
      options: ['Three times the energy for one cation', 'About zero, because the fields cancel at the molecule', 'Half the energy for one cation', 'Infinite'],
      answer: 1,
      why: 'Induction depends on the total field squared, −½α|ΣE|². The three fields cancel, so the energy is about zero, while a sum of pair terms would give three times the one-ion value. This is the many-body nature of polarization.',
    },
    {
      q: 'In a Drude model the polarizability is α = q_D²/k_D. What happens if you make the spring stiffer (larger k_D) but keep α the same?',
      options: ['The Drude charge must be larger, and the Drude particle moves less', 'The Drude charge must be smaller', 'The polarizability goes up', 'Nothing changes at all'],
      answer: 0,
      why: 'To keep q_D²/k_D fixed, q_D grows as √k_D, while the displacement d = q_D E/k_D shrinks. Stiff springs keep the Drude particle close to its core, so the model behaves like a point dipole.',
    },
    {
      q: 'Why do polarizable force fields damp the dipole–dipole interaction at short range (for example with Thole functions)?',
      options: ['To make the code faster', 'Point dipoles that get too close can reinforce each other without limit (polarization catastrophe)', 'To remove the long-range part', 'Because PME requires it'],
      answer: 1,
      why: 'Without damping, two point dipoles closer than about (4α²)^(1/6) have no stable solution and the induced dipoles diverge. Real electron clouds overlap and screen each other; damping mimics that.',
    },
  ],
};
