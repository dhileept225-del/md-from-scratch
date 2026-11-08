import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-molecules',
  lab: true,
  summary: 'So far our atoms only felt Lennard-Jones forces. A molecule needs more: bonds, angles and dihedrals, plus rules for which atoms inside a molecule may also feel LJ. We take one butane molecule, four beads, and compute every term by hand. The numbers show why each term is needed, what the force-field choices (1–4 scaling, combining rules) do to the result, and why constraining bonds lets us take three times longer time steps.',
  minutes: 40,
  sim: 'lab4',
  learn: [
    'Compute bond, angle and dihedral energies and forces from coordinates, with only dot and cross products',
    'See from the numbers why bonded neighbours are excluded from LJ, and how 1–4 scaling and combining rules change results',
    'Iterate SHAKE by hand, and measure how much longer a time step rigid bonds allow',
  ],
  simHint: 'Langevin mode: 40 butane molecules at the chosen temperature. The Newman view looks along the central bond of molecule 1; the bars collect the dihedral of all 40, the amber line is the exact answer from step 4. Switch the 1–4 rule and watch gauche disappear. Time-step mode: pick flexible or rigid bonds and a time step, and read the energy fluctuation; flexible bonds blow up at 10 fs, rigid ones do not.',

  setup: `
<p>Our molecule is <b>n-butane</b>, CH₃–CH₂–CH₂–CH₃, in the <b>united-atom</b> picture: each carbon together with its hydrogens is one bead. Four beads, three bonds, two angles, one dihedral. Parameters are those of the TraPPE force field. From here on we use the units of GROMACS, so the numbers can be compared with real input files.</p>
<div class="calc">beads     CH₃ (15.035)  CH₂ (14.027)  CH₂  CH₃  g/mol
units     length nm, time ps, energy kJ/mol
          force kJ/mol/nm;  F/m comes out in nm/ps²
T = 300 K:  k_B T = 0.0083145 × 300 = <b>2.494 kJ/mol</b>

bond      U = ½ k_b (r − 0.154)²,  k_b = 224 262 kJ/mol/nm²
angle     U = ½ k_θ (θ − 114°)²,   k_θ = 519.66 kJ/mol/rad²
dihedral  U = c₁(1 + cos φ) + c₂(1 − cos 2φ) + c₃(1 + cos 3φ)
          c₁ = 2.9519, c₂ = −0.5670, c₃ = 6.5794 kJ/mol
LJ        CH₃: ε = 0.8148, σ = 0.375
          CH₂: ε = 0.3825, σ = 0.395</div>
<p>(TraPPE keeps the bonds rigid; for the flexible runs we give them a typical C–C spring.) Every energy below is compared with k_B T = 2.494 kJ/mol: anything much larger is almost never visited at 300 K, anything smaller is visited all the time.</p>`,

  steps: [
    {
      title: 'Why a molecule needs bonded terms, and exclusions',
      why: 'If LJ were enough, we could just put four beads together and let them stick. Let us check what LJ says about two carbons at the length of a real C–C bond.',
      work: `<div class="calc">two CH₂ beads, LJ only
LJ minimum at 2^(1/6) σ = 1.1225 × 0.395 = <b>0.4434 nm</b>
a C–C bond is              <b>0.154 nm</b>

LJ at 0.154 nm:  σ/r = 2.565
U = 4 × 0.3825 × (81 081 − 284.7) = <b>123 607 kJ/mol</b>

1–3 pair (two bonds apart), distance
2 × 0.154 × sin(57°) = 0.2583 nm:  U_LJ = <b>230.5 kJ/mol</b></div>`,
      see: 'LJ describes how separate molecules touch, at about 0.4 nm. A chemical bond is a much shorter, much stronger thing: at 0.154 nm LJ gives 123 607 kJ/mol, 50 000 times k_B T, and would blow the molecule apart in one step. Two consequences. First, the shape of a molecule is held by its own terms: bonds, angles, dihedrals. Second, atoms 1 bond apart (1–2) and 2 bonds apart (1–3) must be <b>excluded</b> from LJ, otherwise 230 kJ/mol of fake repulsion fights the angle term. Every MD code keeps such an exclusion list.',
    },
    {
      title: 'The bond: a stiff spring that sets the time step',
      why: 'The bond keeps two beads near 0.154 nm. The simplest description is a spring. We need its force for the dynamics, and its vibration period, because the time step has to follow the fastest motion in the system (Lab 1).',
      work: `<div class="calc">stretched to r = 0.160 nm:
U = ½ × 224 262 × 0.006² = <b>4.037 kJ/mol</b>
F = −224 262 × 0.006     = <b>−1345.6 kJ/mol/nm</b>

thermal stretch √(k_B T / k_b) = √(2.494/224 262)
                              = 0.00334 nm  (2.2 %)

vibration:  reduced mass μ = 14.027 × 14.027 / 28.054
                           = 7.0135 g/mol
ω = √(k_b/μ) = √(224 262 / 7.0135) = 178.8 /ps
period = 2π/ω = <b>35.1 fs</b>   (measured in a run: 37 fs)

all-atom C–H bond, k = 284 512, μ = 0.9300:
ω = 553.1 /ps,  period = <b>11.4 fs</b></div>`,
      see: 'At 300 K the bond only wiggles by 2 %, so for most purposes it is effectively rigid, yet it is by far the fastest motion: one period every 35 fs. With about 15–20 steps per period (Lab 1), this alone limits us to about 2 fs. In an all-atom model the C–H bond vibrates every 11 fs, forcing 0.5–1 fs. We spend most of the computer time following a motion we do not care about. Step 7 shows the fix.',
    },
    {
      title: 'The angle: dot product, energy, force',
      why: 'Bonds alone would let the chain fold flat onto itself. The angle term keeps C–C–C near 114°. To use it we must get the angle from coordinates, then turn the energy into forces on the atoms.',
      work: `<div class="calc">atoms: A (0.154, 0, 0),  B (0, 0, 0),
       C (−0.077, 0.13337, 0)     (B is the middle)
u = A − B,  v = C − B
u·v = 0.154 × (−0.077) = −0.011858
|u||v| = 0.154 × 0.154 = 0.023716
cos θ = −0.5000  ->  θ = <b>120.0°</b>

dθ = 6.0° = 0.10472 rad
U = ½ × 519.66 × 0.10472² = <b>2.849 kJ/mol</b>
dU/dθ = 519.66 × 0.10472 = 54.42 kJ/mol/rad
force on an end atom = 54.42 / 0.154 = <b>353.4 kJ/mol/nm</b>
(perpendicular to its bond, closing the angle)

thermal spread √(k_B T / k_θ) = 0.0693 rad = <b>3.97°</b>
vibration period measured in a run: <b>91 fs</b></div>`,
      see: 'Turning an angle into forces needs nothing beyond the derivative and a lever arm: the torque dU/dθ divided by the bond length. At 300 K the angle wanders by about 4°, much more than the bond (2 %), and it vibrates 2.5 times more slowly (91 fs vs 37 fs). That ratio is exactly what we gain in step 7 when the bonds are frozen and the angle becomes the fastest motion.',
    },
    {
      title: 'The dihedral: where the chemistry is',
      why: 'With bonds and angles fixed, butane can still twist around its central bond. This twist decides whether the chain is straight (trans) or kinked (gauche), and in polymers and lipids it controls the whole shape. The dihedral angle φ needs four atoms and two cross products.',
      work: `<div class="calc">gauche butane (bonds 0.154, angles 114°), nm:
atom 1 (0.00000, 0.00000, 0.00000)
atom 2 (0.15400, 0.00000, 0.00000)
atom 3 (0.21664, 0.14069, 0.00000)
atom 4 (0.17785, 0.22652, 0.12184)

b1 = 2 − 1 = (0.154, 0, 0)
b2 = 3 − 2 = (0.06264, 0.14069, 0)
b3 = 4 − 3 = (−0.03878, 0.08583, 0.12184)
m = b1 × b2 = (0, 0, 0.02167)      plane of 1-2-3
n = b2 × b3 = (0.01714, −0.00763, 0.01083)  plane of 2-3-4
cos φ = m·n / (|m||n|)
      = 0.0002347 / (0.021666 × 0.021666) = 0.5000
φ = <b>60.0°</b>  (gauche)

trans    φ = 180°:  U = <b>0.000</b>
gauche   φ =  60°:  2.9519×1.5 − 0.5670×1.5 = <b>3.577</b>
barrier  φ = 120°:  1.476 − 0.851 + 13.159 = <b>13.784</b>
cis      φ =   0°:  <b>19.063</b> kJ/mol

each gauche / trans = e^(−3.577/2.494) = 0.238
two gauche states:  trans fraction = 1/(1 + 2 × 0.238)
                                   = <b>0.677</b>
whole curve, summed in 1° steps:   <b>0.663</b></div>
<p>How fast does the molecule switch? An escape over a barrier happens at a rate of (vibration frequency) × e^(−barrier/k_B T):</p>
<div class="calc">curvature at trans  c₁ + 4c₂ + 9c₃ = 59.90 kJ/mol/rad²
moment of inertia   I = 0.1488 g/mol nm²
frequency ν = √(59.90/0.1488)/2π = 3.19 /ps
rate = 2 × 3.19 × e^(−13.784/2.494) = 2 × 3.19 × 0.00398
     = 0.0254 /ps  ->  one escape every <b>39 ps</b></div>
<div class="tbl"><table>${head(['10-ps window', '1', '2', '3', '4', '5', '7', '10'])}
${row(['start all trans', '0.905', '0.774', '0.723', '0.717', '0.740', '0.719', '0.675'])}
${row(['start all gauche', '0.159', '0.364', '0.520', '0.572', '0.605', '0.622', '0.651'])}</table></div>
<p>Trans fraction in a Langevin run of 200 molecules at 300 K; one jump every 23 ps per molecule.</p>`,
      see: 'The dihedral needs only two cross products, and its three-term formula gives the whole story of butane: trans is lowest, gauche is 3.6 kJ/mol higher (1.4 k_B T, so still common: about one third of the molecules), cis is almost never visited. The full sum (0.663) differs a little from the three-state count (0.677) because the wells have different widths. Limitation, by numbers: the molecule switches only every few tens of picoseconds. Starting from all-trans or all-gauche, 200 molecules need about 100 ps before both runs agree on 0.65–0.68. A longer chain, or a higher barrier, makes this nanoseconds to microseconds: sampling dihedrals is often the real bottleneck of a simulation (Chapter 7, metadynamics).',
    },
    {
      title: '1–4 pairs: exclude, scale or keep?',
      why: 'Atoms 1 and 4 are three bonds apart. They are not held at a fixed distance by bonds and angles: their distance changes with φ. Should they feel LJ? Different force fields answer differently, and the numbers show this is not a detail.',
      work: `<div class="calc">CH₃ ... CH₃ distance and LJ (ε = 0.8148, σ = 0.375):
trans   r₁₄ = 0.3964 nm   U_LJ =  <b>−0.662</b>
gauche  r₁₄ = 0.3127 nm   U_LJ = <b>+19.135</b>
cis     r₁₄ = 0.2793 nm   U_LJ = <b>+92.870</b>  kJ/mol</div>
<div class="tbl"><table>${head(['1–4 LJ with the TraPPE dihedral', 'gauche − trans', 'trans fraction'])}
${row(['excluded (TraPPE\'s rule)', '3.58 kJ/mol', '<b>0.663</b>'])}
${row(['scaled by 0.5 (OPLS\'s rule)', '13.48 kJ/mol', '<b>0.957</b>'])}
${row(['full', '23.37 kJ/mol', '<b>0.988</b>'])}</table></div>`,
      see: 'In gauche the two CH₃ ends are only 0.313 nm apart, inside σ = 0.375, so their LJ is strongly repulsive (+19 kJ/mol). TraPPE fitted its dihedral with this pair excluded; OPLS fitted its own dihedral with this pair at half strength. Mixing the two, TraPPE\'s dihedral with OPLS\'s rule, turns butane 96 % trans, almost no gauche left: the chains come out far too straight. The lesson: <b>a force field is a package</b>. Dihedral parameters, 1–4 scaling and exclusions were fitted together, and must be used together. In GROMACS this is <code>gen-pairs</code>, <code>fudgeLJ</code> and <code>fudgeQQ</code>; in LAMMPS <code>special_bonds</code>.',
    },
    {
      title: 'Combining rules: CH₃ meets CH₂',
      why: 'We have LJ parameters for CH₃–CH₃ and CH₂–CH₂, but a simulation also needs CH₃–CH₂, and with many atom types there are far too many pairs to fit each one. Combining rules make the cross parameters from the pure ones.',
      work: `<div class="calc">Lorentz–Berthelot (TraPPE, AMBER, CHARMM):
σ = (0.375 + 0.395)/2 = <b>0.3850</b>
ε = √(0.8148 × 0.3825) = <b>0.5582</b>
geometric σ (OPLS):  √(0.375 × 0.395) = <b>0.3849</b>

very different sizes, σ = 0.15 and 0.40 nm (ε = 0.5):
arithmetic σ = 0.2750,  geometric σ = 0.2449  (12 % apart)
at r = 0.25:  U = <b>+2.734</b>  vs  <b>−0.204</b> kJ/mol
contact distance 2^(1/6)σ:  0.3087 vs 0.2749 nm</div>`,
      see: 'For similar atoms like CH₃ and CH₂, the two rules agree to 0.0001 nm and the choice does not matter. For very different sizes, such as a small cation next to a large anion, they disagree by 12 % in σ: at 0.25 nm one rule says strong repulsion, the other mild attraction, and the contact distance moves by 0.034 nm. That directly shifts the radial distribution function and the binding. Like 1–4 scaling, the combining rule belongs to the force field (GROMACS <code>comb-rule</code>, LAMMPS <code>pair_modify mix</code>). Ion–ion and ion–solvent pairs are often fitted separately for exactly this reason.',
    },
    {
      title: 'Constraints: SHAKE by hand, and the time step it buys',
      why: 'The bond is stiff (2 % wiggle) and fast (35 fs). If we freeze it at 0.154 nm, the fastest motion becomes the angle (91 fs) and the time step can grow. SHAKE does this: after each normal step, it moves atoms along their old bond directions until every bond length is right again.',
      work: `<div class="calc">one bond, old direction along x, masses 14.027 each
after an unconstrained step:
atoms (0.001, 0.002), (0.1565, −0.001):  bond 0.1555289
g = (|d|² − b₀²) / (2 (1/m₁ + 1/m₂) d·d_old)
  = 0.0004732 / (2 × 0.14258 × 0.023947) = 0.069302
move: atom 1 by +g/m₁ along d_old, atom 2 by −g/m₂

iteration 0  bond 0.1555289   error +1.5×10⁻³ nm
iteration 1  bond 0.1540075   error +7.5×10⁻⁶
iteration 2  bond 0.1540000   error +1.8×10⁻¹⁰
butane, 3 coupled bonds, to 10⁻¹⁰: <b>about 10 iterations</b>

degrees of freedom: 4 beads × 3 = 12, minus 3 bonds = <b>9</b>
(T computed with 12 would read 25 % low)</div>
<div class="tbl"><table>${head(['time step', 'flexible bonds', 'rigid bonds (SHAKE)'])}
${row(['1 fs', '0.0117', '—'])}
${row(['2 fs', '<b>0.0482</b>', '0.0042'])}
${row(['4 fs', '0.2173', '0.0169'])}
${row(['6 fs', '0.6119', '<b>0.0392</b>'])}
${row(['8 fs', '1.7311', '0.0729'])}
${row(['10–12 fs', 'unstable', '0.1887 (12 fs)'])}</table></div>
<p>Energy fluctuation in NVE (kJ/mol, 20 molecules × 10 ps): smaller means a more accurate trajectory.</p>`,
      see: 'The error falls very fast (1.5×10⁻³ → 7.5×10⁻⁶ → 1.8×10⁻¹⁰ nm), so two or three iterations are enough for one bond; coupled bonds in a chain need about 10, because fixing one bond disturbs its neighbours. The payoff: rigid bonds at <b>6 fs</b> are as accurate as flexible bonds at <b>2 fs</b> (0.039 vs 0.048), three times fewer steps for the same physics, and flexible bonds blow up at 10 fs while rigid ones still run at 12 fs. This matches the period ratio of step 3 (91/37 ≈ 2.5). In all-atom models, constraining only the bonds to H already doubles the time step (1 -> 2 fs). Limitations: the iterations cost time and are hard to parallelise (GROMACS prefers LINCS, and SETTLE for water); constraining angles in long chains converges badly; and every constraint removes a degree of freedom, which must be subtracted when computing the temperature.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['ingredient', 'why (number)', 'limitation / caution (number)'])}
${row(['bonds + 1–2, 1–3 exclusions', 'LJ at 0.154 nm = 123 607 kJ/mol', 'bond period 35 fs (C–H 11 fs) limits dt'])}
${row(['angle', 'keeps 114°, spread 4°', 'period 91 fs: next limit on dt'])}
${row(['dihedral', 'trans 0.66, gauche 0.34', 'one jump per ~30 ps: slow sampling'])}
${row(['1–4 rule', 'must match the dihedral fit', 'wrong rule: trans 0.96 instead of 0.66'])}
${row(['combining rule', 'cross terms for free', 'σ 12 % apart for unlike sizes'])}
${row(['SHAKE / LINCS', 'dt 2 -> 6 fs at equal accuracy', '~10 iterations; N_f 12 -> 9'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the LJ test, the bond and angle by hand, the dihedral from coordinates, the populations and escape rate, the Langevin windows, the 1–4 and combining-rule tables, SHAKE iterations and the time-step test. Run <code>pip install numpy</code>, then <code>python lab_molecules.py</code> (about a minute).',
    src: `# Lab 4.6: a molecule by hand. United-atom n-butane CH3-CH2-CH2-CH3 (4 beads), GROMACS units:
# length nm, mass g/mol, time ps, energy kJ/mol (force kJ/mol/nm; F/m comes out in nm/ps^2).
# Angle, dihedral and LJ parameters: TraPPE-UA. Bond spring: a typical C-C value (TraPPE itself keeps bonds rigid).
import numpy as np, time
kB = 0.0083145; T = 300.0; kT = kB * T                 # kJ/mol/K;  kT = 2.494 kJ/mol at 300 K
M = np.array([15.035, 14.027, 14.027, 15.035])         # CH3, CH2, CH2, CH3
B0, KB = 0.154, 224262.4                               # U = 1/2 KB (r - B0)^2
TH0, KTH = np.radians(114.0), 62500 * kB               # U = 1/2 KTH (theta - TH0)^2
C1, C2, C3 = 355.03 * kB, -68.19 * kB, 791.32 * kB     # U = C1(1 + cos phi) + C2(1 - cos 2phi) + C3(1 + cos 3phi), trans = 180
E3, S3, E2, S2 = 98 * kB, 0.375, 46 * kB, 0.395        # LJ of CH3 and CH2
LJ = lambda r, e, s: 4 * e * ((s / r)**12 - (s / r)**6)
dot = lambda a, b: (a * b).sum(-1)
Udih = lambda p: C1 * (1 + np.cos(p)) + C2 * (1 - np.cos(2 * p)) + C3 * (1 + np.cos(3 * p))

def place(phi, b=B0, th=TH0):                          # 4 positions from bond length, bond angle and dihedral (degrees)
    p = np.radians(phi); r1 = np.zeros(3); r2 = np.array([b, 0, 0]); r3 = r2 + b * np.array([-np.cos(th), np.sin(th), 0])
    bc = (r3 - r2) / b; nv = np.cross(r2 - r1, bc); nv /= np.linalg.norm(nv); mv = np.cross(nv, bc)
    return np.array([r1, r2, r3, r3 + b * (-np.cos(th) * bc + np.sin(th) * np.cos(p) * mv + np.sin(th) * np.sin(p) * nv)])

def energy_forces(x, s14=0.0, bonds=True):             # x has shape (molecules, 4, 3)
    f = np.zeros_like(x); U = np.zeros(x.shape[0])
    if bonds:
        for i in range(3):
            d = x[:, i + 1] - x[:, i]; r = np.linalg.norm(d, axis=-1); U += 0.5 * KB * (r - B0)**2
            g = (-KB * (r - B0) / r)[:, None] * d; f[:, i + 1] += g; f[:, i] -= g
    for j in (1, 2):
        u = x[:, j - 1] - x[:, j]; v = x[:, j + 1] - x[:, j]; lu = np.linalg.norm(u, axis=-1); lv = np.linalg.norm(v, axis=-1)
        c = np.clip(dot(u, v) / lu / lv, -1, 1); th = np.arccos(c); s = np.sqrt(1 - c * c); dU = KTH * (th - TH0); U += 0.5 * KTH * (th - TH0)**2
        fi = (dU / s)[:, None] * (v / (lu * lv)[:, None] - (c / lu**2)[:, None] * u)
        fk = (dU / s)[:, None] * (u / (lu * lv)[:, None] - (c / lv**2)[:, None] * v)
        f[:, j - 1] += fi; f[:, j + 1] += fk; f[:, j] -= fi + fk
    b1, b2, b3 = x[:, 1] - x[:, 0], x[:, 2] - x[:, 1], x[:, 3] - x[:, 2]; m, n = np.cross(b1, b2), np.cross(b2, b3)
    lb2 = np.linalg.norm(b2, axis=-1); phi = np.arctan2(lb2 * dot(b1, n), dot(m, n)); U += Udih(phi)
    dU = -C1 * np.sin(phi) + 2 * C2 * np.sin(2 * phi) - 3 * C3 * np.sin(3 * phi)
    f1 = (dU * lb2 / dot(m, m))[:, None] * m; f4 = -(dU * lb2 / dot(n, n))[:, None] * n
    p = (dot(b1, b2) / lb2**2)[:, None]; q = (dot(b3, b2) / lb2**2)[:, None]
    f[:, 0] += f1; f[:, 3] += f4; f[:, 1] += -f1 - p * f1 + q * f4; f[:, 2] += -f4 + p * f1 - q * f4
    if s14:                                            # LJ between the two CH3 ends (1-4 pair), scaled by s14
        d = x[:, 3] - x[:, 0]; r2 = dot(d, d); sr6 = (S3**2 / r2)**3; U += s14 * 4 * E3 * (sr6**2 - sr6)
        g = (s14 * 24 * E3 * (2 * sr6**2 - sr6) / r2)[:, None] * d; f[:, 3] += g; f[:, 0] -= g
    return U, f, phi

# ---------- Step 1: why bonded terms and exclusions ----------
print(f"Step 1: two CH2 beads with LJ only (eps = {E2:.4f} kJ/mol, sigma = {S2} nm)")
print(f"   LJ minimum at 2^(1/6) sigma = {2**(1/6) * S2:.4f} nm, but a C-C bond is {B0} nm")
print(f"   LJ at the bond length: sigma/r = {S2 / B0:.3f}, U = 4 x {E2:.4f} x ({(S2 / B0)**12:.0f} - {(S2 / B0)**6:.1f}) = {LJ(B0, E2, S2):,.0f} kJ/mol")
r13 = 2 * B0 * np.sin(TH0 / 2)
print(f"   1-3 distance 2 b sin(57) = {r13:.4f} nm: U_LJ = {LJ(r13, E2, S2):.1f} kJ/mol (kT = {kT:.3f})")

# ---------- Step 2: the bond ----------
r = 0.160; mu22 = M[1] * M[2] / (M[1] + M[2]); w = np.sqrt(KB / mu22)
print(f"\\nStep 2: bond, KB = {KB} kJ/mol/nm^2. At r = {r}: U = 1/2 x {KB} x {r - B0:.3f}^2 = {0.5 * KB * (r - B0)**2:.3f} kJ/mol, F = -{KB} x {r - B0:.3f} = {-KB * (r - B0):.1f} kJ/mol/nm")
print(f"   thermal stretch sqrt(kT/KB) = sqrt({kT:.3f}/{KB}) = {np.sqrt(kT / KB):.5f} nm ({100 * np.sqrt(kT / KB) / B0:.1f} % of the bond)")
print(f"   CH2-CH2: reduced mass {mu22:.4f} g/mol, omega = sqrt(KB/mu) = {w:.1f} 1/ps, period 2 pi/omega = {2 * np.pi / w * 1000:.1f} fs")
kch, much = 284512.0, 12.011 * 1.008 / (12.011 + 1.008); wch = np.sqrt(kch / much)       # a typical all-atom C-H bond
print(f"   C-H (all-atom, k = {kch:.0f}): reduced mass {much:.4f}, omega = {wch:.1f} 1/ps, period {2 * np.pi / wch * 1000:.1f} fs")

# ---------- Step 3: the angle ----------
xa = np.array([[0.154, 0, 0], [0, 0, 0], [0.154 * np.cos(np.radians(120)), 0.154 * np.sin(np.radians(120)), 0]])
u, v = xa[0] - xa[1], xa[2] - xa[1]; c = dot(u, v) / np.linalg.norm(u) / np.linalg.norm(v); th = np.arccos(c)
print(f"\\nStep 3: angle, KTH = 62500 K x kB = {KTH:.2f} kJ/mol/rad^2, theta0 = 114 deg")
print(f"   atoms {np.round(xa[0], 5)}, {np.round(xa[1], 5)}, {np.round(xa[2], 5)}: u.v = {dot(u, v):.6f}, |u||v| = {np.linalg.norm(u) * np.linalg.norm(v):.6f}, cos = {c:.4f}, theta = {np.degrees(th):.1f} deg")
d = th - TH0
print(f"   d theta = {np.degrees(d):.1f} deg = {d:.5f} rad;  U = 1/2 x {KTH:.2f} x {d:.5f}^2 = {0.5 * KTH * d**2:.3f} kJ/mol;  dU/dtheta = {KTH * d:.2f};  force on an end atom = {KTH * d:.2f}/{B0} = {KTH * d / B0:.1f} kJ/mol/nm")
print(f"   thermal spread sqrt(kT/KTH) = {np.sqrt(kT / KTH):.4f} rad = {np.degrees(np.sqrt(kT / KTH)):.2f} deg")

# ---------- Step 4: the dihedral ----------
x = place(60.0); b1, b2, b3 = x[1] - x[0], x[2] - x[1], x[3] - x[2]; m, n = np.cross(b1, b2), np.cross(b2, b3)
print(f"\\nStep 4: dihedral of a gauche butane (bonds {B0} nm, angles 114 deg); positions in nm:")
for i in range(4): print(f"   atom {i + 1}: ({x[i, 0]:8.5f}, {x[i, 1]:8.5f}, {x[i, 2]:8.5f})")
print(f"   b1 = {np.round(b1, 5)}, b2 = {np.round(b2, 5)}, b3 = {np.round(b3, 5)}")
print(f"   m = b1 x b2 = {np.round(m, 5)},  n = b2 x b3 = {np.round(n, 5)}")
print(f"   cos phi = m.n/(|m||n|) = {dot(m, n):.7f}/({np.linalg.norm(m):.6f} x {np.linalg.norm(n):.6f}) = {dot(m, n) / np.linalg.norm(m) / np.linalg.norm(n):.4f}  ->  phi = {np.degrees(np.arccos(dot(m, n) / np.linalg.norm(m) / np.linalg.norm(n))):.1f} deg")
for ph, name in ((180, 'trans'), (60, 'gauche'), (120, 'barrier'), (0, 'cis')):
    p = np.radians(ph); print(f"   {name:8s} phi = {ph:3d}: U = {C1:.4f}(1{np.cos(p):+.1f}) + {C2:.4f}(1{-np.cos(2 * p):+.1f}) + {C3:.4f}(1{np.cos(3 * p):+.1f}) = {Udih(p):.3f} kJ/mol")
g = np.exp(-Udih(np.radians(60)) / kT)
print(f"   each gauche vs trans: exp(-{Udih(np.radians(60)):.3f}/{kT:.3f}) = {g:.3f};  two gauche states -> trans fraction 1/(1 + 2 x {g:.3f}) = {1 / (1 + 2 * g):.3f}")
PHG = np.radians(np.arange(-179.5, 180, 1.0)); wgt = np.exp(-Udih(PHG) / kT)
print(f"   whole curve summed in 1-degree steps (|phi| > 120 = trans): trans fraction = {wgt[np.abs(PHG) > np.radians(120)].sum() / wgt.sum():.3f}")
k2 = C1 + 4 * C2 + 9 * C3; rp = B0 * np.sin(np.pi - TH0); I = M[0] * rp**2 / 2; nu = np.sqrt(k2 / I) / (2 * np.pi); bar = Udih(np.radians(120))
print(f"   escape from trans: curvature U'' = C1 + 4 C2 + 9 C3 = {k2:.2f} kJ/mol/rad^2, I = m r_perp^2 / 2 = {I:.4f} g/mol nm^2, nu = {nu:.2f} /ps")
print(f"   rate = 2 x nu x exp(-{bar:.3f}/{kT:.3f}) = 2 x {nu:.2f} x {np.exp(-bar / kT):.5f} = {2 * nu * np.exp(-bar / kT):.4f} /ps  ->  one escape every {1 / (2 * nu * np.exp(-bar / kT)):.0f} ps")

def langevin(start, B=200, ps=100.0, dt=0.002, gamma=5.0, seed=3):
    r = np.random.default_rng(seed); x = np.repeat(place(start)[None], B, 0); v = r.normal(size=(B, 4, 3)) * np.sqrt(kT / M)[None, :, None]
    U, f, phi = energy_forces(x); c = np.exp(-gamma * dt); sn = np.sqrt((1 - c * c) * kT / M)[None, :, None]
    st = np.zeros(B, int); win = []; tr = 0; jumps = 0; nwin = int(10 / dt)
    for s in range(int(ps / dt)):
        v += 0.5 * dt * f / M[None, :, None]; x += dt * v; U, f, phi = energy_forces(x); v += 0.5 * dt * f / M[None, :, None]
        v = c * v + sn * r.normal(size=v.shape)
        ap = np.degrees(np.abs(phi)); new = np.where(ap > 150, 1, np.where((ap > 30) & (ap < 90), 2, st))
        jumps += ((new != st) & (st > 0)).sum(); st = new; tr += (ap > 120).mean()
        if (s + 1) % nwin == 0: win.append(tr / nwin); tr = 0
    return win, jumps, B * ps
t0 = time.time(); wt, jt, tt = langevin(180.0); wg, jg, tg = langevin(60.0, seed=4)
print(f"   Langevin, 200 molecules x 100 ps at 300 K ({time.time() - t0:.0f} s). Trans fraction in 10-ps windows:")
print("   start all trans : " + " ".join(f"{q:.3f}" for q in wt))
print("   start all gauche: " + " ".join(f"{q:.3f}" for q in wg))
print(f"   jumps between trans and gauche: one every {(tt + tg) / (jt + jg):.0f} ps per molecule")

# ---------- Step 5: 1-4 interactions ----------
print(f"\\nStep 5: the two CH3 ends (1-4 pair), LJ eps = {E3:.4f} kJ/mol, sigma = {S3} nm")
for ph, name in ((180, 'trans'), (60, 'gauche'), (0, 'cis')):
    xx = place(ph); r14 = np.linalg.norm(xx[3] - xx[0]); print(f"   {name:7s} r14 = {r14:.4f} nm  U_LJ = {LJ(r14, E3, S3):8.3f} kJ/mol")
r14 = np.array([np.linalg.norm(place(np.degrees(p))[3] - place(np.degrees(p))[0]) for p in PHG])
for s14, name in ((0.0, 'excluded (TraPPE)'), (0.5, 'scaled 0.5 (OPLS rule)'), (1.0, 'full')):
    Ut = Udih(PHG) + s14 * LJ(r14, E3, S3); w14 = np.exp(-(Ut - Ut.min()) / kT)
    dG = (Udih(np.radians(60)) + s14 * LJ(np.linalg.norm(place(60)[3] - place(60)[0]), E3, S3)) - (s14 * LJ(np.linalg.norm(place(180)[3] - place(180)[0]), E3, S3))
    print(f"   1-4 {name:24s}: gauche - trans = {dG:6.2f} kJ/mol, trans fraction = {w14[np.abs(PHG) > np.radians(120)].sum() / w14.sum():.3f}")

# ---------- Step 6: combining rules ----------
print("\\nStep 6: CH3-CH2 cross parameters")
print(f"   Lorentz-Berthelot: sigma = (0.375 + 0.395)/2 = {(S3 + S2) / 2:.4f}, eps = sqrt({E3:.4f} x {E2:.4f}) = {np.sqrt(E3 * E2):.4f}")
print(f"   geometric (OPLS):  sigma = sqrt(0.375 x 0.395) = {np.sqrt(S3 * S2):.4f}")
sa, sb, e = 0.15, 0.40, 0.5
s_ar, s_ge = (sa + sb) / 2, np.sqrt(sa * sb)
print(f"   very different sizes, sigma 0.15 and 0.40 nm: arithmetic {s_ar:.4f}, geometric {s_ge:.4f} ({100 * (s_ar - s_ge) / s_ge:.0f} % apart)")
for rr in (0.25, 0.28):
    print(f"   at r = {rr}: U(arith) = {LJ(rr, e, s_ar):8.3f}, U(geom) = {LJ(rr, e, s_ge):8.3f} kJ/mol (eps = {e})")
print(f"   contact (minimum) distance 2^(1/6) sigma: {2**(1/6) * s_ar:.4f} vs {2**(1/6) * s_ge:.4f} nm")

# ---------- Step 7: constraints (SHAKE / RATTLE) ----------
iM = 1 / M
def shake(x, xo, v, dt, tol=1e-10):
    for it in range(1000):
        err = 0.0
        for i in range(3):
            dd = x[:, i + 1] - x[:, i]; do = xo[:, i + 1] - xo[:, i]; diff = dot(dd, dd) - B0**2
            gg = diff / (2 * (iM[i] + iM[i + 1]) * dot(dd, do)); x[:, i] += (gg * iM[i])[:, None] * do; x[:, i + 1] -= (gg * iM[i + 1])[:, None] * do
            v[:, i] += (gg * iM[i] / dt)[:, None] * do; v[:, i + 1] -= (gg * iM[i + 1] / dt)[:, None] * do; err = max(err, np.abs(diff).max() / (2 * B0**2))
        if err < tol: return it + 1
    return 1000
def rattle_v(x, v, tol=1e-12):
    for it in range(1000):
        err = 0.0
        for i in range(3):
            dd = x[:, i + 1] - x[:, i]; dv = v[:, i + 1] - v[:, i]; k = dot(dd, dv) / (B0**2 * (iM[i] + iM[i + 1]))
            v[:, i] += (k * iM[i])[:, None] * dd; v[:, i + 1] -= (k * iM[i + 1])[:, None] * dd; err = max(err, np.abs(k).max())
        if err < tol: return

print("\\nStep 7: SHAKE on one bond by hand. Old bond along x, length 0.154; after an unconstrained step the atoms are at")
xo = np.array([[0.0, 0, 0], [0.154, 0, 0]]); xn = np.array([[0.0010, 0.0020, 0], [0.1565, -0.0010, 0]]); mi = M[1:3]
for it in range(4):
    dd = xn[1] - xn[0]; diff = dot(dd, dd) - B0**2; do = xo[1] - xo[0]
    print(f"   iteration {it}: atoms ({xn[0, 0]:.6f}, {xn[0, 1]:.6f}), ({xn[1, 0]:.6f}, {xn[1, 1]:.6f}), bond = {np.linalg.norm(dd):.7f}, error = {np.linalg.norm(dd) - B0:+.2e} nm")
    gg = diff / (2 * (1 / mi[0] + 1 / mi[1]) * dot(dd, do))
    if it == 0: print(f"      g = (|d|^2 - b0^2) / (2 (1/m1 + 1/m2) d.d_old) = {diff:.7f} / (2 x {1 / mi[0] + 1 / mi[1]:.5f} x {dot(dd, do):.6f}) = {gg:.6f}")
    xn[0] += gg / mi[0] * do; xn[1] -= gg / mi[1] * do

def nve(dt, ps=10.0, constrained=False, B=20):
    r = np.random.default_rng(0); x = np.repeat(place(180.0)[None], B, 0); v = r.normal(size=(B, 4, 3)) * np.sqrt(kT / M)[None, :, None]
    if constrained: rattle_v(x, v)
    U, f, _ = energy_forces(x, bonds=not constrained); E = []; its = []
    for s in range(int(round(ps / dt))):
        v += 0.5 * dt * f / M[None, :, None]; xo = x.copy(); x = x + dt * v
        if constrained: its.append(shake(x, xo, v, dt))
        U, f, _ = energy_forces(x, bonds=not constrained); v += 0.5 * dt * f / M[None, :, None]
        if constrained: rattle_v(x, v)
        E.append(U + 0.5 * (M[None, :, None] * v**2).sum((1, 2)))
        if not np.isfinite(E[-1]).all() or np.abs(E[-1] - E[0]).max() > 1e3: return None, its
    return np.array(E).std(0).mean(), its
t0 = time.time()
print("   NVE, 20 molecules x 10 ps: energy fluctuation (kJ/mol; smaller = more accurate). Flexible bonds:")
for dt in (1, 2, 4, 6, 8, 10):
    sE, _ = nve(dt / 1000); print(f"     dt = {dt:2d} fs: " + ("unstable (energy explodes)" if sE is None else f"{sE:.4f}"))
print("   Rigid bonds (SHAKE/RATTLE, tolerance 1e-10):")
for dt in (2, 4, 6, 8, 12):
    sE, its = nve(dt / 1000, constrained=True); print(f"     dt = {dt:2d} fs: {sE:.4f}   SHAKE iterations per step: {np.mean(its):.1f}")
print(f"   ({time.time() - t0:.0f} s)")

def record(dt=0.0005, ps=10.0, B=20):                 # vibration periods from the spectrum of a short flexible run
    r = np.random.default_rng(1); x = np.repeat(place(180.0)[None], B, 0); v = r.normal(size=(B, 4, 3)) * np.sqrt(kT / M)[None, :, None]
    U, f, _ = energy_forces(x); bl, an = [], []
    for s in range(int(ps / dt)):
        v += 0.5 * dt * f / M[None, :, None]; x += dt * v; U, f, _ = energy_forces(x); v += 0.5 * dt * f / M[None, :, None]
        u_, w_ = x[:, 0] - x[:, 1], x[:, 2] - x[:, 1]; bl.append(np.linalg.norm(x[:, 2] - x[:, 1], axis=-1))
        an.append(np.arccos(dot(u_, w_) / np.linalg.norm(u_, axis=-1) / np.linalg.norm(w_, axis=-1)))
    out = []
    for s, lo, hi in ((np.array(bl), 15, 100), (np.array(an), 3, 15)):      # look for the bond peak above 15/ps, the angle peak below
        s = s - s.mean(0); P = (np.abs(np.fft.rfft(s, axis=0))**2).mean(1); fr = np.fft.rfftfreq(len(s), dt); k = (fr > lo) & (fr < hi)
        out.append(1000 / fr[k][np.argmax(P[k])])
    return out
pb, pa = record()
print(f"   vibration periods measured from a 10-ps flexible run: bond {pb:.0f} fs, angle {pa:.0f} fs  (ratio {pa / pb:.1f})")
`,
  },

  exercises: [
    { q: 'The C–C bond in our butane is compressed to r = 0.150 nm. What are the bond energy and the force on each atom?', a: 'dr = −0.004 nm. U = ½ × 224 262 × 0.004² = 1.794 kJ/mol. F = −k dr = +897 kJ/mol/nm, pushing the two atoms apart.' },
    { q: 'What is the dihedral energy at φ = 90°, and how often is it visited compared with trans at 300 K?', a: 'cos 90° = 0, cos 180° = −1, cos 270° = 0: U = 2.9519 × 1 + (−0.5670) × 2 + 6.5794 × 1 = 8.397 kJ/mol. Relative probability e^(−8.397/2.494) = e^(−3.367) = 0.035, about 30 times rarer than trans.' },
    { q: 'A box of 1000 rigid water molecules (3 atoms, all distances fixed by SETTLE, total momentum zero). How many degrees of freedom? By how much would T be wrong if you used 3N with N = 3000 atoms?', a: 'Each molecule: 9 coordinates − 3 constraints = 6. N_f = 6000 − 3 = 5997. Using 9000 instead makes T too low by a factor 5997/9000 = 0.67, i.e. 33 % low (300 K would read 200 K).' },
  ],
};
