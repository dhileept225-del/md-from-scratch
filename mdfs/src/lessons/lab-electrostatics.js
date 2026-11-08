import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-electrostatics',
  lab: true,
  summary: 'Charges are where simple MD tricks break down. A cutoff that works perfectly for Lennard-Jones gives nonsense for Coulomb, and even the order in which you add up the ions changes the answer. With two ions in a box we compute the energy by hand: bare, with the minimum image, and with Ewald summation. Then we see what a net charge does, why PME is needed, and how induced dipoles add a second, self-consistent layer.',
  minutes: 40,
  sim: 'lab5',
  learn: [
    'See by numbers why Coulomb cannot be cut off, and why the order of summation matters',
    'Do an Ewald sum by hand for two ions: real space, a few k-vectors, the self term, and the role of α',
    'Compute induced dipoles by iteration, and find the distance where polarizable models blow up',
  ],
  simHint: 'Ewald: move α and the number of k-vectors and watch the three parts change while the total stays on the exact line, until α is too small (real space misses images) or too large (too few k-vectors). Sums: the sphere sums jump around, the cube sums settle on 1.7476. Dipoles: bring the two polarizable atoms closer and watch the iterations converge more slowly, then run away below 0.144 nm.',

  setup: `
<p>We stay in GROMACS units. One new constant: the Coulomb energy of two unit charges at distance r is f/r with</p>
<div class="calc">f = 1/(4πε₀) in these units = <b>138.935 kJ/mol · nm / e²</b>

two ions   Na⁺ (+1) at (0, 0, 0)
           Cl⁻ (−1) at (0.3, 0, 0) nm
box        cubic, periodic, L = 1.0 nm
bare pair  U = −138.935 / 0.3 = <b>−463.12 kJ/mol</b></div>
<p>For comparison, k_B T at 300 K is 2.494 kJ/mol: an ion pair at contact is about 200 k_B T. Electrostatics is the largest energy in an ionic system, so small relative errors in it are large errors in everything else.</p>`,

  steps: [
    {
      title: 'Why Coulomb is different from Lennard-Jones',
      why: 'In Lab 2 we cut LJ off at 2.5σ and lost about 3 % of the energy. Before doing the same with charges, compare how the two fall off with distance, and how many neighbours sit at each distance.',
      work: `<div class="calc">Na⁺–Cl⁻:      Coulomb        LJ (ε = 0.5, σ = 0.30)
r = 0.28 nm   −496.20         +1.5515
r = 0.50      −277.87         −0.0890
r = 1.00      −138.94         −0.0015
r = 2.00       −69.47         −0.0000   kJ/mol

ions in a shell of radius r (thickness 0.1 nm,
30 ions/nm³, a molten salt):
r = 0.5:    9 ions   if charges did not cancel:  2 619
r = 1.0:   38 ions                                5 238
r = 2.0:  151 ions                               10 475 kJ/mol
LJ of the same shells:  0.838,  0.055,  0.003</div>`,
      see: 'LJ falls as 1/r⁶ while the number of neighbours grows as r², so far shells contribute almost nothing: a cutoff is safe. Coulomb falls only as 1/r, so a far shell holds more energy than a near one, 10 000 kJ/mol at 2 nm, unless its positive and negative charges cancel almost perfectly. The whole electrostatic energy is a small difference between huge positive and negative pieces. Any method that treats far shells carelessly, such as a simple cutoff, can be wrong by more than the answer itself.',
    },
    {
      title: 'Truncation and the order of summation',
      why: 'To see the problem cleanly, take crystals whose exact energy is known. A 1D chain + − + − with spacing a has the energy of one ion equal to −(2f/a) × (1 − ½ + ⅓ − ¼ + …) = −(2f/a) ln 2. Add the neighbours one by one, as a cutoff would.',
      work: `<div class="calc">neighbours on each side    sum        error
  1                       1.0000     +0.3069
  2                       0.5000     −0.1931
  3                       0.8333     +0.1402
  4                       0.5833     −0.1098
  5                       0.7833     +0.0902
 10                       0.6456     −0.0475
 11                       0.7365     +0.0434
100                       0.6882     −0.0050
exact  ln 2 =             0.6931</div>
<p>In 3D the rock-salt (NaCl) crystal has Madelung constant 1.747565. Add up all ions around one Na⁺ inside a growing sphere, or inside a growing cube (ions on the cube faces counted ½, edges ¼, corners ⅛):</p>
<div class="tbl"><table>${head(['size (spacings)', '2', '3', '4', '5', '6', '10'])}
${row(['sphere', '0.069', '3.164', '−1.213', '1.890', '−0.937', '−1.708'])}
${row(['cube', '1.7518', '1.7470', '1.7477', '1.7475', '1.7476', '1.7476'])}</table></div>`,
      see: 'In 1D, cutting off after n neighbours leaves an error of about 1/(2n), and its sign flips each time a new ion is added: with 100 neighbours per side it is still 0.7 % off. In 3D it gets worse: a sphere never converges at all (0.07, 3.16, −1.21, 1.89, …) because each new spherical shell is not neutral. The cube, which keeps every layer neutral, settles on 1.7476. The sum depends on the order of addition, which mathematicians call <b>conditionally convergent</b>. The lesson: never cut Coulomb off bluntly. Either make every piece neutral, or use a method that sums the infinite periodic system exactly: Ewald.',
    },
    {
      title: 'Two ions in a periodic box: what the minimum image misses',
      why: 'In a periodic simulation, each ion interacts with the other ion and with all of their images, an infinite lattice of ion pairs. The minimum image (Lab 2) keeps only the nearest copy. For LJ that is enough; let us see what it misses for charges.',
      work: `<div class="tbl"><table>${head(['separation (nm)', 'bare / minimum image', 'exact periodic (Ewald)', 'missed'])}
${row(['0.2', '−694.67', '−707.01', '−12.34'])}
${row(['0.3', '−463.12', '<b>−492.89</b>', '<b>−29.78</b>'])}
${row(['0.4', '−347.34', '−405.60', '−58.26'])}
${row(['0.5', '−277.87', '−380.87', '−103.00'])}</table></div>`,
      see: 'At 0.3 nm the images add −29.8 kJ/mol (6 %, 12 k_B T). The error grows as the ions separate, to −103 kJ/mol at 0.5 nm, where the Cl⁻ sits halfway between the Na⁺ and its image and is pulled by both, but the minimum image sees only one. So the minimum image gives the wrong forces too: it would push the pair together. To get the periodic energy we need all images, and Step 2 showed we cannot just add them up one by one.',
    },
    {
      title: 'Ewald summation by hand',
      why: 'Ewald\'s trick: split 1/r into a short-range part erfc(αr)/r, which dies off fast and can be summed in real space like LJ, and a smooth long-range part erf(αr)/r, which is summed as waves (k-vectors) in reciprocal space, where it converges fast. A self term removes the interaction of each charge with its own smoothed cloud. α decides where the split happens.',
      work: `<div class="calc">α = 5 /nm

real space (only the pair at 0.3 matters):
f × (−1) × erfc(1.5) / 0.3 = −138.935 × 0.033895 / 0.3
                          = <b>−15.697</b>
next image at 0.7: erfc(3.5) = 7×10⁻⁷  (negligible)

self term: −f α/√π × (1² + 1²)
         = −138.935 × 5 / 1.77245 × 2 = <b>−783.857</b>

reciprocal: f (2π/V) Σ_k  e^(−k²/4α²) / k² × |S(k)|²
prefactor f 2π / V = <b>872.954</b>
S(k) = Σ q e^(i k·r) = 1 − e^(i kₓ 0.3)
|S|² = 2 − 2 cos(kₓ × 0.3)

n = (±1,0,0):  k² = (2π)² = 39.478
   e^(−39.478/100) = 0.6738,  |S|² = 2.6180
   2 vectors:  872.954 × 2 × 0.6738/39.478 × 2.618 = <b>78.016</b>
n = (±1,±1,0), (±1,0,±1):  k² = 78.957
   e^(−0.7896) = 0.4540,  8 vectors -> <b>105.139</b>
n = (0,±1,0), (0,0,±1):  kₓ = 0,  |S|² = 0  -> 0</div>
<div class="tbl"><table>${head(['k-vectors (|nᵢ| ≤ k_max)', 'reciprocal', 'total', 'error'])}
${row(['26 (k_max 1)', '230.385', '−569.170', '−76.28'])}
${row(['124 (k_max 2)', '303.127', '−496.427', '−3.54'])}
${row(['342 (k_max 3)', '306.489', '−493.065', '−0.17'])}
${row(['728 (k_max 4)', '306.657', '<b>−492.897</b>', '−0.005'])}
${row(['1330 (k_max 5)', '306.662', '−492.892', '0.000'])}</table></div>
<div class="tbl"><table>${head(['α (k_max 4, minimum image)', 'real', 'reciprocal', 'self', 'total'])}
${row(['2', '−183.461', '11.581', '−313.543', '−485.42'])}
${row(['3', '−94.055', '72.055', '−470.314', '−492.31'])}
${row(['4', '−41.535', '175.743', '−627.085', '<b>−492.88</b>'])}
${row(['5', '−15.697', '306.657', '−783.857', '<b>−492.90</b>'])}
${row(['6', '−5.052', '452.638', '−940.628', '−493.04'])}
${row(['8', '−0.319', '756.167', '−1254.171', '−498.32'])}</table></div>`,
      see: 'Each piece is ordinary arithmetic: one erfc, one square root, and for each k-vector an exponential and a cosine. The three parts change enormously with α (the self term from −314 to −1254), but their sum stays at −492.9 kJ/mol whenever both sums are converged, which is the check that Ewald is right. The limits show up at the edges: with α = 2 the real-space part needs images beyond the minimum image (error +7.5); with α = 8 the k-sum needs more than 728 vectors (error −5.4). Choosing α is balancing the two costs; codes pick it from the cutoff and a tolerance (GROMACS <code>ewald-rtol</code>, LAMMPS <code>kspace_style ... 1e-5</code>).',
    },
    {
      title: 'A charged box: the background and its energy',
      why: 'Ewald needs the box to be neutral; otherwise the infinite lattice of identical charges has infinite energy. Codes then silently add a uniform neutralising background. This is common when one computes the solvation free energy of a single ion, and the background has an energy that depends on the box size.',
      work: `<div class="calc">one ion, charge q, in a cubic box of side L:
E = −ξ f q² / (2L),   ξ = 2.837297 (Wigner constant)

L = 1 nm:  −2.837297 × 138.935 / 2 = <b>−197.10 kJ/mol</b>
L = 2 nm:                            <b>−98.55</b>
L = 3 nm:                            −65.70
L = 5 nm:                            −39.42
(Ewald with background, L = 2: −98.55, same)</div>`,
      see: 'An artificial energy of −99 kJ/mol in a 2 nm box (40 k_B T), and it falls only as 1/L: making the box 5 times bigger still leaves −39 kJ/mol. For a charged system, any energy that is compared between box sizes or between charge states, such as an ion solvation free energy or a redox potential, must be corrected analytically for this term (and its dielectric screening). For a neutral electrolyte it cancels, which is why simulations should be neutral whenever possible.',
    },
    {
      title: 'The cost: why PME',
      why: 'Ewald is exact but its reciprocal sum loops over every k-vector for every charge. How does that scale with the number of atoms N, compared with adding all pairs, or with particle-mesh Ewald (PME), which spreads the charges on a grid and uses the fast Fourier transform?',
      work: `<div class="calc">operations per step (order of magnitude)
N            all pairs N²/2   Ewald ~N^1.5   PME ~N log₂N
1 000        5×10⁵            3×10⁴          1×10⁴
10 000       5×10⁷            1×10⁶          1×10⁵
100 000      5×10⁹            3×10⁷          2×10⁶
1 000 000    5×10¹¹           1×10⁹          2×10⁷</div>`,
      see: 'For 100 000 atoms, PME needs about 20 times fewer operations than optimally tuned Ewald and 3000 times fewer than all pairs (which would not even include the images). That is why every modern code uses PME or a relative (PPPM in LAMMPS). Its price is an interpolation error from the grid, set by the grid spacing (~0.12 nm) and the interpolation order; it is the same α trade-off as in step 4, now with a grid. Limitation: the FFT needs all-to-all communication, which limits scaling on many nodes.',
    },
    {
      title: 'Induced dipoles: polarization by iteration',
      why: 'Real atoms are not fixed charges: an electric field shifts their electron cloud and creates a dipole μ = α E, which in turn changes the field on its neighbours. Polarizable force fields include this, and it matters most next to ions, where fields are strong. We compute it, then find where it fails.',
      work: `<div class="calc">ion +1 at x = 0, polarizable atom at r = 0.25 nm,
α' = α/(4πε₀) = 0.0015 nm³ (1.5 Å³)
field  E = f q / r² = 138.935 / 0.0625 = 2222.96 kJ/mol/(nm e)
dipole μ = α' q / r² = 0.0015 / 0.0625 = <b>0.0240 e nm</b> (1.15 D)
energy U = −½ f α' q² / r⁴ = <b>−26.68 kJ/mol</b>

two atoms on the axis at 0.25 (A) and 0.50 (B):
each also feels the other's dipole field 2μ/d³ (d = 0.25)
iteration   μ_A       μ_B   (e nm)
1           0.02400   0.00600
2           0.02515   0.01061
3           0.02604   0.01083
4           0.02608   0.01100
5           0.02611   0.01101
6           0.02611   0.01101   converged

polarization catastrophe: each round multiplies by 2α'/d³
= 1 at d = (2α')^(1/3) = <b>0.144 nm</b>
d = 0.30   factor 0.111   dipoles × 1.12
d = 0.20   factor 0.375   × 1.60
d = 0.15   factor 0.889   × 9.00
d = 0.14   factor 1.093   diverges</div>`,
      see: 'One induced dipole is worth −27 kJ/mol next to an ion, about 10 k_B T: not a small correction. With two atoms, the mutual induction raises μ_B from 0.006 to 0.011, almost double, and needs about 5 iterations at every time step (or extra degrees of freedom, as in Drude oscillators). Limitation, by numbers: if two polarizable sites come closer than 0.144 nm, each iteration amplifies the dipoles instead of correcting them and the energy goes to −∞. Real polarizable force fields therefore damp the interaction at short range (Thole screening). The cheap alternative many non-polarizable models for ionic liquids use is to scale all ionic charges by about 0.7–0.9, an average version of the same screening.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)'])}
${row(['plain cutoff', 'cheap', '1D chain still 0.7 % off at 100 neighbours; 3D spheres never converge'])}
${row(['minimum image', 'fine for LJ', 'misses −29.8 kJ/mol at 0.3 nm, −103 at 0.5 nm'])}
${row(['Ewald', 'exact: −492.89, independent of α', '~N^1.5; α must balance real and k sums'])}
${row(['PME / PPPM', '~N log N: 20× fewer operations than Ewald at 10⁵ atoms', 'grid error; FFT communication'])}
${row(['net charge', '—', '−ξ f q²/2L: −99 kJ/mol at L = 2 nm, falls only as 1/L'])}
${row(['induced dipoles', '−27 kJ/mol per dipole near an ion', 'iterations each step; diverges below 0.144 nm without damping'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the Coulomb/LJ comparison, the 1D and 3D Madelung sums, the two-ion energies, the Ewald pieces and their convergence, the charged-box energy, the cost table and the induced-dipole iterations. Run <code>pip install numpy</code>, then <code>python lab_electrostatics.py</code> (about a second).',
    src: `# Lab 5.5: electrostatics by hand. GROMACS units: nm, e, kJ/mol; Coulomb constant f = 138.935 kJ mol^-1 nm e^-2.
import numpy as np, itertools, time
from math import erfc, sqrt, pi, exp, log
f = 138.935
def lj(r, e, s): return 4 * e * ((s / r)**12 - (s / r)**6)

# ---------- Step 1: why Coulomb is different ----------
print("Step 1: Na+ - Cl- (LJ eps = 0.5 kJ/mol, sigma = 0.30 nm, for comparison)")
for r in (0.28, 0.5, 1.0, 2.0):
    print(f"   r = {r:4.2f} nm: Coulomb = {-f / r:8.2f}   LJ = {lj(r, 0.5, 0.30):9.4f} kJ/mol")
print("   ions in a spherical shell of radius r, thickness 0.1 nm, ion density 30 /nm^3 (molten-salt-like):")
for r in (0.5, 1.0, 2.0):
    n = 4 * pi * r**2 * 0.1 * 30; print(f"   r = {r}: {n:6.0f} ions;  if their charges did NOT cancel, |shell energy| ~ {n:.0f} x {f:.1f}/{r} = {n * f / r:9.0f} kJ/mol;  LJ shell ~ {n * abs(lj(r, 0.5, 0.30)):.3f}")

# ---------- Step 2: the 1D ionic chain and truncation ----------
print("\\nStep 2: 1D chain + - + - ..., spacing a. Energy of one ion = -(2 f / a) x (1 - 1/2 + 1/3 - 1/4 + ...)")
s = 0.0
for n in range(1, 101):
    s += (-1)**(n + 1) / n
    if n in (1, 2, 3, 4, 5, 10, 11, 100): print(f"   cutoff after {n:3d} neighbours on each side: sum = {s:.4f}   error = {s - log(2):+.4f}")
print(f"   exact: ln 2 = {log(2):.4f}")
print("   3D rock salt (NaCl), Madelung constant = 1.747565. Sum over all ions inside a sphere / a cube around one Na+:")
def madelung(shape, R):
    n = int(R) + 1; g = np.arange(-n, n + 1); X, Y, Z = np.meshgrid(g, g, g, indexing='ij'); r = np.sqrt(X**2 + Y**2 + Z**2)
    q = (-1.0)**(X + Y + Z); inside = (r <= R + 1e-9) if shape == 'sphere' else ((abs(X) <= R) & (abs(Y) <= R) & (abs(Z) <= R))
    inside &= r > 0; w = np.ones_like(r)
    if shape == 'cube':                                       # Evjen: ions on faces count 1/2, edges 1/4, corners 1/8
        for A in (X, Y, Z): w *= np.where(abs(A) == R, 0.5, 1.0)
    return -(q * w / np.where(r > 0, r, 1))[inside].sum()
for R in (2, 3, 4, 5, 6, 8, 10):
    print(f"   radius / half-side {R:2d} spacings: sphere {madelung('sphere', R + 0.5):8.4f}   cube (Evjen) {madelung('cube', R):8.4f}")

# ---------- Ewald for two ions in a periodic box ----------
L = 1.0; V = L**3; pos = np.array([[0.0, 0, 0], [0.3, 0, 0]]); q = np.array([1.0, -1.0])
def ewald(alpha, nreal=2, kmax=6, parts=False):
    Er = 0.0
    for n in itertools.product(range(-nreal, nreal + 1), repeat=3):
        for i in range(2):
            for j in range(2):
                if n == (0, 0, 0) and i == j: continue
                r = np.linalg.norm(pos[j] - pos[i] + L * np.array(n)); Er += 0.5 * f * q[i] * q[j] * erfc(alpha * r) / r
    Ek = 0.0
    for n in itertools.product(range(-kmax, kmax + 1), repeat=3):
        if n == (0, 0, 0): continue
        k = 2 * pi * np.array(n) / L; k2 = k @ k; S = (q * np.exp(1j * pos @ k)).sum()
        Ek += f * 2 * pi / V * exp(-k2 / (4 * alpha**2)) / k2 * abs(S)**2
    Es = -f * alpha / sqrt(pi) * (q**2).sum()
    return (Er, Ek, Es, Er + Ek + Es) if parts else Er + Ek + Es
exact = ewald(5.0, nreal=3, kmax=10)
print(f"\\nStep 3: Na+ at (0,0,0), Cl- at (0.3,0,0), periodic box L = {L} nm")
print(f"   bare pair (= minimum image, cutoff 0.5): -f/0.3 = {-f / 0.3:.2f} kJ/mol")
print(f"   exact periodic energy (Ewald, converged): {exact:.2f} kJ/mol  -> the images add {exact + f / 0.3:+.2f} kJ/mol")
# effect of moving the pair: energy vs separation, bare vs periodic
for d in (0.2, 0.3, 0.4, 0.5):
    pos[1, 0] = d; print(f"   separation {d}: bare {-f / d:8.2f}   periodic {ewald(5.0, nreal=3, kmax=10):8.2f}   difference {ewald(5.0, nreal=3, kmax=10) + f / d:+7.2f}")
pos[1, 0] = 0.3

# ---------- Step 4: Ewald by hand ----------
a = 5.0
print(f"\\nStep 4: Ewald with alpha = {a} /nm")
print(f"   real space, the pair at 0.3: f x (-1) x erfc({a * 0.3}) / 0.3 = -{f:.3f} x {erfc(a * 0.3):.6f} / 0.3 = {-f * erfc(a * 0.3) / 0.3:.3f}")
print(f"   next image of Cl- at 0.7: erfc({a * 0.7}) = {erfc(a * 0.7):.2e};  Na+ image at 1.0: erfc({a * 1.0}) = {erfc(a):.2e}  (negligible)")
print(f"   self term: -f alpha/sqrt(pi) x (1 + 1) = -{f:.3f} x {a}/{sqrt(pi):.5f} x 2 = {-f * a / sqrt(pi) * 2:.3f}")
pref = f * 2 * pi / V
print(f"   reciprocal prefactor f 2 pi / V = {pref:.3f}")
for name, vecs in (("n = (+-1,0,0)", [(1, 0, 0), (-1, 0, 0)]), ("n = (+-1,+-1,0) etc. with nx = +-1", [(sx, sy, sz) for sx in (1, -1) for sy, sz in ((1, 0), (-1, 0), (0, 1), (0, -1))]),
                   ("n = (0,+-1,0), (0,0,+-1)", [(0, 1, 0), (0, -1, 0), (0, 0, 1), (0, 0, -1)])):
    tot = 0
    for n in vecs:
        k = 2 * pi * np.array(n) / L; k2 = k @ k; S = (q * np.exp(1j * pos @ k)).sum(); tot += pref * exp(-k2 / (4 * a**2)) / k2 * abs(S)**2
    n = vecs[0]; k = 2 * pi * np.array(n) / L; k2 = k @ k; S = (q * np.exp(1j * pos @ k)).sum()
    print(f"   {name}: k^2 = {k2:.3f}, exp(-k^2/4a^2) = {exp(-k2 / (4 * a**2)):.4f}, |S|^2 = 2 - 2cos(kx 0.3) = {abs(S)**2:.4f}, {len(vecs)} vectors -> {tot:.3f}")
print("   reciprocal sum vs number of k-vectors used (|n_i| <= kmax):")
for km in (1, 2, 3, 4, 5, 6, 8):
    Er, Ek, Es, E = ewald(a, kmax=km, parts=True); print(f"   kmax = {km}: {(2 * km + 1)**3 - 1:5d} vectors, reciprocal = {Ek:9.3f}, total = {E:9.3f}  (error {E - exact:+.3f})")
print("   alpha decides the split (total must not depend on it). Real space with the minimum image only (cutoff 0.5), kmax = 4:")
for al in (2.0, 3.0, 4.0, 5.0, 6.0, 8.0, 10.0):
    Er, Ek, Es, E = ewald(al, nreal=0, kmax=4, parts=True); print(f"   alpha = {al:4.1f}: real {Er:9.3f}  recip {Ek:9.3f}  self {Es:9.3f}  total {E:9.3f}  (error {E - exact:+.3f})")

# ---------- Step 5: a charged box ----------
xi = 2.837297
print(f"\\nStep 5: one ion alone in a periodic box (neutralised by a uniform background): E = -xi f q^2 / (2L), xi = {xi}")
for LL in (1.0, 2.0, 3.0, 5.0):
    print(f"   L = {LL} nm: E = -{xi} x {f} / (2 x {LL}) = {-xi * f / (2 * LL):8.2f} kJ/mol")
# check with Ewald (single ion, background term)
def ewald1(LL, alpha=3.0, kmax=8, nreal=3):
    Er = 0.0
    for n in itertools.product(range(-nreal, nreal + 1), repeat=3):
        if n == (0, 0, 0): continue
        r = LL * np.linalg.norm(n); Er += 0.5 * f * erfc(alpha * r) / r
    Ek = 0.0
    for n in itertools.product(range(-kmax, kmax + 1), repeat=3):
        if n == (0, 0, 0): continue
        k2 = (2 * pi / LL)**2 * (n[0]**2 + n[1]**2 + n[2]**2); Ek += f * 2 * pi / LL**3 * exp(-k2 / (4 * alpha**2)) / k2
    return Er + Ek - f * alpha / sqrt(pi) - f * pi / (2 * LL**3 * alpha**2)
print(f"   check with Ewald for L = 2: {ewald1(2.0):.2f} kJ/mol")

# ---------- Step 6: cost ----------
print("\\nStep 6: operations per step (order of magnitude)")
for N in (1000, 10_000, 100_000, 1_000_000):
    print(f"   N = {N:>9,d}: all pairs N^2/2 = {N * N / 2:.1e}   Ewald ~ N^1.5 = {N**1.5:.1e}   PME ~ N log2 N = {N * np.log2(N):.1e}")

# ---------- Step 7: induced dipoles ----------
ap = 0.0015                                            # polarizability volume alpha/(4 pi eps0), nm^3 (1.5 A^3)
print(f"\\nStep 7: an ion (+1) at x = 0 and polarizable atoms (alpha' = {ap} nm^3 = 1.5 A^3) on the x axis")
r = 0.25; E0 = f * 1 / r**2                            # field in kJ/mol/(nm e); dipole mu = alpha' E0 / f in e nm
mu = ap * E0 / f
print(f"   one atom at r = {r}: field E = f q/r^2 = {f} / {r**2} = {E0:.2f} kJ/mol/(nm e);  mu = alpha' q / r^2 = {ap} / {r**2} = {mu:.4f} e nm ({mu / 0.020819:.2f} D)")
print(f"   energy U = -1/2 f alpha' q^2 / r^4 = -0.5 x {f} x {ap} / {r**4:.6f} = {-0.5 * f * ap / r**4:.2f} kJ/mol   (ion-pair at 0.25: {-f / r:.1f})")
rA, rB = 0.25, 0.50
print(f"   two atoms at x = {rA} and {rB}: each feels the ion AND the other's dipole (field on the axis 2 mu / d^3)")
EA0, EB0 = 1 / rA**2, 1 / rB**2; d = rB - rA; muA = muB = 0.0
for it in range(6):
    nA = ap * (EA0 + 2 * muB / d**3); nB = ap * (EB0 + 2 * muA / d**3)
    print(f"   iteration {it + 1}: mu_A = {nA:.5f}  mu_B = {nB:.5f} e nm");  muA, muB = nA, nB
print(f"   without mutual induction: mu_A = {ap * EA0:.5f}, mu_B = {ap * EB0:.5f}")
dc = (2 * ap)**(1 / 3)
print(f"   polarization catastrophe: two atoms d apart reinforce each other by a factor 2 alpha'/d^3 per round;  = 1 at d = (2 alpha')^(1/3) = {dc:.4f} nm")
for dd in (0.30, 0.20, 0.16, 0.15, 0.14):
    print(f"   d = {dd:5.3f}: factor 2 alpha'/d^3 = {2 * ap / dd**3:.3f}  -> {'converges' if 2 * ap / dd**3 < 1 else 'diverges'}" + (f", two atoms in the same field: dipoles {1 / (1 - 2 * ap / dd**3):.2f} x the single-atom value" if 2 * ap / dd**3 < 1 else ""))
`,
  },

  exercises: [
    { q: 'Our Na⁺/Cl⁻ pair in a box of L = 2 nm instead of 1 nm: does the minimum-image error at 0.3 nm get bigger or smaller, and why?', a: 'Smaller. The images are now 2 nm away instead of 1 nm, and the pair is a dipole whose field falls as 1/r³, so the image–dipole interaction falls roughly as 1/L³: about 8 times smaller. Ewald gives −3.4 kJ/mol instead of −29.8 (set L = 2.0 in the script to check).' },
    { q: 'Ewald with α = 5 /nm: what is the real-space term of a Na⁺/Cl⁻ pair at 0.25 nm, and how big is the next image at 0.75 nm? (erfc(1.25) = 0.077100, erfc(3.75) = 1.1×10⁻⁷)', a: 'f × (−1) × 0.077100 / 0.25 = −42.85 kJ/mol. The image: −138.935 × 1.1×10⁻⁷ / 0.75 ≈ −2×10⁻⁵ kJ/mol, negligible, so the minimum image is enough in real space once α is large enough.' },
    { q: 'An atom with α\' = 1.0 Å³ = 0.001 nm³ sits 0.30 nm from a −1 ion. What dipole and energy are induced? Where would two such atoms run into the catastrophe?', a: 'μ = 0.001 / 0.09 = 0.0111 e nm (0.53 D), pointing towards the ion. U = −½ × 138.935 × 0.001 / 0.0081 = −8.58 kJ/mol. Catastrophe at d = (2 × 0.001)^(1/3) = 0.126 nm.' },
  ],
};
