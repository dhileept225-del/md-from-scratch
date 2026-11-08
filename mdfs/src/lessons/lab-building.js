import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-building',
  lab: true,
  summary: 'Before the first time step, a simulation needs units, starting positions, starting velocities and a periodic box, plus rules for which pairs to compute. We build a 9-atom system by hand and let the numbers show why each of these choices is made, and what goes wrong without it.',
  minutes: 30,
  sim: 'lab2',
  learn: [
    'Convert argon to reduced units, and place 9 atoms on a lattice instead of at random',
    'Prepare velocities: remove the drift, count degrees of freedom, scale to the target temperature',
    'Use periodic boundaries, the minimum image, a cutoff and a neighbour list, and see what each one fixes',
  ],
  simHint: 'Switch between the lattice and the random start and read the energy. Step through the velocity preparation (raw, drift removed, scaled) and watch T and the drift. Tap an atom: solid lines are its neighbours inside the cutoff, also through the walls. Turn the minimum image off and see neighbours disappear. Press Run to see atoms wrap and the neighbour list being rebuilt.',

  setup: `
<p>To keep every number checkable, we work in <b>two dimensions</b> with <b>9 argon atoms</b> in a square box of side <b>L = 3.6</b> (reduced units, Lesson 2.1). The steps are exactly those of a real 3D simulation; only the counting is smaller.</p>
<div class="calc">atoms        N = 9
box          L = 3.6 × 3.6, periodic in x and y
target       T = 0.72  (86 K for argon)
time step    dt = 0.005
cutoff       r_c = 1.8  (the largest this box allows, see step 6)</div>`,

  steps: [
    {
      title: 'Units: why we do not use joules and metres',
      why: 'Argon quantities in SI units are tiny numbers like 10⁻²¹ J. Writing them out at every step invites mistakes, and the same simulation would have to be repeated for every gas. Reduced units fix both problems.',
      work: `<div class="calc">argon:   σ = 3.405×10⁻¹⁰ m    ε/k_B = 119.8 K    mass = 39.948 g/mol

energy unit   ε = 119.8 × 1.3806×10⁻²³ J/K   = 1.6540×10⁻²¹ J
mass unit     m = 39.948 × 1.6605×10⁻²⁷ kg   = 6.6335×10⁻²⁶ kg
time unit     τ = σ × √(m/ε)
                = 3.405×10⁻¹⁰ × √(4.0106×10⁻⁵)
                = 2.156×10⁻¹² s  = <b>2.156 ps</b>

T*  = 0.72   ->  0.72 × 119.8 = <b>86.3 K</b>
dt* = 0.005  ->  0.005 × 2.156 ps = <b>10.8 fs</b></div>`,
      see: 'In reduced units every number is close to 1, and one simulation describes all Lennard-Jones gases: to get krypton, convert back with krypton\'s σ, ε and mass. Limitation: you must convert back at the end, and a mixture needs one species chosen as the reference.',
    },
    {
      title: 'Starting positions: why a lattice and not random',
      why: 'Atoms must start somewhere. Dropping them at random looks natural, but two atoms almost always land on top of each other, and the LJ energy explodes at short distance.',
      work: `<p><b>Random</b> (9 positions drawn by the script, rounded): atom 1 at (3.40, 1.84) and atom 6 at (3.25, 1.72):</p>
<div class="calc">dx = −0.15,  dy = −0.12
r  = √(0.0225 + 0.0144) = √0.0369 = <b>0.192</b>
U  = 4 × (0.192⁻¹² − 0.192⁻⁶) = <b>1.6 × 10⁹</b>
F  = 9.9 × 10¹⁰
velocity change in ONE step = F × dt = 4.9 × 10⁸
(a typical thermal speed is about 1)</div>
<p><b>Lattice</b>: a 3 × 3 square grid with spacing a = L/3 = 1.2, positions x = (i + ½)·a, y = (j + ½)·a:</p>
<div class="calc">(0.6,0.6) (1.8,0.6) (3.0,0.6) (0.6,1.8) (1.8,1.8) (3.0,1.8) (0.6,3.0) (1.8,3.0) (3.0,3.0)

every atom has 4 neighbours at 1.2 and 4 at 1.2×√2 = 1.697
U(1.2)   = −0.8910
U(1.697): 1.697² = 2.88,  2.88³ = 23.89,  r⁻⁶ = 0.04186,  r⁻¹² = 0.001752
          U = 4 × (0.001752 − 0.04186) = −0.1604
energy per atom = ½ × (4 × −0.8910 + 4 × −0.1604) = <b>−2.1028</b>
(the ½: each pair energy is shared by two atoms)
total = 9 × −2.1028 = −18.925</div>`,
      see: 'The random start puts 1.6 billion ε into one pair; the first step would fire the atoms apart at absurd speeds. The lattice has no overlaps and every atom has the same environment. Limitations: the start is an ordered crystal, so a liquid must first be melted and equilibrated; N must fit the lattice (n² in 2D, 4n³ for fcc in 3D); molecules need packing tools such as Packmol.',
    },
    {
      title: 'Starting velocities: drift, degrees of freedom, scaling',
      why: 'Velocities set the temperature. We draw them from a bell-shaped (Gaussian) distribution, as in a real gas (Maxwell–Boltzmann), but a few random numbers never average exactly to zero, and temperature must be counted correctly.',
      work: `<div class="calc">9 random (vx, vy), two decimals:
(−0.33,−0.88) (−0.66,−0.67) (+0.38,−0.11) (+1.48,−1.83) (−0.00,−0.89)
(+0.78,−2.12) (−0.34,+0.21) (−1.48,+0.99) (+0.18,+1.01)

1. drift (average velocity):
   vx: (−0.33 −0.66 +0.38 +1.48 −0.00 +0.78 −0.34 −1.48 +0.18)/9 = 0.01/9  = +0.0011
   vy: (−0.88 −0.67 −0.11 −1.83 −0.89 −2.12 +0.21 +0.99 +1.01)/9 = −4.29/9 = <b>−0.4767</b>
   kinetic energy of the drift = ½ × 9 × 0.4767² = 1.02  (of 8.87 in total)
   -> subtract the drift from every atom

2. temperature, 2D:  T = Σ m v² / N_f
   N_f = 2N − 2 = 16   (2 directions per atom, minus 2 for the fixed total momentum)
   Σ v² after removing the drift = 15.696
   T = 15.696 / 16 = 0.9810

3. scale every velocity by λ = √(0.72 / 0.9810) = 0.8567   ->  T = <b>0.7200</b>

   wrong count N_f = 2N = 18:  T would read 0.64, 11 % too low</div>`,
      see: 'Without step 1 the whole box would slide downwards, and 11 % of the "thermal" energy would really be motion of the box. With only 9 atoms, using 2N instead of 2N − 2 gives an 11 % error in T; for 10 000 atoms it is 0.01 %. Limitations: 9 random numbers make a very rough bell curve, but any shape relaxes to Maxwell–Boltzmann within a few collisions; what must be right is the total momentum (zero) and T.',
    },
    {
      title: 'Periodic boundaries: why the box has no walls',
      why: 'We want to simulate a piece of a bulk liquid, not a tiny droplet. In a small group of atoms almost all of them sit on the surface, where their environment is completely different.',
      work: `<div class="tbl"><table>${head(['system (no periodic boundaries)', 'atoms on the edge', 'fraction'])}
${row(['3 × 3 (our 2D box)', '8 of 9', '89 %'])}
${row(['10 × 10', '36 of 100', '36 %'])}
${row(['100 × 100', '396 of 10 000', '4 %'])}
${row(['10 × 10 × 10 (3D)', '488 of 1000', '49 %'])}</table></div>
<p>With periodic boundaries the box is surrounded by copies of itself. An atom that leaves on one side comes back on the other:</p>
<div class="calc">x_wrapped = x − L × floor(x / L)
x = 3.75  ->  3.75 − 3.6 × 1 = <b>0.15</b>
x = −0.30 ->  −0.30 − 3.6 × (−1) = <b>3.30</b></div>`,
      see: 'Even 1000 atoms in 3D would be half surface without periodic boundaries. With them, every atom is surrounded by atoms and there is no surface at all. Limitations: the system is artificially periodic, so the box must be larger than any structure or correlation you study, and some properties need size corrections (for example diffusion, Lesson 6.2).',
    },
    {
      title: 'The minimum image: neighbours through the wall',
      why: 'With periodic copies, every atom has infinitely many images. For short-range forces we only want the nearest copy of each other atom. Measuring distances inside the box alone misses neighbours that sit just across a wall.',
      work: `<div class="calc">atom 1 at (0.6, 0.6), atom 3 at (3.0, 0.6)
raw dx = 3.0 − 0.6 = 2.4          (looks far away: beyond the cutoff 1.8)
minimum image: dx − L × round(dx / L) = 2.4 − 3.6 × 1 = <b>−1.2</b>
-> atom 3 is a NEAREST neighbour of atom 1, just across the left wall</div>
<div class="tbl"><table>${head(['', 'neighbours of atom 1', 'energy of atom 1', 'total energy'])}
${row(['with minimum image', '4 at 1.2, 4 at 1.697', '−2.1028', '−18.925'])}
${row(['without (raw distances)', '2 at 1.2, 1 at 1.697', '−0.9712', '−11.975'])}</table></div>`,
      see: 'Forgetting the minimum image throws away 54 % of the energy of a corner atom and 37 % of the total: atoms near the walls would feel almost no neighbours and drift to the middle. The rule needs one line of code and is applied to every pair, every step.',
    },
    {
      title: 'Cutoffs: why we ignore far pairs, and the price',
      why: 'The number of pairs grows as N², but the LJ energy of a far pair is tiny. A cutoff r_c skips pairs farther than r_c. It saves time, but it changes the energy, and the minimum image only works if r_c ≤ L/2.',
      work: `<div class="tbl"><table>${head(['cutoff r_c', 'energy per atom (infinite lattice)', 'error'])}
${row(['1.8 (our box)', '−2.1028', '4.8 %'])}
${row(['2.5', '−2.1445', '2.9 %'])}
${row(['4.0', '−2.2015', '0.3 %'])}
${row(['∞', '−2.2089', '0'])}</table></div>
<div class="calc">jump at the cutoff:  U(1.8) = 4 × (1.8⁻¹² − 1.8⁻⁶) = 4 × (0.00087 − 0.02940) = <b>−0.1141</b>
-> when a pair crosses r = 1.8 the energy jumps by 0.114: energy is not conserved
   fix: shifted potential U(r) − U(r_c), which is exactly 0 at the cutoff

box size rule: r_c ≤ L/2, otherwise an atom could interact with two images of the same atom
   our box: L/2 = 1.8, so r_c = 2.5 is impossible with 9 atoms
   r_c = 2.5 needs L ≥ 5.0:  2D at this density N ≥ (5.0/1.2)² ≈ 17
                              3D liquid (ρ = 0.84):  N ≥ 0.84 × 5³ ≈ 105</div>`,
      see: 'A cutoff of 2.5 already loses about 3 % of the energy in this dense lattice; codes add a "tail correction" for the missing part. The cutoff also creates a small jump in energy unless the potential is shifted. And it is why real simulations need at least ~100 atoms: the box must be twice the cutoff.',
    },
    {
      title: 'Neighbour lists: computing only the pairs that matter',
      why: 'Even with a cutoff, a code that checks all N(N−1)/2 pairs every step to find the close ones wastes most of its time. A neighbour list stores the pairs within r_c + skin, and is rebuilt only when atoms have moved far enough (Lesson 2.6).',
      work: `<div class="tbl"><table>${head(['N', 'all pairs N(N−1)/2'])}
${row(['9', '36'])}
${row(['1 000', '499 500'])}
${row(['100 000', '4 999 950 000'])}</table></div>
<div class="calc">3D liquid, ρ = 0.84, r_c = 2.5, skin = 0.3, list radius 2.8:
   neighbours per atom = (4/3) π × 2.8³ × 0.84 = 77.2
   N = 1000:  1000 × 77.2 / 2 = 38 620 pairs   (13× fewer than 499 500)

when to rebuild: when some atom has moved more than skin/2 = 0.15
   speed ≈ 1.5, dt = 0.005  ->  0.15 / (1.5 × 0.005) ≈ every 20 steps

skin 0.6 instead: 104.8 neighbours (+36 % pairs), rebuild every ~40 steps</div>`,
      see: 'For our 9 atoms the list holds all 36 pairs, so it saves nothing: lists pay off only for large N, where they turn an N² cost into a cost proportional to N. The skin is a trade-off: larger skin, fewer rebuilds but more pairs to check each step.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['choice', 'what it fixes', 'price / limitation'])}
${row(['reduced units', 'awkward 10⁻²¹ numbers; one run for all LJ gases', 'convert back at the end'])}
${row(['lattice start', 'overlaps (U = 1.6×10⁹ for one random pair)', 'must melt and equilibrate'])}
${row(['remove drift, N_f = 2N−2', 'sliding box; 11 % error in T for N = 9', 'none'])}
${row(['periodic boundaries', '89 % surface atoms in a 3×3 group', 'artificial periodicity, size effects'])}
${row(['minimum image', '54 % of a corner atom\'s energy missing', 'requires r_c ≤ L/2'])}
${row(['cutoff (shifted)', 'N² far pairs', '~3 % energy at r_c = 2.5 (tail correction)'])}
${row(['neighbour list', '499 500 pair checks for 1000 atoms', 'rebuilds; memory'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the unit conversion, the lattice and random energies, the velocity preparation, the surface fractions, the minimum-image test, the cutoff table and the neighbour-list counts. Run <code>pip install numpy</code>, then <code>python lab_building.py</code> (under a second).',
    src: `# Lab 2: building a small system by hand. 9 argon atoms in a 2D periodic box, reduced units (sigma = epsilon = m = 1).
import numpy as np

def U(r): return 4 * (r**-12 - r**-6)

# ---------- Step 1: reduced units for argon ----------
kB, amu = 1.380649e-23, 1.66053907e-27
sigma, eps_K, mass = 3.405e-10, 119.8, 39.948                     # m, K (epsilon/kB), g/mol
eps = eps_K * kB; m = mass * amu; tau = sigma * np.sqrt(m / eps)
print("Step 1: argon in SI units")
print(f"   epsilon = 119.8 K x kB = {eps:.4e} J;  m = 39.948 x 1.6605e-27 = {m:.4e} kg")
print(f"   tau = sigma sqrt(m/epsilon) = 3.405e-10 x sqrt({m / eps:.4e}) = {tau:.4e} s = {tau * 1e12:.3f} ps")
print(f"   T* = 1.0 -> {eps_K:.1f} K;  T* = 0.72 -> {0.72 * eps_K:.1f} K;  dt* = 0.005 -> {0.005 * tau * 1e15:.1f} fs")

# ---------- Step 2: starting positions: lattice vs random ----------
n, a = 3, 1.2; L = 3.6; N = n * n                                 # L = n * a
lat = np.array([[(i + 0.5) * a, (j + 0.5) * a] for j in range(n) for i in range(n)])
def pair_list(x, mic=True, rc=1.8):
    out = []
    for i in range(len(x)):
        for j in range(i + 1, len(x)):
            d = x[j] - x[i]
            if mic: d = d - L * np.round(d / L)
            r = np.hypot(*d)
            if r < rc: out.append((i, j, r))
    return out
def energy(x, mic=True, rc=1.8): return sum(U(r) for _, _, r in pair_list(x, mic, rc))
print(f"\\nStep 2: 3 x 3 square lattice, spacing a = {a}, box L = {L}")
print("   positions:", ", ".join(f"({p[0]:.1f},{p[1]:.1f})" for p in lat))
print(f"   U(1.2) = {U(1.2):.4f}, U(1.697) = {U(np.sqrt(2) * 1.2):.4f};  energy per atom = (4 U(1.2) + 4 U(1.697))/2 = {(4 * U(1.2) + 4 * U(np.sqrt(2) * 1.2)) / 2:.4f}")
print(f"   total energy (9 atoms, minimum image, cutoff 1.8) = {energy(lat):.4f} = 9 x {energy(lat) / 9:.4f}")
rng = np.random.default_rng(4); rnd = np.round(rng.random((N, 2)) * L, 2)
pl = pair_list(rnd, rc=10); i, j, rmin = min(pl, key=lambda t: t[2])
print("   random positions:", ", ".join(f"({p[0]:.2f},{p[1]:.2f})" for p in rnd))
print(f"   closest pair: atoms {i + 1} and {j + 1}, r = {rmin:.3f}: U = {U(rmin):.4g}, force = {24 * (2 * rmin**-13 - rmin**-7):.4g}")
print(f"   one step (dt = 0.005) gives a velocity change of F*dt = {24 * (2 * rmin**-13 - rmin**-7) * 0.005:.4g}  (thermal speed ~ 1)")

# ---------- Step 3: velocities from a Gaussian: remove drift, then scale to T ----------
g = np.round(rng.normal(size=(N, 2)), 2)                            # the 'random' numbers, rounded for hand work
print("\\nStep 3: Gaussian random velocities (2 decimals):")
print("   " + ", ".join(f"({v[0]:+.2f},{v[1]:+.2f})" for v in g))
mean = g.mean(0); v = g - mean
print(f"   average velocity (drift) = ({mean[0]:+.4f}, {mean[1]:+.4f});  after subtracting: average = ({v.mean(0)[0]:+.1e}, {v.mean(0)[1]:+.1e})")
KE_drift = 0.5 * N * (mean @ mean)
print(f"   kinetic energy of the drift alone = 1/2 N |v_mean|^2 = {KE_drift:.4f}  (the whole box sliding: not heat)")
s2 = (v**2).sum(); Nf = 2 * N - 2
print(f"   sum of v^2 = {s2:.4f};  degrees of freedom N_f = 2N - 2 = {Nf}")
T_now = s2 / Nf; lam = np.sqrt(0.72 / T_now)
print(f"   T = sum(m v^2) / N_f = {T_now:.4f};  scale factor sqrt(0.72 / T) = {lam:.4f};  new T = {((lam * v)**2).sum() / Nf:.4f}")
print(f"   if you wrongly used N_f = 2N = {2 * N}: T would read {((lam * v)**2).sum() / (2 * N):.4f} instead of 0.72 ({(1 - Nf / (2 * N)) * 100:.1f} % low)")

# ---------- Step 4: periodic boundaries ----------
print("\\nStep 4: wrapping and surface atoms")
x = 3.75; print(f"   an atom moves to x = {x}: wrapped x = x - L*floor(x/L) = {x - L * np.floor(x / L):.2f};  x = -0.30 -> {-0.30 - L * np.floor(-0.30 / L):.2f}")
for k in [3, 10, 100]:
    print(f"   {k} x {k} atoms without periodic boundaries: {k * k - (k - 2) ** 2} of {k * k} on the edge = {100 * (k * k - (k - 2) ** 2) / k**2:.0f} %")
print(f"   10 x 10 x 10 in 3D: {1000 - 512} of 1000 on the surface = {100 * 488 / 1000:.0f} %")

# ---------- Step 5: minimum image ----------
A, C = lat[0], lat[2]
d = C - A; dm = d - L * np.round(d / L)
print(f"\\nStep 5: atom 1 at ({A[0]}, {A[1]}) and atom 3 at ({C[0]}, {C[1]}): raw dx = {d[0]:.1f};  dx - L*round(dx/L) = {d[0]:.1f} - {L:.1f}*{np.round(d[0] / L):.0f} = {dm[0]:.1f}")
e_corner_no = sum(U(r) for i_, j_, r in pair_list(lat, mic=False) if 0 in (i_, j_)) / 2
e_corner = sum(U(r) for i_, j_, r in pair_list(lat, mic=True) if 0 in (i_, j_)) / 2
print(f"   energy of atom 1 (half of its pair energies): with minimum image {e_corner:.4f}, without {e_corner_no:.4f}  ({100 * (1 - e_corner_no / e_corner):.0f} % missing)")
print(f"   total energy: with {energy(lat):.4f}, without {energy(lat, mic=False):.4f}")

# ---------- Step 6: cutoffs ----------
def lattice_sum(rc, R=30):                                           # energy per atom of the infinite square lattice
    e = 0
    for i in range(-R, R + 1):
        for j in range(-R, R + 1):
            r = a * np.hypot(i, j)
            if 0 < r < rc: e += U(r)
    return e / 2
print("\\nStep 6: cutoff")
for rc in [1.8, 2.5, 4.0, 30.0]:
    print(f"   energy per atom of the infinite lattice with cutoff {rc:4.1f}: {lattice_sum(rc):.4f}")
print(f"   U at the cutoff 1.8 = {U(1.8):.4f}: a pair crossing r = 1.8 makes the energy jump by this much (use a shifted potential: U(r) - U(rc))")
print(f"   box must be at least 2 x cutoff: for rc = 2.5 need L >= 5.0 -> 2D at this density N >= (5.0/1.2)^2 = {(5 / 1.2)**2:.0f};  3D liquid (rho 0.84): N >= 0.84 x 5^3 = {0.84 * 125:.0f}")

# ---------- Step 7: neighbour lists ----------
print("\\nStep 7: neighbour lists")
for Nn in [9, 1000, 100000]:
    print(f"   N = {Nn:6d}: all pairs N(N-1)/2 = {Nn * (Nn - 1) // 2:>12,d}")
rl = 2.5 + 0.3; nb = 4 / 3 * np.pi * rl**3 * 0.84
print(f"   3D liquid, list radius rc + skin = {rl}: neighbours per atom = 4/3 pi {rl}^3 x 0.84 = {nb:.1f};  N = 1000 -> {1000 * nb / 2:,.0f} pairs ({1000 * 999 / 2 / (1000 * nb / 2):.0f}x fewer)")
print(f"   rebuild when an atom has moved half the skin, 0.15; at speed ~1.5 and dt = 0.005: every 0.15/(1.5 x 0.005) = {0.15 / (1.5 * 0.005):.0f} steps")
for skin in [0.3, 0.6]:
    rl = 2.5 + skin; nb = 4 / 3 * np.pi * rl**3 * 0.84
    print(f"   skin {skin}: list radius {rl}, {nb:.1f} neighbours per atom, rebuild every ~{skin / 2 / (1.5 * 0.005):.0f} steps")
`,
  },

  exercises: [
    { q: 'Krypton has σ = 3.60 Å, ε/k_B = 171 K and mass 83.80 g/mol. What temperature (in K) and what time step (in fs) correspond to T* = 0.72 and dt* = 0.005?', a: 'T = 0.72 × 171 = 123 K. τ = σ√(m/ε): ε = 171 × 1.3806×10⁻²³ = 2.361×10⁻²¹ J, m = 83.80 × 1.6605×10⁻²⁷ = 1.3915×10⁻²⁵ kg, √(m/ε) = √(5.894×10⁻⁵) = 7.677×10⁻³, τ = 3.60×10⁻¹⁰ × 7.677×10⁻³ = 2.764 ps, so dt = 0.005 × 2.764 ps = 13.8 fs.' },
    { q: 'In a 3D box of N = 4 atoms with zero total momentum, how many degrees of freedom does the kinetic energy have? By how much would T be wrong if you used 3N?', a: 'N_f = 3N − 3 = 9. Using 12 gives T too low by a factor 9/12, i.e. 25 % low.' },
    { q: 'Atom positions in our box: (0.4, 1.0) and (3.3, 1.2). What is their minimum-image distance? Are they inside the cutoff 1.8?', a: 'dx = 2.9 → 2.9 − 3.6 × round(0.806) = 2.9 − 3.6 = −0.7; dy = 0.2. r = √(0.49 + 0.04) = 0.728. Yes, they are inside the cutoff (in fact strongly repelling: U(0.728) ≈ +154).' },
  ],
};
