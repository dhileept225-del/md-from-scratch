import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-structure-factor',
  lab: true,
  summary: 'X-ray and neutron experiments do not measure g(r); they measure the structure factor S(k), how strongly density waves of wavevector k are present. We compute S(k) by hand for four atoms, then for the LJ liquid in two ways: directly from positions and as a Fourier transform of g(r). The numbers show which k values a box allows, what cutting g(r) at L/2 does, and how S(k → 0) is tied to the compressibility.',
  minutes: 30,
  sim: 'lab13',
  learn: [
    'Compute S(k) from positions by hand, and see why a lattice gives a Bragg peak and a liquid a broad one',
    'Get S(k) from g(r) by a Fourier transform, and see the truncation error and the window that fixes it',
    'Link S(k → 0) to the compressibility, and know the smallest k a box can show',
  ],
  simHint: 'A 2D liquid of soft disks evolves by Monte Carlo. Change the density and watch S(k) (dots, computed directly from positions at the k values the box allows) and the transform of g(r) (line). Switch the window on and off: without it the transform wiggles and goes wrong at small k; with it the small-k values agree but the main peak is lowered.',

  setup: `
<p>A density wave e^(ik·r) has wavelength 2π/k. S(k) measures how strongly the particles line up with it:</p>
<div class="calc">S(k) = (1/N) |Σⱼ e^(i k·rⱼ)|²
     = (1/N) [ (Σ cos k·rⱼ)² + (Σ sin k·rⱼ)² ]

periodic box of side L:  only k = (2π/L)(nₓ, n_y, n_z) are allowed
from g(r) (3D):  S(k) = 1 + 4πρ ∫ r² (g(r) − 1) sin(kr)/(kr) dr</div>
<p>System: the 3D LJ liquid of Labs 6.9–6.11 (500 atoms, ρ = 0.8, T = 1.0, L = 8.55), 50 frames.</p>`,

  steps: [
    {
      title: 'S(k) of four atoms by hand',
      why: 'The definition is just two sums of cosines and sines. On four atoms in a 1D box of length 4, we can see exactly when the waves add up and when they cancel.',
      work: `<div class="calc">lattice  x = 0, 1, 2, 3   (L = 4, k = 2πn/4)
n = 1  k = 1.571  Σcos = 0,  Σsin = 0      S = <b>0</b>
n = 2  k = 3.142  0, 0                    S = <b>0</b>
n = 3  k = 4.712  0, 0                    S = <b>0</b>
n = 4  k = 6.283  Σcos = 4,  Σsin = 0      S = 16/4 = <b>4 = N</b>

disordered  x = 0, 1.1, 1.9, 3.2
n = 1  Σcos = 0.165,  Σsin = 0.193   S = 0.0161
n = 2  Σcos = 0.191,  Σsin = −1.206  S = 0.3726
n = 3  Σcos = −0.246, Σsin = 0.151   S = 0.0208
n = 4  Σcos = 2.927,  Σsin = 0.951   S = <b>2.368</b></div>`,
      see: 'On the lattice every wave cancels exactly, except the one whose wavelength equals the spacing (k = 2π/1): there all four atoms are in phase and S = N, a Bragg peak that grows with the number of atoms. Shift the atoms a little and the cancellation is no longer perfect: the peak drops to 2.37 and is shared with the other k values. A liquid is the extreme case: no Bragg peaks, only a broad maximum near k ≈ 2π/(neighbour distance).',
    },
    {
      title: 'The LJ liquid, directly from positions',
      why: 'For a 3D box, every allowed k vector gives one value of S; vectors of the same length are averaged. The smallest k is set by the box: the longest wave that fits is L.',
      work: `<div class="calc">smallest k = 2π/L = 2π/8.550 = <b>0.735</b>
next:  0.735 × √2 = 1.039,  × √3 = 1.273,  × 2 = 1.470</div>
<div class="tbl"><table>${head(['k', '0.735', '1.039', '2.0', '4.0', '6.0', '6.625', '8.0', '10.0'])}
${row(['S(k)', '0.054', '0.050', '0.062', '0.128', '1.226', '<b>2.424</b>', '0.950', '0.702'])}</table></div>
<div class="calc">main peak at k = 6.625  ->  2π/k = 0.948
(about the nearest-neighbour distance, ≈ 1.1)</div>`,
      see: 'S(k) has a sharp main peak near k = 6.6, the signature of the shell of neighbours, damped oscillations around 1 at large k, and very small values (≈ 0.05) at small k: long density waves are strongly suppressed because the liquid is hard to compress. The box decides what we can see: nothing exists below k = 0.735, and the first few k values have only 3–6 vectors each, so they are noisy. To see structure on longer length scales (clusters, nano-segregation in ionic liquids, the "pre-peak") the box must be several times larger than that length.',
    },
    {
      title: 'S(k) from g(r): truncation and the window',
      why: 'S(k) can also be computed from g(r), which most analyses already have. The transform needs g(r) to infinity, but we only know it up to L/2. Cutting the integral there causes errors; a window function that smoothly brings the integrand to zero (here the Lorch window, sin(πr/R)/(πr/R)) is the usual fix.',
      work: `<div class="calc">one term of the integral, r = 1.1008, g = 2.6078, k = 6.75:
r² (g − 1) sin(kr)/(kr) dr = <b>0.005109</b>   (sum 200 such terms)</div>
<div class="tbl"><table>${head(['k', 'direct', 'from g(r)', 'from g(r), window'])}
${row(['0.735', '0.054', '0.046', '0.046'])}
${row(['1.039', '0.050', '<b>−0.085</b>', '0.050'])}
${row(['1.470', '0.053', '0.052', '0.057'])}
${row(['2.0', '0.062', '0.120', '0.059'])}
${row(['4.0', '0.128', '0.075', '0.144'])}
${row(['6.625', '2.424', '2.274', '<b>1.923</b>'])}
${row(['8.0', '0.950', '0.907', '1.009'])}
${row(['10.0', '0.702', '0.717', '0.724'])}</table></div>
<div class="calc">k = 0 from g(r):  <b>0.636</b> (raw),  0.046 (window)</div>`,
      see: 'At large k the transform works (10.0: 0.717 vs 0.702). At small k the raw transform goes wrong, even negative (−0.085 at k = 1.04, impossible for |Σe^(ik·r)|²), and k = 0 comes out 13 times too large (0.636), because the cut at L/2 = 4.3 leaves a step in r²(g − 1) whose transform rings. The window removes the ringing (0.050, 0.046) but smears the main peak (1.92 instead of 2.42). So: compute S(k) directly from positions when the box allows it; use the transform of g(r) for k values between the allowed ones, with a window, and know that it trades small-k accuracy against peak height.',
    },
    {
      title: 'S(0) and the compressibility',
      why: 'As k → 0 the structure factor measures how much the number of particles in a large region fluctuates, which is set by how compressible the liquid is: S(0) = ρk_BTκ_T = k_BT / (∂P/∂ρ). This links structure to thermodynamics and gives a check of both.',
      work: `<div class="calc">pressure at ρ = 0.78, 0.80, 0.82 (256 atoms, 2 runs each):
2.885, 3.294, 3.711
∂P/∂ρ = (3.711 − 2.885)/0.04 = <b>20.64</b>
S(0) = k_BT / (∂P/∂ρ) = 1/20.64 = <b>0.0485</b>
κ_T = 1/(ρ ∂P/∂ρ) = 0.0606       (ideal gas: S(0) = 1)

direct S at the smallest k:  0.0535, 0.0501, 0.0455, 0.0526
g(r) tail offset expected (Lab 6.9):  S(0)/N = 0.0485/500 = 1×10⁻⁴</div>`,
      see: 'The pressure route (0.0485) and the small-k structure factor (0.045–0.054) agree: a dense liquid has density fluctuations 20 times weaker than an ideal gas. This is also why the 1/N offset of g(r) in Lab 6.9 is tiny for a dense liquid (10⁻⁴) but large for a gas (1/N). For mixtures and ionic systems, the same idea gives partial structure factors: the charge–charge structure factor of an electrolyte must go to zero as k² at small k (perfect screening), a strong test of a simulation\'s electrostatics. Limitation: S(0) from a few small-k points needs good statistics and a large box; the pressure route needs runs at several densities.',
    },
    {
      title: 'Comparing with experiment',
      why: 'To compare a simulated S(k) with X-ray or neutron data, the units must match and the atoms must be weighted the way the experiment sees them.',
      work: `<div class="calc">reduced -> real k:  k[Å⁻¹] = k* / σ[Å]
argon (σ = 3.405 Å):  peak k* = 6.625 -> 1.95 Å⁻¹
smallest k of the 500-atom box: 0.735/3.405 = 0.216 Å⁻¹

X-ray:    each atom weighted by its form factor f(k)
          (≈ number of electrons, falls with k)
neutron:  each element weighted by its scattering length b
          (isotope-dependent: H and D differ in sign)
mixture:  S(k) = Σ_αβ wα wβ S_αβ(k) / (Σ wα)²</div>`,
      see: 'For a one-component liquid the simulated S(k) can be compared directly once k is converted (liquid argon\'s main peak is indeed near 2 Å⁻¹). For molecules, mixtures and ionic liquids each pair of species gets its own partial structure factor, weighted by form factors or scattering lengths, and features such as the low-k pre-peak need boxes several nanometres wide. Tools: LAMMPS <code>compute rdf</code> plus a transform, or <code>compute saed/xrd</code>; dedicated codes such as TRAVIS compute weighted S(k) directly.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)'])}
${row(['direct from positions', 'exact at allowed k: peak 2.424', 'only k ≥ 2π/L = 0.735; few vectors at small k'])}
${row(['transform of g(r), raw', 'any k; good at large k (0.717 vs 0.702)', 'small k wrong: −0.085, S(0) = 0.636'])}
${row(['transform with window', 'small k fixed (0.050, 0.046)', 'main peak lowered to 1.92'])}
${row(['S(0) from ∂P/∂ρ', '0.0485, matches small-k S ≈ 0.05', 'runs at several densities'])}
${row(['experiment', 'k* / σ; peak 1.95 Å⁻¹ for argon', 'weights per species; big boxes for pre-peaks'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the four-atom S(k) by hand, the direct S(k) of the 500-atom LJ liquid, g(r) and its transform with and without the window, and the pressure route to S(0). Run <code>pip install numpy</code>, then <code>python lab_structure_factor.py</code> (about 50 s).',
    src: `# Lab: the structure factor S(k) by hand. From positions, from g(r), the smallest k of a box, truncation, S(0) and the
# compressibility. 3D LJ liquid of Labs 6.9-6.11 (rho = 0.8, T = 1.0, cutoff 1.6).
import numpy as np, time

# ---------- Step 1: S(k) of 4 atoms by hand ----------
def Sk1d(x, L, n):
    k = 2 * np.pi * n / L; c, s = np.cos(k * x).sum(), np.sin(k * x).sum(); return k, c, s, (c * c + s * s) / len(x)
print("Step 1: S(k) = |sum_j exp(i k x_j)|^2 / N,  allowed k = 2 pi n / L.  4 atoms, L = 4")
for name, x in (("lattice x = 0, 1, 2, 3", np.array([0.0, 1.0, 2.0, 3.0])), ("disordered x = 0, 1.1, 1.9, 3.2", np.array([0.0, 1.1, 1.9, 3.2]))):
    print(f"   {name}:")
    for n in (1, 2, 3, 4):
        k, c, s, S = Sk1d(x, 4.0, n); print(f"     n = {n}: k = {k:.4f}, sum cos = {c:+.4f}, sum sin = {s:+.4f}, S = ({c:.4f}^2 + {s:.4f}^2)/4 = {S:.4f}")

# ---------- LJ liquid ----------
RC, SK, T0, dt = 1.6, 0.3, 1.0, 0.004
def fcc(nc, a):
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]); g = np.array([[i, j, k] for i in range(nc) for j in range(nc) for k in range(nc)])
    return ((g[:, None, :] + b[None]) * a).reshape(-1, 3)
def pairs(x, L):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); return np.nonzero(np.triu(r2 < (RC + SK)**2, 1))
def forces(x, L, i, j, N):
    d = x[j] - x[i]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); m = r2 < RC * RC
    ir6 = np.where(m, 1 / r2**3, 0); fs = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0); ff = fs[:, None] * d; f = np.zeros((N, 3))
    for k in range(3): f[:, k] = np.bincount(j, ff[:, k], N) - np.bincount(i, ff[:, k], N)
    return f, (fs * r2).sum()
def run(nc, rho, nsteps, seed, every=0, tau=0.5):
    N = 4 * nc**3; L = (N / rho)**(1 / 3); x = fcc(nc, L / nc); r = np.random.default_rng(seed); NF = 3 * N - 3
    v = r.normal(size=(N, 3)); v -= v.mean(0); v *= np.sqrt(NF * T0 / (v**2).sum()); i, j = pairs(x, L); xl = x.copy(); f, w = forces(x, L, i, j, N)
    c = np.exp(-dt / tau); frames = []; Ps = []
    for s in range(-3000, nsteps):
        v += 0.5 * dt * f; x += dt * v
        if np.abs((x - xl + L / 2) % L - L / 2).max() > SK / 2: x %= L; i, j = pairs(x, L); xl = x.copy()
        f, w = forces(x, L, i, j, N); v += 0.5 * dt * f
        K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = r.normal(); S2 = r.chisquare(NF - 1)
        v *= np.sqrt(max(c + (1 - c) * Kt / (NF * K) * (R1**2 + S2) + 2 * R1 * np.sqrt(c * (1 - c) * Kt / (NF * K)), 0))
        if s >= 0:
            Ps.append((N * T0 + w / 3) / L**3)                 # P = rho kT + virial / (3V)
            if every and s % every == 0: frames.append((x % L).copy())
    return np.array(frames), L, np.mean(Ps), np.std(Ps) / np.sqrt(len(Ps) / 100)

t0 = time.time(); X, L, P, _ = run(5, 0.8, 10000, 1, every=200); N = X.shape[1]; rho = N / L**3
print(f"\\nStep 2: 3D LJ liquid, N = {N}, L = {L:.3f}, {len(X)} frames ({time.time() - t0:.0f} s)")
print(f"   smallest wavevector 2 pi / L = 2 pi / {L:.3f} = {2 * np.pi / L:.4f};  next ones x sqrt(2), sqrt(3), 2: {2 * np.pi / L * np.sqrt(2):.3f}, {2 * np.pi / L * np.sqrt(3):.3f}, {4 * np.pi / L:.3f}")
nmax = int(12 * L / (2 * np.pi)); g = np.arange(-nmax, nmax + 1)
nv = np.array([[a, b, c] for a in g for b in g for c in g if (a, b, c) > (0, 0, 0) and a * a + b * b + c * c <= nmax**2], float)
kv = 2 * np.pi * nv / L; km = np.linalg.norm(kv, axis=1); S = np.zeros(len(kv))
for fr in X:
    for c0 in range(0, len(kv), 2000):
        ph = fr @ kv[c0:c0 + 2000].T; S[c0:c0 + 2000] += (np.cos(ph).sum(0)**2 + np.sin(ph).sum(0)**2) / N
S /= len(X)
kb = np.arange(0.5, 12.01, 0.25); kc = 0.5 * (kb[1:] + kb[:-1]); Sb = np.array([S[(km >= a) & (km < b)].mean() if np.any((km >= a) & (km < b)) else np.nan for a, b in zip(kb[:-1], kb[1:])])
k1 = np.argmin(np.abs(km - 2 * np.pi / L))
print(f"   S at the smallest k (6 vectors, |k| = {km[k1]:.3f}): {S[np.abs(km - km[k1]) < 1e-6].mean():.4f}")
kp = np.nanargmax(Sb); print(f"   main peak S = {Sb[kp]:.3f} at k = {kc[kp]:.3f}  ->  2 pi / k = {2 * np.pi / kc[kp]:.3f} (about the neighbour distance)")
def direct(kk, w=0.125): sel = np.abs(km - kk) < (1e-3 if kk < 1.6 else w); return S[sel].mean(), sel.sum()
print("   direct S(k) (k vectors within +-0.125 of the listed value, exact values below 1.6):")
for kk in (0.735, 1.039, 2.0, 4.0, 6.0, 6.625, 8.0, 10.0, 11.0):
    v_, n_ = direct(kk); print(f"     k = {kk:6.3f}: {v_:.3f}   ({n_} vectors)")

# g(r) from the same frames, then the Fourier transform
edges = np.linspace(0, L / 2, 201); mid = 0.5 * (edges[1:] + edges[:-1]); h = np.zeros(200); iu = np.triu_indices(N, 1)
for fr in X:
    d = fr[None] - fr[:, None]; d -= L * np.round(d / L); r = np.sqrt((d**2).sum(-1))[iu]; h += np.histogram(r, edges)[0]
gr = h / len(X) / (rho * 4 / 3 * np.pi * (edges[1:]**3 - edges[:-1]**3) * N / 2); dr = edges[1] - edges[0]; Rm = L / 2
def S_from_g(k, window):
    W = np.sinc(mid / Rm) if window else 1.0                   # Lorch window sin(pi r/R)/(pi r/R)
    return 1 + 4 * np.pi * rho * np.sum(mid**2 * (gr - 1) * np.sinc(k * mid / np.pi) * W) * dr
print(f"\\nStep 3: S(k) = 1 + 4 pi rho integral r^2 (g - 1) sin(kr)/(kr) dr, g(r) known up to L/2 = {Rm:.3f}")
i12 = np.argmin(np.abs(mid - 1.1)); print(f"   one term, r = {mid[i12]:.4f}, g = {gr[i12]:.4f}, k = 6.75: r^2 (g - 1) sin(kr)/(kr) dr = {mid[i12]**2 * (gr[i12] - 1) * np.sin(6.75 * mid[i12]) / (6.75 * mid[i12]) * dr:.6f}")
print("     k      direct   from g(r)   from g(r) with window")
for kk in (0.735, 1.039, 1.470, 2.0, 4.0, 6.625, 8.0, 10.0):
    print(f"   {kk:5.3f}   {direct(kk)[0]:6.3f}   {S_from_g(kk, False):8.3f}   {S_from_g(kk, True):8.3f}")
print(f"   k = 0 from g(r): {S_from_g(1e-6, False):.3f} (raw), {S_from_g(1e-6, True):.3f} (window)")

# ---------- Step 4: S(0) and the compressibility ----------
t0 = time.time(); Pr = {}
for rr in (0.78, 0.80, 0.82):
    Pr[rr] = np.mean([run(4, rr, 10000, 20 + k)[2] for k in range(2)])
dPdrho = (Pr[0.82] - Pr[0.78]) / 0.04
print(f"\\nStep 4: pressure at rho = 0.78, 0.80, 0.82 (256 atoms, 2 runs each, {time.time() - t0:.0f} s): {Pr[0.78]:.3f}, {Pr[0.80]:.3f}, {Pr[0.82]:.3f}")
print(f"   dP/drho = ({Pr[0.82]:.3f} - {Pr[0.78]:.3f}) / 0.04 = {dPdrho:.2f};  S(0) = kT / (dP/drho) = 1 / {dPdrho:.2f} = {T0 / dPdrho:.4f}")
print(f"   compressibility kappa_T = 1/(rho dP/drho) = {1 / (0.8 * dPdrho):.4f};  ideal gas would have S(0) = 1")
small = np.sort(np.unique(np.round(km, 4)))[:4]
print("   direct S at the smallest k values: " + ", ".join(f"k = {kk:.3f}: {S[np.abs(km - kk) < 1e-3].mean():.4f}" for kk in small))
print(f"   g(r) tail offset expected from Lab 6.9: S(0)/N = {T0 / dPdrho:.4f}/{N} = {T0 / dPdrho / N:.1e}")
`,
  },

  exercises: [
    { q: 'Three atoms in a 1D box of length 3 at x = 0, 1, 2. What is S(k) at n = 1, 2, 3?', a: 'k = 2πn/3. n = 1: phases 0, 2π/3, 4π/3; Σcos = 1 − 0.5 − 0.5 = 0, Σsin = 0 + 0.866 − 0.866 = 0, S = 0. n = 2: same cancellation, S = 0. n = 3: all phases multiples of 2π, Σcos = 3, S = 9/3 = 3 = N, the Bragg peak.' },
    { q: 'A simulated ionic liquid shows a pre-peak at k = 3.5 nm⁻¹. What is the smallest cubic box that can show this k value at all, and what size would you use?', a: 'The smallest allowed k is 2π/L, so L ≥ 2π/3.5 = 1.80 nm just to reach it, but then it is a single noisy point. To resolve the peak shape you need several k values below and around it, so a box 3–5 times larger (6–9 nm) is sensible.' },
    { q: 'A liquid has ∂P/∂ρ = 5 k_BT (in reduced units, per unit density). What is S(0)? Is it more or less compressible than our LJ liquid?', a: 'S(0) = k_BT/(∂P/∂ρ) = 1/5 = 0.2, four times our 0.0485: much more compressible (larger density fluctuations), for example a liquid near its critical point, where S(0) grows without limit.' },
  ],
};
