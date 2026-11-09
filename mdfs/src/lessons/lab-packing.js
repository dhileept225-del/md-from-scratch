import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-packing',
  lab: true,
  summary: 'Lab 2.7 started from a neat lattice of 9 atoms. A real project starts from a target density and a few hundred or thousand molecules that must be put in a box without overlapping, relaxed, heated and equilibrated before anything is measured. We do every stage for 256 Lennard-Jones atoms (liquid argon): box size from the density, random insertion and why it fails, steepest descent by hand, NVT and NPT, and three tests that tell you when the system is ready.',
  minutes: 35,
  sim: 'lab16',
  learn: [
    'Turn a target density into a box size, and see why random positions blow up the simulation',
    'Do steepest-descent minimization by hand and see what it fixes and what it does not',
    'Run NVT and NPT from different starts and decide, with numbers, when the system is equilibrated',
  ],
  simHint: 'A 64-atom 2D box at density 0.70. 1 Insert places atoms at random, refusing any closer than d to an existing atom; at d = 0 overlaps are coral, at d ≥ 1 the box jams before all atoms fit. 2 Minimize runs steepest descent until the largest force is below 3.4. 3 NVT heats to 0.72 and marks the detected end of equilibration (t0). Try NVT straight after inserting with d = 0: it blows up.',

  setup: `
<p>The recipe used by every real setup (Packmol, <code>gmx insert-molecules</code>, LAMMPS <code>create_atoms random</code>):</p>
<div class="calc">1. box size from the target density
2. place molecules without bad overlaps
3. minimize the energy (remove the remaining clashes)
4. NVT: give velocities, heat to the target temperature
5. NPT: let the density find its own value (if unknown)
6. decide when equilibrium is reached; only then measure</div>
<p>System: N = 256 LJ atoms, cutoff 2.5σ, k_BT = 0.722 and ρ* = 0.8442 (the liquid near the triple point, the classic test state). Bussi thermostat (τ = 0.5), C-rescale barostat, Δt = 0.005. Argon: σ = 3.405 Å, ε = 0.996 kJ/mol.</p>`,

  steps: [
    {
      title: 'Box size from the target density',
      why: 'Experiments give density in g/cm³; a simulation needs a number of molecules and a box length. The conversion is one line, but getting it wrong by a few percent shifts every property you compute (pressure most of all).',
      work: `<div class="calc">volume per molecule  v = M / (ρ N_A)

argon: M = 39.948 g/mol, ρ = 1.40 g/cm³
  v = 39.948 / (1.40 × 6.022×10²³) = 4.738×10⁻²³ cm³
    = 47.38 Å³   (1 cm³ = 10²⁴ Å³)
  256 atoms: V = 12 130 Å³,  L = V^(1/3) = <b>22.98 Å</b>
  in LJ units: L = 22.98/3.405 = 6.748 σ,  ρ* = <b>0.833</b>

water: M = 18.015, ρ = 0.997
  v = 30.00 Å³;  1000 molecules: L = <b>31.07 Å</b>

1 mol/kg NaCl in 500 waters:
  n = 500 × 0.018015 kg/mol × 1 mol/kg = 9.0 ion pairs
  density of the solution: not known in advance</div>`,
      see: 'Pure liquids with a known density can be built at the right size directly. For a solution, a mixture or a new force field, the density is part of the answer: build at a guess, then let NPT (Step 4) find it. Our test state ρ* = 0.8442 is the same liquid as argon at 1.42 g/cm³.',
    },
    {
      title: 'Random positions: why they fail, and why a tolerance only partly helps',
      why: 'The obvious way to fill a box is to put each atom at a random point. In a dense liquid, that makes many atoms land almost on top of each other, and the LJ repulsion grows as r⁻¹².',
      work: `<div class="calc">expected pairs closer than 0.8σ for random points:
  (N/2) × ρ × (4/3)π(0.8)³ = 128 × 0.8442 × 2.145 = 232
measured (one random box):            231
closest pair:                          r = 0.077σ
U/N = 3.9×10¹¹,  largest force = 1.5×10¹⁶
   (equilibrium liquid: U/N ≈ −5, forces ≈ 10)

NVE, Δt = 0.005, straight from this start:
  step 1:   T = 3.8×10²⁴
  step 50:  T = 2.4×10³⁰       (target 0.722)</div>
<p>Insertion with a tolerance: a new atom is refused if it is closer than d to any placed atom (give up after 20 000 refusals in a row).</p>
<div class="tbl"><table>${head(['d (σ)', 'placed', 'trials', 'ρ* reached', 'packing of d-spheres'])}
${row(['0.8', '256', '1 621', '0.844', '0.226'])}
${row(['0.9', '256', '20 736', '0.844', '0.322'])}
${row(['1.0', '203', '119 624', '0.669', '0.351'])}
${row(['1.1', '144', '67 746', '0.475', '0.331'])}</table></div>
<div class="calc">packing fraction  φ = ρ × (π/6) d³
random insertion jams at φ ≈ 0.38 (no room left for one more)
d = 1.0 at ρ* = 0.8442 needs φ = 0.8442 × 0.524 = 0.442 > 0.38</div>`,
      see: 'Fully random positions give forces 10¹⁵ times too large and the first step sends atoms flying. A tolerance removes the worst clashes, but a tolerance equal to the atom size cannot reach a liquid density: random filling jams long before. That is why Packmol uses a tolerance smaller than the molecules (2 Å) and why a minimization always follows. A lattice start (Lab 2.7) avoids the problem entirely, at the price of having to melt it.',
    },
    {
      title: 'Steepest descent by hand',
      why: 'Minimization moves every atom downhill along its force until no force is large. GROMACS uses a simple, robust rule: move the atom with the largest force by a step h, others in proportion; grow h when the energy falls, cut it hard when it rises.',
      work: `<div class="calc">x_new = x + h F / |F|max
U fell:  accept,  h → 1.2 h
U rose:  reject,  h → 0.2 h

dimer, r = 0.9:  U = 6.6361,  F = 138.66 on each atom
both atoms move h apart, so r grows by 2h
iter  h        new r    U new     decision
1     0.0100   0.9200    4.2826   accept
2     0.0120   0.9440    2.3349   accept
3     0.0144   0.9728    0.8493   accept
4     0.0173   1.0074   −0.1648   accept
5     0.0207   1.0488   −0.7476   accept
(the minimum is at 2^(1/6) = 1.1225)</div>
<p>The random 256-atom box, stopped at the GROMACS default emtol = 10 kJ/mol/nm = 10/(0.996/0.3405) = 3.42 ε/σ:</p>
<div class="tbl"><table>${head(['iteration', 'U/N', '|F|max', 'h'])}
${row(['0', '3.9×10¹¹', '1.6×10¹⁶', '0.010'])}
${row(['5', '2.1×10⁷', '5.7×10¹⁰', '0.025'])}
${row(['20', '134', '1.1×10⁵', '0.38'])}
${row(['50', '−4.97', '81', '0.070'])}
${row(['100', '−5.78', '31', '0.082'])}
${row(['160 (done)', '−6.00', '2.2', ''])}</table></div>
<div class="calc">closest pair after minimization: 1.05σ

NVE from here with velocities at T = 0.722:
  1000 steps: E/N −4.9182 → −4.9177   (stable)
  but T drops to 0.401

near a minimum, kinetic energy is shared equally
with potential energy (equipartition):
  T → 0.722 / 2 = 0.361 for a harmonic system
  0.401 here (the liquid is not harmonic)</div>`,
      see: 'Steepest descent removed every clash in 160 iterations; the start that exploded now conserves energy. But a minimized structure is not an equilibrium liquid: it sits at T = 0, at an energy (−6.00) far below the liquid at T = 0.722 (−5.20), and velocities given at T lose about half their energy to the potential. That is why minimization is always followed by a thermostatted run, not by a measurement.',
    },
    {
      title: 'NVT from two starts, then NPT to find the density',
      why: 'The thermostat brings the temperature to the target within a few τ. Whether the structure follows depends on the start: a lattice has to melt, a minimized random box is already disordered. The FCC order parameter (1 on the lattice, about 1/N in a liquid) shows how far along the melting is.',
      work: `<div class="tbl"><table>${head(['t', 'FCC: U/N', 'T', 'FCC order', 'random: U/N', 'T'])}
${row(['0', '−6.331', '0.717', '0.999', '−5.993', '0.723'])}
${row(['1', '−5.583', '0.552', '0.528', '−5.350', '0.700'])}
${row(['2', '−5.480', '0.629', '0.450', '−5.232', '0.682'])}
${row(['5', '−5.386', '0.732', '0.411', '−5.142', '0.708'])}
${row(['10', '−5.155', '0.737', '0.124', '−5.270', '0.696'])}
${row(['20', '−5.258', '0.703', '0.009', '−5.175', '0.751'])}
${row(['40', '−5.153', '0.749', '0.009', '−5.198', '0.717'])}</table></div>
<div class="calc">FCC: T dips to 0.55 at t = 1: melting absorbs heat
     the crystal survives at order ≈ 0.4 until t ≈ 5–10

NVT pressure (second halves, 4 runs): P* = 0.874

NPT test: start at ρ* = 0.70 (easy to pack),
C-rescale at P* = 0.874, τ_p = 2
  t     ρ*       P
  0     0.700   −4.78
  1     0.787   −0.37
  5     0.832    0.61
  10    0.861    0.66
  50    0.835    0.91
second half: ρ* = <b>0.8421 ± 0.0019</b>   (true 0.8442)</div>`,
      see: 'The temperature is right after t ≈ 1 in both runs, but the lattice start keeps crystal order until t ≈ 10: temperature alone says nothing about whether the structure has relaxed. NPT started 17% too dilute and found the density to 0.25% (1 error bar) within t ≈ 10. This is how the density of the NaCl solution in Step 1 would be found.',
    },
    {
      title: 'When is it equilibrated? Three tests',
      why: 'Everything before equilibrium must be thrown away, and everything after kept. Throwing away too little biases the average; too much wastes data. The tests below replace “it looks flat”.',
      work: `<p>(i) Choose t₀ to maximize the number of independent samples after it: N_eff = (n − n₀)/g, with g the statistical inefficiency of the data after t₀ (Lab 6.8).</p>
<div class="tbl"><table>${head(['series', 't₀', 'g after t₀', 'N_eff'])}
${row(['U/N, FCC start', '8.0', '15.2', '42'])}
${row(['U/N, random start', '2.0', '16.6', '46'])}
${row(['FCC order, FCC start', '10.0', '21.9', '27'])}
${row(['ρ*, NPT from 0.70', '11.3', '18.2', '43'])}</table></div>
<p>(ii) Half-split: compare the means of the first and second halves of what is kept.</p>
<div class="calc">U/N, FCC start, all data:
  −5.302 ± 0.084  vs  −5.190 ± 0.012:  1.3 errors, "passes"
U/N, FCC start, t > 12:
  −5.189 ± 0.016  vs  −5.195 ± 0.008:  0.3 errors
the first test should have failed: the melting trend makes
the first half look very noisy, which inflates its error bar</div>
<p>(iii) Independent replicas (different starts and seeds), t &gt; 12 kept:</p>
<div class="tbl"><table>${head(['start', 'U/N', 'P*'])}
${row(['FCC', '−5.192 ± 0.009', '0.89 ± 0.05'])}
${row(['random + minimized', '−5.200 ± 0.009', '0.87 ± 0.05'])}
${row(['FCC, seed 2', '−5.202 ± 0.016', '0.84 ± 0.09'])}
${row(['random + minimized, seed 2', '−5.204 ± 0.010', '0.83 ± 0.06'])}</table></div>
<div class="calc">spread of the four U/N means: 0.005
typical error bar:                  0.010  → consistent</div>`,
      see: 'Different quantities equilibrate at different times: the energy of the lattice start at t₀ ≈ 8, its structure at t₀ ≈ 10, the random start already at 2. Discard by the slowest quantity you care about. The half-split test can be fooled by the very trend it is meant to detect, so look at the trace as well. Replicas from different starts are the strongest test: if they agree within their error bars, no memory of the start is left.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['stage', 'what goes wrong without it', 'number'])}
${row(['box from density', 'wrong pressure, wrong everything', 'L = 22.98 Å for 256 argon'])}
${row(['careful insertion', 'overlaps: forces 10¹⁶, T = 10²⁴ after 1 step', '231 pairs closer than 0.8σ'])}
${row(['tolerance d', 'd ≥ σ jams below liquid density', 'd = 1.0 stops at ρ* = 0.67'])}
${row(['minimization', 'explosion; it does not equilibrate', '160 iterations, U/N −6.00 (T = 0)'])}
${row(['NVT', 'T right in t ≈ 1, structure in t ≈ 10', 'FCC order 0.999 → 0.009'])}
${row(['NPT', 'unknown density stays wrong', '0.70 → 0.8421 ± 0.0019 (true 0.8442)'])}
${row(['equilibration tests', 'biased averages', 't₀ = 2–11; replicas agree to 0.005'])}</table></div>`,

  code: {
    note: 'The whole lab in one script: box sizes, random and tolerance insertion, steepest descent, NVT from a lattice and from a minimized random box, NPT, and the three equilibration tests. Run <code>pip install numpy</code>, then <code>python lab_build.py</code> (about 3 minutes).',
    src: `# Lab 2.8: packing and equilibrating a real liquid (3D Lennard-Jones, reduced units)
import numpy as np, time
t0 = time.time()
NA = 6.02214076e23

print("Step 1: box size from the target density")
def box(N, M, rho):                      # M in g/mol, rho in g/cm^3 -> volume per molecule and L in angstrom
    v = M / (rho * NA) * 1e24; return v, (N * v) ** (1 / 3)
for name, N, M, rho in [("argon (liquid near 87 K)", 256, 39.948, 1.40), ("water (298 K)", 1000, 18.015, 0.997)]:
    v, L = box(N, M, rho)
    print(f"   {name}: {N} molecules, M = {M} g/mol, rho = {rho} g/cm3 -> {v:.2f} A^3 per molecule, V = {N * v:.0f} A^3, L = {L:.2f} A")
sig = 3.405; v, L = box(256, 39.948, 1.40)
print(f"   argon in LJ units (sigma = {sig} A): L = {L / sig:.3f} sigma, rho* = 256 sigma^3 / V = {256 * sig**3 / (256 * v):.4f}")
print(f"   1 mol/kg NaCl with 500 waters: n_NaCl = 500 x 0.018015 kg/mol x 1 mol/kg = {500 * 0.018015:.2f} -> 9 ion pairs; solution density unknown -> NPT needed")

# ---------- LJ machinery: N = 256, rc = 2.5, energy shifted at rc ----------
RC = 2.5; USH = 4 * (RC**-12 - RC**-6); N = 256; RHO = 0.8442; T0 = 0.722
def forces(x, L):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L); r2 = (d * d).sum(-1)
    np.fill_diagonal(r2, np.inf); m = r2 < RC * RC; ir2 = np.where(m, 1 / r2, 0); ir6 = ir2**3
    U = 0.5 * np.where(m, 4 * (ir6 * ir6 - ir6) - USH, 0).sum()
    ff = np.where(m, 24 * (2 * ir6 * ir6 - ir6) * ir2, 0)          # F/r
    F = (ff[:, :, None] * d).sum(1); W = 0.5 * np.where(m, ff * r2, 0).sum() / 3     # virial  (1/3) sum r.F
    return F, U, W, r2
def fcc(n, L):
    a = L / n; b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    return np.array([(np.array([i, j, k]) + q) * a for i in range(n) for j in range(n) for k in range(n) for q in b])
def order(x, L):                          # FCC order parameter: 1 on the lattice, ~1/N in a liquid
    a = L / 4; ks = 2 * np.pi / a * np.array([[s1, s2, s3] for s1 in (1, -1) for s2 in (1, -1) for s3 in (1, -1)])
    return np.mean(np.abs(np.exp(1j * x @ ks.T).mean(0))**2)
def bussi(v, T, tau, dt, rng):
    nf = 3 * N - 3; K = 0.5 * (v * v).sum(); Kt = 0.5 * nf * T; c = np.exp(-dt / tau); R = rng.normal(); S = (rng.normal(size=nf - 1)**2).sum()
    a2 = c + Kt / (nf * K) * (1 - c) * (R * R + S) + 2 * R * np.sqrt(c * Kt / (nf * K) * (1 - c))
    return v * np.sqrt(a2)
def temp(v): return (v * v).sum() / (3 * N - 3)
def vstart(T, rng):
    v = rng.normal(size=(N, 3)) * np.sqrt(T); v -= v.mean(0); return v * np.sqrt(T / temp(v))
def md(x, v, L, steps, rng, dt=0.005, tau=0.5, every=10, npt=None, flip=False):
    F, U, W, _ = forces(x, L); out = []
    for s in range(steps):
        v += 0.5 * dt * F; x = (x + dt * v) % L; F, U, W, _ = forces(x, L); v += 0.5 * dt * F
        v = bussi(v, T0, tau, dt, rng)
        if npt is not None:                # C-rescale: stochastic volume move, eps = ln V
            P0, taup, beta = npt; V = L**3; P = (N * temp(v) + W) / V
            de = -beta / taup * (P0 - P) * dt + np.sqrt(2 * T0 * beta * dt / (V * taup)) * rng.normal()
            s3 = np.exp(de / 3); x *= s3; L *= s3; v /= s3
        if s % every == 0:
            out.append((s * dt, U / N, temp(v), (N * temp(v) + W) / L**3, N / L**3, order(x, L)))
    return x, v, L, np.array(out)

L0 = (N / RHO) ** (1 / 3)
print(f"\\nStep 2: {N} atoms at rho* = {RHO}: L = {L0:.4f} sigma.  Random positions vs insertion with a tolerance")
rng = np.random.default_rng(1)
xr = rng.random((N, 3)) * L0; F, U, W, r2 = forces(xr, L0)
npair = (r2 < 0.8**2).sum() // 2; exp_pairs = 0.5 * N * RHO * 4 / 3 * np.pi * 0.8**3
print(f"   random: pairs closer than 0.8 sigma = {npair} (expected N/2 x rho x 4/3 pi 0.8^3 = {exp_pairs:.0f});  closest pair r = {np.sqrt(r2.min()):.4f}")
print(f"   U/N = {U / N:.3e},  largest force = {np.abs(F).max():.3e}  (equilibrium liquid: U/N ~ -5, forces ~ 10)")
x = xr.copy(); v = vstart(T0, rng); Fv, U, W, _ = forces(x, L0); E0 = U + 0.5 * (v * v).sum()
for s in range(1, 51):
    v += 0.0025 * Fv; x = (x + 0.005 * v) % L0; Fv, U, W, _ = forces(x, L0); v += 0.0025 * Fv
    if s in (1, 10, 50): print(f"   NVE from random, dt = 0.005, step {s:2d}: T = {temp(v):.3e}  (target 0.722),  E/N = {(U + 0.5 * (v * v).sum()) / N:.3e}")
print("   tolerance insertion (each new atom must be at least d from all placed atoms; give up after 20 000 failures in a row):")
for d in (0.8, 0.9, 1.0, 1.1):
    rg = np.random.default_rng(2); P = np.empty((0, 3)); tries = 0; fails = 0
    while len(P) < N and fails < 20000:
        c = rg.random(3) * L0; tries += 1
        if len(P):
            dd = P - c; dd -= L0 * np.round(dd / L0)
            if ((dd * dd).sum(1) < d * d).any(): fails += 1; continue
        P = np.vstack([P, c]); fails = 0
    phi = len(P) / L0**3 * np.pi / 6 * d**3
    print(f"   d = {d}: placed {len(P)}/{N} after {tries} trials, reached rho* = {len(P) / L0**3:.3f}, packing fraction of d-spheres {phi:.3f}  (random jamming 0.38)")
    if d == 0.9: xpack = P.copy()

print("\\nStep 3: steepest descent, GROMACS rule: x += h F/|F|max; if U fell, accept and h x 1.2, else reject and h x 0.2")
r = 0.9; h = 0.01
def Ud(r): return 4 * (r**-12 - r**-6)
def Fd(r): return 24 * (2 * r**-13 - r**-7)
print(f"   dimer by hand: r = {r}, U = {Ud(r):.4f}, F = {Fd(r):.4f} (minimum at 2^(1/6) = {2**(1/6):.4f})")
for it in range(1, 6):
    rn = r + 2 * h            # both atoms move h along their force, so r grows by 2h
    ok = Ud(rn) < Ud(r)
    print(f"   iter {it}: h = {h:.4f}, try r = {rn:.4f}: U = {Ud(rn):.4f} {'<' if ok else '>='} {Ud(r):.4f} -> {'accept, h x 1.2' if ok else 'reject, h x 0.2'}")
    if ok: r = rn; h *= 1.2
    else: h *= 0.2
def sd(x, L, emtol, maxit=5000, log=()):
    h = 0.01; F, U, _, _ = forces(x, L); it = 0; rows = []
    while it < maxit:
        fm = np.sqrt((F * F).sum(1)).max()
        if it in log: rows.append((it, U / N, fm, h))
        if fm < emtol: break
        xn = (x + h * F / fm) % L; Fn, Un, _, _ = forces(xn, L)
        if Un < U: x, F, U, h = xn, Fn, Un, h * 1.2
        else: h *= 0.2
        it += 1
    return x, U, fm, it, rows
emtol = 10 / (0.996 / 0.3405)   # GROMACS default 10 kJ/mol/nm in argon LJ units
xm, Um, fm, its, rows = sd(xr.copy(), L0, emtol, log=(0, 1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000))
print(f"   random 256-atom start, emtol = 10 kJ/mol/nm = {emtol:.2f} eps/sigma for argon")
for it_, u, f, hh in rows: print(f"      iter {it_:5d}: U/N = {u:12.4e}   |F|max = {f:10.3e}   h = {hh:.2e}")
print(f"   converged after {its} iterations: U/N = {Um / N:.4f}, |F|max = {fm:.3f}")
F, U, W, r2 = forces(xm, L0); print(f"   closest pair now {np.sqrt(r2.min()):.3f} sigma;  U/N = {U / N:.3f} is a T = 0 'inherent structure', below the liquid at T = 0.722")
x = xm.copy(); v = vstart(T0, np.random.default_rng(3)); Fv, U, W, _ = forces(x, L0); E0 = (U + 0.5 * (v * v).sum()) / N
for s in range(1, 1001):
    v += 0.0025 * Fv; x = (x + 0.005 * v) % L0; Fv, U, W, _ = forces(x, L0); v += 0.0025 * Fv
print(f"   NVE from the minimized start, dt = 0.005, 1000 steps: E/N {E0:.4f} -> {(U + 0.5 * (v * v).sum()) / N:.4f} (stable);  T = {temp(v):.3f} (started 0.722)")
print(f"   ({time.time() - t0:.0f} s so far)")

print(f"\\nStep 4: NVT at T* = {T0} (Bussi, tau = 0.5), 8000 steps (t = 40), from two different starts")
runs = {}
for lab, xs, seed in [("FCC lattice", fcc(4, L0), 11), ("random + minimized", xm.copy(), 12), ("FCC lattice, seed 2", fcc(4, L0), 21), ("random + minimized, seed 2", sd(np.random.default_rng(22).random((N, 3)) * L0, L0, emtol)[0], 22)]:
    rg = np.random.default_rng(seed); _, _, _, o = md(xs, vstart(T0, rg), L0, 8000, rg); runs[lab] = o
for lab in list(runs)[:2]:
    o = runs[lab]
    print(f"   {lab}:")
    for t_ in (0, 1, 2, 5, 10, 20, 40):
        k = min(int(t_ / 0.05), len(o) - 1); print(f"      t = {o[k, 0]:5.2f}: U/N = {o[k, 1]:7.4f}  T = {o[k, 2]:.3f}  P = {o[k, 3]:7.4f}  FCC order = {o[k, 5]:.3f}")
Pt = np.mean([runs[k][len(runs[k]) // 2:, 3].mean() for k in runs]); print(f"   NVT average pressure (second halves, 4 runs): P* = {Pt:.3f}")
print(f"   NPT test: start at rho* = 0.70 (easy to pack), C-rescale at P* = {Pt:.3f}, tau_p = 2, 10 000 steps")
Lb = (N / 0.70) ** (1 / 3); rg = np.random.default_rng(5); xb = sd(fcc(4, Lb), Lb, emtol)[0]
_, _, Lend, onpt = md(xb, vstart(T0, rg), Lb, 10000, rg, npt=(Pt, 2.0, 0.1))
for t_ in (0, 1, 2, 5, 10, 20, 50):
    k = min(int(t_ / 0.05), len(onpt) - 1); print(f"      t = {onpt[k, 0]:5.2f}: rho* = {onpt[k, 4]:.4f}  P = {onpt[k, 3]:7.4f}  U/N = {onpt[k, 1]:.4f}")
def ineff(a):
    a = a - a.mean(); n = len(a); f = np.fft.rfft(a, 2 * n); c = np.fft.irfft(f * np.conj(f))[:n] / np.arange(n, 0, -1); c /= c[0]; g = 1
    for k in range(1, n // 3):
        if c[k] < 0: break
        g += 2 * c[k]
    return g
def berr(a): return a.std() * np.sqrt(ineff(a) / len(a))
h2 = onpt[len(onpt) // 2:, 4]; print(f"   NPT density (second half): rho* = {h2.mean():.4f} +- {berr(h2):.4f}   (target {RHO})")
print(f"   ({time.time() - t0:.0f} s so far)")

print("\\nStep 5: when is it equilibrated?  (i) pick t0 that maximises the number of independent samples after it")
def pick_t0(a, dt):
    best = None
    for t in range(0, len(a) * 3 // 4, max(1, len(a) // 40)):
        g = ineff(a[t:]); ne = (len(a) - t) / g
        if best is None or ne > best[2]: best = (t, g, ne)
    return best[0] * dt, best[1], best[2]
for lab, ser, col in [("U/N, FCC start", runs["FCC lattice"], 1), ("U/N, random start", runs["random + minimized"], 1), ("FCC order, FCC start", runs["FCC lattice"], 5), ("density, NPT from 0.70", onpt, 4)]:
    t_, g, ne = pick_t0(ser[:, col], 0.05); print(f"   {lab:24s}: t0 = {t_:5.2f}, g after t0 = {g:5.1f} samples, independent samples = {ne:5.0f}")
print("   (ii) half-split test: mean of first half vs second half of the kept part, in units of their combined error")
for lab, ser, col in [("U/N, FCC start, all", runs["FCC lattice"], 1), ("U/N, FCC start, t > 12", runs["FCC lattice"][240:], 1), ("density NPT, all", onpt, 4), ("density NPT, t > 12", onpt[240:], 4)]:
    y = ser[:, col]; a, b = y[:len(y) // 2], y[len(y) // 2:]; z = (b.mean() - a.mean()) / np.hypot(berr(a), berr(b))
    print(f"   {lab:22s}: halves {a.mean():.4f} +- {berr(a):.4f} | {b.mean():.4f} +- {berr(b):.4f}  -> difference {z:+.1f} errors: {'NOT equilibrated' if abs(z) > 2 else 'consistent'}")
print("   (iii) independent replicas, t > 12 kept:")
for lab, o in runs.items():
    y = o[240:, 1]; print(f"   {lab:28s}: U/N = {y.mean():.4f} +- {berr(y):.4f}   P = {o[240:, 3].mean():.4f} +- {berr(o[240:, 3]):.4f}")
m = np.array([o[240:, 1].mean() for o in runs.values()]); print(f"   spread of the 4 means: {m.std(ddof=1):.4f} (standard deviation); compare with the error bars above")
print(f"\\ntotal time {time.time() - t0:.0f} s")
`,
  },

  exercises: [
    { q: 'You want 500 molecules of ethanol (M = 46.07 g/mol, ρ = 0.789 g/cm³). What box length do you need?', a: 'v = 46.07/(0.789 × 6.022×10²³) = 9.70×10⁻²³ cm³ = 97.0 Å³. V = 500 × 97.0 = 48 490 Å³, L = 48 490^(1/3) = 36.5 Å.' },
    { q: 'Steepest descent reduced |F|max from 10¹⁶ to 2 in 160 iterations, yet minimizers are said to be slow. When is steepest descent slow, and what do codes switch to?', a: 'It is fast at removing large clashes (big forces, simple downhill direction) but slow near the minimum in long, narrow valleys, where it zig-zags: each step follows the local force, which points across the valley more than along it. Codes switch to conjugate gradients or L-BFGS, which remember earlier steps. For MD setup this rarely matters: you only need to remove clashes, not find the exact minimum.' },
    { q: 'Our NVT pressure is P* = 0.87, while the literature equation of state for LJ at this state gives P* ≈ 0. Which is wrong?', a: 'Neither: our potential is cut at 2.5σ and the attractions beyond are missing. The standard tail correction ΔP = (16/3)πρ²[(2/3)r_c⁻⁹ − r_c⁻³] = 16.76 × 0.713 × (0.00017 − 0.064) = −0.76, so the full-potential pressure is ≈ 0.87 − 0.76 ≈ 0.1. This is also why an NPT run with a cut potential gives a slightly different density from experiment unless the tail correction is switched on.' },
  ],
};
