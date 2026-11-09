import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-mc',
  lab: true,
  summary: 'Molecular dynamics is not the only way to sample a liquid. Monte Carlo (MC) makes random trial moves and accepts them with the Boltzmann probability; it needs no forces and no time step, but it also has no time. We do Metropolis moves by hand, check that MC and MD give the same averages for the 9-atom box, tune the MC step size, and see what MC cannot do (dynamics) and what MD cannot do (jumps over barriers).',
  minutes: 30,
  sim: 'lab15',
  learn: [
    'Make Metropolis moves by hand and see why they sample the Boltzmann distribution',
    'Check that MC and MD give the same averages, and tune the MC step size for efficiency',
    'See that MC "diffusion" is meaningless, and that smart MC moves cross barriers MD cannot',
  ],
  simHint: 'Step size: the 9 atoms move by Metropolis; small steps are almost always accepted but barely change anything, large steps are mostly rejected. Watch the acceptance and how fast the energy trace forgets. Barrier: a walker in a double well with Langevin MD, small MC steps, or MC with "flip" moves; only the flips give the 50/50 answer quickly.',

  setup: `
<p>The Metropolis recipe for one trial move:</p>
<div class="calc">1. pick an atom, displace it by a random amount up to ±δ
2. compute the energy change ΔU
3. accept with probability  min(1, e^(−ΔU/k_BT))
   (draw a random number R in [0, 1); accept if R < that)
4. if rejected, the old configuration is counted again

one sweep = N trial moves (each atom tried once on average)</div>
<p>No velocities, no forces, no time step. Systems: 3 atoms on a line, a harmonic well, the 9-atom 2D LJ box of earlier labs (k_BT = 0.72), and the double well of Lab 7.4.</p>`,

  steps: [
    {
      title: 'Metropolis moves by hand',
      why: 'Moves that lower the energy are always accepted; moves that raise it are accepted sometimes, with exactly the Boltzmann factor. Three LJ atoms on a line are enough to see the rule at work.',
      work: `<div class="calc">atoms at x = 0, 1.12, 2.30,  k_BT = 0.72,  U = −1.9595

move atom 2 to 1.20:  U = −1.9012,  ΔU = +0.0583
   p = e^(−0.0583/0.72) = <b>0.9222</b>;  R = 0.40 < p: accept
move atom 2 to 1.05:  U = −1.5580,  ΔU = +0.4015
   p = e^(−0.4015/0.72) = <b>0.5726</b>;  R = 0.40: accept
                                      R = 0.90: <b>reject</b>
   (atom 2 stays at 1.12 and that state is counted again)</div>`,
      see: 'A small uphill move is accepted 92 % of the time, a bigger one 57 %; a move into a much higher energy would almost never be accepted. Counting the old state again after a rejection is essential: it is what gives low-energy states their larger weight. Leaving it out (only recording accepted states) is a classic bug that biases every average.',
    },
    {
      title: 'Does Metropolis really give Boltzmann?',
      why: 'The acceptance rule is built so that, in the long run, each state is visited in proportion to e^(−U/k_BT) (detailed balance: the flow from A to B equals the flow from B to A). Test it on a harmonic well, where every average is known.',
      work: `<div class="calc">U = x²/2, k_BT = 1, trial step uniform in [−1, 1]
200 000 moves, acceptance 0.804

⟨x²⟩ = 1.0071        (exact 1)
⟨x⁴⟩ = 3.0644        (exact 3)
fraction |x| > 2 = 0.0467   (exact 0.0455)</div>`,
      see: 'Even the tails (|x| > 2, visited 4.6 % of the time) come out right. MC needs only energies, never forces, which makes it the natural method for potentials without smooth derivatives (hard spheres, lattice models, discontinuous cutoffs) and for moves that have no physical path at all (step 6).',
    },
    {
      title: 'The same averages from MC and MD',
      why: 'MC and MD sample the same canonical ensemble, so equilibrium averages must agree within error bars. This is a strong check of both codes, and it brings out details that are easy to miss.',
      work: `<div class="tbl"><table>${head(['9-atom LJ box, k_BT = 0.72', '⟨U⟩/N', '⟨P⟩'])}
${row(['MD (Bussi), forces cut at r_c', '−1.8093 ± 0.0066', '0.644 ± 0.046'])}
${row(['MC, energy cut and shifted at r_c', '−1.8130 ± 0.0014', '0.622 ± 0.014'])}
${row(['MC, energy cut without shift (jump −0.114)', '−1.8143 ± 0.0023', '0.598 ± 0.023'])}
${row(['MD with a net momentum left in', '<b>−1.8549 ± 0.0069</b>', '<b>0.311 ± 0.050</b>'])}</table></div>
<div class="calc">MC pressure: kinetic part N k_BT / A = 9 × 0.72 / 12.96 = 0.500
added exactly (MC has no velocities), plus the virial</div>`,
      see: 'MC and MD agree within their error bars. Two details matter. First, MD with truncated forces samples the energy that is shifted to zero at r_c; MC must use the same shifted energy to be the same model (here the difference is within the errors, but for a short cutoff it is not). Second, our first MD run disagreed badly (−1.855, P = 0.31): the initial velocities had a net momentum, the thermostat kept part of the kinetic energy in the drift of the whole box, and the atoms were effectively colder. Comparing with MC exposed the bug. MC gives the kinetic part of the pressure exactly, so its pressure is less noisy (± 0.014 vs ± 0.046).',
    },
    {
      title: 'Step size and efficiency',
      why: 'The only MC parameter is the maximum step δ. Tiny steps are almost always accepted but change little; huge steps change a lot but are almost always rejected. The right measure is not the acceptance but how many sweeps it takes to get an independent sample, the statistical inefficiency g of Lab 6.8.',
      work: `<div class="tbl"><table>${head(['step δ', 'acceptance', 'g (sweeps)', 'independent samples per 1000 sweeps'])}
${row(['0.02', '0.938', '136.5', '7'])}
${row(['0.05', '0.842', '50.7', '20'])}
${row(['0.10', '0.702', '15.6', '64'])}
${row(['0.15', '0.579', '<b>9.4</b>', '<b>107</b>'])}
${row(['0.25', '0.400', '9.7', '103'])}
${row(['0.40', '0.238', '9.9', '101'])}
${row(['0.80', '0.076', '28.3', '35'])}</table></div>
<div class="calc">MD for comparison: g of U = 202 steps (1.0 time unit)
cost: one sweep ≈ one MD step (both need all N(N − 1)/2 pair terms)</div>`,
      see: 'Efficiency is flat and best for acceptances between about 25 % and 60 % (g ≈ 10 sweeps), and drops by a factor of 15 for timid steps (94 % acceptance) and by 3 for reckless ones (8 %). The old rule "aim for 50 %" is a rough guide; measuring g is the real test. For this small atomic system MC even decorrelates the energy with less work than MD (10 sweeps vs 200 steps). For molecules with bonds and for large systems, single-atom moves must be tiny to be accepted, and MD, which moves all atoms together along the forces, usually wins; that is why liquids are mostly simulated with MD.',
    },
    {
      title: 'What MC cannot do: dynamics',
      why: 'An MC sweep is not a unit of time. Its "trajectory" moves however the trial moves let it. Compute a mean squared displacement per sweep anyway and see what it means.',
      work: `<div class="calc">MSD after 200 sweeps -> "D" = MSD / (4 × 200):
step 0.05:  MSD 0.105   "D" = 0.00013 per sweep
step 0.15:  MSD 0.449   "D" = 0.00056 per sweep
step 0.40:  MSD 0.985   "D" = 0.00123 per sweep

MD: MSD after 10 time units = 2.397  ->  D = <b>0.0599</b></div>`,
      see: 'The MC "diffusion coefficient" changes tenfold with the step size, a parameter we chose for efficiency: it says nothing about the liquid. MC gives equilibrium averages (energies, pressures, g(r), free energies) but no diffusion, viscosity, conductivity, lifetimes or spectra. Those need real dynamics, i.e. MD (or kinetic Monte Carlo with physical rates, a different method). The MD value, 0.0599, matches Lab 6.7.',
    },
    {
      title: 'What MD cannot do: jump',
      why: 'MD must follow a continuous path, so a barrier of 5 k_BT is crossed rarely. MC can propose any move whose acceptance it can compute, including a jump straight to the other well. In the double well U(s) = 5(s² − 1)², a "flip" s → −s has ΔU = 0.',
      work: `<div class="calc">exact: half the time on each side
                                   right fraction   crossings
Langevin MD, 400 000 steps         0.532            1
MC, small steps, 400 000 moves     0.509            116
MC + 10 % flips, 400 000 moves     <b>0.501</b>            <b>36 577</b>

flip: U(−s) = U(s)  ->  ΔU = 0  ->  e^0 = 1: always accepted</div>`,
      see: 'MD crossed once in 400 000 steps, so its 0.532 is luck. MC with flip moves crosses tens of thousands of times and nails 0.501. This is where MC shines: moves with no physical path (swapping two molecules, inserting and deleting particles in the grand canonical ensemble, changing the box volume, configurational-bias moves for chains, cluster moves near critical points). The limitation: designing such moves needs insight into the system, and they give no dynamics. Hybrid schemes combine both (MD for dynamics-like sampling plus MC swaps or insertions; hybrid Monte Carlo uses short MD trajectories as trial moves).',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['aspect', 'MC (number)', 'MD (number)'])}
${row(['needs', 'energies only', 'forces and a time step'])}
${row(['equilibrium averages', '⟨U⟩/N −1.8130 ± 0.0014', '−1.8093 ± 0.0066 (agree)'])}
${row(['pressure', 'kinetic part exact: ± 0.014', 'both parts noisy: ± 0.046'])}
${row(['efficiency knob', 'step δ: best at 25–60 % acceptance (g ≈ 10)', 'time step (Lab 1)'])}
${row(['dynamics', '"D" depends on δ (0.00013–0.00123)', 'physical D = 0.0599'])}
${row(['barriers', 'flip moves: 36 577 crossings', '1 crossing in 400 000 steps'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the Metropolis moves by hand, the harmonic-well test, MC vs MD averages (including the net-momentum bug), the step-size table, MC "diffusion" and the double-well comparison. Run <code>pip install numpy</code>, then <code>python lab_mc.py</code> (about 1.5 minutes).',
    src: `# Lab 3.9: Monte Carlo vs molecular dynamics by hand. Metropolis steps, the same averages from MC and MD, step size and
# efficiency, why MC has no dynamics, and moves MD cannot make.
import numpy as np, time
LJ = lambda r: 4 * (r**-12 - r**-6)

# ---------- Step 1: one Metropolis move by hand ----------
kT = 0.72; x = np.array([0.0, 1.12, 2.30])
def U3(x): return LJ(abs(x[1] - x[0])) + LJ(abs(x[2] - x[1])) + LJ(abs(x[2] - x[0]))
print(f"Step 1: three LJ atoms on a line at x = {x.tolist()}, kT = {kT}")
for new, R in ((1.20, 0.40), (1.05, 0.40), (1.05, 0.90)):
    xn = x.copy(); xn[1] = new; dU = U3(xn) - U3(x); p = min(1.0, np.exp(-dU / kT))
    print(f"   move atom 2 to {new}: U {U3(x):.4f} -> {U3(xn):.4f}, dU = {dU:+.4f};  acceptance min(1, exp(-dU/kT)) = {p:.4f};  random number {R}: {'accept' if R < p else 'reject (atom stays, old state counted again)'}")

# ---------- Step 2: Metropolis reproduces Boltzmann (harmonic well) ----------
r = np.random.default_rng(0); s = 0.0; xs = []; acc = 0; d = 1.0
for k in range(200000):
    sn = s + d * (r.random() - 0.5) * 2
    if r.random() < np.exp(-0.5 * (sn**2 - s**2)): s = sn; acc += 1
    xs.append(s)
xs = np.array(xs)
print(f"\\nStep 2: U = x^2/2, kT = 1, trial step uniform in [-1, 1], 200 000 moves: acceptance {acc / len(xs):.3f}")
print(f"   <x^2> = {np.mean(xs**2):.4f} (exact 1);  <x^4> = {np.mean(xs**4):.4f} (exact 3);  fraction |x| > 2: {np.mean(np.abs(xs) > 2):.4f} (exact {0.0455:.4f})")

# ---------- 2D LJ, 9 atoms: MC and MD ----------
L, N, RC, T0 = 3.6, 9, 1.8, 0.72
LAT = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(3) for i in range(3)])
USH = 4 * (RC**-12 - RC**-6)                                    # U(rc) = -0.1141: the jump of a truncated potential
def pair_terms(x, i, shift=True):
    d = x - x[i]; d -= L * np.round(d / L); r2 = (d**2).sum(1); r2[i] = np.inf; m = r2 < RC * RC
    ir6 = np.where(m, r2**-3.0, 0); return (4 * (ir6**2 - ir6) - (USH if shift else 0) * m).sum(), 0.0
def totals(x):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1)[np.triu_indices(N, 1)]; m = r2 < RC * RC
    ir6 = np.where(m, r2**-3.0, 0); return (4 * (ir6**2 - ir6)).sum(), (np.where(m, 24 * (2 * ir6**2 - ir6), 0)).sum()
def mc(nsweep, delta, seed, keep_x=False, shift=True):
    rr = np.random.default_rng(seed); x = LAT.copy(); Us, Ps, acc, X = [], [], 0, []
    for sw in range(-500, nsweep):
        for m in range(N):
            i = rr.integers(N); old = x[i].copy(); u0, _ = pair_terms(x, i, shift)
            x[i] = (x[i] + delta * (rr.random(2) - 0.5) * 2) % L; u1, _ = pair_terms(x, i, shift)
            if rr.random() < np.exp(-(u1 - u0) / T0): acc += sw >= 0
            else: x[i] = old
        if sw >= 0:
            u, w = totals(x); Us.append(u); Ps.append((N * T0 + 0.5 * w) / L**2)
            if keep_x: X.append(x.copy())
    return np.array(Us), np.array(Ps), acc / (nsweep * N), np.array(X)
def md(nsteps, seed, dt=0.005, every=1, keep_x=False, zero_momentum=True):
    rr = np.random.default_rng(seed); x = LAT.copy(); v = rr.normal(size=(N, 2)) * np.sqrt(T0); NF = 16
    if zero_momentum: v -= v.mean(0)
    def F(x):
        d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); mm = r2 < RC * RC
        ir6 = np.where(mm, r2**-3.0, 0); ff = np.where(mm, 24 * (2 * ir6**2 - ir6) / r2, 0); return -(ff[:, :, None] * d).sum(1)
    f = F(x); Us, Ps, X = [], [], []; cb = np.exp(-dt / 1.0); u_ = x.copy()
    for s in range(-2000, nsteps):
        v += 0.5 * dt * f; x += dt * v; u_ += dt * v; x %= L; f = F(x); v += 0.5 * dt * f
        K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = rr.normal(); S2 = rr.chisquare(NF - 1)          # Bussi, tau = 1
        v *= np.sqrt(max(cb + (1 - cb) * Kt / (NF * K) * (R1**2 + S2) + 2 * R1 * np.sqrt(cb * (1 - cb) * Kt / (NF * K)), 0))
        if s >= 0 and s % every == 0:
            u, w = totals(x); Us.append(u); Ps.append((N * T0 + 0.5 * w) / L**2)
            if keep_x: X.append(u_.copy())
    return np.array(Us), np.array(Ps), np.array(X)
def block_err(a, nb=20):
    m = a[:len(a) // nb * nb].reshape(nb, -1).mean(1); return m.std(ddof=1) / np.sqrt(nb)
def g_ineff(a):
    e1 = a.std(ddof=1) / np.sqrt(len(a)); return (block_err(a) / e1)**2

t0 = time.time(); Umd, Pmd, _ = md(100000, 2, every=1); Umc, Pmc, acc, _ = mc(40000, 0.15, 1); Umt, Pmt, _, _ = mc(40000, 0.15, 1, shift=False)
print(f"\\nStep 3: 9-atom 2D LJ box, kT = 0.72 ({time.time() - t0:.0f} s).  MD: 100 000 steps (Bussi);  MC: 40 000 sweeps, step 0.15 (acceptance {acc:.2f})")
print(f"   observable: U = sum of plain LJ over pairs closer than rc, per atom;  P with the kinetic part N kT/A = {N * T0 / L**2:.3f} (exact in MC)")
print(f"   MD (forces cut at rc):                  <U>/N = {Umd.mean() / N:.4f} +- {block_err(Umd) / N:.4f}   <P> = {Pmd.mean():.4f} +- {block_err(Pmd):.4f}")
print(f"   MC, energy cut and shifted at rc:       <U>/N = {Umc.mean() / N:.4f} +- {block_err(Umc) / N:.4f}   <P> = {Pmc.mean():.4f} +- {block_err(Pmc):.4f}")
print(f"   MC, energy cut without shift (jump {USH:.4f} at rc): <U>/N = {Umt.mean() / N:.4f} +- {block_err(Umt) / N:.4f}   <P> = {Pmt.mean():.4f} +- {block_err(Pmt):.4f}")
Ub, Pb, _ = md(100000, 2, every=1, zero_momentum=False)
print(f"   MD started with a net momentum (forgot to remove the drift): <U>/N = {Ub.mean() / N:.4f} +- {block_err(Ub) / N:.4f}   <P> = {Pb.mean():.4f} +- {block_err(Pb):.4f}")

# ---------- Step 4: step size, acceptance, efficiency ----------
print("\\nStep 4: MC step size vs acceptance and statistical inefficiency g of U (in sweeps); efficiency = 1/g per sweep")
t0 = time.time(); eff = {}
for delta in (0.02, 0.05, 0.1, 0.15, 0.25, 0.4, 0.8):
    U_, _, a_, _ = mc(40000, delta, 3); g = g_ineff(U_); eff[delta] = (a_, g)
    print(f"   step {delta:4.2f}: acceptance {a_:.3f}   g = {g:6.1f} sweeps   independent samples per 1000 sweeps = {1000 / g:5.1f}")
print(f"   ({time.time() - t0:.0f} s)")
gmd = g_ineff(Umd); print(f"   MD: g of U = {gmd:.0f} steps = {gmd * 0.005:.2f} time units;  one MD step costs about one sweep (each needs all N(N-1)/2 pair terms)")

# ---------- Step 5: no dynamics in MC ----------
print("\\nStep 5: 'diffusion' from MC sweeps depends on the arbitrary step size; MD gives a physical D")
for delta in (0.05, 0.15, 0.4):
    _, _, a_, X = mc(4000, delta, 5, keep_x=True)
    # unwrap: MC moves are small (< L/2), so unwrap with the minimum image of each sweep's change
    dX = np.diff(X, axis=0); dX -= L * np.round(dX / L); Uw = np.concatenate([X[:1], X[0] + np.cumsum(dX, 0)])
    lag = 200; msd = ((Uw[lag:] - Uw[:-lag])**2).sum(-1).mean(); print(f"   MC step {delta:4.2f}: MSD after {lag} sweeps = {msd:.3f}  ->  'D' = MSD/(4 x {lag}) = {msd / (4 * lag):.5f} per sweep")
_, _, Xmd = md(40000, 6, every=10, keep_x=True)
lag = 200; msd = ((Xmd[lag:] - Xmd[:-lag])**2).sum(-1).mean(); print(f"   MD: MSD after {lag * 10} steps ({lag * 10 * 0.005:.0f} time units) = {msd:.3f}  ->  D = {msd / (4 * lag * 10 * 0.005):.4f} (physical, time units)")

# ---------- Step 6: moves MD cannot make ----------
Uw2 = lambda s: 5.0 * (s**2 - 1)**2; dUw2 = lambda s: 20.0 * s * (s**2 - 1)
print("\\nStep 6: double well U = 5(s^2 - 1)^2, kT = 1: fraction of time on the right (exact 0.5)")
rr = np.random.default_rng(7); s = -1.0; v = 0.0; c = np.exp(-0.005); sn = np.sqrt(1 - c * c); right = 0; cross = 0; side = -1; n = 400000
for k in range(n):
    v += -0.0005 * dUw2(s); s += 0.001 * v; v += -0.0005 * dUw2(s); v = c * v + sn * rr.normal(); right += s > 0
    if s * side < 0 and abs(s) > 0.8: cross += 1; side = -side
print(f"   Langevin MD, {n} steps: right fraction {right / n:.3f}, crossings {cross}")
for name, pflip in (("MC, small steps only", 0.0), ("MC, 10 % flip moves s -> -s", 0.1)):
    rr = np.random.default_rng(8); s = -1.0; right = 0; cross = 0; side = -1; acc_f = [0, 0]
    for k in range(n):
        if rr.random() < pflip: sn_ = -s; acc_f[1] += 1
        else: sn_ = s + 0.3 * (rr.random() - 0.5) * 2
        dU = Uw2(sn_) - Uw2(s)
        if rr.random() < np.exp(-dU): s = sn_; acc_f[0] += pflip > 0 and abs(sn_ + s) < 1e-12
        right += s > 0
        if s * side < 0 and abs(s) > 0.8: cross += 1; side = -side
    print(f"   {name}, {n} moves: right fraction {right / n:.3f}, crossings {cross}")
print("   flip move: U(-s) = U(s), so dU = 0 and exp(-0/kT) = 1: always accepted")
`,
  },

  exercises: [
    { q: 'A trial move gives ΔU = +1.5 k_BT. What is the acceptance probability? What if ΔU = −1.5 k_BT?', a: 'e^(−1.5) = 0.223, so about 1 in 4.5 such moves is accepted. For ΔU = −1.5 the move lowers the energy and is always accepted (probability 1).' },
    { q: 'A program records only accepted configurations. For a two-state system with energies 0 and 1 k_BT, moves always proposing the other state, what fraction of time does it report in the low state, and what is correct?', a: 'From the high state every move (downhill) is accepted; from the low state only a fraction e^(−1) = 0.368. Correct Metropolis counts the low state again after each rejection: P(low) = 1/(1 + e^(−1)) = 0.731. Recording only accepted moves alternates low, high, low, high: it reports 0.5, completely wrong.' },
    { q: 'Why can MC not give the viscosity of a liquid, even with a perfectly tuned step size?', a: 'The viscosity comes from how stress fluctuations decay in time (Green–Kubo, Lab 6.10). MC has no physical time: the order of configurations depends on the trial moves, so time correlation functions, and everything computed from them (D, η, σ, lifetimes), have no meaning.' },
  ],
};
