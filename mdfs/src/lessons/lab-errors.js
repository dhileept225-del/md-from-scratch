import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-errors',
  lab: true,
  summary: 'Every number from a simulation needs an error bar, and the textbook formula s/√N is almost always wrong for MD, because consecutive frames are not independent. We compute error bars by hand on 10 and 16 numbers, then use block averaging and the statistical inefficiency on a process with a known answer and on the pressure of our 9-atom box. The numbers show how many independent samples a run really contains, how errors combine, and how long a run must be.',
  minutes: 30,
  sim: 'lab8',
  learn: [
    'Compute a standard error by hand, then correct it for correlation with the statistical inefficiency g',
    'Do block averaging by hand and read the plateau; check it against independent runs',
    'Propagate errors to ratios and differences, discard equilibration, and estimate how long to run',
  ],
  simHint: 'Choose how strongly each value remembers the previous one (the correlation φ) and press New run. The top plot is the time series, the bottom the block-averaging curve: it rises from the naive error and levels off at the true error (amber line). Press 20 runs: the spread of the run means (cyan bar) matches the plateau, not the naive error.',

  setup: `
<p>An average from a simulation is itself a random number: run again with other random velocities and you get a slightly different value. The error bar says how much. Two facts make MD special:</p>
<div class="calc">1. frames are saved every few steps, but the system
   changes slowly: consecutive frames are correlated
2. the formula  error = s / √N  assumes independent samples

so:  N frames  ≠  N independent samples
the true number is  N_eff = N / g,  g = statistical inefficiency</div>
<p>Notation: s is the standard deviation of the samples, "error" means the standard error of the mean (one σ, about 68 % confidence).</p>`,

  steps: [
    {
      title: 'The textbook error bar',
      why: 'Start with the case where the formula is right: independent measurements. The spread of the samples tells how much one measurement varies; the error of the mean is smaller by √N because random deviations partly cancel.',
      work: `<div class="calc">10 independent pressure readings:
1.2  0.4  2.1  1.7  −0.3  1.1  0.9  2.6  0.2  1.5

mean = 11.4 / 10 = 1.140
Σ (P − mean)² = 7.0640
s = √(7.0640 / 9) = 0.8859        (divide by N − 1)
error = s / √N = 0.8859 / 3.1623 = <b>0.2802</b>

result:  P = <b>1.14 ± 0.28</b></div>`,
      see: 'One reading varies by ±0.89, the mean of ten by ±0.28. To halve the error you need four times more independent samples: errors fall only as 1/√N. This is the formula most analysis scripts use by default, and for MD data it gives error bars that are far too small, as the next steps show.',
    },
    {
      title: 'Correlated data: the statistical inefficiency',
      why: 'In MD, a frame looks much like the one before it. Such samples carry less information than independent ones. The autocorrelation ρ(k) measures how similar values k frames apart are; adding it up tells how many frames count as one independent sample.',
      work: `<div class="calc">16 consecutive readings:
1.0 1.2 1.5 1.4 1.1 0.7 0.5 0.6 0.9 1.3 1.6 1.8 1.5 1.2 0.8 0.9
mean = 1.1250,  s = 0.3786

ρ(k) = ⟨δx(t) δx(t+k)⟩ / ⟨δx²⟩
lag 1:  0.1021 / 0.1344 = <b>0.760</b>
lag 2:  0.0169 / 0.1344 = <b>0.126</b>
lag 3:  −0.556   (negative: stop adding here)

g = 1 + 2 (0.760 + 0.126) = <b>2.77</b>
N_eff = 16 / 2.77 = <b>5.8</b> independent samples

naive error      s/√16      = 0.0946
corrected error  s √(g/N)   = <b>0.1576</b>  (1.66 × larger)</div>`,
      see: 'Sixteen correlated numbers are worth about six independent ones, and the honest error bar is 1.66 times the naive one. For real MD, where frames may be saved every 10 fs while the system needs picoseconds to forget, g can be in the hundreds. The limitation: ρ(k) from a short series is noisy (here it swings to −0.56 at lag 3), so the sum must be cut off sensibly (at the first negative value is a common rule), and the result is only an estimate.',
    },
    {
      title: 'Block averaging by hand',
      why: 'Block averaging avoids computing ρ(k) at all. Cut the series into blocks, average each block, and treat the block means as samples. Once the blocks are longer than the memory of the data, the block means are independent and the naive formula becomes correct for them.',
      work: `<div class="calc">same 16 numbers
block size 1:  16 means          error = 0.0946
block size 2:   8 means
   1.10 1.45 0.90 0.55 1.10 1.70 1.35 0.85
                                 error = 0.1299
block size 4:   4 means
   1.275 0.725 1.400 1.100       error = <b>0.1468</b>
block size 8:   2 means
   1.000 1.250                   error = 0.1250

error = std(block means) / √(number of blocks)</div>`,
      see: 'The error grows with block size (0.095 → 0.130 → 0.147) and stops growing once blocks are longer than the correlation: the plateau, 0.147, agrees with the g-corrected 0.158. Beyond that, too few blocks remain (two blocks of 8 give 0.125, itself very uncertain). In practice: increase the block size until the error levels off, and read the plateau while there are still at least 20–50 blocks. Squaring the ratio plateau/naive gives g directly: (0.147/0.095)² = 2.4.',
    },
    {
      title: 'A test with a known answer',
      why: 'To trust the method, apply it where the exact answer is known: x(t+1) = φ x(t) + noise. Each value keeps a fraction φ of the previous one, the correlation is exactly ρ(k) = φᵏ, and g = (1 + φ)/(1 − φ).',
      work: `<div class="calc">φ = 0.9, 20 000 samples, variance 1
g = (1 + 0.9)/(1 − 0.9) = <b>19</b>
naive error  1/√20 000  = 0.0071
exact error  √(19/20 000) = <b>0.0308</b></div>
<div class="tbl"><table>${head(['block size', '1', '4', '16', '64', '128', '256', '512'])}
${row(['error', '0.0072', '0.0135', '0.0228', '0.0316', '0.0331', '0.0313', '0.0317'])}</table></div>
<div class="calc">200 independent runs: their means scatter by <b>0.0302</b>
runs whose ±1 error bar contains the true mean 0:
naive bars  <b>17 %</b>,  corrected bars  <b>72 %</b>  (should be 68 %)</div>`,
      see: 'The block curve rises from the naive 0.0072 to a plateau of about 0.031 at blocks of 64–512, exactly the exact error 0.0308 and the real scatter of 200 independent runs, 0.0302. With naive error bars, only 17 % of runs would "contain" the true value instead of 68 %: the naive error bars claim 4.4 times more precision than the data have. Here, 20 000 samples are worth about 1000 independent ones.',
    },
    {
      title: 'A real simulation: the pressure of the 9-atom box',
      why: 'Now a quantity from an actual MD run: the instantaneous pressure of our 9-atom LJ box at T = 0.72, recorded every step for 100 000 steps. We compute the naive error, the block-averaging error, and check them against 20 independent runs.',
      work: `<div class="calc">mean = 0.6030,  s = 1.1149
naive error s/√100 000 = <b>0.0035</b></div>
<div class="tbl"><table>${head(['block (steps)', 'blocks', 'error', 'uncertainty of the error'])}
${row(['1', '100 000', '0.0035', '0 %'])}
${row(['10', '10 000', '0.0108', '1 %'])}
${row(['100', '1 000', '0.0237', '2 %'])}
${row(['500', '200', '0.0338', '5 %'])}
${row(['1000', '100', '<b>0.0353</b>', '7 %'])}
${row(['2000', '50', '0.0384', '10 %'])}
${row(['10000', '10', '0.0640', '24 %'])}</table></div>
<div class="calc">plateau (blocks 500–2000) = 0.0358
g = (0.0358 / 0.0035)² = <b>103 steps</b> (0.52 time units)
independent samples in 100 000 steps: <b>969</b>

check: 20 independent runs of 20 000 steps
means scatter by 0.0786  ->  for 100 000 steps 0.0786/√5 = <b>0.0351</b></div>`,
      see: 'The honest result is P = 0.603 ± 0.036, ten times the naive ± 0.0035, and independent runs confirm it (0.0351). The system forgets its pressure only every 103 steps, so saving it every step gives no extra information beyond one sample per ~100 steps. The last rows show the limitation of block averaging: with 10 blocks the error estimate is itself uncertain by ±24 % and drifts up; read the plateau where there are still enough blocks. The uncertainty of an error from n blocks is about 1/√(2(n − 1)).',
    },
    {
      title: 'Errors of derived quantities',
      why: 'Results are often ratios or differences of measured numbers: a Haven ratio, a free-energy difference, a change of D with temperature. Independent errors add in squares: absolute errors for sums and differences, relative errors for products and ratios.',
      work: `<div class="calc">Haven ratio from Lab 6.7 (one run):
D_self = 1.003 ± 0.024 (2.4 %),  D_σ = 0.665 ± 0.171 (25.7 %)
H = 1.003 / 0.665 = 1.508
relative error √(0.024² + 0.257²) = 0.258
H = <b>1.51 ± 0.39</b>
average of 10 independent runs: ± 0.39/√10 = <b>± 0.12</b>

difference of two close numbers:
(5.0 ± 0.2) − (4.0 ± 0.2) = 1.0 ± √(0.2² + 0.2²)
                          = <b>1.0 ± 0.28</b>  (28 % from 4–5 % inputs)</div>`,
      see: 'The error of a ratio is dominated by its noisiest part: the 2.4 % of D_self hardly matters next to the 26 % of the collective term, so the Haven ratio inherits ±26 %. Differences are worse: two numbers known to 4–5 % give a difference known only to 28 %, which is why small free-energy differences or small changes between two force fields need long, careful runs. When the parts are correlated (computed from the same trajectory), these formulas are only approximate; the safest error bar then comes from repeating the whole analysis on independent runs or blocks.',
    },
    {
      title: 'Equilibration and how long to run',
      why: 'Two practical questions remain. If the run starts away from equilibrium (a box that is too large, velocities too hot), the start biases the average. And before running, one wants to know how many steps a target precision needs.',
      work: `<div class="calc">x(t+1) = 0.99 x(t) + noise (true mean 0, g = 199),
started at 30 instead of near 0. Block means of 1000 steps:
3.429  0.408  −0.123  −0.203  −0.177  0.158  −0.034 ...
the start adds 30/(1 − 0.99) = 3000 to the sum: 0.300 to the mean

mean of all 10 000 steps    <b>0.306 ± 0.141</b>  (2σ from 0)
after discarding first 1000 <b>−0.041 ± 0.149</b>

how long to run, P of the LJ box to ± 0.01:
N_eff = (s / target)² = (1.115 / 0.01)² = 12 431
steps = N_eff × g = 12 431 × 103 = <b>1.28 million</b></div>`,
      see: 'Including the equilibration period shifts the mean by 0.30, twice the error bar, a systematic error that no amount of averaging removes. Discarding it costs a little precision (0.141 → 0.149) and removes the bias. Look at block means or running averages, and discard everything before they settle. For planning, the formula N_eff = (s/target)² × g turns a target precision into a number of steps: halving the error costs four times the run. These estimates assume the run samples all relevant states; slow events (Lab 4.6 dihedrals, Lab 7.4 barriers) need far more than the formula suggests.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['method', 'merit (number)', 'limitation (number)'])}
${row(['s/√N', 'correct for independent data (1.14 ± 0.28)', 'MD: 4–10× too small; covers the truth 17 % of the time'])}
${row(['statistical inefficiency g', 'N_eff = N/g (16 -> 5.8; 100 000 -> 969)', 'noisy ρ(k); needs a cut-off rule'])}
${row(['block averaging', 'no ρ(k) needed; plateau 0.0358 = runs 0.0351', 'needs ≥ 20–50 blocks; ±24 % with 10'])}
${row(['independent runs', 'the most honest check', 'costs several equilibrations'])}
${row(['error propagation', 'ratio: relative errors in squares', 'differences: 28 % from 4–5 % inputs'])}
${row(['discard equilibration', 'removes a 2σ bias', 'slightly larger error (0.141 -> 0.149)'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the error bars by hand, the autocorrelation and g, the block averages, the AR(1) test with 200 runs, the pressure of the 9-atom box with block averaging and 20 independent runs, error propagation and the equilibration test. Run <code>pip install numpy</code>, then <code>python lab_errors.py</code> (about 15 s).',
    src: `# Lab: error bars by hand. How many independent samples does a simulation really contain?
import numpy as np, time
rng = np.random.default_rng(0)

# ---------- Step 1: the naive error bar ----------
P = np.array([1.2, 0.4, 2.1, 1.7, -0.3, 1.1, 0.9, 2.6, 0.2, 1.5])
m, s = P.mean(), P.std(ddof=1)
print("Step 1: 10 independent pressure readings:", P)
print(f"   mean = {P.sum():.1f}/10 = {m:.3f};  deviations^2 sum = {((P - m)**2).sum():.4f};  s = sqrt({((P - m)**2).sum():.4f}/9) = {s:.4f}")
print(f"   standard error of the mean = s / sqrt(N) = {s:.4f}/{np.sqrt(10):.4f} = {s / np.sqrt(10):.4f}  ->  P = {m:.2f} +- {s / np.sqrt(10):.2f}")

# ---------- Step 2: correlated data, autocorrelation and statistical inefficiency ----------
X = np.array([1.0, 1.2, 1.5, 1.4, 1.1, 0.7, 0.5, 0.6, 0.9, 1.3, 1.6, 1.8, 1.5, 1.2, 0.8, 0.9])
mx, sx = X.mean(), X.std()                                     # population std for the ACF
print(f"\\nStep 2: 16 consecutive readings (each close to the one before): mean = {mx:.4f}, s = {X.std(ddof=1):.4f}")
dx = X - mx; rho = [1.0]
for k in range(1, 5):
    c = np.mean(dx[:-k] * dx[k:]) / np.mean(dx**2); rho.append(c)
    print(f"   lag {k}: rho = <dx(t) dx(t+{k})> / <dx^2> = {np.mean(dx[:-k] * dx[k:]):.4f} / {np.mean(dx**2):.4f} = {c:.3f}")
g = 1 + 2 * sum(r for r in rho[1:] if r > 0)
print(f"   statistical inefficiency g = 1 + 2 (sum of positive rho) = 1 + 2 x ({' + '.join(f'{r:.3f}' for r in rho[1:] if r > 0)}) = {g:.2f}")
print(f"   effective independent samples N_eff = N / g = 16 / {g:.2f} = {16 / g:.1f}")
print(f"   naive error s/sqrt(16) = {X.std(ddof=1) / 4:.4f};  corrected s sqrt(g/N) = {X.std(ddof=1) * np.sqrt(g / 16):.4f}  ({np.sqrt(g):.2f} x larger)")

# ---------- Step 3: block averaging by hand ----------
print("\\nStep 3: block averaging of the same 16 numbers")
for b in (1, 2, 4, 8):
    nb = 16 // b; means = X.reshape(nb, b).mean(1); err = means.std(ddof=1) / np.sqrt(nb)
    print(f"   block size {b}: {nb:2d} block means {np.round(means, 3)}  ->  error = std/sqrt(n_blocks) = {err:.4f}")

# ---------- Step 4: a process with a known answer (AR(1)) ----------
def ar1(n, phi, rng):
    e = rng.normal(size=n) * np.sqrt(1 - phi**2); x = np.empty(n); x[0] = rng.normal()
    for i in range(1, n): x[i] = phi * x[i - 1] + e[i]
    return x
def block_curve(x, sizes):
    out = []
    for b in sizes:
        nb = len(x) // b
        if nb < 4: break
        mm = x[:nb * b].reshape(nb, b).mean(1); out.append((b, mm.std(ddof=1) / np.sqrt(nb)))
    return out
phi, n = 0.9, 20000
x = ar1(n, phi, rng); gex = (1 + phi) / (1 - phi)
print(f"\\nStep 4: x(t+1) = {phi} x(t) + noise, {n} samples, variance 1. Exact: rho(k) = {phi}^k, g = (1 + {phi})/(1 - {phi}) = {gex:.0f}")
print(f"   naive error = 1/sqrt({n}) = {1 / np.sqrt(n):.4f};  exact error = sqrt(g/N) = {np.sqrt(gex / n):.4f}")
sizes = [1, 2, 4, 8, 16, 32, 64, 128, 256, 512, 1024]
print("   block size:  " + "  ".join(f"{b:6d}" for b, _ in block_curve(x, sizes)))
print("   error:       " + "  ".join(f"{e:6.4f}" for _, e in block_curve(x, sizes)))
means = np.array([ar1(n, phi, np.random.default_rng(100 + k)).mean() for k in range(200)])
print(f"   200 independent runs: the means scatter by {means.std(ddof=1):.4f}  (exact {np.sqrt(gex / n):.4f})")
cover = np.mean(np.abs(means) < 1 / np.sqrt(n)); cover2 = np.mean(np.abs(means) < np.sqrt(gex / n))
print(f"   fraction of runs whose +-1 error bar contains the true mean 0: naive {cover:.2f}, corrected {cover2:.2f} (should be 0.68)")

# ---------- Step 5: a real simulation: pressure of the 9-atom LJ box ----------
L, N, RC, dt, T0 = 3.6, 9, 1.8, 0.005, 0.72
LAT = np.array([[(i + 0.5) * 1.2, (j + 0.5) * 1.2] for j in range(3) for i in range(3)])
def fv(x):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m_ = r2 < RC * RC
    ir6 = np.where(m_, r2**-3.0, 0); ff = np.where(m_, 24 * (2 * ir6**2 - ir6) / r2, 0)
    return -(ff[:, :, None] * d).sum(1), 0.5 * (ff * np.where(m_, r2, 0)).sum()
def pressure_run(nsteps, seed, start_T=None):
    r = np.random.default_rng(seed); x = LAT.copy(); v = r.normal(size=(N, 2)) * np.sqrt(start_T or T0); f, w = fv(x)
    c = np.exp(-dt); sn = np.sqrt((1 - c * c) * T0); out = np.empty(nsteps)
    for k in range(nsteps):
        v += 0.5 * dt * f; x += dt * v; x %= L; f, w = fv(x); v += 0.5 * dt * f; v = c * v + sn * r.normal(size=v.shape)
        out[k] = ((v**2).sum() + w) / (2 * L * L)
    return out
t0 = time.time(); Pt = pressure_run(102000, 1)[2000:]
print(f"\\nStep 5: instantaneous pressure of the 9-atom box (Langevin, T = 0.72), 100 000 steps ({time.time() - t0:.0f} s)")
print(f"   mean = {Pt.mean():.4f}, s = {Pt.std(ddof=1):.4f};  naive error s/sqrt(N) = {Pt.std(ddof=1) / np.sqrt(len(Pt)):.4f}")
bc = block_curve(Pt, [1, 10, 100, 500, 1000, 2000, 5000, 10000])
for b, e in bc: print(f"   block {b:6d} steps: {len(Pt) // b:6d} blocks, error {e:.4f}  (+- {100 / np.sqrt(2 * (len(Pt) // b - 1)):.0f} % uncertainty of the error itself)")
plat = np.mean([e for b, e in bc if 500 <= b <= 2000]); gP = (plat / bc[0][1])**2
print(f"   plateau (blocks 500-2000) = {plat:.4f};  plateau / naive = {plat / bc[0][1]:.1f}  ->  g = {gP:.0f} steps ({gP * dt:.2f} time units): only {len(Pt) / gP:.0f} independent samples in 100 000")
t0 = time.time(); runs = np.array([pressure_run(22000, 10 + k)[2000:].mean() for k in range(20)])
print(f"   check: 20 independent runs of 20 000 steps ({time.time() - t0:.0f} s): means scatter by {runs.std(ddof=1):.4f}  ->  for 100 000 steps {runs.std(ddof=1) / np.sqrt(5):.4f}  (plateau {plat:.4f}, naive {bc[0][1]:.4f})")

# ---------- Step 6: errors of derived quantities ----------
print("\\nStep 6: propagating errors")
Ds, dDs, Dq, dDq = 1.003, 0.024, 0.665, 0.171
H = Ds / Dq; dH = H * np.sqrt((dDs / Ds)**2 + (dDq / Dq)**2)
print(f"   Haven = D_self / D_sigma = {Ds}/{Dq} = {H:.3f};  relative errors {100 * dDs / Ds:.1f} % and {100 * dDq / Dq:.1f} %")
print(f"   relative error of a ratio: sqrt({dDs / Ds:.3f}^2 + {dDq / Dq:.3f}^2) = {np.sqrt((dDs / Ds)**2 + (dDq / Dq)**2):.3f}  ->  Haven = {H:.2f} +- {dH:.2f}")
print(f"   with 10 independent runs averaged: error / sqrt(10) -> +- {dH / np.sqrt(10):.2f}")
a, da, b_, db = 5.0, 0.2, 4.0, 0.2
print(f"   a difference of two close numbers: ({a} +- {da}) - ({b_} +- {db}) = {a - b_} +- {np.sqrt(da**2 + db**2):.2f}  ({100 * np.sqrt(da**2 + db**2) / (a - b_):.0f} % error from 4-5 % inputs)")

# ---------- Step 7: equilibration and how long to run ----------
y = ar1(10000, 0.99, np.random.default_rng(3)); y0 = 30.0
for i in range(10000): y[i] += y0 * 0.99**i                      # same noise, but started far from equilibrium (like a box that is too big)
print("\\nStep 7: x(t+1) = 0.99 x(t) + noise (true mean 0, g = 199), started at x = 30 instead of near 0. Block means of 1000 steps:")
print("   " + " ".join(f"{q:.3f}" for q in y.reshape(10, 1000).mean(1)))
ee = np.sqrt(199 / 10000); ee9 = np.sqrt(199 / 9000)
print(f"   start term adds sum 30 x 0.99^k = 30 / 0.01 = 3000 to the sum, i.e. {3000 / 10000:.3f} to the mean")
print(f"   mean of all 10 000: {y.mean():.3f} +- {ee:.3f};  after discarding the first 1000: {y[1000:].mean():.3f} +- {ee9:.3f}")
target = 0.01; Neff = (Pt.std(ddof=1) / target)**2
print(f"   how long to run: P of the LJ box to +-{target}: N_eff = (s/target)^2 = ({Pt.std(ddof=1):.3f}/{target})^2 = {Neff:,.0f} independent samples = {Neff * gP:,.0f} steps")
`,
  },

  exercises: [
    { q: 'A run gives D = 2.40×10⁻⁹ m²/s with a naive error of 0.01×10⁻⁹, but block averaging levels off at 0.08×10⁻⁹. What is g, and how many independent samples did the 50 000 frames contain?', a: 'g = (0.08/0.01)² = 64. N_eff = 50 000/64 ≈ 780. Report D = (2.40 ± 0.08)×10⁻⁹ m²/s.' },
    { q: 'For x(t+1) = φ x(t) + noise with φ = 0.5, what is g? How many samples do you need for an error of 0.01 if the variance is 1?', a: 'g = 1.5/0.5 = 3. N_eff = (1/0.01)² = 10 000, so N = 3 × 10 000 = 30 000 samples.' },
    { q: 'Two force fields give densities 1.052 ± 0.004 and 1.046 ± 0.003 g/cm³. Is the difference significant?', a: 'Difference 0.006 ± √(0.004² + 0.003²) = 0.006 ± 0.005. It is only 1.2 σ: not significant. Longer runs (or more independent runs) are needed before claiming one force field gives a higher density.' },
  ],
};
