import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-lifetimes',
  lab: true,
  summary: 'How long does an ion pair last? How long does a solvent molecule stay in a shell? These "lifetimes" are among the most quoted numbers in simulation papers, and among the most definition-dependent. We compute continuous and intermittent correlation functions by hand on a 12-frame on/off sequence, then measure first-shell lifetimes in the LJ liquid and change the cutoff, the saving interval and the tolerance for brief excursions, to see how much each choice moves the answer.',
  minutes: 30,
  sim: 'lab12',
  learn: [
    'Compute continuous and intermittent correlation functions by hand, and read a lifetime from their area',
    'See how much the cutoff, the saving interval and a tolerance t* change a lifetime, with numbers',
    'Compare lifetimes with diffusion times, and know what to report with any lifetime',
  ],
  simHint: 'A model pair rattles in a well and now and then escapes. Move the cutoff, the saving interval and the tolerance t*, and watch the continuous lifetime (cyan) change a lot while the intermittent one (amber) changes much less. The strip at the top shows one pair: bonded frames are lit.',

  setup: `
<p>For each pair (an ion pair, an ion and a solvent molecule, two neighbours) define h(t) = 1 if the pair is "bonded" (closer than a cutoff) at frame t, and 0 otherwise. Two correlation functions are used:</p>
<div class="calc">intermittent  C_I(t) = ⟨h(0) h(t)⟩ / ⟨h⟩
              still (or again) bonded after t; breaks allowed
continuous    C_C(t) = ⟨h(0) H(t)⟩ / ⟨h⟩
              H(t) = 1 only if bonded at every frame from 0 to t

lifetime τ = area under C(t)   (= decay time if C is exponential)</div>
<p>Data: a 12-frame sequence by hand, and the 3D LJ liquid of Labs 6.9–6.10 (108 atoms, ρ = 0.8, T = 1.0), with "neighbours" defined by the first minimum of g(r) at r = 1.45.</p>`,

  steps: [
    {
      title: 'Continuous and intermittent by hand',
      why: 'The two functions answer different questions. The continuous one asks "has this pair stayed together without interruption?", the intermittent one asks "is it together again (or still) after time t?". On a short sequence we can count every case.',
      work: `<div class="calc">one pair, frames every 1 ps:   1 1 1 0 1 1 1 1 0 0 0 1
⟨h⟩ = 8/12 = 0.667;  unbroken stretches: 3, 4, 1

lag   intermittent                  continuous
0     8/8 = 1.000                   8/8 = 1.000
1     5/7 = 0.714                   5/8 = 0.625
2     4/7 = 0.571                   3/8 = 0.375
3     3/7 = 0.429                   1/8 = 0.125
4     4/7 = 0.571                   0/8 = 0.000

intermittent: bonded at t and at t + lag / bonded at t
continuous:   origins that stay bonded lag frames
              without a break / all bonded origins
              (a stretch of length L gives L − lag origins)

lifetime = area (lags 0–4, first term halved, × 1 ps):
continuous <b>1.62 ps</b>        intermittent <b>2.79 ps</b>
check: continuous area = Σ L² / (2 Σ L) = 26/16 = 1.62 ps</div>`,
      see: 'The continuous function dies as soon as each stretch ends; the intermittent one even rises again at lag 4, because the pair re-forms. The intermittent lifetime is therefore always the longer one (2.79 vs 1.62 ps here). Neither is "right": the continuous lifetime describes how long an uninterrupted contact lasts, the intermittent one how long until the partners have really separated and forgotten each other. The formula Σ L²/(2ΣL) shows that the continuous area depends only on the stretch lengths, which is why anything that cuts a stretch short (next steps) changes it strongly.',
    },
    {
      title: 'Forgiving short breaks: the tolerance t*',
      why: 'Thermal motion makes a bonded pair cross the cutoff briefly without really separating. Many analyses therefore forgive breaks shorter than a tolerance t*. It sounds harmless; let us count what it does.',
      work: `<div class="calc">same sequence, breaks of ≤ 1 frame forgiven (t* = 1 ps):
1 1 1 1 1 1 1 1 0 0 0 1
stretches: 8, 1
continuous lifetime = Σ L²/(2 Σ L) = 65/18 = <b>3.61 ps</b>
(was 1.62 ps:  2.2 times longer from one forgiven frame)</div>`,
      see: 'Joining two stretches of 3 and 4 into one of 8 more than doubles the continuous lifetime, because the area grows with L². A tolerance is reasonable physics (a brief excursion is not a real break), but its value must be stated: with different t* two papers can report lifetimes that differ by a factor of 2 for the same simulation. The intermittent function already ignores short breaks, so it hardly changes.',
    },
    {
      title: 'First-shell neighbours in the LJ liquid',
      why: 'Now the real case: every pair of the 108-atom liquid, 2000 frames saved every 0.02 time units, neighbours defined by the first minimum of g(r) (r < 1.45). We compute both functions and remove the long-time baseline of the intermittent one: two random atoms are neighbours a fraction ⟨h⟩ of the time anyway.',
      work: `<div class="calc">⟨h⟩ = 0.0991  (10.60 neighbours per atom)</div>
<div class="tbl"><table>${head(['t', '0.02', '0.1', '0.2', '0.5', '1.0', '2.0', '5.0'])}
${row(['continuous', '0.984', '0.923', '0.856', '0.693', '0.496', '0.257', '0.035'])}
${row(['intermittent', '0.984', '0.928', '0.878', '0.784', '0.683', '0.544', '0.316'])}
${row(['intermittent − baseline, normalized', '0.983', '0.920', '0.864', '0.760', '0.648', '0.493', '0.241'])}</table></div>
<div class="calc">baseline correction: (C_I − ⟨h⟩) / (1 − ⟨h⟩)
lifetimes:  continuous <b>1.46</b>,  intermittent <b>3.25</b> time units
(intermittent integrated until it falls below 0.01, at t = 12)</div>`,
      see: 'Within the first 0.2 time units the two functions agree; after that the continuous one falls much faster, because neighbours that step out briefly and come back are lost to it. The intermittent function decays to the baseline ⟨h⟩ = 0.099, not to zero: without subtracting it the lifetime would come out too long (and grow with the integration limit). Limitation of the intermittent lifetime: its tail is long and noisy (here integrated to t = 12), so it needs long trajectories, exactly like the collective quantities of Lab 6.8.',
    },
    {
      title: 'How much the definition moves the answer',
      why: 'Three choices enter every lifetime: the cutoff, how often frames are saved, and the tolerance t*. We change each one on the same trajectory.',
      work: `<div class="tbl"><table>${head(['cutoff', 'neighbours', 'continuous', 'intermittent'])}
${row(['1.30', '8.12', '0.723', '2.660'])}
${row(['1.45 (g(r) minimum)', '10.60', '1.460', '3.251'])}
${row(['1.60', '12.87', '2.172', '3.519'])}</table></div>
<div class="tbl"><table>${head(['saving interval', 'continuous lifetime', 'mean unbroken stretch'])}
${row(['0.02', '1.460', '1.241'])}
${row(['0.10', '1.526', '1.356'])}
${row(['0.50', '<b>2.064</b>', '2.210'])}</table></div>
<div class="tbl"><table>${head(['tolerance t*', '0', '0.04', '0.10', '0.20'])}
${row(['continuous lifetime', '1.460', '1.509', '1.684', '2.025'])}</table></div>`,
      see: 'The continuous lifetime triples between cutoffs of 1.30 and 1.60 (0.72 → 2.17); the intermittent one changes by 30 % (2.66 → 3.52). Saving frames less often hides short excursions and makes the continuous lifetime longer (1.46 → 2.06 at 25 times larger intervals); a tolerance does the same on purpose (2.03 with t* = 0.2). So the continuous lifetime is mostly a property of the definition, the intermittent one mostly a property of the liquid. Practical rule: put the cutoff at the minimum of g(r), prefer the intermittent function (with the baseline removed) for "how long do partners stay together", and always state cutoff, saving interval and t*.',
    },
    {
      title: 'Compare with diffusion',
      why: 'A lifetime is only "long" or "short" compared with something. The natural yardstick is how long the atoms take to diffuse across the shell or across one diameter, using D from Lab 6.10.',
      work: `<div class="calc">D = 0.0688 (Lab 6.10),  t = Δr² / (6D)
shell width 0.35:   0.1225 / 0.4128 = <b>0.30</b>
one diameter 1.0:   1.0000 / 0.4128 = <b>2.42</b>

intermittent lifetime 3.25  ≈  time to diffuse one diameter
continuous lifetime   1.46  ≈  5 times the shell-crossing time</div>`,
      see: 'In this simple liquid, partners stay together about as long as it takes to diffuse one diameter: nothing special binds them, they just drift apart. In an ionic liquid or electrolyte the interesting question is whether an ion pair lives much longer than this diffusion time (a real, long-lived pair that moves as one neutral unit and lowers the conductivity, Lab 6.7) or not (a transient contact). Comparing τ with σ²/(6D) turns a bare lifetime into that physical statement. Limitation: lifetimes from a single cutoff ignore the shape of the free-energy barrier between bonded and free states; for barrier-controlled processes the PMF of Lab 7.5 and a rate theory give a better-defined number.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['choice', 'effect (number)', 'advice'])}
${row(['continuous vs intermittent', '1.46 vs 3.25 (same data)', 'say which; intermittent for "partners"'])}
${row(['cutoff 1.30 -> 1.60', 'continuous ×3, intermittent ×1.3', 'cutoff at the g(r) minimum'])}
${row(['saving interval ×25', 'continuous 1.46 -> 2.06', 'save often; report the interval'])}
${row(['tolerance t* = 0.2', 'continuous 1.46 -> 2.03', 'report t*'])}
${row(['intermittent baseline', 'decays to ⟨h⟩ = 0.099, not 0', 'subtract and normalize'])}
${row(['yardstick', 'σ²/(6D) = 2.42 ≈ intermittent 3.25', 'compare τ with diffusion times'])}</table></div>`,

  code: {
    note: 'This script prints every number on this page: the hand counts for the 12-frame sequence, the tolerance example, the LJ liquid correlation functions and lifetimes, and the cutoff, saving-interval and tolerance tables. Run <code>pip install numpy</code>, then <code>python lab_lifetimes.py</code> (about 30 s).',
    src: `# Lab: lifetimes and residence times by hand. Continuous vs intermittent correlation functions, tolerance, cutoff and
# sampling interval, measured for first-shell neighbours in the 3D LJ liquid of Labs 6.9-6.10.
import numpy as np, time

# ---------- Steps 1-3: one pair by hand ----------
h = np.array([1, 1, 1, 0, 1, 1, 1, 1, 0, 0, 0, 1])
print("Step 1: one pair, 1 = inside the cutoff, frames every 1 ps:", "".join(map(str, h)))
print(f"   <h> = {h.sum()}/{len(h)} = {h.mean():.4f}")
def intermittent(h, k): return np.mean(h[:len(h) - k] * h[k:]) / np.mean(h[:len(h) - k]) if k else 1.0
def runs(h):
    out, n = [], 0
    for x in list(h) + [0]:
        if x: n += 1
        elif n: out.append(n); n = 0
    return out
def continuous(h, k):                                # fraction of bonded origins that stay bonded without a break for k frames
    rl = runs(h); return sum(max(0, L - k) for L in rl) / sum(rl)
rl = runs(h); print(f"   unbroken stretches: {rl}")
for k in range(5):
    num = int(np.sum(h[:len(h) - k] * h[k:])); den = int(h[:len(h) - k].sum()); cnum = sum(max(0, L - k) for L in rl)
    print(f"   lag {k}: intermittent = pairs bonded at t and t+{k} / bonded at t = {num}/{den} = {intermittent(h, k):.3f};   continuous = {cnum}/{sum(rl)} = {continuous(h, k):.3f}")
ti = 0.5 + sum(intermittent(h, k) for k in range(1, 5)); tc = 0.5 + sum(continuous(h, k) for k in range(1, 5))
print(f"   lifetimes (sum, lags 0-4, x 1 ps, first term halved): continuous {tc:.2f} ps, intermittent {ti:.2f} ps")
print(f"   check: area under the continuous C = sum L^2 / (2 sum L) = {sum(L * L for L in rl)} / {2 * sum(rl)} = {sum(L * L for L in rl) / (2 * sum(rl)):.2f} ps")
def tolerate(h, tstar):                              # fill breaks of <= tstar frames that are bonded on both sides
    g = h.copy(); i = 0
    while i < len(g):
        if g[i] == 0:
            j = i
            while j < len(g) and g[j] == 0: j += 1
            if 0 < i and j < len(g) and j - i <= tstar: g[i:j] = 1
            i = j
        else: i += 1
    return g
g1 = tolerate(h, 1)
r1 = runs(g1); print(f"\\nStep 3: forgive breaks of 1 frame (t* = 1 ps): {''.join(map(str, g1))}  ->  stretches {r1},  continuous lifetime sum L^2/(2 sum L) = {sum(L * L for L in r1)}/{2 * sum(r1)} = {sum(L * L for L in r1) / (2 * sum(r1)):.2f} ps")

# ---------- Steps 4-6: the LJ liquid ----------
RC, SK, rho, T0, dt = 1.6, 0.3, 0.8, 1.0, 0.004
def fcc(nc, a):
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]]); g = np.array([[i, j, k] for i in range(nc) for j in range(nc) for k in range(nc)])
    return ((g[:, None, :] + b[None]) * a).reshape(-1, 3)
def pairs(x, L):
    d = x[None] - x[:, None]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); return np.nonzero(np.triu(r2 < (RC + SK)**2, 1))
def forces(x, L, i, j, N):
    d = x[j] - x[i]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); m = r2 < RC * RC
    ir6 = np.where(m, 1 / r2**3, 0); ff = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0)[:, None] * d; f = np.zeros((N, 3))
    for k in range(3): f[:, k] = np.bincount(j, ff[:, k], N) - np.bincount(i, ff[:, k], N)
    return f
def trajectory(nc=3, nsteps=10000, every=5, seed=1, tau=1.0):
    N = 4 * nc**3; L = (N / rho)**(1 / 3); x = fcc(nc, L / nc); r = np.random.default_rng(seed); NF = 3 * N - 3
    v = r.normal(size=(N, 3)); v -= v.mean(0); v *= np.sqrt(NF * T0 / (v**2).sum()); i, j = pairs(x, L); xl = x.copy(); f = forces(x, L, i, j, N)
    c = np.exp(-dt / tau); out = []
    for s in range(-3000, nsteps):
        v += 0.5 * dt * f; x += dt * v
        if np.abs((x - xl + L / 2) % L - L / 2).max() > SK / 2: x %= L; i, j = pairs(x, L); xl = x.copy()
        f = forces(x, L, i, j, N); v += 0.5 * dt * f
        K = 0.5 * (v**2).sum(); Kt = 0.5 * NF * T0; R1 = r.normal(); S2 = r.chisquare(NF - 1)
        v *= np.sqrt(max(c + (1 - c) * Kt / (NF * K) * (R1**2 + S2) + 2 * R1 * np.sqrt(c * (1 - c) * Kt / (NF * K)), 0))
        if s >= 0 and s % every == 0: out.append((x % L).copy())
    return np.array(out), L
t0 = time.time(); X, L = trajectory(); N = X.shape[1]; iu = np.triu_indices(N, 1)
D = X[:, None, :, :] - X[:, :, None, :]; D -= L * np.round(D / L); R = np.sqrt((D**2).sum(-1))[:, iu[0], iu[1]]; del D
fdt = 5 * dt
print(f"\\nStep 4: 3D LJ liquid (108 atoms, rho = 0.8, T = 1.0), {len(X)} frames every {fdt} time units ({time.time() - t0:.0f} s); first-shell cutoff at the g(r) minimum")
def runlengths(Hm, val=True):                        # lengths of unbroken stretches of \`val\` in every column (vectorised)
    B = (Hm == val).astype(np.int8); P = np.zeros((B.shape[0] + 2, B.shape[1]), np.int8); P[1:-1] = B
    dif = np.diff(P, axis=0); st = np.argwhere(dif.T == 1); en = np.argwhere(dif.T == -1)
    return en[:, 1] - st[:, 1], st, en                  # lengths, (column, start), (column, end)
def C_all(Hm, kmax):
    Hf = Hm.astype(float); avg = Hf.mean(); n = Hf.shape[0]
    F = np.fft.rfft(Hf, 2 * n, axis=0); num = np.fft.irfft(np.abs(F)**2, axis=0)[:kmax].sum(1)      # sum over pairs and origins of h(t)h(t+k)
    cs = np.concatenate([[0], np.cumsum(Hf.sum(1))]); den = np.array([cs[n - k] for k in range(kmax)])  # bonded origins usable at lag k
    ci = num / den
    rl = runlengths(Hm)[0]; cc = np.array([np.maximum(rl - k, 0).sum() / rl.sum() for k in range(kmax)])
    return ci, cc, avg, rl
def tolerate_all(Hm, ts):                              # fill breaks of <= ts frames that have bonds on both sides
    G = Hm.copy(); L_, st, en = runlengths(Hm, False); nF = Hm.shape[0]
    for (col, a), (_, b), n in zip(st, en, L_):
        if n <= ts and a > 0 and b < nF: G[a:b, col] = True
    return G
for rc in (1.45,):
    Hm = R < rc; ci, cc, avg, rl = C_all(Hm, 600)
    print(f"   cutoff {rc}: <h> = {avg:.4f} (= {avg * (N - 1):.2f} neighbours per atom)")
    for tt in (0.02, 0.1, 0.2, 0.5, 1.0, 2.0, 5.0):
        k = int(round(tt / fdt)); print(f"     t = {tt:4.2f}: continuous {cc[k]:.3f}   intermittent {ci[k]:.3f}   intermittent minus baseline, normalised {(ci[k] - avg) / (1 - avg):.3f}")
    tc = fdt * (cc.sum() - 0.5); 
    cin = (ci - avg) / (1 - avg); kz = np.argmax(cin < 0.01) if np.any(cin < 0.01) else len(cin); ti = fdt * (cin[:kz].sum() - 0.5)
    print(f"   lifetimes: continuous {tc:.3f}, intermittent (baseline removed, integrated to where it falls below 0.01, t = {kz * fdt:.1f}) {ti:.3f} time units")
print("\\nStep 5: dependence on the definition")
print("   cutoff   neighbours   continuous   intermittent (time units)")
for rc in (1.30, 1.45, 1.60):
    ci, cc, avg, rl = C_all(R < rc, 600); cin = (ci - avg) / (1 - avg); kz = np.argmax(cin < 0.01) if np.any(cin < 0.01) else len(cin)
    print(f"   {rc:5.2f}    {avg * (N - 1):6.2f}       {fdt * (cc.sum() - 0.5):6.3f}       {fdt * (cin[:kz].sum() - 0.5):6.3f}")
print("   saving interval (cutoff 1.45):")
for sub in (1, 5, 25):
    ci, cc, avg, rl = C_all((R < 1.45)[::sub], 600 // sub); dtt = fdt * sub
    print(f"   every {dtt:5.2f}: continuous lifetime {dtt * (cc.sum() - 0.5):.3f}   (mean unbroken stretch {dtt * rl.mean():.3f})")
print("   tolerance t* (cutoff 1.45, frames every 0.02):")
H45 = (R < 1.45)
for tstar_t in (0.0, 0.04, 0.1, 0.2):
    ts = int(round(tstar_t / fdt)); Ht = tolerate_all(H45, ts) if ts else H45
    ci, cc, avg, rl = C_all(Ht, 600); print(f"   t* = {tstar_t:4.2f}: continuous lifetime {fdt * (cc.sum() - 0.5):.3f}")

# ---------- Step 6: compare with diffusion ----------
Dinf = 0.0688
print(f"\\nStep 6: time to diffuse the width of the first shell and one diameter (D = {Dinf} from Lab 6.10)")
for dr in (0.35, 1.0):
    print(f"   distance {dr}: t = dr^2 / (6 D) = {dr**2:.4f} / {6 * Dinf:.4f} = {dr**2 / (6 * Dinf):.3f}")
`,
  },

  exercises: [
    { q: 'A pair has the sequence 1 1 0 1 1 1 1 0 1 1 (frames 1 ps apart). Compute the continuous lifetime with t* = 0 and with t* = 1 ps.', a: 't* = 0: stretches 2, 4, 2; Σ L² = 4 + 16 + 4 = 24, Σ L = 8, lifetime 24/16 = 1.5 ps. t* = 1: all breaks are single frames between bonded frames, so one stretch of 10; lifetime 100/20 = 5.0 ps.' },
    { q: 'An intermittent function decays to 0.12 at long times instead of 0, and the average number of partners is 0.12 of all possible pairs. What is the normalized value at a time where C_I = 0.40?', a: '(0.40 − 0.12)/(1 − 0.12) = 0.28/0.88 = 0.318.' },
    { q: 'Li⁺–anion pairs have an intermittent lifetime of 2 ns. The ions diffuse with D ≈ 2×10⁻¹¹ m²/s and the pair distance is about 0.5 nm. Are these long-lived pairs?', a: 'Time to diffuse 0.5 nm: (0.5×10⁻⁹)²/(6 × 2×10⁻¹¹) = 2.5×10⁻¹⁹/1.2×10⁻¹⁰ = 2.1×10⁻⁹ s ≈ 2 ns. The lifetime equals the diffusion time, so the pairs are not especially long-lived: partners separate about as fast as diffusion alone would separate them.' },
  ],
};
