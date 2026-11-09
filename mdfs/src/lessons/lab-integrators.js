import { tex } from '../tex.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'lab-integrators',
  lab: true,
  summary: 'Velocity Verlet keeps the energy steady for billions of steps, while “better” fourth-order methods slowly leak it away. The reason is a hidden conserved quantity, the shadow Hamiltonian. We compute it by hand, check it in an LJ cluster, see how single precision breaks energy conservation (and why it depends on where the atoms are), and speed up a stiff system with multiple time steps (RESPA), including the resonance that limits it.',
  minutes: 30,
  sim: 'lab17',
  learn: [
    'Compute the shadow Hamiltonian that Verlet conserves, and see why Verlet does not drift but RK4 does',
    'See how float32 rounding causes energy drift, and why it grows with the distance from the origin',
    'Split fast and slow forces with RESPA, count the savings, and find the resonance that limits the outer step',
  ],
  simHint: 'Shadow: Verlet points on a harmonic oscillator stay on the dashed shadow ellipse forever; switch to RK4 and the points spiral inward; push h past 2 and Verlet explodes. Precision: the same 7-atom cluster in float32 and float64; move it away from the origin and the float32 energy starts to wander. RESPA: drag the outer step H near H × omega / pi = 1 and the energy blows up, even though every inner step is tiny.',

  setup: `
<p>Velocity Verlet (Lab 1.5), one step of size h, mass 1:</p>
<div class="calc">p ← p + (h/2) F(x)
x ← x + h p
p ← p + (h/2) F(x)</div>
<p>Systems: a harmonic oscillator (ω = 1), a 13-atom LJ cluster (the icosahedron, no cutoff, k_BT ≈ 0.25), a fast-plus-slow spring, and a liquid of 32 stiff dimers. All in LJ units.</p>`,

  steps: [
    {
      title: 'The shadow Hamiltonian, by hand',
      why: 'Verlet does not conserve the true energy E exactly: it wobbles by about h². But for a harmonic oscillator it conserves a slightly different quantity exactly, forever. That hidden quantity is why the energy never drifts away.',
      work: `<div class="calc">E  = p²/2 + x²/2
H̃ = p²/2 + (x²/2)(1 − h²/4)     (conserved by Verlet)

h = 0.5,  start x = 1, p = 0
step  x          p          E          H̃
0     +1.000000  +0.000000  0.500000   0.468750
1     +0.875000  −0.468750  0.492676   0.468750
2     +0.531250  −0.820312  0.477570   0.468750
3     +0.054688  −0.966797  0.468843   0.468750
4     −0.435547  −0.871582  0.474678   0.468750

step 1 in detail:
  p = 0 − 0.25 × 1          = −0.25
  x = 1 + 0.5 × (−0.25)     = 0.875
  p = −0.25 − 0.25 × 0.875  = −0.46875
  H̃ = 0.46875²/2 + (0.875²/2)(0.9375) = 0.468750</div>
<div class="calc">shadow frequency (2/h) arcsin(h/2) = 1.0107: 1.07% fast
stable only while h ω &lt; 2
  h = 2.1: 1 − h²/4 = −0.10 &lt; 0, no closed orbit
  after 20 steps |x| = 1.5×10⁵</div>`,
      see: 'E moves between 0.469 and 0.5 but H̃ stays at 0.468750 to every digit. Verlet follows the exact orbit of a slightly different system (a spring 1 − h²/4 times as stiff). Its energy error is bounded, and its only lasting error is a small phase shift (1% at h = 0.5). Above h ω = 2 no shadow orbit exists and the motion blows up: the stability limit of Lesson 1.4.',
    },
    {
      title: 'The shadow energy in an LJ cluster, and why RK4 drifts',
      why: 'For any smooth force, Verlet conserves a shadow energy that differs from E by terms of order h². Adding the first correction should remove most of the wobble, and what is left should shrink as h⁴. A non-symplectic method has no shadow energy at all, however accurate it is per step.',
      work: `<div class="calc">H̃ = E + h² ( p·∇²U·p / 12 − |F|²/24 )
(reduces to the exact H̃ of Step 1 for a spring, to order h²)

13 LJ atoms, icosahedron (U = −44.3268), k_BT = 0.25, t = 50</div>
<div class="tbl"><table>${head(['h', 'std of E', 'std of H̃', 'ratio'])}
${row(['0.008', '1.16×10⁻³', '7.5×10⁻⁶', '154'])}
${row(['0.004', '2.99×10⁻⁴', '5.0×10⁻⁷', '594'])}
${row(['0.002', '7.70×10⁻⁵', '3.2×10⁻⁸', '2429'])}</table></div>
<div class="calc">halving h:  E wobble ÷ 3.88   (h²: ÷ 4)
            H̃ wobble ÷ 15–16 (h⁴: ÷ 16)</div>
<p>Verlet vs RK4 (fourth order) on the oscillator, h = 0.2, 100 000 steps:</p>
<div class="tbl"><table>${head(['method', 'E after 10 steps', 'E range over run', 'E at the end'])}
${row(['Verlet', '0.49675', '0.49500 – 0.50000', '0.49700'])}
${row(['RK4', '0.4999951', '0.45768 – 0.50000', '0.45768'])}</table></div>
<div class="calc">RK4 energy factor per step: 1 − h⁶/72 = 0.999999111
after 100 000 steps: 0.999999111^100000 = 0.915</div>`,
      see: 'One correction term makes the conserved quantity 150–2400 times steadier, with the h⁴ scaling expected. RK4 is a thousand times more accurate per step than Verlet, yet it loses 8.5% of the energy over the run, steadily and without limit; Verlet stays within 1%. Over long MD runs, having a conserved shadow energy matters more than the order of the method.',
    },
    {
      title: 'Single vs double precision',
      why: 'GPUs are much faster in single precision (float32, about 7 digits). Every position update is then rounded. If the rounding errors were truly random they would cancel, but a coordinate far from the origin has fewer digits left for the small step it takes each time.',
      work: `<div class="calc">float32 keeps 24 binary digits: spacing = 2⁻²³ × 2^(exponent)
  at 1:     1.19×10⁻⁷       at 100:   7.63×10⁻⁶
  at 10:    9.54×10⁻⁷       at 1000:  6.10×10⁻⁵
float64 at 1000: 1.1×10⁻¹³

step size: v h ≈ 0.45 × 0.004 = 0.0018 per step
relative rounding of that step:
  near 0:     1×10⁻⁷ / 0.0018 = 7×10⁻⁵
  at 1000:   6×10⁻⁵ / 0.0018 = 3%</div>
<div class="tbl"><table>${head(['precision', 'cluster at', 'drift of E over t = 200', 'rms of E − E₀'])}
${row(['float64', '0', '−1.9×10⁻⁶', '1.21×10⁻³'])}
${row(['float32', '0', '−3.7×10⁻⁶', '1.20×10⁻³'])}
${row(['float32', '10', '−2.8×10⁻³', '0.97×10⁻³'])}
${row(['float32', '100', '−1.4×10⁻²', '6.1×10⁻³'])}
${row(['float32', '1000', '+5.3×10⁻²', '5.7×10⁻²'])}</table></div>`,
      see: 'Near the origin, float32 is as good as float64: the wobble is the h² shadow fluctuation, not rounding. At distance 1000 the drift is 27 000 times larger. That is why GROMACS mixed precision stores coordinates relative to the box, keeps the box small (nanometres), and accumulates some sums in double precision, and why long NVE runs for transport properties are often checked in double precision.',
    },
    {
      title: 'RESPA by hand: fast and slow forces on different clocks',
      why: 'The time step is set by the fastest motion (bond vibrations), but the most expensive forces (long-range electrostatics, non-bonded interactions) change slowly. RESPA (reversible reference system propagator) applies the slow force once per long step H and the fast force every short step h = H/n. It is still symplectic and time reversible.',
      work: `<div class="calc">outer step H:
  p ← p + (H/2) F_slow
  n inner Verlet steps of h = H/n with F_fast
  p ← p + (H/2) F_slow

x'' = −(ω_f² + ω_s²) x,  ω_f = 10, ω_s = 1, ω = 10.05
H = 0.1, n = 2 (h = 0.05), start x = 1, p = 0, E = 50.5
slow half kick:   p = 0 − 0.05 × 1 × 1        = −0.05000
inner 1: kick p = −0.05 − 0.025 × 100 × 1      = −2.55000
         drift x = 1 + 0.05 × (−2.55)          = +0.87250
         kick p = −2.55 − 0.025 × 100 × 0.8725 = −4.73125
inner 2: kick −6.91250, drift +0.52687, kick −8.22969
slow half kick:   p = −8.25603,  E = 48.10
(slow force: 1 evaluation; fast force: 2)</div>
<p>Resonance scan, inner step ≈ 0.005, t = 400:</p>
<div class="tbl"><table>${head(['H', 'H ω / π', 'max |E/E₀ − 1|'])}
${row(['0.05', '0.16', '1.5×10⁻³'])}
${row(['0.20', '0.64', '2.0×10⁻²'])}
${row(['0.30', '0.96', '0.22'])}
${row(['0.3126', '1.00', '<b>10⁶ (blows up)</b>'])}
${row(['0.33', '1.06', '0.19'])}
${row(['0.40', '1.28', '2.3×10⁻²'])}
${row(['0.625', '2.00', '<b>10⁶ (blows up)</b>'])}</table></div>`,
      see: 'The inner steps are tiny, yet the run explodes when the outer step equals half a period of the fast motion (H = π/ω) or a multiple of it: the slow kicks then always hit the vibration at the same phase and pump energy into it. In practice the outer step must stay well below half the fastest period, about H &lt; 0.3 × period.',
    },
    {
      title: 'RESPA in a liquid of stiff molecules',
      why: 'Now a real case: 32 dimers with stiff bonds (k = 2000) in a 3D LJ liquid. The bonds are cheap and fast; the LJ forces are expensive and slower. Counting LJ evaluations per unit time measures the cost.',
      work: `<div class="calc">bond: ω = √(k/μ) = √(2000/0.5) = 63.2,  period 0.0993
Verlet limit h &lt; 2/ω = 0.032;  resonance H = π/ω = 0.050
64 atoms, ρ = 0.6, T = 1.5, NVE for t = 10</div>
<div class="tbl"><table>${head(['method', 'LJ evals / time', 'bond evals / time', 'std of E/N', 'drift E/N'])}
${row(['Verlet h = 0.001', '1000', '1000', '8.0×10⁻⁵', '+2×10⁻⁶'])}
${row(['Verlet h = 0.004', '250', '250', '1.3×10⁻³', '+3.8×10⁻⁴'])}
${row(['Verlet h = 0.008', '125', '125', '5.9×10⁻³', '−1.1×10⁻³'])}
${row(['RESPA H = 0.004, n = 4', '250', '1000', '3.6×10⁻⁴', '+1.2×10⁻⁴'])}
${row(['RESPA H = 0.008, n = 8', '125', '1000', '1.4×10⁻³', '−6.8×10⁻⁴'])}
${row(['RESPA H = 0.02, n = 20', '50', '1000', '5.5×10⁻²', '+0.18'])}
${row(['RESPA H = 0.04, n = 40', '25', '1000', 'blew up', ''])}</table></div>`,
      see: 'RESPA with H = 0.008 matches Verlet with h = 0.004 at half the LJ cost, and at equal LJ cost it is 4 times steadier. The gain stops near H ≈ 0.01: beyond that the LJ forces themselves are too fast for the long step (H = 0.04 blows up, below the bond resonance at 0.05). In real codes the slow part is usually the Ewald reciprocal-space force (LAMMPS <code>run_style respa</code>, OpenMM <code>MTSIntegrator</code>); constraints (SHAKE, Lesson 4.4) are the other way to remove the fastest motions.',
    },
  ],

  summaryTable: `<div class="tbl"><table>${head(['topic', 'key number', 'lesson'])}
${row(['shadow Hamiltonian', 'H̃ = 0.468750 at every step', 'Verlet follows the exact orbit of a nearby system'])}
${row(['stability', 'h ω &lt; 2', 'no shadow orbit beyond it'])}
${row(['shadow correction', 'E wobble ∝ h², H̃ wobble ∝ h⁴', 'the conserved quantity exists for real forces'])}
${row(['RK4 vs Verlet', '−8.5% vs bounded, 10⁵ steps', 'symplectic beats high order over long runs'])}
${row(['float32', 'drift ×27 000 at distance 1000', 'keep coordinates small; check in double'])}
${row(['RESPA resonance', 'blows up at H = π/ω', 'outer step well below half the fast period'])}
${row(['RESPA savings', 'same accuracy, ½ the LJ cost', 'limited by the slow force’s own time scale'])}</table></div>`,

  code: {
    note: 'The whole lab in one script: the shadow Hamiltonian table, the LJ13 shadow energy, Verlet vs RK4, float32 vs float64, RESPA by hand with the resonance scan, and the dimer liquid. Run <code>pip install numpy</code>, then <code>python lab_integrators.py</code> (about 15 seconds).',
    src: `# Lab 1.8: integrator extras -- shadow Hamiltonian, Verlet vs RK4, single vs double precision, RESPA
import numpy as np, time
t0 = time.time()

print("Step 1: the shadow Hamiltonian of a harmonic oscillator (omega = 1, m = 1), velocity Verlet")
def vv_h(x, p, h, w=1.0):
    p = p - 0.5 * h * w * w * x; x = x + h * p; p = p - 0.5 * h * w * w * x; return x, p
h = 0.5; x, p = 1.0, 0.0
E = lambda x, p: 0.5 * p * p + 0.5 * x * x
Hs = lambda x, p, h: 0.5 * p * p + 0.5 * x * x * (1 - h * h / 4)
print(f"   h = {h}: shadow energy H~ = p^2/2 + (x^2/2)(1 - h^2/4)")
for n in range(5):
    print(f"   step {n}: x = {x:+.6f}  p = {p:+.6f}   E = {E(x, p):.6f}   H~ = {Hs(x, p, h):.6f}")
    x, p = vv_h(x, p, h)
print(f"   shadow frequency (2/h) arcsin(h/2) = {2 / h * np.arcsin(h / 2):.6f} (true 1): the Verlet orbit turns {100 * (2 / h * np.arcsin(h / 2) - 1):.2f}% fast")
print(f"   stable only while h*omega < 2: h = 2.1 -> 1 - h^2/4 = {1 - 2.1**2 / 4:.4f} < 0: H~ is no longer a closed ellipse, the orbit grows")
x, p = 1.0, 0.0
for n in range(20): x, p = vv_h(x, p, 2.1)
print(f"   check, h = 2.1, after 20 steps: |x| = {abs(x):.3e}")

print("\\nStep 2: shadow energy of a 13-atom LJ cluster;  H~ = H + h^2 (p.Hess.p/12 - |F|^2/24)")
def lj(x):
    d = x[:, None, :] - x[None, :, :]; r2 = (d * d).sum(-1); np.fill_diagonal(r2, 1.0); ir6 = 1 / r2**3
    U = np.triu(4 * (ir6 * ir6 - ir6), 1).sum(); f = 24 * (2 * ir6 * ir6 - ir6) / r2; np.fill_diagonal(f, 0)
    return (f[:, :, None] * d).sum(1), U
def pHp(x, p):
    d = x[:, None, :] - x[None, :, :]; r2 = (d * d).sum(-1); np.fill_diagonal(r2, 1.0); r = np.sqrt(r2)
    w = p[:, None, :] - p[None, :, :]; rw = (d * w).sum(-1) / r; ww = (w * w).sum(-1)
    U1 = 4 * (-12 * r**-13 + 6 * r**-7); U2 = 4 * (156 * r**-14 - 42 * r**-8)
    c = U2 * rw**2 + U1 / r * (ww - rw**2); np.fill_diagonal(c, 0); return np.triu(c, 1).sum()
g = (1 + 5**0.5) / 2
ico = np.array([[0, s1, s2 * g] for s1 in (1, -1) for s2 in (1, -1)] + [[s1, s2 * g, 0] for s1 in (1, -1) for s2 in (1, -1)] + [[s2 * g, 0, s1] for s1 in (1, -1) for s2 in (1, -1)], float)
x13 = np.vstack([[0, 0, 0], ico / np.linalg.norm(ico[0]) * 1.1])
for _ in range(2000):
    F, U = lj(x13); x13 += 0.002 * F
rng = np.random.default_rng(1); p13 = rng.normal(size=(13, 3)) * np.sqrt(0.25); p13 -= p13.mean(0)
print(f"   icosahedron minimum U = {lj(x13)[1]:.4f} (literature -44.3268);  start at kT = 0.25")
def run13(h, T, dtype=np.float64, offset=0.0, every=10):
    x = (x13 + offset).astype(dtype); p = p13.astype(dtype); hh = dtype(h); F, U = lj(x); out = []
    for s in range(int(round(T / h))):
        p = p + dtype(0.5) * hh * F; x = x + hh * p; F, U = lj(x); p = p + dtype(0.5) * hh * F
        if s % every == 0:
            x64, p64 = x.astype(np.float64), p.astype(np.float64); F64, U64 = lj(x64)
            Ek = 0.5 * (p64 * p64).sum(); Et = U64 + Ek
            out.append((s * h, Et, Et + h * h * (pHp(x64, p64) / 12 - (F64 * F64).sum() / 24)))
    return np.array(out)
rows = []
for h in (0.008, 0.004, 0.002):
    o = run13(h, 50, every=max(1, int(0.02 / h))); rows.append((h, o[:, 1].std(), o[:, 2].std()))
    print(f"   h = {h}: std of E = {o[:, 1].std():.3e}   std of H~ = {o[:, 2].std():.3e}   ratio {o[:, 1].std() / o[:, 2].std():.0f}")
print(f"   halving h: E fluctuation / {rows[0][1] / rows[1][1]:.2f}, {rows[1][1] / rows[2][1]:.2f} (h^2 -> 4);   H~ fluctuation / {rows[0][2] / rows[1][2]:.1f}, {rows[1][2] / rows[2][2]:.1f} (h^4 -> 16)")
print("   Verlet vs RK4 (4th order, not symplectic), harmonic oscillator, h = 0.2, 100 000 steps:")
def rk4(x, p, h):
    f = lambda x, p: (p, -x)
    k1 = f(x, p); k2 = f(x + h / 2 * k1[0], p + h / 2 * k1[1]); k3 = f(x + h / 2 * k2[0], p + h / 2 * k2[1]); k4 = f(x + h * k3[0], p + h * k3[1])
    return x + h / 6 * (k1[0] + 2 * k2[0] + 2 * k3[0] + k4[0]), p + h / 6 * (k1[1] + 2 * k2[1] + 2 * k3[1] + k4[1])
for name, stepf in (("Verlet", lambda x, p: vv_h(x, p, 0.2)), ("RK4", lambda x, p: rk4(x, p, 0.2))):
    x, p = 1.0, 0.0; Es = []
    for n in range(100000):
        x, p = stepf(x, p)
        if n % 10 == 0: Es.append(E(x, p))
    Es = np.array(Es); print(f"   {name:6s}: E after 10 steps {Es[1]:.7f}, min/max over run {Es.min():.5f}/{Es.max():.5f}, final {Es[-1]:.5f}  (start 0.5)")
print(f"   RK4 per-step energy factor 1 - h^6/72 = {1 - 0.2**6 / 72:.9f};  ^100000 = {(1 - 0.2**6 / 72)**100000:.4f}")

print(f"\\nStep 3: single vs double precision ({time.time() - t0:.0f} s)")
for X in (1, 10, 100, 1000):
    print(f"   float32 spacing at {X:5d}: {np.spacing(np.float32(X)):.3e}   float64: {np.spacing(np.float64(X)):.3e}")
print("   LJ13, h = 0.004, t = 200 (50 000 steps), cluster placed at a distance 'offset' from the origin")
for dt_, off in ((np.float64, 0), (np.float32, 0), (np.float32, 10), (np.float32, 100), (np.float32, 1000)):
    o = run13(0.004, 200, dt_, off, every=50); tt = o[:, 0]; e = o[:, 1] - o[0, 1]; s = np.polyfit(tt, o[:, 1], 1)[0]
    print(f"   {dt_.__name__:8s} offset {off:5d}: drift of E = {s * 200:+.3e} over t = 200,  rms of E - E0 = {np.sqrt((e * e).mean()):.3e}")

print(f"\\nStep 4: RESPA by hand: x'' = -(w_f^2 + w_s^2) x with w_f = 10 (fast spring) and w_s = 1 (slow spring) ({time.time() - t0:.0f} s)")
wf, ws = 10.0, 1.0
def respa(x, p, H, n):
    p -= 0.5 * H * ws * ws * x; h = H / n
    for _ in range(n): p -= 0.5 * h * wf * wf * x; x += h * p; p -= 0.5 * h * wf * wf * x
    p -= 0.5 * H * ws * ws * x; return x, p
Ef = lambda x, p: 0.5 * p * p + 0.5 * (wf * wf + ws * ws) * x * x
x, p = 1.0, 0.0; H = 0.1; hn = 0.05
print(f"   one outer step H = {H}, two inner steps h = {hn}:  start x = {x}, p = {p}, E = {Ef(x, p):.4f}")
p -= 0.5 * H * ws * ws * x; print(f"   slow half kick:  p -= (H/2) ws^2 x = 0.05 x 1 x {x} -> p = {p:+.5f}")
for k in range(2):
    p -= 0.5 * hn * wf * wf * x; pa = p; x += hn * p; xb = x; p -= 0.5 * hn * wf * wf * x
    print(f"   inner {k + 1}: fast kick p = {pa:+.5f}; drift x = {xb:+.5f}; fast kick p = {p:+.5f}")
p -= 0.5 * H * ws * ws * x; print(f"   slow half kick: p = {p:+.5f};  E = {Ef(x, p):.4f}   (slow force computed once, fast force twice)")
print("   resonance scan: inner step fixed near 0.005, outer step H varied, t = 400:")
wt = np.sqrt(wf * wf + ws * ws)
for H in (0.05, 0.1, 0.2, 0.28, 0.30, 0.3126, 0.33, 0.4, 0.6, 0.625, 0.65):
    n = int(np.ceil(H / 0.005)); x, p = 1.0, 0.0; E0 = Ef(x, p); mx = 0
    for _ in range(int(400 / H)):
        x, p = respa(x, p, H, n); mx = max(mx, abs(Ef(x, p) / E0 - 1))
        if mx > 1e6: break
    print(f"   H = {H:.4f} (H x omega/pi = {H * wt / np.pi:.3f}), n = {n:3d}: max |E/E0 - 1| = {mx:.3e}")

print(f"\\nStep 5: RESPA in a liquid of 32 stiff dimers (bond k = 2000, LJ between molecules), 3D, rho = 0.6, T = 1.5 ({time.time() - t0:.0f} s)")
N = 64; L = (N / 0.6) ** (1 / 3); RC = 2.3; USH = 4 * (RC**-12 - RC**-6); KB = 2000.0
partner = np.arange(N) ^ 1
def lj_pbc(x):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L); r2 = (d * d).sum(-1)
    m = (r2 < RC * RC); np.fill_diagonal(m, False); m[np.arange(N), partner] = False; r2 = np.where(m, r2, 1.0); ir6 = 1 / r2**3
    U = 0.5 * np.where(m, 4 * (ir6 * ir6 - ir6) - USH, 0).sum(); f = np.where(m, 24 * (2 * ir6 * ir6 - ir6) / r2, 0); return (f[:, :, None] * d).sum(1), U
def bond(x):
    a, b = x[0::2], x[1::2]; d = a - b; d -= L * np.round(d / L); r = np.sqrt((d * d).sum(1)); f = (-KB * (r - 1.0) / r)[:, None] * d
    F = np.zeros_like(x); F[0::2] = f; F[1::2] = -f; return F, 0.5 * KB * ((r - 1.0)**2).sum()
cs = np.array([[(i + .5) * L / 4, (j + .5) * L / 4, (k + .5) * L / 2] for i in range(4) for j in range(4) for k in range(2)])
x = np.repeat(cs, 2, 0); x[0::2, 2] -= 0.5; x[1::2, 2] += 0.5
rng = np.random.default_rng(3); v = rng.normal(size=(N, 3)) * np.sqrt(1.5)
F = lj_pbc(x)[0] + bond(x)[0]
for s in range(4000):   # Langevin equilibration with a small step
    v += 0.0005 * F; x = (x + 0.001 * v) % L; F = lj_pbc(x)[0] + bond(x)[0]; v += 0.0005 * F
    c = np.exp(-0.001); v = c * v + np.sqrt((1 - c * c) * 1.5) * rng.normal(size=v.shape)
v -= v.mean(0); xs, vs = x.copy(), v.copy()
mu = 0.5; w = np.sqrt(KB / mu); print(f"   bond frequency sqrt(k/mu) = sqrt({KB:.0f}/0.5) = {w:.1f}, period {2 * np.pi / w:.4f};  Verlet limit h < 2/omega = {2 / w:.4f};  resonance H = pi/omega = {np.pi / w:.4f}")
def energy(x, v): return lj_pbc(x)[1] + bond(x)[1] + 0.5 * (v * v).sum()
def run(H, n, T=10.0):
    x, v = xs.copy(), vs.copy(); Fs = lj_pbc(x)[0]; Ff = bond(x)[0]; Es = []; h = H / n; nlj = 0; nb = 0
    for s in range(int(round(T / H))):
        v += 0.5 * H * Fs
        for _ in range(n): v += 0.5 * h * Ff; x = (x + h * v) % L; Ff = bond(x)[0]; v += 0.5 * h * Ff; nb += 1
        Fs = lj_pbc(x)[0]; nlj += 1; v += 0.5 * H * Fs
        if not np.isfinite(v).all() or (v * v).sum() > 1e8: return None, 1 / H, n / H
        if s % max(1, int(0.01 / H)) == 0: Es.append(energy(x, v) / N)
    return np.array(Es), nlj / T, nb / T
print("   method                    LJ evals per time   bond evals per time   std E/N      drift E/N over t = 10")
for name, H, n in (("Verlet h = 0.001", 0.001, 1), ("Verlet h = 0.004", 0.004, 1), ("Verlet h = 0.008", 0.008, 1), ("RESPA H = 0.004, n = 4", 0.004, 4), ("RESPA H = 0.008, n = 8", 0.008, 8), ("RESPA H = 0.02, n = 20", 0.02, 20), ("RESPA H = 0.04, n = 40", 0.04, 40), ("RESPA H = 0.05, n = 50", 0.05, 50)):
    Es, a, b = run(H, n)
    if Es is None: print(f"   {name:24s}  {a:8.0f}            {b:8.0f}              blew up"); continue
    tt = np.linspace(0, 10, len(Es)); print(f"   {name:24s}  {a:8.0f}            {b:8.0f}              {Es.std():.2e}   {np.polyfit(tt, Es, 1)[0] * 10:+.2e}")
print(f"\\ntotal time {time.time() - t0:.0f} s")
`,
  },

  exercises: [
    { q: 'For h = 1 on the oscillator, what is the shadow energy for x = 1, p = 0, and what range does the true energy E cover?', a: 'H̃ = 0 + (1/2)(1 − 1/4) = 0.375. Along the orbit E = H̃ + h²x²/8 = 0.375 + x²/8. x² is largest when p = 0: x² = 2 × 0.375/0.75 = 1, giving E = 0.5; it is smallest (0) when x = 0, giving E = 0.375. So E swings between 0.375 and 0.5 (25%) but never drifts.' },
    { q: 'Why is the float32 drift at the origin no worse than float64, when float32 has 9 fewer digits?', a: 'Near the origin the rounding error of a step (≈ 10⁻⁷) is far smaller than the Verlet shadow wobble (≈ 10⁻³ in energy). The rounding is also nearly unbiased there. Precision only matters once the rounding of each displacement becomes a sizeable fraction of the displacement itself, which happens far from the origin, or for very small steps.' },
    { q: 'A water model has O–H stretches with a period of about 9 fs. Using RESPA with bonds on the inner loop, what outer step causes resonance, and why do most water simulations use SHAKE instead?', a: 'Resonance at H = half the period = 4.5 fs (and 9 fs, …); the practical limit is about 0.3 × 9 ≈ 3 fs. With rigid water (SETTLE/SHAKE) the stretch is removed entirely; the next fastest motions (librations, about 20–30 fs period) allow 2 fs steps with plain Verlet, at no extra cost and without the resonance problem.' },
  ],
};
