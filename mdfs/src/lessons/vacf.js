import { tex } from '../tex.js';

export default {
  id: 'vacf',
  summary: 'The velocity autocorrelation function (VACF) asks how long an atom remembers its velocity. Its shape tells you whether atoms fly freely, rattle in a cage or vibrate in a lattice; its integral gives the diffusion coefficient (Green–Kubo), and its Fourier transform gives the vibrational spectrum.',
  minutes: 14,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/vacf.mp4', poster: 'videos/vacf.jpg' },
  sim: 'vacf',
  learn: [
    'What the VACF is and how to compute it from saved velocities (all time origins)',
    'Typical shapes for gas, liquid and solid, and what the negative dip means',
    'Green–Kubo diffusion and the vibrational density of states from the same function',
  ],
  simHint: 'Watch the liquid VACF drop through zero after about 0.15 τ and come back: the atom bounced off its neighbours. The Green–Kubo D(t) overshoots and settles near 0.03, the same as the MSD widget. Switch to Solid: C(t) oscillates and D(t) goes to zero. Gas: no dip at all, just a slow decay. Open "Spectrum" to compare the liquid and solid vibrations.',

  idea: `
<p>Take an atom's velocity at some moment, <b>v</b>(0), and again a time t later, <b>v</b>(t). Their dot product is positive if the atom is still moving the same way, negative if it has turned around, and zero if there is no relation. Averaged over all atoms and all starting times this is the <b>velocity autocorrelation function</b>:</p>
<p style="text-align:center">C(t) = ⟨ <b>v</b>(t) · <b>v</b>(0) ⟩, &nbsp; C(0) = ⟨v²⟩ = 3kT/m</p>
<p><b>Three shapes</b> (code tab, Parts 1 and 2):</p>
<ul>
<li><b>Dilute gas</b>: atoms fly straight until a rare collision. C(t) decays slowly and almost exponentially, never negative.</li>
<li><b>Dense liquid</b>: an atom quickly hits the "cage" of its neighbours and bounces back. C(t) falls through zero after about 0.15 τ (0.3 ps in argon), dips to about −0.13 and then decays to zero.</li>
<li><b>Solid</b>: atoms vibrate about lattice sites, so C(t) oscillates and its integral is zero.</li>
</ul>
<p><b>Green–Kubo.</b> The diffusion coefficient is the area under the VACF: D = ⅓∫C(t)dt. It is exactly equivalent to the MSD slope (they are two views of the same motion). For the LJ liquid both give D* ≈ 0.031 from the same run. The running integral overshoots, because the positive early part comes first, then settles once the negative dip has been added.</p>
<p><b>Spectrum.</b> The Fourier transform of C(t) is the <b>vibrational density of states</b> (VDOS): which frequencies the atoms move at (Part 3). A solid has peaks and S(0) = 0; a liquid has a broad band and a nonzero value at ω = 0, proportional to D. For molecules the VDOS shows the bond-stretch and bend frequencies and can be compared with IR, Raman and inelastic neutron spectra (with the right weighting).</p>
<p>Long-time tails: in liquids C(t) decays as t<sup>−3/2</sup> at long times (hydrodynamics), so the Green–Kubo integral converges slowly; noise in the tail is the main source of error.</p>
<div class="callout"><b>Key idea.</b> C(t) = ⟨v(t)·v(0)⟩ measures velocity memory. Its area gives D (Green–Kubo, same as the MSD), the negative dip shows caging, and its Fourier transform gives the vibrational spectrum.</div>`,

  algorithm: `
<pre class="pseudo">save velocities every n steps (short interval: a few fs for molecules, ~0.01 τ here)
for each lag k = 0 ... K:
    C[k] = mean over atoms i and origins t0 of v_i(t0 + k) · v_i(t0)
D(t)  = (1/3) * cumulative trapezoid integral of C            # running Green-Kubo integral
D     = plateau value of D(t)
S(ω)  = (2 / C[0]) Σ_k C[k] w[k] cos(ω k Δt) Δt                 # VDOS, w = smooth window</pre>
<ol class="steps">
  <li>Run in NVE (or with a very weak thermostat) for dynamical properties; strong thermostats change the VACF.</li>
  <li>The FFT gives C(t) for all lags at once (Wiener–Khinchin): pad to twice the length, |FFT|², inverse FFT.</li>
  <li>For molecules, use atom (or centre-of-mass) velocities depending on what you want: COM for diffusion, atoms for spectra; mass-weight the VDOS for the total.</li>
  <li>Tools: MDAnalysis, TRAVIS (power spectra, IR/Raman with dipoles), GROMACS <code>gmx velacc</code>, LAMMPS <code>compute vacf</code> (single origin).</li>
</ol>`,

  code: {
    note: 'Part 1 computes the liquid VACF from a 30 τ NVE run and compares Green–Kubo with the MSD slope from the same trajectory. Part 2 shows the gas and solid. Part 3 Fourier-transforms the VACF into the vibrational density of states. Run <code>pip install numpy</code>, then <code>python vacf.py</code> (about 30 s).',
    src: `import numpy as np, time
rng = np.random.default_rng(3)

def fcc(n, rho):
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + b[None]).reshape(-1, 3)) * a, n * a

def forces(x, L, rc=2.5):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)
    r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
    ir6 = np.where(m, 1 / r2**3, 0.0); f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
    return (f[:, :, None] * d).sum(1)

def run(ncell, rho, T, nequil, nprod, every, dt=0.005):
    """Equilibrate with velocity rescaling, then NVE; return velocities and unwrapped positions every \`every\` steps."""
    x, L = fcc(ncell, rho); N = len(x)
    v = rng.normal(0, np.sqrt(T), x.shape); v -= v.mean(0); f = forces(x, L)
    for s in range(nequil):
        v += 0.5 * dt * f; x = (x + dt * v) % L; f = forces(x, L); v += 0.5 * dt * f
        v *= np.sqrt(T * (3 * N - 3) / (v**2).sum())
    u = x.copy(); V, U = [], []
    for s in range(nprod):
        v += 0.5 * dt * f; x = (x + dt * v) % L; u += dt * v; f = forces(x, L); v += 0.5 * dt * f
        if s % every == 0: V.append(v.copy()); U.append(u.copy())
    return np.array(V), np.array(U), dt * every

def vacf(V, max_lag):
    """C(t) = <v(t0 + t) . v(t0)>, averaged over atoms and all time origins."""
    return np.array([(V[k:] * V[:len(V) - k]).sum(-1).mean() for k in range(max_lag + 1)])

def integrate(y, h):                     # trapezoid rule, cumulative
    return np.concatenate([[0], np.cumsum(0.5 * (y[1:] + y[:-1]) * h)])

# ---------- Part 1: liquid VACF and Green-Kubo diffusion ----------
t0 = time.time(); T = 0.722
V, U, h = run(4, 0.8442, T, 1500, 6000, 2)
C = vacf(V, 300); t = np.arange(301) * h; D_gk = integrate(C, h) / 3
print(f"Part 1: LJ liquid (rho* 0.8442, T* {T}), {len(V)} frames, {time.time() - t0:.0f} s")
print(f"   C(0) = {C[0]:.3f}  (should be 3 kT/m = {3 * T:.3f})")
k0 = np.argmax(C < 0); kmin = np.argmin(C)
print(f"   C(t) first crosses zero at t = {t[k0]:.2f} tau, minimum C/C(0) = {C[kmin] / C[0]:.3f} at t = {t[kmin]:.2f} tau (atoms bounce back off their cage)")
for tt in [0.5, 1.0, 1.5, 3.0]:
    print(f"   Green-Kubo running integral D(t) = 1/3 int C: t = {tt:.1f}: {D_gk[int(round(tt / h))]:.4f}")
msd = np.array([((U[k:] - U[:len(U) - k])**2).sum(-1).mean() for k in range(100, 301, 50)])
D_msd = np.polyfit(t[100:301:50], msd, 1)[0] / 6
print(f"   Einstein (MSD slope / 6, 1-3 tau) from the same run: D = {D_msd:.4f}")

# ---------- Part 2: gas and solid ----------
print("\\nPart 2: gas vs solid (108 atoms)")
Vg, _, hg = run(3, 0.05, 2.0, 1000, 6000, 10)
Cg = vacf(Vg, 200); tg = np.arange(201) * hg; tau_c = integrate(Cg / Cg[0], hg)[-1]
print(f"   gas (rho 0.05, T 2): C(t) never negative: min C/C(0) = {Cg.min() / Cg[0]:.3f};  decay time int C/C(0) = {tau_c:.2f} tau;  D = {integrate(Cg, hg)[-1] / 3:.2f}")
Vs, _, hs = run(3, 1.1, 0.3, 500, 3000, 2)
Cs = vacf(Vs, 300); ts = np.arange(301) * hs
print(f"   solid (rho 1.1, T 0.3): C/C(0) oscillates: min {Cs.min() / Cs[0]:.3f} at t = {ts[np.argmin(Cs)]:.2f} tau;  D = 1/3 int C up to 3 tau = {integrate(Cs, hs)[-1] / 3:.4f} (zero)")

# ---------- Part 3: vibrational density of states (Fourier transform of the VACF) ----------
def vdos(C, h, nfreq=400):
    w = np.linspace(0, 40, nfreq)                                    # angular frequency, 1/tau
    win = np.hanning(2 * len(C))[len(C):]                            # smooth cut-off at the end of C
    return w, np.array([2 * h * np.sum(C * win * np.cos(wi * np.arange(len(C)) * h)) for wi in w]) / C[0]
print("\\nPart 3: vibrational density of states S(omega) = (2 / C(0)) int C(t) cos(omega t) dt")
for name, CC, hh, D in [("liquid", C, h, D_gk[-1]), ("solid", Cs, hs, 0.0)]:
    w, S = vdos(CC, hh)
    print(f"   {name:6s}: peak at omega = {w[np.argmax(S)]:4.1f} / tau;  S(0) = {S[0]:.3f}  (theory 2 D m / kT = {2 * D / (T if name == 'liquid' else 0.3):.3f}: nonzero only if atoms diffuse)")
`,
  },

  maths: () => `
<h4>1. Definition and normalization</h4>
${tex(String.raw`C(t) = \big\langle \mathbf v_i(t_0 + t)\cdot\mathbf v_i(t_0)\big\rangle_{i,\,t_0}, \qquad C(0) = \langle v^2\rangle = \frac{3k_BT}{m}`, true)}
<h4>2. Green–Kubo relation for diffusion</h4>
${tex(String.raw`D = \frac13\int_0^\infty C(t)\,dt`, true)}
<p>It follows from the MSD: r(t) − r(0) = ∫v dt, so MSD(t) = 2∫₀ᵗ(t − s) C(s) ds, whose slope at long times is 2∫C = 6D.</p>
<h4>3. Simple model: Langevin (exponential) decay</h4>
${tex(String.raw`C(t) = \frac{3k_BT}{m}\,e^{-\gamma t} \;\Rightarrow\; D = \frac{k_BT}{m\gamma}`, true)}
<p>A good description of a dilute gas or a Brownian particle; it cannot produce the negative dip of a liquid.</p>
<h4>4. Vibrational density of states</h4>
${tex(String.raw`S(\omega) = \frac{2}{C(0)}\int_0^\infty C(t)\cos(\omega t)\,dt, \qquad S(0) = \frac{6D}{C(0)} = \frac{2Dm}{k_BT}`, true)}
<h4>5. Long-time tail (Alder and Wainwright 1970)</h4>
${tex(String.raw`C(t) \propto t^{-3/2} \quad (t\to\infty, \text{ 3D fluid})`, true)}`,

  quiz: [
    {
      q: 'The VACF of a dense liquid becomes negative after a short time. What does that mean physically?',
      options: ['The temperature is negative', 'Atoms bounce back off the cage of their neighbours, reversing their velocity', 'The simulation is unstable', 'The atoms are diffusing faster'],
      answer: 1,
      why: 'Negative correlation means v(t) points on average opposite to v(0): the atom has collided with its neighbours and rebounded. This caging is the hallmark of dense liquids.',
    },
    {
      q: 'How do you get the diffusion coefficient from the VACF in 3D?',
      options: ['D = C(0) / 3', 'D = (1/3) ∫ C(t) dt', 'D = slope of C(t) / 6', 'D = the time where C crosses zero'],
      answer: 1,
      why: 'The Green–Kubo relation D = (1/3)∫₀^∞ C(t) dt. It gives the same D as the long-time MSD slope divided by 6.',
    },
    {
      q: 'The vibrational density of states of a crystal at ω = 0 is zero, but for a liquid it is not. Why?',
      options: ['Crystals have no vibrations', 'S(0) is proportional to the diffusion coefficient, which is zero in a crystal', 'Liquids are hotter', 'It is a numerical artefact'],
      answer: 1,
      why: 'S(0) = 2∫C/C(0) = 2Dm/kT. Atoms in a crystal do not diffuse (D = 0), so there is no zero-frequency motion.',
    },
  ],
};
