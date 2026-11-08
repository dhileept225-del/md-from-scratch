import { tex } from '../tex.js';

export default {
  id: 'conductivity',
  summary: 'Ionic conductivity measures how much current ions carry in an electric field. From an equilibrium MD run it follows from the fluctuations of the total charge displacement (Einstein–Helfand) or of the electric current (Green–Kubo). Comparing it with the Nernst–Einstein estimate from self-diffusion reveals ion pairing and correlated motion.',
  minutes: 16,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/conductivity.mp4', poster: 'videos/conductivity.jpg' },
  sim: 'conductivity',
  learn: [
    'The Einstein–Helfand and Green–Kubo formulas for conductivity from an equilibrium run',
    'The Nernst–Einstein approximation, ionicity and the Haven ratio, and why ion pairs lower σ',
    'Why conductivity is so noisy, and how to convert to S/m and mS/cm',
  ],
  simHint: 'With no pairing (rate 0) the red and dashed blue curves lie on top of each other: σ = σ_NE. Raise the pairing rate: more ions pair up, D drops a little but the red curve (real conductivity) drops much more. Make pairs live longer or diffuse faster and watch the ionicity. Note how long you must wait before the red curve is smooth: collective properties need long runs.',

  idea: `
<p>Apply a small electric field <b>E</b> and ions drift: cations one way, anions the other. The current density is J = σE, where σ is the <b>ionic conductivity</b>. You could simulate exactly that (non-equilibrium MD), but linear-response theory gives σ from an ordinary equilibrium run, by watching how the ions' charge spreads out spontaneously.</p>
<p><b>Einstein–Helfand.</b> Add up every ion's displacement weighted by its charge, M(t) = Σ q<sub>i</sub>Δ<b>r</b><sub>i</sub>(t). Its mean square grows linearly in time, and the slope gives σ, exactly as the MSD gives D:</p>
<p style="text-align:center">σ = (1 / 6VkT) · d/dt ⟨|Σ q<sub>i</sub> Δ<b>r</b><sub>i</sub>(t)|²⟩</p>
<p>The equivalent Green–Kubo form integrates the autocorrelation of the current J = Σq<sub>i</sub><b>v</b><sub>i</sub>.</p>
<p><b>Nernst–Einstein.</b> Expanding the square gives self terms (each ion with itself, i.e. its own MSD) and cross terms (pairs of different ions). Keeping only the self terms gives the <b>Nernst–Einstein</b> estimate from the self-diffusion coefficients (Lesson 6.2):</p>
<p style="text-align:center">σ<sub>NE</sub> = (e² / VkT) (N<sub>+</sub>q<sub>+</sub>²D<sub>+</sub> + N<sub>−</sub>q<sub>−</sub>²D<sub>−</sub>)</p>
<p>The cross terms describe correlated motion. If a cation and an anion move together as a neutral pair, both still diffuse (they count in D), but together they carry no current. Then σ &lt; σ<sub>NE</sub> (code tab, Part 1).</p>
<ul>
<li><b>Ionicity</b> (or degree of dissociation) = σ/σ<sub>NE</sub>; its inverse is the <b>Haven ratio</b> H<sub>R</sub> = σ<sub>NE</sub>/σ.</li>
<li>For many ionic liquids and concentrated electrolytes, measured ionicities are roughly 0.5–0.8. Cation–cation and anion–anion correlations matter too, not only ion pairs: the full picture needs the Onsager coefficients (Lesson 6.6).</li>
</ul>
<p><b>Statistics.</b> D averages over N ions, but M(t) is one collective number per frame, so σ converges about N times more slowly (Part 2). Use long runs, many independent replicas, and block averaging, and always quote an error bar. In molecular systems remember that σ is also affected by the force field (for example, polarizability or charge scaling, Lesson 5.4).</p>
<div class="callout"><b>Key idea.</b> σ = slope of ⟨|Σq<sub>i</sub>Δr<sub>i</sub>|²⟩ / 6VkT. Keeping only self terms gives Nernst–Einstein. Their ratio (ionicity) measures correlated ion motion. Conductivity needs much more sampling than diffusion.</div>`,

  algorithm: `
<pre class="pseudo">save UNWRAPPED positions every n steps
M(t) = Σ_i q_i u_i(t)                                   # one 3-vector per frame
for each lag k:  MSD_M[k] = < |M(t0 + k) - M(t0)|² >_t0       # all time origins
σ    = slope of MSD_M (linear region) / (6 V kB T)
σ_NE = (N+ q+² D+ + N- q-² D-) / (V kB T)                    # D from the per-species MSD
ionicity = σ / σ_NE;   Haven ratio = σ_NE / σ
errors: blocks or independent replicas (≥ 5);  check that the slope is linear
SI units: σ [S/m] = σ_reduced · e² / (kB T · length · time) with your MD units</pre>
<ol class="steps">
  <li>Use the same lag-time window for σ and for the D values in σ_NE.</li>
  <li>Check that the total charge is exactly zero. Then M(t) does not change when all ions shift together, so any centre-of-mass drift cancels.</li>
  <li>Green–Kubo alternative: save velocities often, compute ⟨J(t)·J(0)⟩, integrate, and divide by 3VkT.</li>
  <li>Non-equilibrium check: apply a weak field E, measure the steady drift current J, and use σ = J/E (field small enough to be linear).</li>
  <li>Tools: TRAVIS, MDAnalysis plus your own script (as here), GROMACS <code>gmx current</code>, LAMMPS via <code>compute</code> of the charge-weighted displacement.</li>
</ol>`,

  code: {
    note: 'A fast model of ion pairing: 500 cations and 500 anions that diffuse freely or as neutral pairs, with pairs forming and breaking at random. Part 1 compares the Einstein–Helfand conductivity with Nernst–Einstein (and with the exact answer for this model), Part 2 shows how noisy σ is compared with D, and Part 3 converts to S/m and mS/cm. Run <code>pip install numpy</code>, then <code>python conductivity.py</code> (about 25 s).',
    src: `import numpy as np
rng = np.random.default_rng(5)
# Reduced units (kT = 1, charge e = 1). A model of ions that diffuse freely or as neutral +- pairs (ion pairing).
# Free ions: diffusion D_free. A bound pair moves as one body with D_pair: both ions get the SAME displacement.
# Pairs form and break at random (rates k_on, k_off), so every ion is sometimes free and sometimes paired.

def simulate(n_pairs_frac, n=500, L=20.0, D_free=1.0, D_pair=0.5, k_off=1.0, dt=0.01, nsteps=30000, every=10):
    """n cations + n anions; returns unwrapped positions every \`every\` steps, charges, volume, frame interval."""
    p_bound = n_pairs_frac                                         # target fraction of ions in pairs
    k_on = k_off * p_bound / max(1e-12, 1 - p_bound)              # detailed balance: bound/free = k_on/k_off
    r = rng.random((2 * n, 3)) * L; partner = -np.ones(n, int)     # partner[i] = anion index bound to cation i
    q = np.r_[np.ones(n), -np.ones(n)]; traj = []; fb = []
    free_a = lambda: np.setdiff1d(np.arange(n), partner[partner >= 0])
    for s in range(nsteps):
        bound = partner >= 0
        brk = bound & (rng.random(n) < k_off * dt); partner[brk] = -1                         # pairs break
        for i in np.where((partner < 0) & (rng.random(n) < k_on * dt))[0]:                  # free cations pair up
            fa = free_a()
            if len(fa): partner[i] = rng.choice(fa)
        step = rng.normal(0, np.sqrt(2 * D_free * dt), r.shape)                              # free motion
        b = np.where(partner >= 0)[0]
        common = rng.normal(0, np.sqrt(2 * D_pair * dt), (len(b), 3))
        step[b] = common; step[n + partner[b]] = common                                      # pair moves together
        r += step
        if s % every == 0: traj.append(r.copy()); fb.append((partner >= 0).mean())
    return np.array(traj), q, L**3, dt * every, np.mean(fb)            # time-averaged fraction of ions in pairs

def slope(t, y): return np.polyfit(t, y, 1)[0]
lags = np.arange(10, 51, 5)                                       # lag times 1 to 5 (pair lifetime is 1)

def analyse(traj, q, V, h):
    t = lags * h
    msd = lambda X: np.array([((X[k:] - X[:-k])**2).sum(-1).mean() for k in lags])
    Dp, Dm = slope(t, msd(traj[:, q > 0])) / 6, slope(t, msd(traj[:, q < 0])) / 6
    sig_NE = ((q > 0).sum() * Dp + (q < 0).sum() * Dm) / V                # Nernst-Einstein: ions independent
    M = (q[None, :, None] * traj).sum(1)                                  # collective charge displacement
    msdM = np.array([((M[k:] - M[:-k])**2).sum(-1).mean() for k in lags])
    sig = slope(t, msdM) / (6 * V)                                        # Einstein-Helfand conductivity
    bl = [slope(t, np.array([((Mb[k:] - Mb[:-k])**2).sum(-1).mean() for k in lags])) / (6 * V) for Mb in np.array_split(M, 4)]
    bD = [slope(t, msd(Tb[:, q > 0])) / 6 for Tb in np.array_split(traj, 4)]
    return Dp, Dm, sig_NE, sig, np.std(bl, ddof=1) / 2, np.std(bD, ddof=1) / 2

# ---------- Part 1: free ions vs paired ions ----------
print("Part 1: 500 cations + 500 anions, 300 time units each (D_free = 1, D_pair = 0.5)")
print("   paired   D+      D-      sigma_NE   sigma (Einstein-Helfand)   sigma/sigma_NE   exact (1-f)D_free / D_self")
for frac in [0.0, 0.3, 0.6, 0.9]:
    traj, q, V, h, fb = simulate(frac)
    Dp, Dm, sNE, s, ds, dD = analyse(traj, q, V, h)
    exact = (1 - fb) * 1.0 / ((1 - fb) * 1.0 + fb * 0.5)
    print(f"   {fb:5.2f}    {Dp:.3f}   {Dm:.3f}   {sNE:.4f}     {s:.4f} +- {ds:.4f}            {s / sNE:.2f} +- {ds / sNE:.2f}   {exact:.2f}")
    if frac == 0.0: rel = (ds / s, dD / Dp)
print("   paired ions still diffuse (D stays large) but a neutral pair carries no current: sigma falls below sigma_NE")
print("   measured and exact agree to within about two standard errors: even 1000 ions and 300 time units give noisy sigma")

# ---------- Part 2: why conductivity is so noisy ----------
print(f"\\nPart 2: relative statistical error from the same run (free ions): sigma {100 * rel[0]:.0f} %, D+ {100 * rel[1]:.1f} %")
print("   D averages over 500 ions; the collective charge displacement is ONE number per frame -> ~N times longer runs needed")

# ---------- Part 3: from D to S/m ----------
print("\\nPart 3: Nernst-Einstein in SI units for a 1 M 1:1 salt with D+ = 1.0e-9 and D- = 1.5e-9 m^2/s at 298 K")
e, kB, NA = 1.602176634e-19, 1.380649e-23, 6.02214076e23
n = 1.0 * 1000 * NA                                                         # ions of each sign per m^3
sigma_NE = e**2 / (kB * 298) * (n * 1.0e-9 + n * 1.5e-9)
print(f"   sigma_NE = e^2 (n+ D+ + n- D-) / kT = {sigma_NE:.2f} S/m = {sigma_NE * 10:.1f} mS/cm")
print(f"   with an ionicity sigma/sigma_NE = 0.6 (typical of ionic liquids and concentrated electrolytes): {sigma_NE * 0.6 * 10:.1f} mS/cm")
print(f"   Haven ratio H = sigma_NE / sigma = {1 / 0.6:.2f}")
`,
  },

  maths: () => `
<h4>1. Einstein–Helfand and Green–Kubo</h4>
${tex(String.raw`\sigma = \lim_{t\to\infty}\frac{1}{6Vk_BT}\frac{d}{dt}\Big\langle \Big|\sum_i q_i\big[\mathbf r_i(t) - \mathbf r_i(0)\big]\Big|^2\Big\rangle = \frac{1}{3Vk_BT}\int_0^\infty \langle \mathbf J(t)\cdot\mathbf J(0)\rangle\,dt, \qquad \mathbf J = \sum_i q_i\mathbf v_i`, true)}
<h4>2. Self and distinct terms</h4>
${tex(String.raw`\Big\langle\big|\textstyle\sum_i q_i\Delta\mathbf r_i\big|^2\Big\rangle = \underbrace{\textstyle\sum_i q_i^2\langle|\Delta\mathbf r_i|^2\rangle}_{\text{self} \;\to\; \sigma_{NE}} + \underbrace{\textstyle\sum_{i\ne j} q_iq_j\langle\Delta\mathbf r_i\cdot\Delta\mathbf r_j\rangle}_{\text{distinct (correlations)}}`, true)}
<h4>3. Nernst–Einstein, ionicity, Haven ratio</h4>
${tex(String.raw`\sigma_{NE} = \frac{1}{Vk_BT}\big(N_+q_+^2D_+ + N_-q_-^2D_-\big), \qquad \text{ionicity} = \frac{\sigma}{\sigma_{NE}}, \qquad H_R = \frac{\sigma_{NE}}{\sigma}`, true)}
<h4>4. Pair model used in the code tab</h4>
${tex(String.raw`D_{\text{self}} = (1 - f)D_{\text{free}} + fD_{\text{pair}}, \qquad \frac{\sigma}{\sigma_{NE}} = \frac{(1 - f)D_{\text{free}}}{(1 - f)D_{\text{free}} + fD_{\text{pair}}}`, true)}
<p>Here f is the fraction of time an ion spends in a pair (valid at times longer than the pair lifetime).</p>
<h4>5. Molar conductivity</h4>
${tex(String.raw`\Lambda = \frac{\sigma}{c}, \qquad \Lambda_{NE} = \frac{F^2}{RT}\big(\nu_+z_+^2D_+ + \nu_-z_-^2D_-\big)`, true)}`,

  quiz: [
    {
      q: 'In a simulation, many Li⁺ ions travel together with an anion as neutral pairs. How does the true conductivity σ compare with the Nernst–Einstein estimate σ_NE?',
      options: ['σ > σ_NE', 'σ < σ_NE', 'σ = σ_NE always', 'σ = 0'],
      answer: 1,
      why: 'Paired ions contribute to self-diffusion (so to σ_NE) but a neutral pair carries no net current. The cross terms are negative, so σ < σ_NE (ionicity < 1, Haven ratio > 1).',
    },
    {
      q: 'Why does the conductivity converge much more slowly than the self-diffusion coefficient in the same run?',
      options: ['Conductivity needs a smaller time step', 'D averages over all N ions, while the charge displacement M(t) is a single collective quantity', 'Conductivity needs wrapped coordinates', 'Because of the thermostat'],
      answer: 1,
      why: 'Each frame gives N independent samples for D but only one value of M(t) for σ, so σ has roughly √N larger relative error for the same run length.',
    },
    {
      q: 'σ_NE = 9.4 S/m for a 1 M salt solution and the measured conductivity is 5.6 S/m. What is the Haven ratio?',
      options: ['0.6', '1.7', '3.8', '15'],
      answer: 1,
      why: 'H_R = σ_NE / σ = 9.4 / 5.6 ≈ 1.7, i.e. an ionicity of 0.6.',
    },
  ],
};
