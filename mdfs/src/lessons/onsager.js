import { tex } from '../tex.js';

export default {
  id: 'onsager',
  summary: 'Onsager transport coefficients L_ij describe how the motion of one species is coupled to another. From them follow the conductivity and the transference number t+, the fraction of current carried by the cation. They show why ion clusters can make t+ small or even negative, and why the t+ estimated from self-diffusion (as in NMR) can mislead.',
  minutes: 18,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/onsager.mp4', poster: 'videos/onsager.jpg' },
  sim: 'onsager',
  learn: [
    'Onsager coefficients from collective displacements, and their self and distinct parts',
    'Conductivity and transference numbers from L++, L−− and L+−, and how clusters change them',
    'Why the reference frame matters, and why a binary molten salt has t+ fixed by the masses',
  ],
  simHint: 'Start with 20 % pairs: t+ (Onsager) stays near the free-ion value while σ/σ_NE drops. Then add [LiA₂]⁻ triplets and make free Li⁺ slow and triplets fast: t+ falls below zero while the NMR-like t+ stays around 0.4–0.5. Watch the amber Li–A curve: positive cross-correlation means cations and anions move together.',

  idea: `
<p>Conductivity (Lesson 6.5) tells you how much current flows. In a battery you also want to know <i>which ions</i> carry it: ideally Li⁺, because anion current only builds up concentration gradients. The <b>transference number</b> t₊ is the fraction of the current carried by the cations.</p>
<p><b>Onsager coefficients.</b> For each species i, add up the displacements of all its ions: ΔR<sub>i</sub>(t) = Σ Δ<b>r</b> over ions of i. The coefficients</p>
<p style="text-align:center">L<sup>ij</sup> = (1/6VkT) · d/dt ⟨ΔR<sub>i</sub>(t) · ΔR<sub>j</sub>(t)⟩</p>
<p>describe how strongly species i and j move together. For a 1:1 salt there are three: L<sup>++</sup>, L<sup>−−</sup> and L<sup>+−</sup>. They connect forces (gradients of electrochemical potential) to fluxes, like a matrix version of Fick's law. Each diagonal L<sup>ii</sup> splits into a <b>self</b> part (each ion's own MSD, i.e. N<sub>i</sub>D<sub>i</sub>) and a <b>distinct</b> part (different ions of the same species moving together); L<sup>+−</sup> is purely distinct.</p>
<p><b>Conductivity and transference</b> follow directly (z<sub>±</sub> = ±1):</p>
<p style="text-align:center">σ = e²(L<sup>++</sup> + L<sup>−−</sup> − 2L<sup>+−</sup>), &nbsp;&nbsp; t₊ = (L<sup>++</sup> − L<sup>+−</sup>) / (L<sup>++</sup> + L<sup>−−</sup> − 2L<sup>+−</sup>)</p>
<p><b>What clusters do</b> (code tab, Part 1):</p>
<ul>
<li><b>Neutral pairs</b> (LiA) raise L<sup>+−</sup>: they lower σ but leave t₊ almost unchanged, because they carry no charge.</li>
<li><b>Negative clusters</b> such as [LiA₂]⁻ carry Li⁺ <i>towards the positive electrode</i>. If they are mobile and free Li⁺ is slow, t₊ can become negative, as reported for some concentrated electrolytes and polymer electrolytes.</li>
<li>The "apparent" t₊ from self-diffusion, D₊/(D₊ + D₋) (what pulsed-field-gradient NMR gives), ignores all cross terms. Here it stays around 0.45 while the true t₊ is −0.24.</li>
</ul>
<p><b>The reference frame matters.</b> Fluxes must be measured relative to something: the solvent (solvent-fixed frame), the centre of mass (barycentric frame, what MD naturally conserves) or the centre of volume. In a pure binary molten salt or ionic liquid with no solvent, momentum conservation forces m₊ΔR₊ + m₋ΔR₋ = 0. Then all three L's are tied together and the barycentric t₊ = m₋/(m₊ + m₋), whatever the dynamics (Part 2). Experimental transference numbers are usually defined relative to the solvent or the anion; convert the frame before you compare.</p>
<div class="callout"><b>Key idea.</b> L<sup>ij</sup> = slope of ⟨ΔR<sub>i</sub>·ΔR<sub>j</sub>⟩/6VkT. σ and t₊ need all of them, including cross terms. Clusters make t₊ differ from D₊/(D₊ + D₋), even in sign. Always state the reference frame.</div>`,

  algorithm: `
<pre class="pseudo">save UNWRAPPED positions; subtract the reference-frame motion (e.g. solvent centre of mass, or total COM)
R_i(t) = Σ_{ions of species i} u(t)                   # one vector per species per frame
for each pair (i, j) and lag k:
    C_ij[k] = < (R_i(t0+k) - R_i(t0)) · (R_j(t0+k) - R_j(t0)) >_t0
L_ij = slope of C_ij in the linear region / (6 V kB T)
self part of L_ii = N_i D_i / (V kB T);   distinct part = L_ii - self
σ   = e² Σ_ij z_i z_j L_ij
t_+ = (z_+² L_++ + z_+ z_- L_+-) · e² / σ</pre>
<ol class="steps">
  <li>Statistics are collective (one R per species), so the same advice as for σ applies: long runs, several replicas, error bars. Off-diagonal L's are the noisiest.</li>
  <li>For an electrolyte with solvent, work in the solvent-fixed frame to compare with Hittorf or concentration-cell experiments; Bruce–Vincent (polarization) measurements give a related but different quantity in concentrated systems.</li>
  <li>Check: L<sup>+−</sup> > 0 means cations and anions move together (pairing); negative values mean they move oppositely (as in the barycentric frame of a molten salt).</li>
  <li>Further reading: Wheeler and Newman (2004) for the concentrated-solution framework and reference frames; Fong et al. (2020, 2021) for Onsager coefficients from MD of polymer electrolytes. Most groups compute the L's with their own scripts, as here.</li>
</ol>`,

  code: {
    note: 'Part 1 uses a cluster model (free Li⁺ and A⁻, neutral LiA pairs, negative [LiA₂]⁻ triplets) to compute L++, L−−, L+−, σ and t+ and compares them with exact values and with the apparent t+ from self-diffusion. Part 2 runs real momentum-conserving MD of a binary mixture to show that the barycentric t+ is fixed by the masses. Run <code>pip install numpy</code>, then <code>python onsager.py</code> (a few seconds).',
    src: `import numpy as np, time
rng = np.random.default_rng(6)
# Reduced units: kT = 1, e = 1, volume V = 1 (we report V kT L_ij, so V and kT drop out).
# ---------- the Onsager estimator ----------
def onsager(R_plus, R_minus, h, lags):
    """L_ij = (1/6) d/dt < dR_i(t) . dR_j(t) >, where R_i(t) = sum of unwrapped positions of all ions of species i."""
    t = lags * h; out = {}
    for name, A, B in [("++", R_plus, R_plus), ("--", R_minus, R_minus), ("+-", R_plus, R_minus)]:
        c = np.array([((A[k:] - A[:-k]) * (B[k:] - B[:-k])).sum(-1).mean() for k in lags])
        out[name] = np.polyfit(t, c, 1)[0] / 6
    return out

def transference(L, zp=1, zm=-1):
    """Cation transference number from Onsager coefficients (frame of the reference species, here the solvent)."""
    num = zp**2 * L["++"] + zp * zm * L["+-"]
    den = zp**2 * L["++"] + zm**2 * L["--"] + 2 * zp * zm * L["+-"]
    return num / den, den                                              # den = sigma (in units e^2 / V kT)

# ---------- Part 1: a cluster model of a salt in a solvent ----------
# cluster type: (number of clusters, Li per cluster, anions per cluster, diffusion coefficient)
# Li+ (fast), A- (slower), neutral pair LiA, negative triplet LiA2-, positive triplet Li2A+.  Clusters do not exchange here.
def cluster_run(spec, T_run=5000.0, dt=0.1):
    """R+(t), R-(t) for non-exchanging Brownian clusters. The sum of n independent Brownian paths with diffusion D
    is itself a Brownian path with diffusion n*D, so we generate one summed path per cluster type (exact, and fast)."""
    steps = int(T_run / dt); Rp = np.zeros((steps, 3)); Rm = np.zeros((steps, 3))
    for (n, a, b, D) in spec:
        S = rng.normal(0, np.sqrt(2 * n * D * dt), (steps, 3)).cumsum(0)   # sum over all clusters of this type
        Rp += a * S; Rm += b * S
    return Rp, Rm, dt

def exact(spec):
    """For non-exchanging Brownian clusters: V kT L_ij = sum over cluster types of n * (count_i) * (count_j) * D."""
    L = {"++": sum(n * a * a * D for n, a, b, D in spec), "--": sum(n * b * b * D for n, a, b, D in spec), "+-": sum(n * a * b * D for n, a, b, D in spec)}
    Np, Nm = sum(n * a for n, a, b, D in spec), sum(n * b for n, a, b, D in spec)
    Dp = sum(n * a * D for n, a, b, D in spec) / Np; Dm = sum(n * b * D for n, a, b, D in spec) / Nm   # self diffusion of each species
    return L, Dp, Dm, Np, Nm

cases = {
    "all free":                [(300, 1, 0, 1.2), (300, 0, 1, 0.8)],
    "40% in pairs":            [(180, 1, 0, 1.2), (180, 0, 1, 0.8), (120, 1, 1, 0.5)],
    "slow Li+, mobile LiA2-":  [(200, 1, 0, 0.3), (50, 1, 1, 0.4), (200, 1, 2, 0.5)],   # electroneutral: 200 Li+ = 200 LiA2-
}
print("Part 1: Li salt in a solvent, cluster model (clusters do not exchange; run of 5000 time units; L in units of 1/(V kT))")
print("   case                     L++      L--      L+-     t+ (Onsager)  exact   apparent D+/(D+ + D-)   sigma/sigma_NE")
lags = np.arange(10, 51, 5)                                             # lag times 1 to 5
for name, spec in cases.items():
    Rp, Rm, h = cluster_run(spec)
    L = onsager(Rp, Rm, h, lags); tp, sig = transference(L)
    Le, Dp, Dm, Np, Nm = exact(spec); tpe, sige = transference(Le)
    t_app = Np * Dp / (Np * Dp + Nm * Dm)                               # what pulsed-field-gradient NMR suggests
    sig_ne = Np * Dp + Nm * Dm
    print(f"   {name:23s} {L['++']:7.1f}  {L['--']:7.1f}  {L['+-']:7.1f}    {tp:6.2f}      {tpe:5.2f}       {t_app:5.2f}                  {sig / sig_ne:5.2f}")
print("   neutral pairs lower sigma but leave t+ unchanged; Li+ inside mobile negative clusters moves the WRONG way in a field -> t+ < 0")

# ---------- Part 2: a pure molten salt: momentum conservation fixes t+ in the barycentric frame ----------
print("\\nPart 2: pure binary 'molten salt' with masses m+ = 1, m- = 3 (momentum-conserving MD, LJ mixture)")
def lj_mixture(N=108, rho=0.8, T=1.5, mp=1.0, mm=3.0, nequil=500, nprod=4000, dt=0.004, every=10):
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(3) for j in range(3) for k in range(3)])
    x = ((g[:, None, :] + b[None]).reshape(-1, 3)) * a; L = 3 * a
    m = np.where(np.arange(N) % 2 == 0, mp, mm)[:, None]; plus = np.arange(N) % 2 == 0
    v = rng.normal(0, 1, x.shape) * np.sqrt(T / m); v -= (m * v).sum(0) / m.sum()     # zero total momentum
    def forces(x):
        d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L); r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf)
        msk = r2 < 6.25; ir6 = np.where(msk, 1 / r2**3, 0.0); f = np.where(msk, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
        return (f[:, :, None] * d).sum(1)
    f = forces(x); u = x.copy(); Rp, Rm = [], []
    for s in range(nequil + nprod):
        v += 0.5 * dt * f / m; x = (x + dt * v) % L; u += dt * v; f = forces(x); v += 0.5 * dt * f / m
        if s >= nequil and s % every == 0: Rp.append(u[plus].sum(0)); Rm.append(u[~plus].sum(0))
    return np.array(Rp), np.array(Rm), dt * every
t0 = time.time(); Rp, Rm, h = lj_mixture()
L = onsager(Rp, Rm, h, np.arange(5, 31, 5)); tp, _ = transference(L)
print(f"   m+ R+ + m- R- stays constant (total momentum is zero):  |m+ dR+ + m- dR-| at the end = {np.abs(1 * (Rp[-1] - Rp[0]) + 3 * (Rm[-1] - Rm[0])).max():.1e}")
print(f"   L++ = {L['++']:.3f}, L-- = {L['--']:.3f}, L+- = {L['+-']:.3f};  m+^2 L++ : m-^2 L-- : -m+ m- L+- = "
      f"{L['++']:.3f} : {9 * L['--']:.3f} : {-3 * L['+-']:.3f}  (all equal)")
print(f"   t+ = {tp:.3f}  vs  m- / (m+ + m-) = {3 / 4:.3f}: in the barycentric frame a binary salt's t+ is fixed by the masses ({time.time() - t0:.0f} s)")
`,
  },

  maths: () => `
<h4>1. Onsager coefficients from equilibrium fluctuations</h4>
${tex(String.raw`L^{ij} = \lim_{t\to\infty}\frac{1}{6Vk_BT}\frac{d}{dt}\big\langle \Delta\mathbf R_i(t)\cdot\Delta\mathbf R_j(t)\big\rangle, \qquad \Delta\mathbf R_i(t) = \sum_{\alpha\in i}\big[\mathbf r_\alpha(t) - \mathbf r_\alpha(0)\big]`, true)}
<p>They relate fluxes to driving forces: J<sub>i</sub> = −Σ<sub>j</sub> L<sup>ij</sup>∇μ̃<sub>j</sub> (μ̃: electrochemical potential), and are symmetric, L<sup>ij</sup> = L<sup>ji</sup> (Onsager reciprocity).</p>
<h4>2. Self and distinct parts</h4>
${tex(String.raw`L^{ii} = \underbrace{\frac{N_iD_i}{Vk_BT}}_{\text{self}} + L^{ii}_{\text{distinct}}, \qquad L^{+-} = L^{+-}_{\text{distinct}}`, true)}
<h4>3. Conductivity and transference number</h4>
${tex(String.raw`\sigma = e^2\sum_{ij} z_iz_jL^{ij}, \qquad t_+ = \frac{e^2\big(z_+^2L^{++} + z_+z_-L^{+-}\big)}{\sigma}`, true)}
<h4>4. Apparent (Nernst–Einstein) transference number</h4>
${tex(String.raw`t_+^{\text{app}} = \frac{N_+D_+}{N_+D_+ + N_-D_-} \quad\big(= D_+/(D_+ + D_-) \text{ for a 1:1 salt}\big)`, true)}
<h4>5. Binary molten salt in the barycentric frame</h4>
${tex(String.raw`m_+\Delta\mathbf R_+ + m_-\Delta\mathbf R_- = 0 \;\Rightarrow\; m_+^2L^{++} = m_-^2L^{--} = -m_+m_-L^{+-} \;\Rightarrow\; t_+ = \frac{m_-}{m_+ + m_-}`, true)}`,

  quiz: [
    {
      q: 'Adding many neutral LiA ion pairs to an electrolyte mainly…',
      options: ['raises the conductivity', 'lowers the conductivity but hardly changes t+', 'makes t+ negative', 'has no effect on any transport property'],
      answer: 1,
      why: 'Neutral pairs diffuse but carry no charge, so they reduce σ (and raise L+−) without changing which charged species carry the current.',
    },
    {
      q: 'How can the Li⁺ transference number become negative?',
      options: ['If Li⁺ is heavier than the anion', 'If much of the Li⁺ moves inside mobile negatively charged clusters such as [LiA₂]⁻', 'If the temperature is too low', 'It cannot; t+ is always between 0 and 1'],
      answer: 1,
      why: 'In a field, negative clusters drift towards the positive electrode and take their Li⁺ with them. If this outweighs the Li⁺ carried by positive species, the net Li⁺ flux points backwards and t+ < 0.',
    },
    {
      q: 'In a pure binary molten salt with cation mass 7 and anion mass 280 (in g/mol), what is the barycentric transference number of the cation?',
      options: ['0.024', '0.5', '0.976', 'It depends on the diffusion coefficients'],
      answer: 2,
      why: 'Momentum conservation gives t+ = m−/(m+ + m−) = 280/287 ≈ 0.976 in the barycentric frame, independent of the dynamics. This is why the reference frame must be stated.',
    },
  ],
};
