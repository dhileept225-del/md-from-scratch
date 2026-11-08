import { tex } from '../tex.js';

export default {
  id: 'barostats',
  summary: 'Pressure in MD comes from the virial: atoms hitting the walls plus the forces between them. It is very noisy. A barostat lets the box breathe so the pressure is right on average. Berendsen gets the density right but the fluctuations wrong; C-rescale and Parrinello–Rahman/MTTK get both right.',
  minutes: 15,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/barostats.mp4', poster: 'videos/barostats.jpg' },
  sim: 'barostats',
  learn: [
    'How pressure is computed from the virial, and why the instantaneous value is so noisy',
    'Berendsen vs stochastic cell rescaling (C-rescale) vs Parrinello–Rahman / MTTK',
    'Choosing τ<sub>P</sub>, the compressibility, and isotropic vs semi-isotropic vs anisotropic coupling',
  ],
  simHint: 'Start with None (NVT) and look at how much the instantaneous pressure jumps. Pick C-rescale and change the target pressure: the box adjusts. Now compare the spread of box area for Berendsen at τ = 0.2 and τ = 2, then for C-rescale. Finally try Piston with τ = 2 and change the target: the box overshoots and oscillates.',

  idea: `
<p><b>Where pressure comes from.</b> Pressure has two parts. The <i>kinetic</i> part is the momentum carried by atoms moving through any surface (the ideal-gas term Nk<sub>B</sub>T/V). The <i>virial</i> part comes from forces acting across that surface: repulsion pushes the box outwards, attraction pulls it in. MD computes both every step from the velocities and from Σ <b>r</b>·<b>F</b> over pairs, a quantity we already have from the force loop.</p>
<p><b>It is noisy.</b> The instantaneous pressure fluctuates strongly, much more than the temperature. For a few thousand water molecules it swings by hundreds of bar from step to step, even when the average is exactly 1 bar. Only long averages mean anything. And like the energy (Lesson 2.5), a cutoff needs a <b>tail correction</b> for the pressure.</p>
<p><b>Barostats.</b> To simulate at constant pressure (NPT, Lesson 3.2) the box volume becomes a variable. All positions are scaled with the box, so the atoms are not squeezed against one wall.</p>
<ul>
<li><b>Berendsen:</b> rescales the box so the pressure relaxes to P₀ with time constant τ<sub>P</sub>. Robust and fast for equilibration, and the average density is fine, but the volume fluctuations are wrong (they depend on τ<sub>P</sub>).</li>
<li><b>C-rescale</b> (stochastic cell rescaling, Bernetti and Bussi, 2020): Berendsen plus a correctly chosen random term, exactly as V-rescale is for temperature. Correct NPT fluctuations, still first-order, so no oscillations. A very good default.</li>
<li><b>Parrinello–Rahman / MTTK:</b> the box has a "mass" and its own equation of motion, like a piston. Correct ensemble and smooth dynamics, but it oscillates if started far from equilibrium; equilibrate with Berendsen or C-rescale first.</li>
</ul>
<p><b>How the box may deform.</b> <i>Isotropic</i>: all sides scale together (liquids). <i>Semi-isotropic</i>: the xy plane and z are coupled separately (membranes, interfaces). <i>Anisotropic / fully flexible</i>: each box vector can change, including shape (crystals, phase transitions).</p>
<div class="callout"><b>Key idea.</b> P = kinetic + virial, averaged over a long run. Equilibrate with Berendsen or C-rescale; produce with C-rescale or Parrinello–Rahman/MTTK. Choose the coupling (isotropic, semi-isotropic, anisotropic) to match the physics.</div>`,

  algorithm: `
<p>One NPT step (thermostat inside the MD step, barostat after it):</p>
<pre class="pseudo">MD step with a thermostat
P = (2 K / 3 + W / 3) / V              # W = Σ r·F over pairs

# Berendsen
de = -(beta / tau_P) * (P0 - P) * dt

# C-rescale (adds the right noise; also scale v)
de = -(beta / tau_P) * (P0 - P - kB T / V) * dt
     + sqrt(2 kB T beta dt / (V tau_P)) * gaussian()

mu = exp(de / 3)                       # de = change of ln V
positions *= mu;  box *= mu
velocities /= mu                       # C-rescale only</pre>
<ol class="steps">
  <li>β is the (estimated) isothermal compressibility; for water use 4.5 × 10⁻⁵ bar⁻¹. A wrong β only changes the effective coupling speed.</li>
  <li>Typical τ<sub>P</sub>: 1–2 ps for Berendsen and C-rescale, 2–5 ps for Parrinello–Rahman. τ<sub>P</sub> should be several times τ<sub>T</sub>.</li>
  <li>Discard the part of the run where the volume is still drifting. Then average.</li>
  <li>GROMACS: <code>pcoupl = C-rescale</code> or <code>Parrinello-Rahman</code>, <code>pcoupltype = isotropic / semiisotropic / anisotropic</code>, <code>tau-p</code>, <code>compressibility</code>, <code>ref-p</code>. LAMMPS: <code>fix npt temp ... iso / aniso / tri ...</code> (MTTK) or <code>fix press/berendsen</code>.</li>
</ol>`,

  code: {
    note: 'Part 1: the instantaneous virial pressure in an NVT run. Part 2: NPT with a Langevin thermostat and either Berendsen or C-rescale, each at two coupling times. All give the same average pressure, but only C-rescale gives volume fluctuations that do not depend on τ<sub>P</sub>. Run <code>pip install numpy</code>, then <code>python barostats.py</code> (about 25 s).',
    src: `import numpy as np
rng = np.random.default_rng(4)

def fcc(n, rho):
    N = 4 * n**3; L = (N / rho) ** (1 / 3); a = L / n
    b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    c = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return (c[:, None] + b[None]).reshape(-1, 3) * a, L

def forces(pos, L, rc=2.5):
    """LJ forces and virial W = sum over pairs of r.F (truncated at rc, no tail correction)."""
    d = pos[:, None] - pos[None]; d -= L * np.round(d / L)
    r2 = (d ** 2).sum(-1); np.fill_diagonal(r2, np.inf)
    inv6 = np.where(r2 < rc * rc, 1 / r2**3, 0.0)
    rF = 24 * (2 * inv6**2 - inv6)
    return ((rF / r2)[..., None] * d).sum(1), 0.5 * rF.sum()

def run(barostat, P0, T0=1.5, steps=10000, dt=0.004, gamma=1.0, beta=0.5, tau_p=0.25, rho=0.6):
    """Langevin (BAOAB) thermostat + an isotropic barostat applied once per step."""
    pos, L = fcc(3, rho); N = len(pos); V = L**3
    v = rng.normal(0, np.sqrt(T0), pos.shape); F, W = forces(pos, L)
    c1 = np.exp(-gamma * dt); out = []
    for s in range(steps):
        v += 0.5 * dt * F; pos += 0.5 * dt * v
        v = c1 * v + np.sqrt((1 - c1**2) * T0) * rng.normal(size=v.shape)
        pos += 0.5 * dt * v; pos %= L
        F, W = forces(pos, L); v += 0.5 * dt * F
        K = 0.5 * (v**2).sum(); P = (2 * K / 3 + W / 3) / V          # virial pressure
        if barostat == "berendsen":                                    # deterministic: relax P to P0
            de = -beta / tau_p * (P0 - P) * dt
        elif barostat == "c-rescale":                                  # Bernetti-Bussi: + noise and kT/V term
            de = (-beta / tau_p * (P0 - P - T0 / V) * dt
                  + np.sqrt(2 * T0 * beta / (V * tau_p) * dt) * rng.normal())
        else:
            de = 0.0
        mu = np.exp(de / 3); pos *= mu; L *= mu; V = L**3
        if barostat == "c-rescale": v /= mu                            # velocities scaled inversely
        out.append((P, V))
    return np.array(out)

# ---- Part 1: pressure in NVT is a very noisy quantity ----
o = run(None, 0, steps=3000)[1000:]
print("Part 1: NVT, 108 LJ atoms, rho = 0.6, T = 1.5")
print(f"  <P> = {o[:, 0].mean():.3f}   std(P) = {o[:, 0].std():.3f}   (the instantaneous P swings by about 25% of its mean)")

# ---- Part 2: NPT with two barostats and two coupling times ----
print("\\nPart 2: NPT at P0 = 1.0, T = 1.5 (10000 steps each, first 1500 discarded)")
print("  barostat    tau_p   <P>     <V>     std(V)")
for b in ["berendsen", "c-rescale"]:
    for tp in [0.2, 1.0]:
        o = run(b, 1.0, steps=10000, tau_p=tp)[1500:]
        print(f"  {b:10s}  {tp:4.1f}   {o[:, 0].mean():.3f}   {o[:, 1].mean():6.1f}   {o[:, 1].std():5.2f}")
print("All runs agree on <P> and (within noise) on <V>. A physical fluctuation cannot depend on tau_p:")
print("C-rescale's std(V) stays the same (within noise), Berendsen's shrinks as the coupling gets weaker.")
`,
  },

  maths: () => `
<h4>1. Virial pressure</h4>
${tex(String.raw`P = \frac{N k_B T}{V} + \frac{1}{3V}\Big\langle \sum_{i<j} \mathbf r_{ij}\cdot\mathbf F_{ij} \Big\rangle \qquad (\text{2D: } \tfrac{1}{2A})`, true)}
<p>For a pair potential, r·F = −r u′(r). Use the kinetic temperature from the same step for the first term.</p>
<h4>2. Tail correction (LJ, cutoff r<sub>c</sub>)</h4>
${tex(String.raw`P_{\text{tail}} = \frac{16}{3}\pi\rho^2\varepsilon\sigma^3\left[\frac{2}{3}\Big(\frac{\sigma}{r_c}\Big)^9 - \Big(\frac{\sigma}{r_c}\Big)^3\right]`, true)}
<h4>3. Volume fluctuations and compressibility</h4>
${tex(String.raw`\kappa_T = -\frac{1}{V}\left(\frac{\partial V}{\partial P}\right)_T = \frac{\langle \delta V^2\rangle}{k_B T\,\langle V\rangle}`, true)}
<p>This only works with a barostat that samples NPT correctly (C-rescale, Parrinello–Rahman, MTTK), not with Berendsen.</p>
<h4>4. Berendsen and C-rescale (isotropic, ε = ln V)</h4>
${tex(String.raw`d\varepsilon = -\frac{\beta_T}{\tau_P}\Big(P_0 - P - \frac{k_BT}{V}\Big)dt + \sqrt{\frac{2k_BT\beta_T}{V\tau_P}}\,dW`, true)}
<p>Drop the k<sub>B</sub>T/V term and the noise and you get Berendsen.</p>
<h4>5. Parrinello–Rahman (box as a dynamical variable)</h4>
${tex(String.raw`W\,\ddot{\mathbf h} = V\,(\mathbf P - P_0\mathbf 1)\,(\mathbf h^{T})^{-1}`, true)}
<p><b>h</b> is the matrix of box vectors and W the "piston mass". Being second order, the box oscillates with a period set by W (or τ<sub>P</sub>).</p>`,

  quiz: [
    {
      q: 'The instantaneous pressure in your NPT simulation of water jumps between −300 and +300 bar. The target is 1 bar. What does this mean?',
      options: ['The barostat is broken', 'Normal: instantaneous pressure is very noisy; check the long-time average', 'The time step is too large', 'The system is boiling'],
      answer: 1,
      why: 'Pressure fluctuations are large, especially in small systems. Only the average over a long run should match the target.',
    },
    {
      q: 'You want the compressibility from volume fluctuations. Which barostat must you NOT use?',
      options: ['C-rescale', 'Parrinello–Rahman', 'MTTK', 'Berendsen'],
      answer: 3,
      why: 'Berendsen gives correct average densities but wrong volume fluctuations, which even depend on τ_P.',
    },
    {
      q: 'Which pressure coupling suits a lipid bilayer in water?',
      options: ['Isotropic', 'Semi-isotropic (xy together, z separately)', 'No barostat', 'Only Berendsen'],
      answer: 1,
      why: 'A membrane is different in the plane and along the normal. Semi-isotropic coupling lets the area and the thickness adjust independently.',
    },
  ],
};
