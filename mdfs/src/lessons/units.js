import { tex } from '../tex.js';

export default {
  id: 'units',
  summary: 'Atoms are tiny and fast, so MD codes use special units: ångströms, femtoseconds, kcal/mol. Different codes pick different ones, which makes unit mix-ups one of the most common bugs. Reduced units avoid the problem entirely.',
  minutes: 12,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/units.mp4', poster: 'videos/units.jpg' },
  sim: 'units',
  learn: [
    'The unit systems used by LAMMPS, GROMACS and others',
    'The conversion traps: ×4.184, ×10, ×1000',
    'Reduced Lennard-Jones units, and how to convert results back to real units',
  ],
  simHint: 'Switch on "Show in reduced units" and watch four different gases fall onto one curve. Then pick argon and move T* and ρ* to find liquid argon. Switch to xenon: the same reduced state point is a different real temperature and density.',

  idea: `
<p>In SI units, an atom is 10<sup>−10</sup> m across, vibrates every 10<sup>−14</sup> s, and has an interaction energy of 10<sup>−21</sup> J. Numbers like that are awkward, so MD codes use units that fit atoms.</p>
<ul>
  <li><b>Length:</b> ångström (1 Å = 0.1 nm) or nanometre.</li>
  <li><b>Time:</b> femtosecond (1 fs = 10<sup>−15</sup> s) or picosecond (1 ps = 1000 fs).</li>
  <li><b>Energy:</b> kcal/mol, kJ/mol (1 kcal = 4.184 kJ) or eV.</li>
  <li><b>Mass:</b> g/mol (atomic mass units).</li>
</ul>
<p>Each code chooses its own. LAMMPS <code>units real</code> uses Å, fs and kcal/mol. <code>units metal</code> uses Å, ps and eV. GROMACS uses nm, ps and kJ/mol. Copy a parameter from one world to the other without converting it, and your forces are off by 4.184 or your distances by 10. The code will not complain. It will simply simulate the wrong system.</p>
<p><b>Reduced units</b> sidestep all of this. For a Lennard-Jones system we set σ = ε = m = 1 and measure everything relative to them. Then argon, krypton and xenon all become the <i>same</i> simulation. Only at the end do we convert back with each gas's own σ, ε and m. This idea is called the <b>law of corresponding states</b>.</p>
<div class="callout"><b>Key idea.</b> Write down the unit of every input before you run. In your own code, use reduced units and convert only when you report results.</div>`,

  algorithm: `
<p>Converting a reduced-unit simulation to a real substance:</p>
<ol class="steps">
  <li>Look up σ (Å), ε/k<sub>B</sub> (K) and m (g/mol) for your substance.</li>
  <li>Compute the time unit τ = σ√(m/ε). For argon this is about 2.16 ps.</li>
  <li>Multiply each reduced result by its unit:
    <ul><li>T = T*·(ε/k<sub>B</sub>)</li><li>t = t*·τ</li><li>ρ = ρ*·m/σ³</li><li>P = P*·ε/σ³</li><li>E = E*·ε</li></ul></li>
</ol>
<pre class="pseudo">sigma, eps_K, mass = 3.405, 119.8, 39.948      # argon
tau = sigma*1e-10 * sqrt(mass*amu / (eps_K*kB))
T_real  = T_star * eps_K                       # kelvin
dt_real = dt_star * tau                        # seconds</pre>
<p><b>Unit check before every run:</b> list each input (cutoff, ε, σ, charges, Δt, temperature) next to the unit your code expects. Most "my simulation explodes" questions end here.</p>`,

  code: {
    note: 'Computes the real-unit value of "1" in reduced units for four noble gases, converts a liquid state point to argon, and shows the kcal/kJ trap. Run <code>pip install numpy</code>, then <code>python units.py</code>.',
    src: `import numpy as np

kB, amu, NA = 1.380649e-23, 1.66053907e-27, 6.02214076e23

# typical Lennard-Jones parameters: sigma (A), epsilon/kB (K), mass (amu)
gases = {"Ne": (2.75, 36.0, 20.18), "Ar": (3.405, 119.8, 39.948),
         "Kr": (3.60, 171.0, 83.80), "Xe": (4.10, 221.0, 131.29)}

def lj_units(sigma_A, eps_K, mass_amu):
    """Real-unit value of '1' for each reduced quantity."""
    sigma, eps, m = sigma_A * 1e-10, eps_K * kB, mass_amu * amu
    tau = sigma * np.sqrt(m / eps)                 # time unit
    return {"tau_ps": tau * 1e12,
            "T_K": eps_K,                          # T* = 1  ->  eps/kB
            "rho_gcc": m / sigma**3 / 1000,        # rho* = 1 -> m/sigma^3
            "P_bar": eps / sigma**3 / 1e5,
            "E_kJmol": eps * NA / 1000}

print(f"{'gas':4s}{'tau (ps)':>10s}{'T (K)':>9s}{'rho (g/cm3)':>13s}{'P (bar)':>10s}")
for name, p in gases.items():
    u = lj_units(*p)
    print(f"{name:4s}{u['tau_ps']:10.3f}{u['T_K']:9.1f}{u['rho_gcc']:13.3f}{u['P_bar']:10.0f}")

# a liquid state point in reduced units, converted to argon
Tstar, rhostar, dtstar = 0.72, 0.84, 0.005
u = lj_units(*gases["Ar"])
print(f"\\nargon at T*={Tstar}, rho*={rhostar}:  T = {Tstar*u['T_K']:.0f} K, "
      f"density = {rhostar*u['rho_gcc']:.2f} g/cm3, dt = {dtstar*u['tau_ps']*1000:.1f} fs")

# the classic trap: the same number in kcal/mol vs kJ/mol
eps_kcal = 0.238                                   # argon epsilon in kcal/mol
print(f"\\nargon epsilon: {eps_kcal} kcal/mol = {eps_kcal*4.184:.3f} kJ/mol "
      f"(typing 0.238 into a kJ/mol code makes it 4.184x too weak)")`,
  },

  maths: () => `
<h4>1. Reduced (dimensionless) variables</h4>
${tex(String.raw`r^* = \frac{r}{\sigma},\quad E^* = \frac{E}{\varepsilon},\quad t^* = \frac{t}{\tau},\quad T^* = \frac{k_B T}{\varepsilon},\quad \rho^* = \rho\,\sigma^3,\quad P^* = \frac{P\sigma^3}{\varepsilon},\quad F^* = \frac{F\sigma}{\varepsilon}`, true)}
<h4>2. Where τ comes from</h4>
<p>The only combination of σ (length), ε (energy) and m (mass) with units of time is</p>
${tex(String.raw`\tau = \sigma\sqrt{\frac{m}{\varepsilon}} \qquad \text{(argon: } \tau \approx 2.16\ \text{ps)}`, true)}
<h4>3. Corresponding states</h4>
<p>In reduced units, Newton's equations for any LJ system become identical:</p>
${tex(String.raw`\frac{d^2\mathbf{r}_i^*}{dt^{*2}} = \sum_{j} 24\left[2\,(r_{ij}^*)^{-14} - (r_{ij}^*)^{-8}\right]\mathbf{r}_{ij}^*`, true)}
<p>So one simulation at (T*, ρ*) describes every noble gas at its own real temperature and density.</p>
<h4>4. Conversion factors to keep nearby</h4>
${tex(String.raw`1\ \tfrac{\text{kcal}}{\text{mol}} = 4.184\ \tfrac{\text{kJ}}{\text{mol}} = 0.04336\ \text{eV} = 503.2\ \text{K}\cdot k_B, \qquad 1\ \text{eV} = 96.485\ \tfrac{\text{kJ}}{\text{mol}}`, true)}
${tex(String.raw`1\ \text{nm} = 10\ \text{Å}, \qquad 1\ \text{ps} = 1000\ \text{fs}, \qquad k_B = 0.0019872\ \tfrac{\text{kcal}}{\text{mol K}} = 0.0083145\ \tfrac{\text{kJ}}{\text{mol K}}`, true)}
<p class="muted">Reduced units only work this neatly when there is a single σ, ε and m. Real force fields with many atom types use real units, which is exactly why the unit check matters.</p>`,

  quiz: [
    {
      q: 'Which units does GROMACS use for length, time and energy?',
      options: ['Å, fs, kcal/mol', 'nm, ps, kJ/mol', 'm, s, J', 'Å, ps, eV'],
      answer: 1,
      why: 'GROMACS uses nm, ps and kJ/mol. LAMMPS "real" uses Å, fs and kcal/mol, and LAMMPS "metal" uses Å, ps and eV.',
    },
    {
      q: 'An argon simulation in reduced units runs at T* = 1.0. What is the real temperature? (ε/kB = 119.8 K)',
      options: ['1 K', 'About 120 K', 'About 300 K', '273 K'],
      answer: 1,
      why: 'T = T*·ε/kB = 1.0 × 119.8 K ≈ 120 K.',
    },
    {
      q: 'You type an ε given in kcal/mol into a code that expects kJ/mol. What happens?',
      options: ['The code stops with an error', 'The interaction is 4.184 times too weak', 'The interaction is 4.184 times too strong', 'Nothing: units do not matter'],
      answer: 1,
      why: 'The number is read as kJ/mol, and 1 kJ is smaller than 1 kcal. So the interaction is 4.184 times too weak, and the code runs without any warning.',
    },
  ],
};
