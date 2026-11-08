import { tex } from '../tex.js';

export default {
  id: 'real-codes',
  summary: 'Everything in your toy code has a one-line counterpart in LAMMPS and GROMACS: the lattice, the pair potential, the neighbour list, the integrator and the thermostat. This lesson runs the same Lennard-Jones liquid in our NumPy code and in LAMMPS and gets identical numbers, then shows how to translate between unit systems and how a real simulation workflow is organized.',
  minutes: 15,
  // VIDEO: paste your YouTube video ID here (the part after "watch?v="). Empty = local MP4.
  video: { youtube: '', mp4: 'videos/real-codes.mp4', poster: 'videos/real-codes.jpg' },
  sim: 'realcodes',
  learn: [
    'How each part of a toy MD code maps onto LAMMPS and GROMACS input commands',
    'Unit systems (LAMMPS lj, real, metal; GROMACS nm, ps, kJ/mol) and how to convert between them',
    'The standard workflow (build, minimize, NVT, NPT, production, analysis) and basic sanity checks',
  ],
  simHint: 'Pick a gas and a code, then move the sliders: the input file and the conversions update live. Compare the same state point in "LAMMPS lj" and "real" units, and check that T*, ρ* and Δt* become sensible kelvin, g/cm³ and femtoseconds. Copy an input and run it: lmp -in in.file (LAMMPS) or gmx grompp + gmx mdrun (GROMACS). (The LAMMPS inputs were tested: lj, real and metal give the same reduced energy. Note that NVE started from a lattice cools to about half of T*, because half the kinetic energy turns into potential energy as the crystal melts; equilibrate with NVT first.)',

  idea: `
<p>Research codes such as <b>LAMMPS</b> and <b>GROMACS</b> do exactly what our toy code does, only much faster and with many more options. Every piece maps onto an input command:</p>
<ul>
<li><b>Starting structure</b>: <code>lattice</code> + <code>create_atoms</code> or <code>read_data</code> (LAMMPS); a <code>.gro</code>/<code>.pdb</code> file (GROMACS).</li>
<li><b>Force field</b>: <code>pair_style lj/cut 2.5</code> + <code>pair_coeff</code> (LAMMPS); <code>[ atomtypes ]</code> in a <code>.top</code> file plus cutoffs in the <code>.mdp</code> (GROMACS).</li>
<li><b>Neighbour list</b> (Lesson 2.6): <code>neighbor 0.3 bin</code>; GROMACS sets its Verlet buffer automatically.</li>
<li><b>Integrator</b>: <code>fix nve</code> (velocity Verlet); GROMACS uses leap-frog (<code>integrator = md</code>), equivalent for most purposes.</li>
<li><b>Thermostat and barostat</b> (Chapter 3): <code>fix nvt</code>/<code>fix npt</code> (Nosé–Hoover); <code>tcoupl = v-rescale</code>, <code>pcoupl = C-rescale</code>.</li>
</ul>
<p><b>Same model, same numbers.</b> If you give LAMMPS exactly the same positions, velocities, potential, cutoff and time step as our NumPy code, you get the same energies to 10 significant figures (code tab, Part 2). That is a powerful way to test your own code, or to understand what a command really does.</p>
<p><b>Units are the classic trap.</b> LAMMPS has several unit systems: <code>lj</code> (reduced: σ = ε = m = 1), <code>real</code> (Å, fs, kcal/mol, atm) and <code>metal</code> (Å, ps, eV, bar), among others. GROMACS always uses nm, ps, kJ/mol and bar. For argon (σ = 3.405 Å, ε/k<sub>B</sub> = 119.8 K), τ = 2.156 ps, so Δt* = 0.005 is 10.8 fs, T* = 0.722 is 86.5 K and ρ* = 0.8442 is 1.419 g/cm³ (Part 3).</p>
<p><b>A real workflow</b> has stages: build the system (Packmol, <code>gmx insert-molecules</code>, moltemplate, fftool), minimize the energy, equilibrate in NVT, then NPT to find the density, then a production run, then analysis (Chapter 6). At each stage check the basics: energy conservation in NVE, temperature and pressure fluctuating around their targets, a sensible density, and a sensible g(r).</p>
<p><b>Speed.</b> Real codes split the box into domains, one per CPU core (domain decomposition), use neighbour lists and PME (Lesson 5.3), and run on GPUs. That is how they handle millions of atoms where our NumPy code handles hundreds.</p>
<div class="callout"><b>Key idea.</b> A research code is your toy code with better engineering. Map each concept to its command, keep the units straight, start from a working example, and validate against something you understand.</div>`,

  algorithm: `
<pre class="pseudo">1. build      : coordinates + topology/force field   (Packmol, gmx insert-molecules, moltemplate)
2. minimize   : LAMMPS minimize / GROMACS integrator = steep   -> remove overlaps
3. NVT        : heat to the target T (short, e.g. 100 ps)
4. NPT        : equilibrate density (check that the box volume has levelled off)
5. production : NVT or NPT, save coordinates (and velocities) at the interval your analysis needs
6. analysis   : g(r), MSD, conductivity ... (Chapter 6) with error bars from blocks or replicas
sanity checks : NVE drift, <T>, <P>, density vs experiment, g(r) shape, no atoms lost</pre>
<ol class="steps">
  <li>LAMMPS: one input script (<code>in.*</code>) runs everything; <code>lmp -in in.lj</code> or <code>mpirun -np 8 lmp -in in.lj</code>. GPU and KOKKOS packages accelerate it.</li>
  <li>GROMACS: <code>gmx grompp</code> combines <code>.mdp</code> + <code>.gro</code> + <code>.top</code> into a <code>.tpr</code>; <code>gmx mdrun</code> runs it, using GPUs automatically when available.</li>
  <li>Keep the input files, code version and random seeds with every data set, so that a run can be reproduced.</li>
</ol>`,

  code: {
    note: 'Part 1 runs 500 NVE steps of our NumPy LJ code. Part 2 writes the same system as a LAMMPS data file and input and runs it through the LAMMPS Python module, which gives identical energies. Part 3 converts to argon in real units and runs LAMMPS again. Part 4 writes GROMACS .gro/.top/.mdp files for the same system (GROMACS itself is not run here). Run <code>pip install numpy lammps</code> (the LAMMPS wheel may also need <code>pip install mpich</code>), then <code>python real_codes.py</code> (a few seconds).',
    src: `import numpy as np, os, time
rng = np.random.default_rng(7)
# The same Lennard-Jones liquid three ways: our NumPy code, LAMMPS (reduced and real units), and GROMACS input files.

# ---------- Part 1: our toy code ----------
def fcc(n, rho):
    a = (4 / rho) ** (1 / 3); b = np.array([[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]])
    g = np.array([[i, j, k] for i in range(n) for j in range(n) for k in range(n)])
    return ((g[:, None, :] + b[None]).reshape(-1, 3) + 0.25) * a, n * a      # +0.25: keep atoms off the box faces

def energy_forces(x, L, rc=2.5):
    d = x[:, None, :] - x[None, :, :]; d -= L * np.round(d / L)
    r2 = (d**2).sum(-1); np.fill_diagonal(r2, np.inf); m = r2 < rc * rc
    ir6 = np.where(m, 1 / r2**3, 0.0)
    pe = 0.5 * np.sum(4 * (ir6**2 - ir6))                                 # plain cut-off, no shift (same as LAMMPS lj/cut)
    f = np.where(m, 24 * (2 * ir6**2 - ir6) / r2, 0.0)
    return pe, (f[:, :, None] * d).sum(1)

x0, L = fcc(4, 0.8442); N = len(x0); T0 = 0.722
v0 = rng.normal(0, np.sqrt(T0), x0.shape); v0 -= v0.mean(0); v0 *= np.sqrt(T0 * (3 * N - 3) / (v0**2).sum())
dt, nsteps = 0.005, 500
x, v = x0.copy(), v0.copy(); pe, f = energy_forces(x, L); toy = {0: pe / N}
t0 = time.time()
for s in range(1, nsteps + 1):
    v += 0.5 * dt * f; x = (x + dt * v) % L; pe, f = energy_forces(x, L); v += 0.5 * dt * f
    if s % 100 == 0: toy[s] = pe / N
print(f"Part 1: our NumPy code, {N} atoms, {nsteps} steps of NVE ({time.time() - t0:.1f} s)")
print("   step:  " + "  ".join(f"{s:>8d}" for s in toy) + "\\n   PE/N:  " + "  ".join(f"{e:8.5f}" for e in toy.values()))

# ---------- Part 2: the same run in LAMMPS (units lj) ----------
def write_data(fname, pos, vel, box, mass, comment):
    with open(fname, "w") as fh:
        fh.write(f"{comment}\\n\\n{len(pos)} atoms\\n1 atom types\\n\\n0 {box:.10f} xlo xhi\\n0 {box:.10f} ylo yhi\\n0 {box:.10f} zlo zhi\\n\\n")
        fh.write(f"Masses\\n\\n1 {mass}\\n\\nAtoms # atomic\\n\\n")
        for i, p in enumerate(pos): fh.write(f"{i + 1} 1 {p[0]:.12f} {p[1]:.12f} {p[2]:.12f}\\n")
        fh.write("\\nVelocities\\n\\n")
        for i, w in enumerate(vel): fh.write(f"{i + 1} {w[0]:.12f} {w[1]:.12f} {w[2]:.12f}\\n")

write_data("lj.data", x0, v0, L, 1.0, "LJ liquid, reduced units")
lammps_lj = """# LJ liquid in reduced units: the same model as our NumPy code
units           lj
atom_style      atomic
boundary        p p p
read_data       lj.data
pair_style      lj/cut 2.5          # plain cut-off at 2.5 sigma, no shift
pair_coeff      1 1 1.0 1.0         # epsilon sigma
neighbor        0.3 bin             # Verlet list with a 0.3 sigma skin (Lesson 2.6)
neigh_modify    every 1 delay 0 check yes
timestep        0.005
fix             1 all nve           # velocity Verlet
thermo_style    custom step temp pe ke etotal press
thermo          100
run             500
"""
open("in.lj", "w").write(lammps_lj)

def run_lammps(script_lines, n_chunks, chunk, key="pe"):
    from lammps import lammps
    lmp = lammps(cmdargs=["-log", "none", "-screen", "none"])
    for line in script_lines: lmp.command(line)
    out = {}; lmp.command("run 0"); out[0] = lmp.get_thermo(key)
    for c in range(1, n_chunks + 1): lmp.command(f"run {chunk} pre no post no"); out[c * chunk] = lmp.get_thermo(key)
    lmp.close(); return out

try:
    import lammps  # noqa: F401
    lines = [l for l in lammps_lj.splitlines() if l.strip() and not l.startswith("#") and not l.startswith("run")]
    lj = run_lammps(lines, 5, 100)                                     # in lj units, thermo pe is per atom
    print("\\nPart 2: LAMMPS, units lj (in.lj + lj.data written)")
    print("   PE/N:  " + "  ".join(f"{e:8.5f}" for e in lj.values()))
    print(f"   largest difference from our code: {max(abs(lj[s] - toy[s]) for s in toy):.1e}  (same model, same start -> same numbers)")
    have_lammps = True
except Exception as err:
    print(f"\\nPart 2: LAMMPS not available ({type(err).__name__}); in.lj and lj.data were written. Run:  lmp -in in.lj"); have_lammps = False

# ---------- Part 3: argon in real units (LAMMPS units real) ----------
sigma, eps_kcal, mass = 3.405, 0.2381, 39.948                         # Angstrom, kcal/mol, g/mol (argon)
eps_J = eps_kcal * 4184 / 6.02214076e23; m_kg = mass * 1e-3 / 6.02214076e23; kB = 1.380649e-23
v_unit = np.sqrt(eps_J / m_kg) * 1e10 / 1e15                           # LJ velocity unit in Angstrom/fs
tau_fs = sigma / v_unit
print(f"\\nPart 3: reduced -> real units for argon (sigma = {sigma} A, epsilon = {eps_kcal} kcal/mol, m = {mass} g/mol)")
print(f"   time unit tau = {tau_fs / 1000:.3f} ps, so dt = 0.005 tau = {0.005 * tau_fs:.2f} fs")
print(f"   T* = {T0} -> {T0 * eps_J / kB:.1f} K;  rho* = 0.8442 -> {0.8442 * mass / 6.02214076e23 / (sigma * 1e-8)**3:.3f} g/cm^3")
write_data("argon.data", x0 * sigma, v0 * v_unit, L * sigma, mass, "liquid argon, real units")
lammps_real = f"""# the same liquid as argon in real units (Angstrom, fs, kcal/mol, K)
units           real
atom_style      atomic
read_data       argon.data
pair_style      lj/cut {2.5 * sigma:.4f}
pair_coeff      1 1 {eps_kcal} {sigma}
neighbor        {0.3 * sigma:.3f} bin
neigh_modify    every 1 delay 0 check yes
timestep        {0.005 * tau_fs:.4f}
fix             1 all nve
thermo_style    custom step temp pe etotal press
thermo          100
run             500
"""
open("in.argon", "w").write(lammps_real)
if have_lammps:
    lines = [l for l in lammps_real.splitlines() if l.strip() and not l.startswith("#") and not l.startswith("run")]
    real = run_lammps(lines, 5, 100)                                   # real units: total PE in kcal/mol
    print("   PE/N in kcal/mol from LAMMPS real units, divided by epsilon:")
    print("   PE/N:  " + "  ".join(f"{e / N / eps_kcal:8.5f}" for e in real.values()))
    print(f"   largest difference from reduced units: {max(abs(real[s] / N / eps_kcal - toy[s]) for s in toy):.1e}")

# ---------- Part 4: GROMACS input files for the same system ----------
nm, kJ = sigma / 10, eps_kcal * 4.184
with open("argon.gro", "w") as fh:
    fh.write("liquid argon\\n" + f"{N:5d}\\n")
    for i, p in enumerate(x0 * nm): fh.write(f"{1 + i:5d}AR      AR{1 + i:5d}{p[0]:8.3f}{p[1]:8.3f}{p[2]:8.3f}\\n")
    fh.write(f"{L * nm:10.5f}{L * nm:10.5f}{L * nm:10.5f}\\n")
open("argon.top", "w").write(f"""[ defaults ]
; nbfunc  comb-rule  gen-pairs
  1       2          no            ; comb-rule 2: sigma/epsilon given directly
[ atomtypes ]
; name  mass     charge  ptype  sigma (nm)  epsilon (kJ/mol)
  AR    {mass}  0.0     A      {nm:.4f}      {kJ:.4f}
[ moleculetype ]
  AR  1
[ atoms ]
  1  AR  1  AR  AR  1  0.0
[ system ]
liquid argon
[ molecules ]
AR  {N}
""")
open("argon.mdp", "w").write(f"""integrator      = md              ; leap-frog (equivalent to velocity Verlet for NVE)
dt              = {0.005 * tau_fs / 1000:.5f}           ; ps
nsteps          = 500
cutoff-scheme   = Verlet
verlet-buffer-tolerance = -1      ; NVE: set the pair-list buffer by hand
rlist           = {2.5 * nm + 0.1:.4f}            ; nm, rvdw + 0.1 nm skin
nstlist         = 10
rvdw            = {2.5 * nm:.4f}            ; nm
vdw-modifier    = None            ; plain cut-off, like lj/cut
DispCorr        = no
coulombtype     = Cut-off
rcoulomb        = {2.5 * nm:.4f}
tcoupl          = no              ; NVE
pcoupl          = no
gen-vel         = yes
gen-temp        = {T0 * eps_J / kB:.1f}
nstenergy       = 100
""")
print(f"\\nPart 4: GROMACS files written (argon.gro/top/mdp): sigma = {nm:.4f} nm, epsilon = {kJ:.4f} kJ/mol, dt = {0.005 * tau_fs / 1000:.5f} ps")
print("   run with:  gmx grompp -f argon.mdp -c argon.gro -p argon.top -o argon.tpr && gmx mdrun -deffnm argon")
print("   (GROMACS draws new velocities here, so energies agree on average, not step by step; GROMACS was not run in this script)")
`,
  },

  maths: () => `
<h4>1. Reduced (LJ) units</h4>
${tex(String.raw`r^* = \frac{r}{\sigma}, \quad E^* = \frac{E}{\varepsilon}, \quad T^* = \frac{k_BT}{\varepsilon}, \quad \rho^* = \rho\sigma^3, \quad P^* = \frac{P\sigma^3}{\varepsilon}, \quad t^* = \frac{t}{\tau}, \;\; \tau = \sigma\sqrt{\frac{m}{\varepsilon}}`, true)}
<h4>2. Argon numbers</h4>
${tex(String.raw`\tau = 3.405\,\text{Å}\sqrt{\frac{39.948\ \text{g/mol}}{0.2381\ \text{kcal/mol}}} = 2.156\ \text{ps}, \qquad T = 0.722 \times 119.8\ \text{K} = 86.5\ \text{K}`, true)}
<h4>3. Useful conversions</h4>
${tex(String.raw`1\ \text{kcal/mol} = 4.184\ \text{kJ/mol} = 0.04336\ \text{eV}, \qquad 1\ \text{atm} = 1.01325\ \text{bar}, \qquad k_B = 0.0019872\ \text{kcal/(mol K)}`, true)}
<h4>4. Lattice constant from density (fcc, 4 atoms per cell)</h4>
${tex(String.raw`a = \left(\frac{4}{\rho}\right)^{1/3}, \qquad N = 4n^3 \text{ atoms in } n\times n\times n \text{ cells}`, true)}`,

  quiz: [
    {
      q: 'In LAMMPS with "units real", what are the units of the time step?',
      options: ['picoseconds', 'femtoseconds', 'reduced τ', 'nanoseconds'],
      answer: 1,
      why: 'units real uses Å, fs, kcal/mol and atm. units metal uses ps and eV; units lj is reduced. GROMACS always uses ps.',
    },
    {
      q: 'Our NumPy code and LAMMPS give identical energies for 500 steps. What had to be the same?',
      options: ['Only the potential', 'Positions, velocities, potential (including the cutoff treatment), time step and integrator', 'Only the number of atoms', 'Nothing; MD is always reproducible across codes'],
      answer: 1,
      why: 'Step-by-step agreement needs the same starting state and exactly the same model and integrator. A shifted potential or different cutoff would change the energies.',
    },
    {
      q: 'What is the usual order of stages in a production MD study?',
      options: ['NPT, minimize, build, NVT', 'Build, minimize, NVT, NPT, production, analysis', 'Production, then minimize', 'Analysis first'],
      answer: 1,
      why: 'Build the system, remove overlaps by minimization, heat at constant volume, equilibrate the density at constant pressure, then run production and analyse.',
    },
  ],
};
