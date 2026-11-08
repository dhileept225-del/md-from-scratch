// Input-file builder: one Lennard-Jones fluid described in reduced units, written out as a LAMMPS input
// (units lj, real or metal) or as GROMACS topology + mdp, with every unit conversion shown.
const kB = 1.380649e-23, NA = 6.02214076e23;
const GASES = { Ne: { s: 2.75, e: 36.0, m: 20.18 }, Ar: { s: 3.405, e: 119.8, m: 39.948 }, Kr: { s: 3.60, e: 171.0, m: 83.80 }, Xe: { s: 4.10, e: 221.0, m: 131.29 } };   // sigma (A), epsilon/kB (K), mass (g/mol); same as Lesson 2.1
const SIZES = { 108: 3, 256: 4, 500: 5, 864: 6, 2048: 8, 4000: 10 };

export function mountRealcodes(root) {
  const st = { gas: 'Ar', code: 'lj', ens: 'nvt', T: 0.722, rho: 0.8442, P: 0.0, dt: 0.005, rc: 2.5, N: 500 };
  root.classList.add('sim');
  root.innerHTML = `
  <div class="sim-view" style="aspect-ratio:auto;min-height:560px;min-width:0;display:flex;flex-direction:column;gap:10px;padding:12px">
    <div class="code" style="max-width:none;flex:1;min-width:0"><button class="btn sm ghost copy" data-k="copy">Copy</button><pre data-o="out" style="margin:0;max-height:520px;overflow:auto;font-size:13px;padding-top:46px"></pre></div>
  </div>
  <div class="sim-ctrl">
    <div class="seg">${Object.keys(GASES).map(g => `<button data-g="${g}" class="${g === st.gas ? 'on' : ''}">${g}</button>`).join('')}</div>
    <div class="seg"><button data-c="lj" class="on">LMP lj</button><button data-c="real">LMP real</button><button data-c="metal">LMP metal</button><button data-c="gmx">GROMACS</button></div>
    <div class="seg"><button data-e="nve">NVE</button><button data-e="nvt" class="on">NVT</button><button data-e="npt">NPT</button></div>
    <label><span>Temperature <i class="nc">T</i>*</span><output data-o="T"></output><input type="range" data-k="T" min="0.5" max="3" step="0.01" value="0.722"></label>
    <label><span>Density <i class="nc">ρ</i>*</span><output data-o="rho"></output><input type="range" data-k="rho" min="0.05" max="1.1" step="0.01" value="0.8442"></label>
    <label data-np hidden><span>Pressure <i class="nc">P</i>*</span><output data-o="P"></output><input type="range" data-k="P" min="0" max="5" step="0.1" value="0"></label>
    <label><span>Time step <i class="nc">Δt</i>*</span><output data-o="dt"></output><input type="range" data-k="dt" min="0.001" max="0.01" step="0.001" value="0.005"></label>
    <div class="seg">${Object.keys(SIZES).map(n => `<button data-n="${n}" class="${+n === st.N ? 'on' : ''}">${n}</button>`).join('')}</div>
    <div class="readouts">
      <div><span>Temperature</span><b data-o="TK"></b></div>
      <div><span>Density</span><b data-o="rg"></b></div>
      <div><span>Time unit <i class="nc">τ</i> / time step</span><b data-o="tau"></b></div>
      <div><span>Box length / cutoff</span><b data-o="box"></b></div>
    </div>
  </div>`;
  const $ = s => root.querySelector(s), o = {}; root.querySelectorAll('[data-o]').forEach(e => o[e.dataset.o] = e);
  const seg = (attr, f) => root.querySelectorAll(`[data-${attr}]`).forEach(b => b.onclick = () => { root.querySelectorAll(`[data-${attr}]`).forEach(z => z.classList.toggle('on', z === b)); f(b.dataset[attr]); render(); });
  seg('g', v => st.gas = v); seg('c', v => st.code = v); seg('e', v => { st.ens = v; $('[data-np]').hidden = v !== 'npt'; }); seg('n', v => st.N = +v);
  ['T', 'rho', 'P', 'dt'].forEach(k => $(`[data-k=${k}]`).oninput = e => { st[k] = +e.target.value; render(); });
  $('[data-k=copy]').onclick = () => { navigator.clipboard?.writeText(o.out.textContent); $('[data-k=copy]').textContent = 'Copied'; setTimeout(() => $('[data-k=copy]').textContent = 'Copy', 1200); };
  const f = (x, d = 4) => (+x.toPrecision(6)).toFixed(d).replace(/\.?0+$/, '') || '0';

  function units() {
    const g = GASES[st.gas], epsJ = g.e * kB, m = g.m * 1e-3 / NA, sig = g.s * 1e-10;
    const tau = sig / Math.sqrt(epsJ / m);                               // s
    return { g, epsJ, tau, T: st.T * g.e, rho: st.rho * g.m / NA / (g.s * 1e-8) ** 3, P: st.P * epsJ / sig ** 3,   // K, g/cm3, Pa
      n: SIZES[st.N], a: Math.cbrt(4 / st.rho) };                         // cells per side, lattice constant (sigma)
  }
  function render() {
    const u = units(), g = u.g, n = u.n, Lr = n * u.a, dtfs = st.dt * u.tau * 1e15, ens = st.ens;
    let txt = '';
    if (st.code === 'lj') {
      txt = `# ${st.N} LJ atoms (${st.gas}) in reduced units: T* = ${st.T}, rho* = ${st.rho}
units           lj
atom_style      atomic
lattice         fcc ${st.rho}                 # in lj units the lattice argument is the reduced density
region          box block 0 ${n} 0 ${n} 0 ${n}
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create ${st.T} 12345 mom yes rot yes dist gaussian
pair_style      lj/cut ${st.rc}
pair_coeff      1 1 1.0 1.0
pair_modify     tail yes                       # long-range correction for energy and pressure
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
timestep        ${st.dt}
${ens === 'nve' ? 'fix             1 all nve' : ens === 'nvt' ? `fix             1 all nvt temp ${st.T} ${st.T} ${f(100 * st.dt)}       # damping time = 100 steps` : `fix             1 all npt temp ${st.T} ${st.T} ${f(100 * st.dt)} iso ${st.P} ${st.P} ${f(1000 * st.dt)}`}
thermo_style    custom step temp pe etotal press density
thermo          1000
run             100000`;
    } else if (st.code === 'real' || st.code === 'metal') {
      const real = st.code === 'real', eps = real ? g.e * kB * NA / 4184 : g.e * 8.617333e-5, dt = real ? dtfs : dtfs / 1000, P = real ? u.P / 101325 : u.P / 1e5;
      txt = `# ${st.N} ${st.gas} atoms in ${st.code} units (${real ? 'Angstrom, fs, kcal/mol, atm' : 'Angstrom, ps, eV, bar'}): T = ${f(u.T, 1)} K, rho = ${f(u.rho, 3)} g/cm^3
units           ${st.code}
atom_style      atomic
lattice         fcc ${f(u.a * g.s, 4)}               # lattice constant in Angstrom
region          box block 0 ${n} 0 ${n} 0 ${n}
create_box      1 box
create_atoms    1 box
mass            1 ${g.m}
velocity        all create ${f(u.T, 2)} 12345 mom yes rot yes dist gaussian
pair_style      lj/cut ${f(st.rc * g.s, 4)}
pair_coeff      1 1 ${f(eps, 6)} ${g.s}           # epsilon (${real ? 'kcal/mol' : 'eV'}), sigma (Angstrom)
pair_modify     tail yes
neighbor        ${f(0.3 * g.s, 3)} bin
neigh_modify    every 1 delay 0 check yes
timestep        ${f(dt, 5)}                    # ${real ? 'fs' : 'ps'}
${ens === 'nve' ? 'fix             1 all nve' : ens === 'nvt' ? `fix             1 all nvt temp ${f(u.T, 2)} ${f(u.T, 2)} ${f(100 * dt, 4)}` : `fix             1 all npt temp ${f(u.T, 2)} ${f(u.T, 2)} ${f(100 * dt, 4)} iso ${f(P, 2)} ${f(P, 2)} ${f(1000 * dt, 3)}   # pressure in ${real ? 'atm' : 'bar'}`}
thermo_style    custom step temp pe etotal press density
thermo          1000
run             100000`;
    } else {
      const nm = g.s / 10, kJ = g.e * kB * NA / 1000, dtps = dtfs / 1000, Lnm = Lr * nm, bar = u.P / 1e5;
      txt = `; ===== ${st.gas.toLowerCase()}.top =====
[ defaults ]
; nbfunc comb-rule gen-pairs
  1      2         no
[ atomtypes ]
; name  mass      charge ptype sigma(nm)  epsilon(kJ/mol)
  ${st.gas.toUpperCase().padEnd(4)}  ${String(g.m).padEnd(8)}  0.0    A     ${f(nm, 5)}   ${f(kJ, 5)}
[ moleculetype ]
  ${st.gas.toUpperCase()}  1
[ atoms ]
  1  ${st.gas.toUpperCase()}  1  ${st.gas.toUpperCase()}  ${st.gas.toUpperCase()}  1  0.0
[ system ]
${st.N} ${st.gas} atoms
[ molecules ]
${st.gas.toUpperCase()}  ${st.N}

; ===== md.mdp =====
integrator      = md
dt              = ${f(dtps, 5)}          ; ps
nsteps          = 100000
cutoff-scheme   = Verlet
rvdw            = ${f(st.rc * nm, 4)}           ; nm
vdw-modifier    = Potential-shift
DispCorr        = EnerPres           ; long-range correction
coulombtype     = Cut-off
rcoulomb        = ${f(st.rc * nm, 4)}
${ens === 'nve' ? 'tcoupl          = no\nverlet-buffer-tolerance = -1\nrlist           = ' + f(st.rc * nm + 0.1, 4) : `tcoupl          = v-rescale\ntc-grps         = System\ntau-t           = 0.1               ; ps\nref-t           = ${f(u.T, 1)}`}
${ens === 'npt' ? `pcoupl          = C-rescale\ntau-p           = 2.0               ; ps\nref-p           = ${f(bar, 1)}             ; bar\ncompressibility = 1e-4              ; 1/bar: set for your liquid` : 'pcoupl          = no'}
gen-vel         = yes
gen-temp        = ${f(u.T, 1)}
nstenergy       = 1000

; ===== starting box: ${st.N} atoms in a cube of ${f(Lnm, 3)} nm =====
; e.g. gmx insert-molecules -ci ${st.gas.toLowerCase()}.gro -nmol ${st.N} -box ${f(Lnm, 3)} -o start.gro
;      (then energy-minimize before MD)`;
    }
    o.out.textContent = txt;
    o.T.textContent = st.T.toFixed(3); o.rho.textContent = st.rho.toFixed(3); o.P.textContent = st.P.toFixed(1); o.dt.textContent = st.dt.toFixed(3);
    o.TK.textContent = `${f(u.T, 1)} K`; o.rg.textContent = `${f(u.rho, 3)} g/cm³`;
    o.tau.textContent = `${f(u.tau * 1e12, 3)} ps / ${f(dtfs, 2)} fs`;
    o.box.textContent = `${f(Lr * g.s / 10, 3)} nm / ${f(st.rc * g.s / 10, 3)} nm`;
  }
  render();
  return () => {};
}
