// Research · Start here: run real LAMMPS in the browser, change one line, check the result by hand.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

const MELT = `units           lj
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0

velocity        all create 3.0 87287 loop geom

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5

neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no

fix             1 all nve

thermo          50
run             250
`;

// ---------- the hand calculation, for any FCC LJ script of this family ----------
function readScript(s) {
  const num = re => { const m = s.match(re); return m ? parseFloat(m[1]) : NaN; };
  const block = s.match(/region\s+\S+\s+block\s+0\s+(\d+)\s+0\s+(\d+)\s+0\s+(\d+)/);
  const pc = s.match(/pair_coeff\s+1\s+1\s+([\d.eE+-]+)\s+([\d.eE+-]+)(?:\s+([\d.eE+-]+))?/);
  const p = {
    rho: num(/lattice\s+fcc\s+([\d.eE+-]+)/), T: num(/velocity\s+all\s+create\s+([\d.eE+-]+)/),
    rcGlobal: num(/pair_style\s+lj\/cut\s+([\d.eE+-]+)/),
    eps: pc ? parseFloat(pc[1]) : NaN, sig: pc ? parseFloat(pc[2]) : NaN, rc: pc && pc[3] ? parseFloat(pc[3]) : NaN,
    cells: block ? [+block[1], +block[2], +block[3]] : null,
    shifted: /pair_modify[^\n]*shift\s+yes/.test(s), tail: /pair_modify[^\n]*tail\s+yes/.test(s),
    units: (s.match(/^\s*units\s+(\S+)/m) || [])[1],
  };
  if (!isFinite(p.rc)) p.rc = p.rcGlobal;
  return p;
}
export function handStep0(s) {
  const p = readScript(s);
  if (p.units !== 'lj' || ![p.rho, p.T, p.eps, p.sig, p.rc].every(isFinite) || !p.cells) return { ok: false, why: 'The hand calculator only follows scripts like the one above: units lj, an fcc lattice, a block region starting at 0, one lj/cut pair_coeff line and velocity all create.' };
  if (p.shifted || p.tail) return { ok: false, why: 'pair_modify shift or tail changes the energy formula; the hand calculator does not include it yet.' };
  const N = 4 * p.cells[0] * p.cells[1] * p.cells[2], a = Math.cbrt(4 / p.rho);
  if (p.rc + 0.3 > 0.5 * a * Math.min(...p.cells)) return { ok: false, why: 'The box is too small for this cutoff.' };
  const base = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]], shells = new Map();
  let E = 0, W = 0;
  for (let i = -5; i <= 5; i++) for (let j = -5; j <= 5; j++) for (let k = -5; k <= 5; k++) for (const b of base) {
    const r = a * Math.hypot(i + b[0], j + b[1], k + b[2]);
    if (r > 1e-9 && r < p.rc) {
      const key = r.toFixed(6); shells.set(key, (shells.get(key) || 0) + 1);
      const s6 = (p.sig / r) ** 6; E += 0.5 * 4 * p.eps * (s6 * s6 - s6); W += 0.5 * 24 * p.eps * (2 * s6 * s6 - s6);
    }
  }
  const f = (N - 1) / N;
  return { ok: true, N, a, shells: [...shells].map(([k, n]) => [+k, n]).sort((x, y) => x[0] - y[0]),
    T: p.T, E, TotEng: E + 1.5 * p.T * f, Press: p.rho * p.T * f + p.rho * W / 3 };
}

// ---------- the four planned changes ----------
const EXP = [
  { key: 'cut', title: 'Cut the force shorter: 2.5 → 2.2', edit: s => s.replace(/lj\/cut 2\.5/, 'lj/cut 2.2').replace(/1 1 1\.0 1\.0 2\.5/, '1 1 1.0 1.0 2.2'),
    lines: 'pair_style      lj/cut 2.2\npair_coeff      1 1 1.0 1.0 2.2',
    calc: `the 4th shell (r = 2.3753) is now beyond 2.2: drop it
S3 = 12 + 6/8 + 24/27              = 13.6388889
S6 = 12 + 6/64 + 24/729            = 12.1266718
E_pair = 2(0.12697593 × 12.1266718 − 0.35633682 × 13.6388889)
       = 2(1.53979542 − 4.86003829)  = <b>−6.6404857</b>
TotEng = −6.64048574 + 4.49887500   = <b>−2.1416107</b>
Press  = 2.53196685 + 3.3768 × (3.07959084 − 4.86003829)
       = 2.53196685 − 6.01221495    = <b>−3.4802481</b>` },
  { key: 'eps', title: 'Weaker attraction: ε = 1.0 → 0.5', edit: s => s.replace(/1 1 1\.0 1\.0/, '1 1 0.5 1.0'),
    lines: 'pair_coeff      1 1 0.5 1.0 2.5',
    calc: `every energy and every force is multiplied by ε
E_pair = 0.5 × (−6.77336805)        = <b>−3.386684</b>
TotEng = −3.38668403 + 4.49887500   = <b>1.112191</b>
Press  = 2.53196685 + 0.5 × (−6.23531727)
       = 2.53196685 − 3.11765864    = <b>−0.58569179</b>
(the kinetic parts do not depend on ε)` },
  { key: 'rho', title: 'Squeeze the crystal: ρ = 0.8442 → 0.9', edit: s => s.replace(/fcc 0\.8442/, 'fcc 0.9'),
    lines: 'lattice         fcc 0.9',
    calc: `x = ρ²/2 = 0.81/2 = 0.405     x² = 0.164025
(the 5th shell, 1.6441 × √(5/2) = 2.5996, is still outside 2.5)
E_pair = 2(0.164025 × 12.1296015 − 0.405 × 13.8263889)
       = 2(1.98955789 − 5.59968750)  = <b>−7.2202592</b>
TotEng = −7.22025922 + 4.49887500   = <b>−2.7213842</b>
Press  = 0.9 × 3 × 0.99975 + 3.6 × (3.97911577 − 5.59968750)
       = 2.69932500 − 5.83405822    = <b>−3.1347332</b>` },
  { key: 'T', title: 'Start colder: T = 3.0 → 1.0', edit: s => s.replace(/create 3\.0/, 'create 1.0'),
    lines: 'velocity        all create 1.0 87287 loop geom',
    calc: `step 0: atoms still on the lattice, only velocities changed
E_pair = <b>−6.7733681</b>   (does not depend on velocities)
KE/N   = 1.5 × 1.0 × 0.99975        = 1.49962500
TotEng = −6.77336805 + 1.49962500   = <b>−5.2737431</b>
Press  = 0.8442 × 1.0 × 0.99975 + (−6.23531727)
       = 0.84398895 − 6.23531727    = <b>−5.3913283</b>` },
];

const PY = `# Step 0 of the LAMMPS melt example, by hand (no LAMMPS needed)
from fractions import Fraction

def step0(rho=0.8442, T=3.0, eps=1.0, rc=2.5, N=4000):
    a = (4 / rho) ** (1 / 3)                     # fcc cube edge
    shells = {1: 12, 2: 6, 3: 24, 4: 12, 5: 24}  # k: atoms at r = a*sqrt(k/2)
    x = rho**2 / 2                               # r_k^-6 = x / k^3
    S3 = sum(n / k**3 for k, n in shells.items() if a * (k / 2) ** 0.5 < rc)
    S6 = sum(n / k**6 for k, n in shells.items() if a * (k / 2) ** 0.5 < rc)
    E = 2 * eps * (x * x * S6 - x * S3)          # energy per atom
    f = (N - 1) / N                              # 3 of 3N degrees of freedom removed
    TotEng = E + 1.5 * T * f
    Press = rho * T * f + 4 * rho * eps * (2 * x * x * S6 - x * S3)
    return E, TotEng, Press

for change, kw in [("original", {}), ("cutoff 2.2", {"rc": 2.2}), ("epsilon 0.5", {"eps": 0.5}),
                   ("density 0.9", {"rho": 0.9}), ("T 1.0", {"T": 1.0})]:
    E, TE, P = step0(**kw)
    print(f"{change:12s} E_pair {E:.8g}  TotEng {TE:.8g}  Press {P:.8g}")
`;

export default {
  id: 'lammps-start',
  title: 'Run LAMMPS, change it, check it by hand',
  minutes: 25,
  summary: 'LAMMPS is the program most researchers use for molecular dynamics. Here it runs inside your browser. You run its classic first example, change one line at a time, and get every number LAMMPS prints at step 0 with a calculator, to the last digit.',
  learn: [
    'What LAMMPS is, and how an input script tells it what to do',
    'Run the official melt example live and read its output',
    'Predict the effect of four changes by hand, then check them against LAMMPS',
  ],

  html: () => `
  <section class="block prose">
    <h2>What LAMMPS is</h2>
    <p>LAMMPS (Large-scale Atomic/Molecular Massively Parallel Simulator) is a free molecular dynamics program from Sandia National Laboratories. Inside, it does what you built in Chapters 1 to 3: it computes forces, moves atoms with velocity Verlet and controls temperature. It does this for anything from 4 atoms to billions, on one laptop or thousands of processors.</p>
    <p>You do not write a program to use LAMMPS. You write an <b>input script</b>: a plain text list of commands. LAMMPS reads it from top to bottom, carries out each command, and prints a <b>log</b> with the results.</p>
    <p>The LAMMPS on this page is the real program, compiled to WebAssembly so that it runs inside your browser. Nothing is sent to a server. Its version is printed in the log when it runs.</p>
  </section>

  <section class="block">
    <h2>How every page in this section works</h2>
    <div class="how cycle">
      <div class="card"><div class="ico">1</div><h3>Understand</h3><p>Read the input or the code, one line at a time.</p></div>
      <div class="card"><div class="ico">2</div><h3>Modify</h3><p>Change one thing, and predict the new result by hand first.</p></div>
      <div class="card"><div class="ico">3</div><h3>Observe</h3><p>Run LAMMPS and see what changed.</p></div>
      <div class="card"><div class="ico">4</div><h3>Match</h3><p>Compare LAMMPS with your hand calculation, digit by digit.</p></div>
    </div>
  </section>

  <section class="block prose">
    <h2>1 · Understand: the melt example, line by line</h2>
    <p>This script ships with LAMMPS (<code>examples/melt/in.melt</code>). It builds a Lennard-Jones crystal and gives it so much kinetic energy that it melts.</p>
    <div class="tbl wide"><table>${head(['Line', 'What it does', 'Number you can check'])}
${row(['<code>units lj</code>', 'Lennard-Jones units: σ = ε = m = 1 (Lesson 2.1)', 'energies in ε, lengths in σ'])}
${row(['<code>atom_style atomic</code>', 'atoms are points: position, velocity, type', ''])}
${row(['<code>lattice fcc 0.8442</code>', 'an FCC crystal with 0.8442 atoms per σ³', 'cube edge a = (4/0.8442)<sup>1/3</sup> = 1.6796'])}
${row(['<code>region box block 0 10 0 10 0 10</code>', 'a box 10 cells long in each direction', 'L = 10a = 16.796'])}
${row(['<code>create_box 1 box</code>', 'one atom type, in that box', ''])}
${row(['<code>create_atoms 1 box</code>', 'fill the box: 4 atoms per FCC cell', 'N = 4 × 10³ = 4000'])}
${row(['<code>mass 1 1.0</code>', 'mass of type 1', 'm = 1'])}
${row(['<code>velocity all create 3.0 87287 loop geom</code>', 'random velocities, scaled to exactly T = 3.0; 87287 is the random seed', 'T = 3 at step 0'])}
${row(['<code>pair_style lj/cut 2.5</code>', 'Lennard-Jones force, cut at 2.5σ (Lesson 2.5)', ''])}
${row(['<code>pair_coeff 1 1 1.0 1.0 2.5</code>', 'ε = 1, σ = 1, cutoff 2.5 for type 1 with type 1', ''])}
${row(['<code>neighbor 0.3 bin</code>', 'neighbour lists reach 2.5 + 0.3 = 2.8 (Lesson 2.6)', ''])}
${row(['<code>neigh_modify every 20 delay 0 check no</code>', 'rebuild the lists every 20 steps', '250/20 → 12 rebuilds'])}
${row(['<code>fix 1 all nve</code>', 'move all atoms with velocity Verlet (Lesson 1.5)', ''])}
${row(['<code>thermo 50</code>', 'print a line of results every 50 steps', 'steps 0, 50, …, 250'])}
${row(['<code>run 250</code>', 'do 250 steps', ''])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Run it in your browser</h2>
    <p class="muted">Press Run. The first time, your browser downloads LAMMPS (about 3 MB). You can edit the script before running.</p>
    <div class="lmp" id="runner">
      <textarea class="lmp-in" spellcheck="false" wrap="off" aria-label="LAMMPS input script"></textarea>
      <div class="row"><button class="btn sm primary" data-k="run">▶ Run LAMMPS</button><button class="btn sm ghost" data-k="reset">Reset script</button><span class="lmp-status muted"></span></div>
      <div class="lmp-thermo tbl"></div>
      <details class="lmp-log"><summary>Full log</summary><pre></pre></details>
    </div>
  </section>

  <section class="block prose">
    <h2>3 · Match: step 0 by hand</h2>
    <p>At step 0 every atom still sits on the perfect crystal, so every atom has exactly the same neighbours. Within the cutoff 2.5 they form four shells, at distances r = a√(k/2) for k = 1, 2, 3, 4.</p>
    <div class="tbl"><table>${head(['shell k', 'atoms n', 'distance r', 'r⁻⁶', 'energy per atom ½·n·4(r⁻¹² − r⁻⁶)'])}
${row(['1', '12', '1.18765', '0.356337', '−5.5046614'])}
${row(['2', '6', '1.67960', '0.044542', '−0.5106972'])}
${row(['3', '24', '2.05708', '0.013198', '−0.6251271'])}
${row(['4', '12', '2.37531', '0.005568', '−0.1328823'])}
${row(['5', '24', '2.65569', '', 'beyond 2.5: not counted'])}
${row(['<b>sum</b>', '54', '', '', '<b>−6.7733681</b>'])}</table></div>
    <p>A shortcut makes this a two-line calculation. Because a³ = 4/ρ, every r⁻⁶ is a simple fraction of one number, x = ρ²/2:</p>
    <div class="calc">r_k⁶ = (a² k/2)³ = 2k³/ρ²    so   r_k⁻⁶ = x/k³,  r_k⁻¹² = x²/k⁶
x  = 0.8442²/2 = <b>0.35633682</b>      x² = 0.12697593
S3 = 12/1 + 6/8 + 24/27 + 12/64       = 13.8263889
S6 = 12/1 + 6/64 + 24/729 + 12/4096   = 12.1296015

E_pair = ½ Σ n·4(r⁻¹² − r⁻⁶) = 2(x² S6 − x S3)
       = 2(1.54016742 − 4.92685145)   = <b>−6.7733681</b>

KE per atom = (3/2) T (N−1)/N           (3 degrees of freedom of
            = 1.5 × 3 × 3999/4000       the whole crystal's motion
            = 4.4988750                  are removed)
TotEng = −6.77336805 + 4.49887500     = <b>−2.2744931</b>

Press = ρT(N−1)/N + (ρ/3)·½ Σ n·r·F(r),   r·F = 24(2r⁻¹² − r⁻⁶)
      = 0.8442 × 3 × 0.99975 + 4ρ(2x² S6 − x S3)
      = 2.53196685 + 3.3768 × (3.08033484 − 4.92685145)
      = 2.53196685 − 6.23531727      = <b>−3.7033504</b></div>
    <div class="lmp-compare" data-cmp="base"><p class="muted">Run LAMMPS above to fill in its numbers.</p></div>
  </section>

  <section class="block prose">
    <h2>4 · Modify, predict, observe, match</h2>
    <p>Each change below edits one or two lines. Do the prediction by hand first (it is worked out for you), then press the button: it changes the script, runs LAMMPS for step 0 only, and compares.</p>
    ${EXP.map((e, i) => `
    <div class="exp" data-exp="${i}">
      <h3>${String.fromCharCode(65 + i)}. ${e.title}</h3>
      <div class="calc">${e.lines}</div>
      <p><b>Predict by hand:</b></p>
      <div class="calc">${e.calc}</div>
      <button class="btn sm primary" data-run-exp="${i}">Run this change</button>
      <div class="lmp-compare" data-cmp="${e.key}"></div>
    </div>`).join('')}
    <h3>Your own change</h3>
    <p>Edit the script in step 2 yourself (for example <code>lattice fcc 0.8</code> or a cutoff of 1.5) and press Run: the page works out the hand prediction for your script and compares it with LAMMPS.</p>
  </section>

  <section class="block prose">
    <h2>What the four changes teach</h2>
    <div class="tbl wide"><table>${head(['Change', 'What moved', 'Why'])}
${row(['A. cutoff 2.2', 'E_pair rises by exactly 0.1328823', 'the 12 atoms of the 4th shell no longer count: the reported energy depends on the cutoff, so papers must state it'])}
${row(['B. ε = 0.5', 'E_pair halves, pressure changes, kinetic parts do not', 'ε multiplies every pair energy and force, never the velocities'])}
${row(['C. ρ = 0.9', 'every shell moves closer; only x changes', 'all distances scale with a = (4/ρ)<sup>1/3</sup>'])}
${row(['D. T = 1.0', 'TotEng and Press move; E_pair does not', 'at step 0 the positions are the same; only velocities changed'])}</table></div>
  </section>

  <section class="block prose">
    <h2>Check every number in Python</h2>
    <div class="code"><button class="btn sm ghost copy">Copy</button><pre class="py"></pre></div>
  </section>

  <section class="block prose">
    <h2>Coming next: changing the C++ code itself</h2>
    <p>Here you changed the <b>input</b>. Later pages change the <b>source code</b> of LAMMPS with the same four steps: read the C++ lines that compute the force, change them (for example the 12 in r⁻¹² to a 9), predict the new energy by hand, and compare.</p>
    <p>A browser cannot compile LAMMPS itself (that takes about 10 minutes on a desktop computer). So each code-change page will come with the modified LAMMPS already compiled: you will see the exact lines changed, run the original and the modified program side by side here, and match the difference by hand. You can also rebuild it on your own computer with the steps in Part 0.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Change the density to 0.8. Is the 5th shell inside the cutoff now? Predict E_pair.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>a = (4/0.8)<sup>1/3</sup> = 1.70998. The 4th shell is at 1.70998 × √2 = 2.4183 (inside); the 5th at 1.70998 × √2.5 = 2.7037 (outside). x = 0.32, x² = 0.1024: E_pair = 2(0.1024 × 12.1296015 − 0.32 × 13.8263889) = <b>−6.3647465</b>. LAMMPS prints −6.3647465.</div></div>
    <div class="ex"><p><b>2.</b> Set the cutoff to 1.5 in both lines. Which shells remain, and what is E_pair?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only the first shell (1.18765 &lt; 1.5 &lt; 1.67960). E_pair = 2 × 12 × (x² − x) = 24 × (0.12697593 − 0.35633682) = <b>−5.5046614</b>, exactly the first row of the shell table.</div></div>
    <div class="ex"><p><b>3.</b> In change B the total energy became positive. Does that mean the crystal flies apart?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Not by itself. The kinetic energy (4.499 per atom) is now larger than the binding energy (3.387 per atom), so the system is far above its melting point and the pressure is high. But the atoms are in a periodic box: they cannot fly away, they form a hot, dense fluid. A positive total energy only means the system would not stay bound if the box walls were removed.</div></div>
  </section>`,

  mount(root) {
    const R = root.querySelector('#runner'), ta = R.querySelector('.lmp-in'), status = R.querySelector('.lmp-status');
    const thermoBox = R.querySelector('.lmp-thermo'), logBox = R.querySelector('.lmp-log pre');
    ta.value = MELT;
    root.querySelector('pre.py').textContent = PY;
    let busy = false;

    const cmpTable = (hand, line, title) => {
      if (!line) return `<p class="lmp-err">LAMMPS printed no step-0 line.</p>`;
      const L = { E: line.E_pair, TotEng: line.TotEng, Press: line.Press };
      const rows = [['E_pair', hand.E, L.E], ['TotEng', hand.TotEng, L.TotEng], ['Press', hand.Press, L.Press]].map(([k, h, l]) => {
        const m = matchMarkup(g8(h), l);
        return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${l}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`;
      }).join('');
      return `<div class="tbl"><table><tr><th>${title}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${rows}</table></div>`;
    };
    const step0 = res => {
      const { cols, rows } = res.thermo; const r = rows.find(w => w[0] === '0');
      if (!r) return null; const o = {}; cols.forEach((c, i) => { o[c] = r[i]; }); return o;
    };
    const showThermo = res => {
      const { cols, rows } = res.thermo;
      thermoBox.innerHTML = rows.length ? `<table><tr>${cols.map(c => `<th>${c}</th>`).join('')}</tr>${rows.map(w => `<tr>${w.map(x => `<td>${x}</td>`).join('')}</tr>`).join('')}</table>` : '';
      logBox.textContent = res.lines.join('\n');
    };

    async function run(script, target) {
      if (busy) return; busy = true;
      root.querySelectorAll('[data-k=run],[data-run-exp]').forEach(b => { b.disabled = true; });
      try {
        status.textContent = 'Loading LAMMPS and running… (the first run downloads about 3 MB)';
        const res = await runLammps(script);
        showThermo(res);
        const ver = (res.version.match(/\(([^)]+)\)/) || [])[1];
        if (res.error) { status.innerHTML = `<span class="lmp-err">${res.error}</span>`; }
        else status.textContent = `Done in ${res.seconds.toFixed(1)} s${ver ? ` · LAMMPS ${ver}` : ''}`;
        const s0 = step0(res), hand = handStep0(script);
        const box = root.querySelector(`[data-cmp="${target}"]`);
        if (box) box.innerHTML = !res.error && hand.ok ? cmpTable(hand, s0, 'step 0') : res.error ? '' : `<p class="muted">${hand.why}</p>`;
      } catch (e) {
        status.innerHTML = `<span class="lmp-err">Could not start LAMMPS in this browser: ${String(e.message || e)}</span>`;
      }
      busy = false;
      root.querySelectorAll('[data-k=run],[data-run-exp]').forEach(b => { b.disabled = false; });
    }

    R.querySelector('[data-k=run]').onclick = () => {
      const s = ta.value, isBase = s.trim() === MELT.trim();
      run(s, isBase ? 'base' : 'own');
    };
    R.querySelector('[data-k=reset]').onclick = () => { ta.value = MELT; };
    // a comparison box for free edits, placed under the runner
    const own = document.createElement('div'); own.className = 'lmp-compare'; own.dataset.cmp = 'own'; R.appendChild(own);
    root.querySelectorAll('[data-run-exp]').forEach(b => b.onclick = () => {
      const e = EXP[+b.dataset.runExp];
      run(e.edit(MELT).replace(/run\s+250/, 'run             0'), e.key);
    });
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    const copy = root.querySelector('.copy');
    copy.onclick = async () => { try { await navigator.clipboard.writeText(PY); copy.textContent = 'Copied ✓'; } catch { copy.textContent = 'Select the code to copy'; } setTimeout(() => copy.textContent = 'Copy', 2000); };
    return () => {};
  },
};
