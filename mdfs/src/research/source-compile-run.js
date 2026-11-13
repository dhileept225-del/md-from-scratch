// Research · Part 0.1: source code, compiling, running.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/src/pair_lj_cut.cpp';

const TWO = r => `units           lj
atom_style      atomic
atom_modify     map yes
region          box block -5 5 -5 5 -5 5
create_box      1 box
create_atoms    1 single 0.0 0.0 0.0
create_atoms    1 single ${r} 0.0 0.0
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
variable        fx1 equal fx[1]
thermo_style    custom step pe v_fx1
thermo_modify   norm no format float %.8g
run             0
`;

function hand(r) {
  if (r >= 2.5) return { E: 0, F: 0, lines: `r = ${r} is at or beyond the cutoff 2.5:
the pair is not counted, so E = <b>0</b> and F = <b>0</b>` };
  const s6 = r ** -6, s12 = s6 * s6, E = 4 * (s12 - s6), Fr = 24 * (2 * s12 - s6) / r, F1 = -Fr;
  return { E, F: F1, lines: `r      = ${r}
r⁻⁶    = 1/${r}⁶            = ${s6.toFixed(8)}
r⁻¹²   = (r⁻⁶)²             = ${s12.toFixed(8)}
E      = 4(r⁻¹² − r⁻⁶)      = 4 × (${(s12 - s6).toFixed(8)}) = <b>${g8(E)}</b>
F(r)   = 24(2r⁻¹² − r⁻⁶)/r  = ${g8(Fr)}   (${Fr > 0 ? 'positive: the atoms push apart' : 'negative: the atoms pull together'})
atom 1 sits at x = 0, atom 2 at x = ${r}, so the force on atom 1 points
${Fr > 0 ? 'away from atom 2 (towards −x)' : 'towards atom 2 (towards +x)'}:  fx[1] = −F(r) = <b>${g8(F1)}</b>` };
}

const PY = `# Two Lennard-Jones atoms: energy and force by hand (epsilon = sigma = 1, cutoff 2.5)
def lj(r, rc=2.5):
    if r >= rc:
        return 0.0, 0.0
    s6 = r ** -6                         # r^-6
    E = 4 * (s6 * s6 - s6)               # 4(r^-12 - r^-6)
    F = 24 * (2 * s6 * s6 - s6) / r      # F(r) = -dE/dr, positive = repulsive
    return E, -F                         # force on atom 1 (at x = 0) along x

for r in [0.9, 1.0, 1.1, 2 ** (1 / 6), 1.5, 2.4, 2.6]:
    E, fx1 = lj(r)
    print(f"r = {r:.6f}   PotEng {E:.8g}   fx[1] {fx1:.8g}")
`;

export default {
  id: 'source-compile-run',
  title: 'Source code, compiling, running',
  minutes: 15,
  summary: 'LAMMPS is text written by people (the source code) that a compiler turns into a program your computer can run. You find the exact lines where LAMMPS computes the Lennard-Jones energy, see why changing the input needs no compiling but changing the formula does, and check a two-atom run by hand.',
  learn: [
    'What source code, a compiler and a program are, using LAMMPS as the example',
    'Where the Lennard-Jones formula lives in the LAMMPS source, line by line',
    'Which changes need compiling and which do not, checked with a live two-atom run',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: three things, one path</h2>
    <p>Every program you use, LAMMPS included, starts as <b>source code</b>: ordinary text files written by people. A computer cannot run that text directly. A second program, the <b>compiler</b>, translates it once into machine instructions and saves the result as the <b>program</b> (also called the executable). For LAMMPS the program is a single file named <code>lmp</code>.</p>
    <svg class="flow" viewBox="0 0 760 120" role="img" aria-label="source code, compiler, program, input script, log">
      <defs><marker id="fa" viewBox="0 0 10 10" refX="9" refY="5" markerWidth="6" markerHeight="6" orient="auto"><path d="M0 0L10 5L0 10z" fill="#8FA3C7"/></marker></defs>
      <g fill="none" stroke="#8FA3C7" stroke-width="1.4"><path d="M138 50H168" marker-end="url(#fa)"/><path d="M300 50H330" marker-end="url(#fa)"/><path d="M462 50H492" marker-end="url(#fa)"/><path d="M462 50" /><path d="M624 50H654" marker-end="url(#fa)"/></g>
      <g font-family="Space Grotesk, sans-serif" text-anchor="middle">
        <rect x="8" y="22" width="130" height="56" rx="10" fill="rgba(56,225,255,.08)" stroke="#38E1FF"/><text x="73" y="46" fill="#F2F5FA" font-size="14" font-weight="700">source code</text><text x="73" y="64" fill="#8FA3C7" font-size="11">pair_lj_cut.cpp, …</text>
        <rect x="170" y="22" width="130" height="56" rx="10" fill="none" stroke="#8FA3C7"/><text x="235" y="46" fill="#F2F5FA" font-size="14" font-weight="700">compiler</text><text x="235" y="64" fill="#8FA3C7" font-size="11">once, minutes</text>
        <rect x="332" y="22" width="130" height="56" rx="10" fill="rgba(255,194,75,.1)" stroke="#FFC24B"/><text x="397" y="46" fill="#F2F5FA" font-size="14" font-weight="700">program lmp</text><text x="397" y="64" fill="#8FA3C7" font-size="11">runs many times</text>
        <rect x="494" y="22" width="130" height="56" rx="10" fill="none" stroke="#8FA3C7"/><text x="559" y="46" fill="#F2F5FA" font-size="14" font-weight="700">input script</text><text x="559" y="64" fill="#8FA3C7" font-size="11">read while running</text>
        <rect x="656" y="22" width="96" height="56" rx="10" fill="none" stroke="#8FA3C7"/><text x="704" y="46" fill="#F2F5FA" font-size="14" font-weight="700">log</text><text x="704" y="64" fill="#8FA3C7" font-size="11">results</text>
        <text x="380" y="106" fill="#8FA3C7" font-size="12">change the source → compile again · change the input → just run again</text>
      </g>
    </svg>
    <p>A kitchen is a good picture. The source code is a recipe book; the compiler is the training that turns a new cook into someone who knows every recipe by heart; the program is that trained cook; the input script is today's order; the log is the meal. Changing today's order is instant. Changing a recipe means training the cook again.</p>
  </section>

  <section class="block prose">
    <h2>2 · LAMMPS in numbers</h2>
    <p>We downloaded LAMMPS (stable release 22 Jul 2025, update 6) and compiled it on a small 2-processor cloud computer. These are the real counts:</p>
    <div class="tbl wide"><table>${head(['What', 'Number'])}
${row(['source files in all of LAMMPS (.cpp)', '1,860'])}
${row(['lines of C++ in all of LAMMPS (.cpp and .h)', '1.19 million'])}
${row(['optional packages', '94'])}
${row(['core only (no packages): files and lines', '363 .cpp + 382 .h, 235,476 lines'])}
${row(['our build: core + MOLECULE, KSPACE, RIGID', '451 files compiled'])}
${row(['size of the program lmp', '11.6 MB'])}
${row(['time to compile everything', '266 s (about 4½ minutes) on 2 processors'])}
${row(['time to recompile after changing one file', '3.8 s: one file compiled, then the pieces joined again (13–16 s in a git download, where a second file that records the git state is rebuilt too: page 0.6)'])}</table></div>
    <p>The last two rows matter most. Compiling everything takes minutes, but after changing one file the compiler only redoes that file and joins the pieces again, which takes seconds. That is why changing LAMMPS is practical.</p>
  </section>

  <section class="block prose">
    <h2>3 · Where the Lennard-Jones formula lives</h2>
    <p>The force between Lennard-Jones atoms is computed in one file, <code>src/pair_lj_cut.cpp</code>. You do not need to know C++ to read the lines that matter: they are the formula of Lesson 1.2 written with a few abbreviations.</p>
    <p><b>First, four numbers are prepared once</b> (<a href="${GH}#L509-L512" target="_blank" rel="noopener">lines 509–512</a>):</p>
    <div class="calc">lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);</div>
    <div class="tbl wide"><table>${head(['Code', 'Plain words', 'Maths'])}
${row(['<code>*</code>', 'multiply', '×'])}
${row(['<code>pow(sigma, 12.0)</code>', 'σ to the power 12', 'σ¹²'])}
${row(['<code>[i][j]</code>', 'for atom types i and j (here only 1 and 1)', ''])}
${row(['<code>lj1 … lj4</code>', 'four numbers', '48εσ¹², 24εσ⁶, 4εσ¹², 4εσ⁶'])}</table></div>
    <p><b>Then, for every pair of atoms closer than the cutoff, every step</b> (<a href="${GH}#L114-L130" target="_blank" rel="noopener">lines 114–130</a>):</p>
    <div class="calc">r2inv = 1.0 / rsq;
r6inv = r2inv * r2inv * r2inv;
forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
fpair = factor_lj * forcelj * r2inv;
…
evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];</div>
    <div class="tbl wide"><table>${head(['Line', 'Plain words', 'Maths'])}
${row(['<code>rsq</code>', 'the squared distance, computed just before: dx² + dy² + dz²', 'r²'])}
${row(['<code>r2inv = 1.0 / rsq</code>', 'one divided by r²', '1/r²'])}
${row(['<code>r6inv = r2inv * r2inv * r2inv</code>', '(1/r²) three times', '1/r⁶'])}
${row(['<code>forcelj = …</code>', 'the force times the distance', 'r·F(r)'])}
${row(['<code>fpair = … * r2inv</code>', 'divide by r² once more (factor_lj is 1 here)', 'F(r)/r'])}
${row(['<code>evdwl = …</code>', 'the pair energy (offset is 0 unless shifted)', 'E(r)'])}</table></div>
    <p>Put the four numbers in, with ε = σ = 1:</p>
    <div class="calc">evdwl   = (1/r⁶) × (4 × (1/r⁶) − 4)    = 4(r⁻¹² − r⁻⁶)     the LJ energy
forcelj = (1/r⁶) × (48 × (1/r⁶) − 24)  = 48r⁻¹² − 24r⁻⁶    = r × F(r)
fpair   = forcelj / r²                  = F(r)/r
x-force on atom i = delx × fpair = F(r) × (delx/r)       (delx = x_i − x_j)</div>
    <p>LAMMPS never takes a square root here: it works with r² and 1/r² because that is faster. <code>offset</code> is zero unless the energy is shifted (<code>pair_modify shift yes</code>), and <code>factor_lj</code> is 1 except for bonded neighbours in molecules.</p>
  </section>

  <section class="block prose">
    <h2>4 · Two kinds of change</h2>
    <div class="tbl wide"><table>${head(['You want to change…', 'Where', 'Compile again?', 'Time'])}
${row(['the distance between atoms, ε, σ, the cutoff, temperature, number of steps', 'the input script', 'no', 'run again, instant'])}
${row(['what is printed, how often, which thermostat', 'the input script', 'no', 'run again, instant'])}
${row(['the formula itself, for example r⁻¹² → r⁻⁹', 'pair_lj_cut.cpp', 'yes', 'about 4 s for one file'])}
${row(['a new kind of force or thermostat', 'a new .cpp file', 'yes', 'seconds to minutes'])}</table></div>
    <p>The program reads the input script while it runs, so input changes need nothing else. The formula, however, is fixed inside the program by the compiler: to change it you change the source and compile again. Part E of this section does exactly that.</p>
  </section>

  <section class="block">
    <h2>5 · Try it: two atoms, live</h2>
    <p class="muted">LAMMPS runs in your browser. Change the distance, predict by hand, run, compare. This is an input change, so no compiling is needed.</p>
    <div class="lmp" id="two">
      <div class="row"><label class="rlab">Distance r <input type="number" class="r-in" value="1.1" min="0.8" max="3" step="0.01"></label>
        <button class="btn sm primary" data-k="run">▶ Run LAMMPS</button><span class="lmp-status muted"></span></div>
      <div class="row quick"><span class="muted small">Try:</span>${['0.9', '1.0', '1.1', '1.122462', '1.5', '2.4', '2.6'].map(v => `<button class="btn sm ghost" data-r="${v}">${v}</button>`).join('')}</div>
      <p><b>Predict by hand</b></p>
      <div class="calc hand"></div>
      <div class="lmp-compare"></div>
      <details class="lmp-log"><summary>Input script and full log</summary><pre class="inp"></pre><pre class="out"></pre></details>
    </div>
  </section>

  <section class="block prose">
    <h2>What you just saw</h2>
    <ul>
      <li>Every distance you tried gave the same numbers by hand and from LAMMPS: the four lines of section 3 are the whole calculation.</li>
      <li>At r = 1.0 the energy is exactly zero and the push is 24; at r = 2<sup>1/6</sup> = 1.122462 the force is almost zero (the bottom of the well); beyond 2.5 LAMMPS prints zero because of the cutoff.</li>
      <li>None of this needed compiling. Changing the 12 in r⁻¹² would.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Check every number in Python</h2>
    <div class="code"><button class="btn sm ghost copy">Copy</button><pre class="py"></pre></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Which of these need LAMMPS to be compiled again: (a) a smaller time step, (b) printing the pressure every 10 steps, (c) a new way of computing temperature, (d) adding the MANYBODY package?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>(a) and (b) are input changes: no. (c) changes source code: yes. (d) adds source files to the program: yes, LAMMPS must be compiled with the package switched on.</div></div>
    <div class="ex"><p><b>2.</b> With ε = 2 and σ = 1, what are lj3 and lj4, and what energy does LAMMPS print at r = 1.1?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>lj3 = 4 × 2 × 1 = 8 and lj4 = 8. E = 8(r⁻¹² − r⁻⁶) = 2 × (−0.98337245) = <b>−1.9667449</b>. Change <code>pair_coeff 1 1 1.0 1.0 2.5</code> to <code>pair_coeff 1 1 2.0 1.0 2.5</code> in the log's input to check.</div></div>
    <div class="ex"><p><b>3.</b> Why does LAMMPS compute <code>r2inv</code> and never r itself?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The LJ formula only needs even powers of r (r⁻⁶, r⁻¹², and F/r = forcelj × r⁻²). These come from r² = dx² + dy² + dz² by multiplying, with no square root. A square root is several times slower than a multiplication, and this line runs billions of times in a long simulation.</div></div>
  </section>`,

  mount(root) {
    const box = root.querySelector('#two'), rin = box.querySelector('.r-in'), status = box.querySelector('.lmp-status');
    const handBox = box.querySelector('.hand'), cmp = box.querySelector('.lmp-compare');
    root.querySelector('pre.py').textContent = PY;
    let busy = false;
    const getR = () => { const v = parseFloat(rin.value); return isFinite(v) && v > 0.5 && v < 5 ? v : 1.1; };
    const showHand = () => { const h = hand(getR()); handBox.innerHTML = h.lines; cmp.innerHTML = ''; return h; };
    showHand();
    rin.oninput = showHand;
    async function run() {
      if (busy) return; busy = true;
      const r = getR(), h = showHand(), script = TWO(r);
      box.querySelectorAll('button').forEach(b => { b.disabled = true; });
      status.textContent = 'Loading LAMMPS and running… (the first run downloads about 3 MB)';
      const res = await runLammps(script);
      const ver = (res.version.match(/\(([^)]+)\)/) || [])[1];
      box.querySelector('pre.inp').textContent = script; box.querySelector('pre.out').textContent = res.lines.join('\n');
      if (res.error) status.innerHTML = `<span class="lmp-err">${res.error}</span>`;
      else {
        status.textContent = `Done in ${res.seconds.toFixed(2)} s${ver ? ` · LAMMPS ${ver}` : ''}`;
        const r0 = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols;
        if (r0) {
          const L = { PotEng: r0[c.indexOf('PotEng')], fx: r0[c.indexOf('v_fx1')] };
          const line = (k, hv, lv) => { const m = matchMarkup(g8(hv), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; };
          cmp.innerHTML = `<div class="tbl"><table><tr><th>r = ${r}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${line('energy E', h.E, L.PotEng)}${line('fx on atom 1', h.F, L.fx)}</table></div>`;
        }
      }
      box.querySelectorAll('button').forEach(b => { b.disabled = false; });
      busy = false;
    }
    box.querySelector('[data-k=run]').onclick = run;
    box.querySelectorAll('[data-r]').forEach(b => b.onclick = () => { rin.value = b.dataset.r; run(); });
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    const copy = root.querySelector('.copy');
    copy.onclick = async () => { try { await navigator.clipboard.writeText(PY); copy.textContent = 'Copied ✓'; } catch { copy.textContent = 'Select the code to copy'; } setTimeout(() => copy.textContent = 'Copy', 2000); };
    return () => {};
  },
};
