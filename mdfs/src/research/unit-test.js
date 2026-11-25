// Research · E.11: test it like the developers. LAMMPS's own unit-test machinery (googletest + YAML reference files)
// applied to pair lj/mine from page E.3: a reference generated with LAMMPS's mie/cut, a pass at 5e-14, and two
// deliberate bugs, each caught by a different test; the failure number reproduced by hand.
import { NATIVE } from './unittest-native.js';

const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const SETUP = `cd lammps/build
cmake -D ENABLE_TESTING=on .     # googletest is downloaded (here: from a local copy, the sandbox blocks the URL)
make -j2 test_pair_style          # 15 s: the library is already built

# make a reference: run an implementation you trust on the test system, write all numbers to YAML
./test_pair_style mol-pair-mie_ref.yaml -g mol-pair-lj_mine.yaml
# (then edit prerequisites, pair_style and pair_coeff for lj/mine; the numbers stay)

# test
./test_pair_style mol-pair-lj_mine.yaml -s      # -s prints error statistics`;

export default {
  id: 'unit-test',
  title: 'Test it like the developers',
  minutes: 18,
  summary: 'LAMMPS checks every force style against a YAML file of reference energies, forces and stresses for a 29-atom molecular system, and repeats the comparison after newton off, a restart, a data file and through single(). You build that test program, write a reference for pair lj/mine (page E.3) with a different implementation, LAMMPS\'s own mie/cut, see lj/mine pass at a tolerance of 5 × 10⁻¹⁴, then plant two bugs: forgetting the exponents in the restart file fails 186 comparisons in the restart tests, a typo in single() fails only the single test. The failure number for the first, 3.0167, you reproduce by hand from two YAML files.',
  learn: [
    'The force-style test: in.fourmol, init and run values, newton on/off, restart, data, single, extract',
    'Generating a reference (-g) from an independent implementation',
    'Relative error with the smaller value as divisor, and what 5e-14 allows',
    'Each copy of the formula needs its own test',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: what the developers' test checks</h2>
    <div class="prose"><p>${gh('unittest/force-styles/test_pair_style.cpp')} reads a YAML file that names a pair style, its coefficients and an input (in.fourmol: 29 atoms in four small molecules, 5 atom types, bonds and special-bond exclusions). It builds the system and compares with the numbers stored in the file:</p>
<ul><li><b>plain</b>: energy, all 87 force components and the 6 virial components at step 0 and after 4 MD steps, with newton on and off; then again after <code>write_restart</code>/<code>read_restart</code>, without the f·r virial shortcut, and after <code>write_data</code>/<code>read_data</code>.</li>
<li><b>single</b>: energy and force of a few pairs through single().</li>
<li><b>extract</b>: the arrays named in the file are reachable (fix adapt needs them, page E.4).</li>
<li>omp, kokkos, gpu, intel, opt: the accelerated versions of the style, skipped when they do not exist.</li></ul>
<p>Every comparison is a relative error, |a − b| / min(|a|, |b|), against the file's <code>epsilon</code>, here 5 × 10⁻¹⁴: a few units in the last place, enough for the summation order to change between newton on and off, not enough for anything else.</p></div>
  </section>

  <section class="block">
    <h2>2 · A reference from someone else's code</h2>
    <pre class="lmp-mini">${esc(SETUP)}</pre>
    <div class="prose"><p>A reference made by lj/mine itself would only prove that lj/mine agrees with lj/mine. So it is generated with mie/cut (as mie/ref, page E.3), every pair 14–6, the coefficients of LAMMPS's own lj/cut test:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.refyamlhead)}</pre>
    <div class="prose"><p>The resulting file for lj/mine (the numbers, all from mie/cut):</p></div>
    <pre class="lmp-mini ut-yaml"></pre>
  </section>

  <section class="block">
    <h2>3 · Match: lj/mine passes</h2>
    <pre class="lmp-mini">${esc(NATIVE.good)}</pre>
    <div class="prose"><p>The statistics (-s) show where the differences are:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.stats)}</pre>
    <div class="prose"><p>With newton on, forces and stress are identical to mie/cut (error 0): same arithmetic, same order. With newton off each pair is computed on both sides and sums come in a different order: errors up to 2 × 10⁻¹⁴ in the stress, still inside 5 × 10⁻¹⁴. Restart and data-file values are identical.</p></div>
  </section>

  <section class="block">
    <h2>4 · Observe: two planted bugs</h2>
    <div class="prose"><p><b>Bug 1</b>: the restart lines of page E.3 left out.</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.bug1)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.bug1out)}</pre>
    <div class="prose ut-b1"></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Reproduce the failure value</button></div>
    <div class="prose ut-a"></div>
    <div class="prose"><p><b>Bug 2</b>: a one-letter typo in single() (rn where rm belongs); compute() is untouched.</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.bug2)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.bug2out)}\n${esc(NATIVE.bug2fail)}</pre>
    <div class="prose"><p>Now plain passes, every MD number is right, and only the single test fails. Nothing in a normal run would show this bug; pair_write, compute pair/local and some fixes would quietly report wrong energies.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> How do you make your YAML file part of LAMMPS's test suite?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Put it in unittest/force-styles/tests/ with a name starting mol-pair-; CMake collects those files and registers one ctest test each, so <code>ctest -R lj_mine</code> (or plain <code>ctest</code>) runs it. A contribution of a new pair style to LAMMPS is expected to come with such a file.</div></div>
    <div class="ex"><p><b>2.</b> Why 5 × 10⁻¹⁴ and not 0?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Newton off, OpenMP threads, GPUs and different compilers all change the order of floating-point sums, and section 3 showed 2 × 10⁻¹⁴ from newton off alone. 5 × 10⁻¹⁴ is a few ulps; a real formula error, like bug 1, is off by 10¹⁴ times more.</div></div>
    <div class="ex"><p><b>3.</b> Bug 1 failed 87 + 6 comparisons twice. Why twice?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Two tests read the restart file: the plain restart check and the "nofdotr" check, which restarts again and recomputes the virial pair by pair instead of with the f·r shortcut. Both see Lennard-Jones forces where the reference has 14–6.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('.ut-yaml').textContent = NATIVE.yaml.split('\n').slice(0, 40).join('\n') + '\n…';
    const c = NATIVE.bug1counts;
    $('.ut-b1').innerHTML = `<p>${Object.entries(c).map(([k, v]) => `${v} × ${esc(k.replace('EXPECT_', ''))}`).join(', ')} and the restart energy:</p><pre class="lmp-mini">${esc(NATIVE.bug1energy)}</pre><p>After read_restart the style has its default exponents, 12 and 6, so the system is computed with Lennard-Jones. The data-file test passes, because it re-issues the pair_style line from the YAML file, exponents included: only a restart shows the bug.</p>`;
    $('[data-k=a]').onclick = () => {
      const mie = Number(NATIVE.mieE), lj = Number(NATIVE.ljcutE), err = Math.abs(lj - mie) / Math.min(Math.abs(lj), Math.abs(mie)), shown = 3.0167068520609837;
      $('.ut-a').innerHTML = `<p>If the restarted style is plain 12–6 LJ, its energy must be the one in LAMMPS's own lj/cut test file, mol-pair-lj_cut.yaml, which uses the same system and coefficients: init_vdwl = ${NATIVE.ljcutE}. The reference (14–6, from mie/cut) is ${NATIVE.mieE}. |a − b| / min(|a|, |b|) = ${s17(err)}, against the printed ${s17(shown)} ${err === shown ? OK : BAD}: the restarted lj/mine is exactly LAMMPS's lj/cut, so the diagnosis is certain. A failing test that you can reproduce from first principles tells you not only that something is wrong but what.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
