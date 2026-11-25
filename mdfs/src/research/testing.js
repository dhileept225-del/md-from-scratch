// Research · F.3: regression testing. LAMMPS's own tester (tools/regression-tests/run_tests.py) run on twelve example
// inputs, then fooled twice: a harmless refactor it flags everywhere, and a real bug it misses almost everywhere.
// Every pass/fail decision, the ulp statistics, the Lyapunov exponent and the bug's energy shift are redone by hand.
import { NATIVE } from './testing-native.js';
import { crpow } from './crpow.js';
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const RT = 'tools/regression-tests/run_tests.py';

const CHECK = `# run_tests.py, for each run in the log (lines 622–660, shortened)
thermo_step = nthermo_steps - 1                  # only the LAST thermo line of each run
for each column:
    abs_diff = abs(float(val) - float(ref))
    if abs(float(ref)) > EPSILON:                # epsilon: 1e-16
        rel_diff = abs_diff / abs(float(ref))
    else:
        rel_diff = abs_diff / abs(float(ref) + nugget)   # nugget: 1.0
    if quantity in config['tolerance']:
        if abs_diff > abs_tol: FAILED             # both checks are counted
        if rel_diff > rel_tol: FAILED`;
const CONFIG = `# tools/regression-tests/config_serial.yaml (part)
  args: "-cite none"
  tolerance:
    PotEng: { abs: 1e-4, rel: 1e-7 }
    TotEng: { abs: 1e-4, rel: 1e-7 }
    Press:  { abs: 1e-4, rel: 1e-7 }
    Temp:   { abs: 1e-4, rel: 1e-7 }
    E_vdwl: { abs: 1e-3, rel: 1e-7 }
  timeout: 180`;
const CMD = `cd lammps
python3 tools/regression-tests/run_tests.py \\
    --lmp-bin=build/lmp \\
    --config-file=tools/regression-tests/config_serial.yaml \\
    --example-folders="examples/melt;examples/min;examples/crack;examples/friction;\\
examples/indent;examples/shear;examples/flow;examples/obstacle" --verbose
# writes failure.yaml, progress.yaml, output.xml (JUnit) and run.log`;
const POWDIFF = `--- a/src/pair_lj_cut.cpp        (PairLJCut::compute, line 115)
-        r6inv = r2inv * r2inv * r2inv;
+        r6inv = pow(r2inv, 3.0);`;
const OFFDIFF = `--- a/src/pair_lj_cut.cpp        (PairLJCut::init_one, line 516)
-    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));
+    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) + pow(ratio, 6.0));`;
const INFO = {
  'in.melt': ['lj/cut 2.5', '250 steps, T = 3 liquid'],
  'in.min': ['lj/cut 2.5, shift yes', '1000 steps at T ≈ 3.3, then minimize'],
  'in.min.box': ['lj/cut 2.5, shift yes', 'same, then 6 box minimizations'],
  'in.crack': ['lj/cut 2.5', '5000 steps, T = 0.01 solid'],
  'in.friction': ['lj/cut 2.5', '20 000 steps, T = 0.1'],
  'in.indent': ['lj/cut 2.5', '2 × 30 000 steps'],
  'in.indent.min': ['lj/cut 2.5', '7 minimizations'],
  'in.shear': ['eam', '3100 steps'],
  'in.shear.void': ['eam', '3100 steps'],
  'in.flow.couette': ['lj/cut 1.12246', '10 000 steps'],
  'in.flow.pois': ['lj/cut 1.12246', '10 000 steps'],
  'in.obstacle': ['lj/cut 1.12246', '25 000 steps'],
};
const verdict = x => (x.note || x.rows.some(r => r[3] === 'FAILED' || r[4] === 'FAILED'));
const cnt = k => NATIVE.runs[k].filter(verdict).length;

const CHAOS = d => `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 5 0 5 0 5
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 3.0 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
neighbor 0.3 bin
neigh_modify every 20 delay 0 check no
fix 1 all nve
${d ? `group one id 1\ndisplace_atoms one move ${d} 0 0\n` : ''}dump d all custom 50 b.dump id x y z
dump_modify d sort id format float %.17g
thermo 1200
run 1200`;
const MELT = `# 3d Lennard-Jones melt  (examples/melt/in.melt, unchanged)
units           lj
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
run             250`;

export default {
  id: 'testing',
  title: 'Unit tests and regression',
  minutes: 16,
  summary: 'Page E.11 covered the unit tests: one style, a few atoms, tight tolerances. This page is about the other half, LAMMPS\'s regression tester, which reruns the inputs in examples/ and compares the last thermo line with the log files shipped next to them. It runs cleanly on twelve inputs, then gets fooled twice: a mathematically harmless refactor fails eight of them, and a real bug passes ten. Both results are explained to the last digit.',
  learn: [
    'What run_tests.py compares: last thermo line, absolute and relative tolerance',
    'Why a 1-ulp change fails a chaotic run (and λ ≈ 5.4 per τ for the LJ liquid)',
    'Why a real bug can pass: coverage, not tolerance',
    'Where unit tests and regression tests each belong',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: what the tester compares</h2>
    <div class="prose"><p>Most folders in examples/ contain an input <code>in.name</code> and one or more reference logs, <code>log.DATE.name.g++.N</code>, produced by a developer when the example was last updated (N is the number of MPI ranks). ${gh(RT)} runs each input with the LAMMPS binary you give it, reads the thermo table of every run in the new log and in the reference log, and compares only <b>the last thermo line of each run</b>. In serial it picks the <code>.g++.1</code> log (${gh(RT, 285)}). A run fails on a column if the difference is larger than the absolute <i>or</i> the relative tolerance:</p></div>
    <pre class="lmp-mini">${esc(CHECK)}</pre>
    <pre class="lmp-mini">${esc(CONFIG)}</pre>
    <div class="prose"><p>Two more ways to fail: a different number of thermo lines (a minimizer that stopped at another step) and a different number of columns. With relative tolerance 1e-7, and reference values printed with 8 significant digits, the check effectively asks for agreement in the last printed digit.</p></div>
  </section>

  <section class="block">
    <h2>2 · Run it</h2>
    <pre class="lmp-mini">${esc(CMD)}</pre>
    <div class="tbl ts-base"></div>
    <div class="prose"><p>All ${NATIVE.runs.pass.length} inputs and ${NATIVE.runs.pass.reduce((s, x) => s + 2 * x.rows.length, 0)} checks pass, in ${NATIVE.regtime} s. The min logs date from November 2018 and melt's from April 2021: after 1000 steps of a hot liquid plus a minimization, a 2025 LAMMPS compiled with a different compiler still reproduces them in every printed digit. Section 3 shows how demanding that is.</p></div>
  </section>

  <section class="block">
    <h2>3 · Modify (1): a refactor that changes nothing</h2>
    <div class="prose"><p>Replace the three multiplications that make r⁻⁶ by a call to pow. Mathematically this is the same number; a reviewer would approve it on sight. Rebuild and rerun the twelve inputs:</p></div>
    <pre class="lmp-mini">${esc(POWDIFF)}</pre>
    <div class="tbl ts-pow"></div>
    <div class="prose"><p><b>${cnt('pow')} of 12 fail.</b> The two formulas round differently: <code>r2inv*r2inv*r2inv</code> rounds twice, <code>pow</code> once (glibc's pow gives the correctly rounded result). In ${(100 * NATIVE.powDiff / NATIVE.powN).toFixed(1)} % of ${NATIVE.powN.toExponential(0).replace('e+', 'e')} random pair distances between 0.8 and 2.5 σ they differ in the last bit (C program, native). One bit in a quarter of the pairs, at the first step, is all it takes.</p>
    <p>An MD trajectory is chaotic: two trajectories that start 10⁻¹⁵ apart separate like e<sup>λt</sup>. The same melt run with both builds for 1500 steps, positions dumped with 17 digits, RMS distance between the two trajectories:</p></div>
    <div class="tbl ts-div"></div>
    <div class="prose"><p>A factor of about 15 every 100 steps, from 10⁻¹⁴ to 10⁻² in 1000 steps (5 τ); after that the two runs are simply two different liquids. Hence the pattern in the table: melt passes because it stops at step 250, when the difference is still below 10⁻¹², far below the 8th digit; crack passes because a solid at T = 0.01 is almost harmonic and the difference barely grows; shear and shear.void pass because they use EAM and never execute the changed line. Every long lj/cut run fails.</p></div>
  </section>

  <section class="block">
    <h2>4 · Modify (2): a real bug that passes</h2>
    <div class="prose"><p>Restore the line, and plant a sign error in the energy shift that <code>pair_modify shift yes</code> subtracts from every pair, so that E(r<sub>c</sub>) = 0:</p></div>
    <pre class="lmp-mini">${esc(OFFDIFF)}</pre>
    <div class="tbl ts-off"></div>
    <div class="prose"><p><b>${12 - cnt('off')} of 12 pass.</b> The offset is computed only when shift is switched on, and only in.min and in.min.box switch it on; the other ten never execute the wrong line. Even in the two that do, Temp and Press pass: the offset is a constant per pair and does not enter the force, so the MD part of the trajectory is unchanged. Only the energy is wrong, and through the energy the minimizer, which stopped after ${NATIVE.off.minBug} instead of ${NATIVE.off.minOk} iterations ("mismatched num steps").</p>
    <p>The lesson is not about tolerances. A regression test can only see code that its inputs run, and the 800-odd inputs in examples/ were written to demonstrate features, not to cover branches.</p></div>
  </section>

  <section class="block">
    <h2>5 · Match by hand</h2>
    <div class="prose"><p><b>a.</b> All ${NATIVE.runs.pass.length * 3} test runs above (${['pass', 'pow', 'off'].reduce((s, k) => s + NATIVE.runs[k].reduce((t, x) => t + x.rows.length, 0), 0)} printed rows), every PASSED/FAILED label recomputed from the two printed numbers and the tolerances:</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Recompute the tester's decisions</button></div>
    <div class="prose ts-asum"></div>
    <div class="prose"><p><b>b.</b> The 1-ulp claim. The native program also printed 300 (r2inv, pow(r2inv, 3)) pairs as bit patterns. Check that glibc's pow is the correctly rounded cube, count how often the product of three differs, and draw 20 000 fresh distances here:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Compare pow and r·r·r</button></div>
    <div class="prose ts-bsum"></div>
    <div class="prose"><p><b>c.</b> λ from the native series, a straight-line fit of ln(RMS) against t = step × 0.005 over steps 100–1100. Then measure it yourself: two 500-atom melts in your browser, the second with atom 1 moved by δ before the first step.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Fit the native λ</button>
      <label class="lc-in">δ (σ) <select class="lc-u ts-d"><option>1e-15</option><option selected>1e-12</option><option>1e-9</option><option>1e-6</option></select></label>
      <button class="btn sm" data-k="c2">▶ Two runs in the browser</button><span class="ts-stc"></span></div>
    <div class="prose ts-csum"></div>
    <div class="tbl ts-ctab"></div>
    <div class="prose"><p><b>d.</b> The planted bug, to the last digit. Each pair inside the cutoff gets an extra −8ε(σ/r<sub>c</sub>)⁶ = −8 × 0.4⁶ = −0.032768, so PotEng per atom must drop by 0.032768 × N<sub>pairs</sub>/N. in.min rerun with <code>thermo_modify format float %.17g</code>, both builds:</p></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Count the pairs from the energy</button></div>
    <div class="prose ts-dsum"></div>
    <div class="prose"><p><b>e.</b> Run the regression test on your browser's LAMMPS: in.melt verbatim, compared with log.8Apr21.melt.g++.1 by the rule of section 1.</p></div>
    <pre class="lmp-mini">${esc(MELT)}</pre>
    <div class="row"><button class="btn sm primary" data-k="e">▶ Run in.melt and test it</button><span class="ts-ste"></span></div>
    <div class="prose ts-esum"></div>
  </section>

  <section class="block prose">
    <h2>6 · What each kind of test is for</h2>
    <ul><li><b>Unit tests</b> (page E.11) check that one style computes the right numbers: a few atoms, every function (compute, single, restart, newton on/off), tolerances near 1e-13. They are the test of <i>correctness</i>, and they run on every pull request on GitHub (unittest-linux.yml and friends).</li>
    <li><b>The regression tester</b> is a <i>change detector</i>. On any chaotic input, "passed" means "the same floating-point operations in the same order as when the log was made"; "failed" means "something in the executed code changed", harmless or not. It cannot tell the refactor of section 3 from a bug, and it cannot see code its inputs do not run.</li>
    <li>For a change you made on purpose, a failed regression test is the start of the work: confirm the change is intended, check the physics with something independent (a unit test, a conserved quantity, an average over a long run), then regenerate the reference logs (<code>--gen-ref</code>). <code>--quick</code> picks only the inputs that use styles your branch changed.</li>
    <li>For your own research code, the same two layers: a few tiny systems compared to many digits at step 0 and a few steps on (where chaos has not set in), and a handful of production-like inputs compared on averages with error bars, not on the last line.</li></ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> in.indent.min is a minimization with lj/cut, like in.min. Why did the offset bug not change it at all, not even the number of minimizer steps?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>It does not use <code>pair_modify shift yes</code>, so offset_flag is 0 and init_one takes the <code>else offset[i][j] = 0.0</code> branch. The wrong line is never executed; the binary computes exactly what the correct one does.</div></div>
    <div class="ex"><p><b>2.</b> in.melt passed the refactor test only because it stops at step 250. With the same input but <code>run 1000</code>, would TotEng pass? The native 17-digit TotEng of the two builds: step 800, −2.2848842478978826 and −2.2848843419787919; step 1000, −2.2861530727022896 and −2.2861482026300548.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. At step 800 the difference is 9.4e-8, relative 4.1e-8: still a pass (and invisible in 8 printed digits). At step 1000 it is 4.87e-6, relative 2.1e-6 &gt; 1e-7: fail. Between steps 800 and 1000 the harmless refactor becomes a "regression".</div></div>
    <div class="ex"><p><b>3.</b> Why do in.min and in.min.box report 3 and 7 checked rows on the refactor test, when the passing run printed 6 and 18?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The tester stops comparing at the first run whose number of thermo lines differs from the reference ("mismatched num steps") and marks the input failed. The minimizer started from a different state, converged after a different number of iterations, and its run table had a different length, so later runs are never compared.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    const res = k => NATIVE.runs[k].map(x => [x.inp, ...INFO[x.inp], verdict(x) ? `${BAD} ${x.note ? 'steps differ' : x.rows.filter(r => r[3] === 'FAILED' || r[4] === 'FAILED').length + ' rows'}` : `${OK} ${2 * x.rows.length} checks`]);
    $('.ts-base').innerHTML = `<table>${head(['input', 'pair style', 'what it runs', 'result'])}${res('pass').map(row).join('')}</table>`;
    $('.ts-pow').innerHTML = `<table>${head(['input', 'pair style', 'what it runs', 'pow refactor'])}${res('pow').map(row).join('')}</table>`;
    $('.ts-off').innerHTML = `<table>${head(['input', 'pair style', 'what it runs', 'offset bug'])}${res('off').map(row).join('')}</table>`;
    $('.ts-div').innerHTML = `<table>${head(['step', 't (τ)', 'RMS distance (σ)', 'growth per 100 steps'])}${NATIVE.rms.map((v, i) => row([(i + 1) * 100, ((i + 1) * 0.5).toFixed(1), v.toExponential(2), i ? (v / NATIVE.rms[i - 1]).toFixed(1) : ''])).join('')}</table>`;

    $('[data-k=a]').onclick = () => {
      const T = NATIVE.tol, EPS = 1e-16, nug = 1.0; let n = 0, ok = 0; const out = [];
      for (const k of ['pass', 'pow', 'off']) {
        let kn = 0, kok = 0;
        for (const x of NATIVE.runs[k]) for (const [q, val, ref, aS, rS] of x.rows) {
          const v = Number(val), r = Number(ref), abs = Math.abs(v - r), rel = Math.abs(r) > EPS ? abs / Math.abs(r) : abs / Math.abs(r + nug);
          const at = q === 'E_vdwl' ? T.eabs : T.abs, a = abs > at ? 'FAILED' : 'PASSED', b = rel > T.rel ? 'FAILED' : 'PASSED';
          kn++; n++; if (a === aS && b === rS) { kok++; ok++; }
        }
        out.push(`${k === 'pass' ? 'clean build' : k === 'pow' ? 'pow refactor' : 'offset bug'}: ${kok}/${kn}`);
      }
      const ip = NATIVE.runs.pow.find(x => x.inp === 'in.indent.min').rows.find(r => r[3] === 'PASSED' && r[4] === 'FAILED');
      const v = Number(ip[1]), r = Number(ip[2]);
      $('.ts-asum').innerHTML = `<p>${out.join(', ')}: ${ok} of ${n} labels identical ${ok === n ? OK : BAD}. One row shows why both checks exist: in in.indent.min on the refactor test, ${ip[0]} ${ip[1]} against ${ip[2]} differs by ${Math.abs(v - r).toExponential(2)}, well inside abs 1e-4, but relative ${(Math.abs(v - r) / Math.abs(r)).toExponential(2)} &gt; 1e-7: absolute PASSED, relative FAILED, and the input fails. For values near 1 the relative check is the stricter one; for values near zero (a Press of 1e-9) the absolute one saves the test from dividing by almost nothing.</p>`;
    };

    $('[data-k=b]').onclick = () => {
      const h2d = h => { const b = new DataView(new ArrayBuffer(8)); b.setBigUint64(0, BigInt('0x' + h)); return b.getFloat64(0); };
      const S = NATIVE.samp.split(' '); let cr = 0, dif = 0;
      for (const w of S) { const a = h2d(w.slice(0, 16)), p = h2d(w.slice(16)); if (crpow(a, 3) === p) cr++; if (a * a * a !== p) dif++; }
      let s = 0x9e3779b9 >>> 0; const rnd = () => { s ^= s << 13; s >>>= 0; s ^= s >>> 17; s ^= s << 5; s >>>= 0; return s / 4294967296; };
      let m = 0, up = 0; const M = 20000;
      for (let k = 0; k < M; k++) { const r = 0.8 + 1.7 * rnd(), x = 1.0 / (r * r), c = crpow(x, 3), q = x * x * x; if (c !== q) { m++; if (c > q) up++; } }
      $('.ts-bsum').innerHTML = `<p>glibc pow(r2inv, 3.0) equal to the correctly rounded cube: ${cr}/${S.length} ${cr === S.length ? OK : BAD}. The product of three differs from it in ${dif} of the ${S.length}. In your browser, ${M} new distances: ${m} differ (${(100 * m / M).toFixed(1)} %; native ${(100 * NATIVE.powDiff / NATIVE.powN).toFixed(1)} % of ${NATIVE.powN.toExponential(0).replace('e+', 'e')}), half of them up and half down (${up} up). Always by one unit in the last place, about 1e-16 relative: the smallest change a double can make.</p>`;
    };

    const fit = (t, y) => { const n = t.length; let sx = 0, sy = 0, sxx = 0, sxy = 0; for (let i = 0; i < n; i++) { sx += t[i]; sy += y[i]; sxx += t[i] * t[i]; sxy += t[i] * y[i]; } return (n * sxy - sx * sy) / (n * sxx - sx * sx); };
    $('[data-k=c]').onclick = () => {
      const r = NATIVE.rms.slice(0, 11), t = r.map((_, i) => (i + 1) * 100 * 0.005), lam = fit(t, r.map(Math.log));
      $('.ts-csum').innerHTML = `<p>Native λ = ${lam.toFixed(2)} τ⁻¹: the difference grows by e<sup>λ·0.5</sup> = ${Math.exp(lam * 0.5).toFixed(1)} per 100 steps, a factor of 10 every ${(Math.log(10) / lam).toFixed(2)} τ = ${Math.round(Math.log(10) / lam / 0.005)} steps. From one ulp (10⁻¹⁶ relative) to the 8th printed digit is 8 decades, about ${Math.round(8 * Math.log(10) / lam / 0.005)} steps; no tolerance a little looser would change the verdict on a long run, it would only move the step where it flips.</p>`;
    };
    const frames = txt => { const L = txt.split('\n'), F = []; for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const st = +L[i + 1], n = +L[i + 3], lo = L[i + 5].trim().split(/\s+/).map(Number), X = new Float64Array(3 * n); for (let a = 0; a < n; a++) { const w = L[i + 9 + a].trim().split(/\s+/); for (let c = 0; c < 3; c++) X[3 * a + c] = Number(w[1 + c]); } F.push({ st, X, Lb: lo[1] - lo[0] }); i += 8 + n; } return F; };
    $('[data-k=c2]').onclick = async e => {
      const b = e.currentTarget, st = $('.ts-stc'), d = $('.ts-d').value; b.disabled = true; st.textContent = ' Running 1 of 2…';
      const A = await runLammps(CHAOS(null), { read: ['b.dump'] }); if (!alive) return; st.textContent = ' Running 2 of 2…';
      const B = await runLammps(CHAOS(d), { read: ['b.dump'] }); if (!alive) return; b.disabled = false;
      if (A.error || B.error || !A.files['b.dump'] || !B.files['b.dump']) { st.innerHTML = `<span class="lmp-err">${esc(A.error || B.error || 'no dump')}</span>`; return; } st.textContent = '';
      const fa = frames(A.files['b.dump']), fb = frames(B.files['b.dump']), rows = [], tt = [], yy = [];
      for (let k = 0; k < fa.length; k++) {
        const X = fa[k].X, Y = fb[k].X, Lb = fa[k].Lb, n = X.length / 3; let s2 = 0;
        for (let j = 0; j < X.length; j++) { let dd = X[j] - Y[j]; dd -= Lb * Math.round(dd / Lb); s2 += dd * dd; }
        const rms = Math.sqrt(s2 / n); if (fa[k].st % 100 === 0) rows.push([fa[k].st, rms ? rms.toExponential(2) : '0']);
        if (fa[k].st >= 100 && rms > 0 && rms < 1e-3) { tt.push(fa[k].st * 0.005); yy.push(Math.log(rms)); }
      }
      const lam = tt.length > 2 ? fit(tt, yy) : NaN;
      $('.ts-ctab').innerHTML = `<table>${head(['step', `RMS distance, δ = ${d}`])}${rows.map(row).join('')}</table>`;
      $('.ts-csum').innerHTML = `<p>Your browser, 500 atoms, δ = ${d}: λ = ${lam.toFixed(2)} τ⁻¹ from ${tt.length} frames below 10⁻³ σ (native, 4000 atoms, a 1-ulp change in a quarter of the pairs: 5.43). λ is a property of the liquid, not of the size or the kind of perturbation. Try another δ: each factor of 1000 moves the moment the runs part by ln(1000)/λ ≈ ${Math.round(Math.log(1000) / (lam || 5.43) / 0.005)} steps, and nothing else changes.</p>`;
    };

    $('[data-k=d]').onclick = () => {
      const o = NATIVE.off, N = 800, c = 8 * crpow(0.4, 6);
      const p0 = (Number(o.ok['0']) - Number(o.bug['0'])) * N / c, p1 = (Number(o.ok['1000']) - Number(o.bug['1000'])) * N / c;
      // hand count on the sq2 lattice: 20 × 20 cells, basis (0,0) and (½,½), a = (2/ρ)^(1/2)
      const a = Math.sqrt(2 / 0.8442), P = []; for (let i = 0; i < 20; i++) for (let j = 0; j < 20; j++) { P.push([i * a, j * a]); P.push([(i + 0.5) * a, (j + 0.5) * a]); }
      const Lb = 20 * a; let np = 0, dmin = 9;
      for (let i = 0; i < P.length; i++) for (let j = i + 1; j < P.length; j++) { let dx = P[i][0] - P[j][0], dy = P[i][1] - P[j][1]; dx -= Lb * Math.round(dx / Lb); dy -= Lb * Math.round(dy / Lb); const r = Math.sqrt(dx * dx + dy * dy); if (r < 2.5) np++; else dmin = Math.min(dmin, r); }
      $('.ts-dsum').innerHTML = `<p>8 × 0.4⁶ = ${s17(c)} in doubles. Step 0: PotEng ${o.ok['0']} correct, ${o.bug['0']} with the bug; ΔE·N/0.032768 = ${p0.toFixed(9)} pairs. Step 1000: ${o.ok['1000']} and ${o.bug['1000']}, ${p1.toFixed(9)} pairs. Integers to within 10⁻⁸, as they must be if the bug does nothing but subtract a constant per pair inside the cutoff.</p>
      <p>Counted independently: on the initial sq2 lattice by hand here, ${np} pairs closer than 2.5 σ ${np === Math.round(p0) ? OK : BAD} (20 neighbours per atom; the next shell is at ${dmin.toFixed(3)} σ, safely outside). Natively, <code>compute pair/local dist</code> listed ${o.np0} pairs at step 0 ${o.np0 === Math.round(p0) ? OK : BAD} and ${o.np1} at step 1000 ${o.np1 === Math.round(p1) ? OK : BAD}, identical in both builds because the forces, and so the trajectory, are unchanged.</p>`;
    };

    $('[data-k=e]').onclick = async e => {
      const b = e.currentTarget, st = $('.ts-ste'); b.disabled = true; st.textContent = ' Running 4000 atoms × 250 steps…';
      const R = await runLammps(MELT); if (!alive) return; b.disabled = false;
      if (R.error || !R.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no thermo')}</span>`; return; } st.textContent = '';
      const cols = R.thermo.cols, last = R.thermo.rows[R.thermo.rows.length - 1], ref = NATIVE.runs.pass.find(x => x.inp === 'in.melt').rows;
      const T = NATIVE.tol; let allp = true, same = 0;
      const tr = ref.map(([q, val, rv]) => { const bv = String(last[cols.indexOf(q)]).trim(), v = Number(bv), r = Number(rv), abs = Math.abs(v - r), rel = abs / Math.abs(r); const p = abs <= T.abs && rel <= T.rel; allp = allp && p; if (bv === val) same++; return [q, bv, rv, abs.toExponential(1), rel.toExponential(1), p ? `${OK} PASSED` : `${BAD} FAILED`]; });
      $('.ts-esum').innerHTML = `<table>${head(['column', 'your browser', 'log.8Apr21.melt.g++.1', 'abs', 'rel', 'verdict'])}${tr.map(row).join('')}</table><p>${allp ? 'Passed' : 'Failed'}: your browser's WebAssembly LAMMPS ${allp ? 'passes' : 'fails'} the same regression test as the native build, and prints the same strings as the native run in ${same} of ${tr.length} columns. A 2021 reference, a 2025 source, a different compiler and a different instruction set: for this input the arithmetic, and its order, has not changed.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
