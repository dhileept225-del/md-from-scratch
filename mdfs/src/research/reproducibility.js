// Research · F.10: reproducibility. Three levels: bitwise (same binary, same input, same ranks), lost to rank count
// and compiler flags (measured divergence), and statistical (seed ensembles compared with a t-test). The browser
// reproduces the native runs bit for bit, reruns the ensembles, and shows a seed that changes the phase.
import { NATIVE } from './reproducibility-native.js';
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const mean = a => a.reduce((s, x) => s + x, 0) / a.length;
const sd = a => { const m = mean(a); return Math.sqrt(a.reduce((s, x) => s + (x - m) ** 2, 0) / (a.length - 1)); };
const avgFrom = (rows, col, from = 1000) => mean(rows.filter(r => Number(r[0]) >= from).map(r => Number(r[col])));

const BUILDS = `# the two builds compared below (same source, same machine)
cmake ../cmake -D CMAKE_BUILD_TYPE=Release                       # -O3 -DNDEBUG, the default
cmake ../cmake -D CMAKE_BUILD_TYPE=Release -D BUILD_MPI=no \\
      -D CMAKE_CXX_FLAGS_RELEASE="-O3 -ffast-math -march=native"   # lets the compiler reorder and fuse arithmetic`;
const CHECK = `# record with every result you publish
lmp -h                      # version, git commit, compiler, packages, sizes of tagint/bigint
cat build/CMakeCache.txt    # the compiler FLAGS: lmp -h does not show -O3 or -ffast-math
mpirun -np 4 lmp ...        # the rank count, and the processor grid in the log
# plus: the input, every data/potential file it reads (with checksums), and the seeds`;

export default {
  id: 'reproducibility',
  title: 'Reproducibility',
  minutes: 13,
  summary: 'What does it mean to reproduce an MD result? Three answers, each measured here. Bitwise: the same binary, input and rank count give the same bits, on this machine and in your browser. Lost: one more MPI rank or one compiler flag changes the last bit at step 0, and chaos makes it a different trajectory by step 1500. Statistical: averages over seeds agree within error bars, and a t-test says whether two builds differ. Plus a warning: near a phase boundary the seed can decide the phase.',
  learn: [
    'Determinism: same binary + input + ranks = same bits',
    'Why rank count and -ffast-math break it, and how fast',
    'Comparing builds or codes properly: ensembles, standard errors, a t-test',
    'What to record so that someone else can reproduce your numbers',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three levels</h2>
    <div class="prose"><p>LAMMPS has no hidden randomness: every random number comes from a seed in the input, and the arithmetic is the same every time the same binary runs the same input on the same number of ranks. That is <b>bitwise</b> reproducibility. Change the order of any floating-point sum (another rank count, another compiler flag, another code) and the last bit changes. Pages F.3 and F.5 showed that one bit grows until the trajectories have nothing in common. What survives is <b>statistical</b> reproducibility: the averages, within their error bars.</p></div>
    <pre class="lmp-mini">${esc(NATIVE.inp)}</pre>
    <div class="prose"><p>256 LJ atoms, NVE, started at T = 3 (so a liquid near T ≈ 1.64), 5000 steps. Run natively twice, with 1, 2 and 4 MPI ranks, and with a second build:</p></div>
    <pre class="lmp-mini">${esc(BUILDS)}</pre>
  </section>

  <section class="block">
    <h2>2 · Observe: where bits are lost</h2>
    <div class="tbl rp-div"></div>
    <div class="prose rp-divsum"></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run the same input in your browser</button><span class="rp-sta"></span></div>
    <div class="prose rp-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Match: statistical reproducibility</h2>
    <div class="prose"><p>The right comparison between two builds (or two codes) is not one trajectory against another but an ensemble against an ensemble. Ten seeds per build; for each run, the average potential energy over steps 1000–5000; then the mean over seeds, its standard error s/√10, and Welch's t for the difference.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Compare the two builds</button></div>
    <div class="prose rp-bsum"></div>
    <div class="tbl rp-btab"></div>
    <div class="row"><button class="btn sm" data-k="c">▶ Rerun the ten seeds in your browser</button><span class="rp-stc"></span></div>
    <div class="prose rp-csum"></div>
  </section>

  <section class="block">
    <h2>4 · Modify: when the seed decides the phase</h2>
    <div class="prose"><p>The same input started at T = 1.2 instead of 3.0 settles near T ≈ 0.65, close to the melting line at this density. The ten seeds natively:</p></div>
    <div class="tbl rp-bi"></div>
    <div class="prose rp-bisum"></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Run seeds 11 and 99 here</button><span class="rp-std"></span></div>
    <div class="prose rp-dsum"></div>
  </section>

  <section class="block">
    <h2>5 · What to record</h2>
    <pre class="lmp-mini">${esc(CHECK)}</pre>
    <div class="prose"><p>What <code>info config</code> printed for the default build (the -ffast-math build printed the same, apart from "OpenMP not enabled"):</p></div>
    <pre class="lmp-mini">${esc(NATIVE.info)}</pre>
    <div class="prose"><p>Nothing there distinguishes -O3 from -O3 -ffast-math -march=native. If bitwise reproduction matters to you (a regression test, a debugging session), keep the CMakeCache.txt with the binary. For a paper, report averages with error bars from independent runs, and the version and settings that produced them.</p><p>One trap found while building this page: when one LAMMPS instance runs several inputs in turn (the Python interface, or the browser here), <code>clear</code> does not delete variables. A second <code>variable seed index 11</code> is silently ignored, and every "new" seed reruns the first one. Use equal-style variables or <code>variable seed delete</code> before redefining.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> With 2 ranks the potential energy already differs at step 0, before any atom has moved. Where does that difference come from?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>From the final sum. Each rank adds the energies of its own pairs, then MPI_Allreduce adds the rank totals. The pairs are the same, the forces identical, but the sum is grouped differently: (a + b) + c is not always a + (b + c) in floating point. The difference is 10⁻¹³ in a total of −1734, about one bit. The forces differ too, by the same kind of rounding where ghost forces are summed back, and that is what then grows.</div></div>
    <div class="ex"><p><b>2.</b> The fast-math build was 14 % faster on this test. Would you use it for production?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Possibly, after checking. -ffast-math allows reassociation and fused multiply-adds (harmless here: the averages agree), but also assumes there are no NaN or infinity, which can break code that checks for them (page F.4's NaN would silently become something else). Validate the observables you care about with an ensemble comparison like section 3, and keep the default build for debugging and regression tests.</div></div>
    <div class="ex"><p><b>3.</b> A colleague's code gives PE = −4.735 for this system and yours −4.751. Different physics?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Not from one run each: single runs here spread from −4.728 to −4.769, a standard deviation of about 0.01. With ten seeds each the standard error is about 0.003, and a difference of 0.016 would then be about 4 standard errors, significant. One run against one run can only tell you that the codes are not wildly different.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    // equal-style, not index: the browser keeps one LAMMPS instance, and index variables survive 'clear' (a second definition is ignored)
    const SEED = (inp, s, T0) => inp.replace('variable seed index 87287', `variable seed equal ${s}`).replace('variable T0 index 3.0', `variable T0 equal ${T0}`);
    const A = NATIVE.a1, rows = A.map((r, i) => [r[0], NATIVE.np2[i][2], NATIVE.np4[i][2], NATIVE.fm[i][2]].map(Number));
    const d = (k, i) => Math.abs(rows[i][k] - Number(A[i][2]));
    const pick = [0, 1, 5, 10, 15, 20, 30, 50];
    $('.rp-div').innerHTML = `<table>${head(['step', '|ΔPE| total: 2 ranks', '4 ranks', '-ffast-math build'])}${pick.map(i => row([A[i][0], d(1, i).toExponential(1), d(2, i).toExponential(1), d(3, i).toExponential(1)])).join('')}</table>`;
    $('.rp-divsum').innerHTML = `<p>Run twice with the same binary and 1 rank: every thermo value identical (the logs differ only in the timings). With 2 or 4 ranks, or the other build, the total potential energy differs by ~10⁻¹³ already at step 0, at the level of the last bit of −1734. It is 10⁻¹¹ at step 500, 10⁻⁵ at step 1000, and of order 0.1, as large as the natural fluctuations, by step 1500. After that the runs are independent samples of the same liquid.</p>`;

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.rp-sta'); b.disabled = true; st.textContent = ' Running…';
      const R = await runLammps(SEED(NATIVE.inp, 87287, '3.0')); if (!alive) return; b.disabled = false; st.textContent = '';
      if (R.error) { st.innerHTML = `<span class="lmp-err">${esc(R.error)}</span>`; return; }
      let same = 0, n = 0, sf = 0; R.thermo.rows.forEach((r, i) => r.forEach((v, c) => { n++; if (String(v) === A[i][c]) same++; if (String(v) === NATIVE.fm[i][c]) sf++; }));
      $('.rp-asum').innerHTML = `<p>WebAssembly in your browser against native x86-64 with 1 rank: ${same} of ${n} thermo strings identical ${same === n ? OK : BAD}, over all 5000 steps. Against the -ffast-math build: ${sf} of ${n}. Bitwise reproducibility needs the same arithmetic in the same order, not the same machine. This browser build happens to be compiled from the same source without reordering flags, and IEEE-754 doubles round the same way everywhere.</p>`;
    };

    const stats = arr => { const v = arr.map(x => x[0]); return { m: mean(v), s: sd(v), se: sd(v) / Math.sqrt(v.length) }; };
    $('[data-k=b]').onclick = () => {
      const a = stats(NATIVE.n), f = stats(NATIVE.f), t = (a.m - f.m) / Math.sqrt(a.se ** 2 + f.se ** 2);
      const nu = (a.se ** 2 + f.se ** 2) ** 2 / (a.se ** 4 / 9 + f.se ** 4 / 9);
      $('.rp-btab').innerHTML = `<table>${head(['seed', 'default build', '-ffast-math build'])}${NATIVE.seeds.map((s, i) => row([s, NATIVE.n[i][0].toFixed(5), NATIVE.f[i][0].toFixed(5)])).join('')}${row(['<b>mean ± SE</b>', `${a.m.toFixed(5)} ± ${a.se.toFixed(5)}`, `${f.m.toFixed(5)} ± ${f.se.toFixed(5)}`])}</table>`;
      $('.rp-bsum').innerHTML = `<p>Difference ${(a.m - f.m).toFixed(5)} per atom; Welch t = ${t.toFixed(2)} with ${nu.toFixed(1)} degrees of freedom. |t| &lt; 2.1 (the 5 % level for ~17 dof) ${Math.abs(t) < 2.1 ? OK : BAD}: no evidence that the builds differ. A single pair of runs (seed 11) differs by ${Math.abs(NATIVE.n[0][0] - NATIVE.f[0][0]).toFixed(4)}, ${(Math.abs(NATIVE.n[0][0] - NATIVE.f[0][0]) / a.se).toFixed(1)} standard errors of the mean, which would look like a real difference if you had only those two numbers. Single runs scatter by s = ${a.s.toFixed(4)}; ten of them pin the mean to ±${a.se.toFixed(4)}.</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.rp-stc'); b.disabled = true; const got = [];
      for (const [i, s] of NATIVE.seeds.entries()) {
        st.textContent = ` seed ${s} (${i + 1}/10)…`;
        const R = await runLammps(SEED(NATIVE.inp, s, '3.0')); if (!alive) return;
        if (R.error) { st.innerHTML = `<span class="lmp-err">${esc(R.error)}</span>`; b.disabled = false; return; }
        got.push(avgFrom(R.thermo.rows, 2));
      }
      b.disabled = false; st.textContent = '';
      const same = got.filter((x, i) => x === NATIVE.n[i][0] || Math.abs(x - NATIVE.n[i][0]) < 1e-9).length;
      $('.rp-csum').innerHTML = `<p>Your browser's ten averages: ${got.map(x => x.toFixed(5)).join(', ')}. ${same} of 10 equal to the native default build's ${same === 10 ? OK : BAD}: the same ensemble, bit for bit, so the statistics above are your statistics too.</p>`;
    };

    const bi = NATIVE.bi;
    $('.rp-bi').innerHTML = `<table>${head(['seed', '⟨PE⟩ per atom', '⟨T⟩', 'state'])}${NATIVE.seeds.map((s, i) => row([s, bi[i][0].toFixed(4), bi[i][1].toFixed(3), bi[i][0] > -5.9 ? 'liquid' : 'crystal (so far)'])).join('')}</table>`;
    const nl = bi.filter(x => x[0] > -5.9).length;
    $('.rp-bisum').innerHTML = `<p>${10 - nl} seeds are still crystalline at the end (PE ≈ −5.96, T ≈ 0.65; seed 11's mean-square displacement after 5000 steps is 0.056 σ², vibration only), ${nl} have melted (PE ≈ −5.86, T ≈ 0.58: melting takes latent heat out of the kinetic energy at constant total energy; seed 99's MSD is 2.1 σ² and growing). The crystal is metastable here, and given more steps more seeds might melt. Same input, same total energy, two phases. Averaging over all ten would give a number that describes neither. With a seed ensemble you see the bimodality and can ask the right question (metastability, nucleation, a too-small box). With one run you would report whichever phase your seed happened to give.</p>`;

    $('[data-k=d]').onclick = async e => {
      const b = e.currentTarget, st = $('.rp-std'); b.disabled = true; const out = [];
      for (const s of [11, 99]) {
        st.textContent = ` seed ${s}…`;
        const R = await runLammps(SEED(NATIVE.inp, s, '1.2')); if (!alive) return;
        if (R.error) { st.innerHTML = `<span class="lmp-err">${esc(R.error)}</span>`; b.disabled = false; return; }
        out.push([s, avgFrom(R.thermo.rows, 2), avgFrom(R.thermo.rows, 1)]);
      }
      b.disabled = false; st.textContent = '';
      const nat = s => bi[NATIVE.seeds.indexOf(s)];
      $('.rp-dsum').innerHTML = `<p>${out.map(([s, p, t]) => `Seed ${s}: ⟨PE⟩ = ${p.toFixed(6)}, ⟨T⟩ = ${t.toFixed(4)} (native ${nat(s)[0].toFixed(6)}, ${nat(s)[1].toFixed(4)}) ${Math.abs(p - nat(s)[0]) < 1e-9 ? OK : BAD}`).join('. ')}. One is still a crystal after 5000 steps, the other has melted, in your browser exactly as natively.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
