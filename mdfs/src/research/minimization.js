// Research · D.14: minimization. MinCG + linemin_quadratic reproduced bit for bit (with LAMMPS's neighbour-list pair order),
// a 4-atom LJ cluster where the default minimizer reports "energy tolerance" at |F| = 0.21 or freezes for 1000 iterations,
// the mechanism (a secant through two overshoot points, accepted by the EMACH slack), 1000 random starts, six minimizers,
// and a one-line guard tried in the hand code.
import { runLammps } from './lammps-runner.js';
import { minimizeCG, randomStarts, lmpScript, STOP_TEXT } from './min-hand.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const e3 = x => (x === 0 ? '0' : x.toExponential(3));

const P0 = [[5.0, 5.0, 5.0], [6.3, 5.1, 5.0], [5.5, 6.2, 5.1], [5.6, 5.5, 6.3]];
const NSTART = 1000;

// LAMMPS output → per-minimize summary
function summary(L) {
  const i = L.findIndex(x => /^\s*Step\s+PotEng/.test(x)), rows = [];
  if (i >= 0) for (let k = i + 1; k < L.length && !/^Loop/.test(L[k]); k++) { const t = L[k].trim().split(/\s+/); if (t.length === 3 && /^\d+$/.test(t[0])) rows.push(t.map(Number)); }
  const get = re => (L.find(x => re.test(x)) || '').split('= ')[1] || '';
  return { rows, stop: get(/Stopping criterion/), it: get(/Iterations, force/), alpha: get(/Final line search alpha/), builds: get(/Neighbor list builds/) };
}
// one worker call per chunk: clear + script + marker, repeated
async function batch(starts, opt, onProgress) {
  const out = [];
  for (let c = 0; c < starts.length; c += 100) {
    const part = starts.slice(c, c + 100);
    const sc = part.map(P => `clear\n${lmpScript(P, { thermo: 1000000, ...opt })}\nprint ZZEND`).join('\n');
    const res = await runLammps(sc);
    if (res.error) return { error: res.error };
    let cur = [];
    for (const l of res.lines) { if (l.trim() === 'ZZEND') { out.push(summary(cur)); cur = []; } else cur.push(l); }
    if (onProgress) onProgress(Math.min(c + 100, starts.length));
  }
  return { out };
}
// plain LJ along a line (all pairs), for the scan of the last line search
function eAlong(x, h, a) {
  const y = x.map((v, i) => v + a * h[i]), f = new Array(12).fill(0); let e = 0;
  for (let i = 0; i < 4; i++) for (let j = i + 1; j < 4; j++) {
    const d = [0, 1, 2].map(k => y[3 * i + k] - y[3 * j + k]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    if (rsq < 6.25) { const r2 = 1 / rsq, r6 = r2 * r2 * r2, fp = r6 * (48 * r6 - 24) * r2; for (let k = 0; k < 3; k++) { f[3 * i + k] += d[k] * fp; f[3 * j + k] -= d[k] * fp; } e += r6 * (4 * r6 - 4); }
  }
  let fh = 0; for (let i = 0; i < 12; i++) fh += f[i] * h[i];
  return { e, fh };
}

export default {
  id: 'minimization',
  title: 'Minimization',
  minutes: 25,
  summary: '"Stopping criterion = energy tolerance" does not mean "at a minimum". You rebuild the default minimizer, Polak–Ribière CG with the quadratic line search, down to the neighbour-list pair order, and match a 4-atom LJ cluster bit for bit. Then you watch it stop with "energy tolerance" at |F| = 0.21, or sit still for 1000 iterations, and find why: a secant through two overshoot points gives a step of 10⁻¹⁰, and the EMACH slack accepts it. 1000 random starts against live LAMMPS, six minimizers compared, and a one-line guard tried in the hand code.',
  learn: [
    'The stopping tests of minimize and what each message means',
    'MinCG: Polak–Ribière β, restart every 3N iterations',
    'linemin_quadratic: backtracking, the secant step α₀ and the EMACH slack',
    'Reading the summary: final Fnorm, and what "max atom move" really is',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: one outer loop, one line search</h2>
    <div class="prose"><p><code>minimize etol ftol maxiter maxeval</code> runs ${gh('src/min_cg.cpp', 38)} (the default <code>min_style cg</code>). Each iteration is one line search along h, then a new h:</p>
<pre class="lmp-mini">linemin(ecurrent, alpha_final);                       // move x along h
if (|E − E_prev| &lt; etol·½(|E| + |E_prev| + 1e-8))  return ETOL;   // "energy tolerance"
if (f·f &lt; ftol²)                                       return FTOL;   // "force tolerance" (2-norm)
beta = MAX(0, (f·f − f·g)/(g·g));                     // Polak–Ribière
if ((niter+1) % ndof == 0) beta = 0;                  // restart: ndof = 3N
g = f;   h = g + beta·h;   if (g·h &lt;= 0) h = g;</pre>
<p>The line search, ${gh('src/min_linesearch.cpp', 325)}, starts at α = min(1, dmax/max|h|) (dmax = 0.1) and halves α until the energy drops enough (backtracking). At every trial it also fits the slope f·h linearly between this α and the previous one and jumps to its zero:</p>
<pre class="lmp-mini">relerr = |1 − (½(α−α_prev)(fh + fh_prev) + E)/E_prev|;      // line 453
alpha0 = α − (α−α_prev)·fh/(fh − fh_prev);                   // line 454
if (relerr &lt;= 0.1 &amp;&amp; 0 &lt; alpha0 &lt; alphamax) {
  E = energy at alpha0;
  if (E − E_original &lt; EMACH) return success;               // line 458, EMACH = 1e-8
}
if (E − E_original &lt;= −0.4·α·(f·h)) return success;          // backtracking test
α_prev = α;  α *= 0.5;   if (0.4·α·(f·h) &lt;= EMACH) return ZEROALPHA;</pre>
<p>Two details matter below. The α₀ step is accepted when the energy has not risen by more than 10⁻⁸: it does not have to go down. And relerr is measured against the total energy E_prev, not against the change along the line.</p>
<p>The test system: four LJ atoms, cutoff 2.5σ. By hand, the minimum is a regular tetrahedron with every pair at r = 2<sup>1/6</sup>σ, where r⁻⁶ = ½ and each pair contributes 4(¼ − ½) = −1: <b>E_min = −6</b>, F = 0.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: the default minimizer, by hand</h2>
    <div class="prose"><p>The hand code (min-hand.js, 120 lines) does MinCG::iterate and linemin_quadratic line by line, and builds LAMMPS's neighbour list the same way (bins of 0.5·(2.5+0.3)σ, the half stencil, the sort at setup, a rebuild when an atom moves more than skin/2), because the order in which pair forces are added changes the last bit, and later iterations amplify it.</p></div>
    <pre class="lmp-mini">${esc(lmpScript(P0, { etol: '1.0e-10' }))}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run LAMMPS and the hand code</button><span class="lmp-status muted mn-sta"></span></div>
    <div class="tbl mn-a"></div>
    <div class="prose mn-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: a fixed point of the line search</h2>
    <div class="prose"><p>The last line search (iteration 17), from the hand code's log. E₀ and f·h(0) are the values at the start of the line:</p></div>
    <div class="tbl mn-ls"></div>
    <div class="prose mn-lssum"></div>
    <div class="prose"><p>With etol = 0 only ftol can stop it. Run 1000 iterations:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run with etol = 0</button><span class="lmp-status muted mn-stb"></span></div>
    <div class="tbl mn-b"></div>
    <div class="prose mn-bsum"></div>
  </section>

  <section class="block">
    <h2>4 · Modify: 1000 starts, six minimizers, one guard</h2>
    <div class="prose"><p>${NSTART} random starts: four atoms within ±0.8σ of the box centre, no pair closer than 0.9σ, coordinates on a 0.001σ grid (seeded, so every reader gets the same set). Each one: <code>minimize 0.0 1.0e-10 1000 10000</code> in LAMMPS and in the hand code.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run ${NSTART} starts</button><span class="lmp-status muted mn-stc"></span></div>
    <div class="tbl mn-c"></div>
    <div class="prose mn-csum"></div>
    <div class="row"><button class="btn sm primary" data-k="d" disabled>▶ The stalled starts with five other minimizers</button><span class="lmp-status muted mn-std"></span></div>
    <div class="tbl mn-d"></div>
    <div class="prose mn-dsum"></div>
    <div class="prose"><p>The hand code can do what the browser build cannot: change one line and rerun all ${NSTART}. Two candidate fixes for line 458: <b>strict</b> accepts α₀ only if E(α₀) &lt; E₀ (what the documentation says: "provided the energy change is downhill"); <b>bracket</b> keeps the EMACH slack but skips the α₀ step when both slopes are negative, so the secant is never extrapolated from two overshoot points. These are my changes to a copy, not LAMMPS.</p></div>
    <div class="row"><button class="btn sm primary" data-k="e">▶ Rerun all starts with each variant</button></div>
    <div class="tbl mn-e"></div>
    <div class="prose mn-esum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Iteration 17 changed E by 4.75 × 10⁻¹⁰. Show that the energy test with etol = 10⁻¹⁰ passes.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>etol·½(|E| + |E_prev| + 10⁻⁸) = 10⁻¹⁰ × ½ × 12.0 = 6.0 × 10⁻¹⁰ &gt; 4.75 × 10⁻¹⁰. The test compares the change with the size of E, not with the force: at |F| = 0.21 it already passes.</div></div>
    <div class="ex"><p><b>2.</b> Four atoms have ndof = 12. At which iterations does CG restart with h = f?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>β = 0 when (niter + 1) % 12 = 0, i.e. at the end of iterations 11, 23, 35 …, so iterations 12, 24, 36 … search along the force. Once stalled, β is 0 anyway: f·f − f·g = 0 because the force no longer changes.</div></div>
    <div class="ex"><p><b>3.</b> A log says "Stopping criterion = energy tolerance" and "Final line search alpha, max atom move = 0.042 0.0062". Converged?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>You cannot tell from these lines. "Max atom move" is the last <i>trial</i> α times the largest force component (finish.cpp:220), not the step taken. Look at "Force two-norm initial, final" against what you need, and if it is large, rerun with <code>min_modify line forcezero</code> or <code>min_style fire</code>.</div></div>
    <div class="ex"><p><b>4.</b> Why not simply use the strict test (E must go down)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Near the minimum a step lowers E by about |F|²α. At |F| = 10⁻⁷ that is 10⁻¹⁵ or less, the last bit of E = −6 (8.9 × 10⁻¹⁶). Energy can no longer show progress, so a strict test stops with "linesearch alpha is zero" around |F| ≈ 10⁻⁷. The EMACH slack is what lets the quadratic search converge forces past the resolution of the energy; the guard keeps that and only refuses the extrapolation.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };
    let H17 = null, starts = null, stalled = null;

    const showLineSearch = H => {
      const L = H.log[H.log.length - 1], E0 = H.rows[H.rows.length - 2][1];
      $('.mn-ls').innerHTML = `<table>${head(['trial', 'α', 'E(α) − E₀', 'f·h(α)', 'relerr', 'α₀', 'E(α₀) − E₀', ''])}${L.trials.map((t, k) => row([k + 1, e3(t.alpha), e3(t.de), e3(t.fh), e3(t.relerr), e3(t.alpha0), t.de0 === undefined ? '–' : e3(t.de0), t.took === 'alpha0' ? '<b>accepted (&lt; EMACH)</b>' : t.de0 !== undefined ? 'rejected, halve' : 'halve'])).join('')}</table>`;
      // scan the same line by hand: the real minimum, and how many more halvings backtracking needed
      const x0 = L.x0, h = L.h, t = L.trials, last = t[t.length - 1], base = eAlong(x0, h, 0);
      let lo = 0, hi = last.alpha;
      for (let k = 0; k < 200; k++) { const m = 0.5 * (lo + hi); if (eAlong(x0, h, m).fh > 0) lo = m; else hi = m; }
      const star = eAlong(x0, h, lo); let a = last.alpha, more = 0;
      do { a *= 0.5; more++; } while (eAlong(x0, h, a).e - base.e > -0.4 * a * L.fdothall && more < 30);
      const r = last.fh / last.fhprev;
      $('.mn-lssum').innerHTML = `<p>E₀ = ${s17(E0)}, f·h(0) = ${e3(L.fdothall)}, α_max = ${e3(L.alphamax)}. All five trials overshoot: E is higher than at the start and the slope f·h has turned negative. From trial 2 on, the secant runs through two such points, α and α_prev = 2α, both beyond the minimum. With r = fh/fh_prev, line 454 becomes α₀ = α(2r − 1)/(r − 1), and here r = ${r.toFixed(9)}: the line through the two points crosses zero almost at α = 0, α₀ = ${e3(last.alpha0)}. relerr = ${e3(last.relerr)} passes easily, because it is divided by the total energy (−6), not by the change along the line. The tiny step changes E by ${e3(last.de0)}, "less than EMACH", so the line search reports success.</p>
<p>By hand, bisecting f·h = 0 along the same line: the minimum is at α* = ${e3(lo)}, ${e3(base.e - star.e)} lower in energy. Backtracking alone would have accepted after ${more} more halving${more > 1 ? 's' : ''} (α = ${e3(a)}).</p>
<p>Why the steps shrink tenfold per iteration (the α column in section 2, from iteration 13): α₀ is proportional to 2r − 1, so each tiny step moves the atoms towards the place where r is exactly ½, and there α₀ = 0. The iterations converge, but to a fixed point of the line search, not of the energy.</p>`;
    };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.mn-sta'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(lmpScript(P0, { etol: '1.0e-10' })); b.disabled = false;
      if (res.error) return err(st, res.error); st.textContent = '';
      const L = summary(res.lines), H = minimizeCG(P0, { etol: 1e-10, keepLog: true }); H17 = H;
      let ok = 0;
      $('.mn-a').innerHTML = `<table>${head(['iter', 'PotEng (LAMMPS)', 'Fnorm (LAMMPS)', 'trials', 'step taken', 'α of that step', 'hand'])}${L.rows.map((r, k) => {
        const h = H.rows[k], same = h && h[0] === r[0] && h[1] === r[1] && h[2] === r[2]; ok += same; const lg = k > 0 ? H.log[k - 1] : null;
        return row([r[0], `<span class="mono">${s17(r[1])}</span>`, `<span class="mono">${s17(r[2])}</span>`, lg ? lg.trials.length : '', lg ? (lg.trials[lg.trials.length - 1].took === 'alpha0' ? 'α₀' : 'backtrack') : '', lg ? e3(lg.trials[lg.trials.length - 1].took === 'alpha0' ? lg.trials[lg.trials.length - 1].alpha0 : lg.trials[lg.trials.length - 1].alpha) : '', same ? OK : BAD]); }).join('')}</table>`;
      const sameStop = L.stop === STOP_TEXT[H.stop], sameIt = L.it === `${H.niter} ${H.neval}`, sameB = +L.builds === H.builds - 1;
      $('.mn-asum').innerHTML = `<p>${ok} of ${L.rows.length} rows identical to the last bit (PotEng and the thermo Fnorm). LAMMPS: "${esc(L.stop)}", iterations and force evaluations ${esc(L.it)}, ${esc(L.builds)} neighbour-list builds; hand: "${STOP_TEXT[H.stop]}" ${sameStop ? OK : BAD}, ${H.niter} ${H.neval} ${sameIt ? OK : BAD}, ${H.builds - 1} builds after setup ${sameB ? OK : BAD}.</p>
<p>So LAMMPS says it stopped on energy tolerance. But E = ${s17(H.e)} is ${e3(H.e + 6)} above the tetrahedron's −6, and |F| = ${H.fnorm2.toFixed(4)}: no minimum. From iteration 13 on, every line search needs five trials and ends with an α₀ step about ten times shorter than the one before; the energy changes shrink the same way (6 × 10⁻⁷, 6 × 10⁻⁸, 5 × 10⁻⁹, 5 × 10⁻¹⁰) until one falls below etol.</p>`;
      showLineSearch(H);
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.mn-stb'); b.disabled = true; st.textContent = 'Running 1000 iterations…';
      const res = await runLammps(lmpScript(P0, {})); b.disabled = false;
      if (res.error) return err(st, res.error); st.textContent = '';
      const L = summary(res.lines), H = minimizeCG(P0, { keepLog: true });
      let ok = 0; L.rows.forEach((r, k) => { const h = H.rows[k]; if (h && h[0] === r[0] && h[1] === r[1] && h[2] === r[2]) ok++; });
      const last = H.rows[H.rows.length - 1], frozen = H.rows.findIndex(r => r[1] === last[1] && r[2] === last[2]);
      const lg = H.log[H.log.length - 1], tr = lg.trials[lg.trials.length - 1], x = H.S.x, f = H.S.f;
      let moved = 0, maxmove = 0; for (let i = 0; i < 12; i++) { if (x[i] + tr.alpha0 * f[i] !== x[i]) moved++; maxmove = Math.max(maxmove, Math.abs(tr.alpha0 * f[i])); }
      let fmax = 0; for (const v of f) fmax = Math.max(fmax, Math.abs(v));
      $('.mn-b').innerHTML = `<table>${head(['iter', 'PotEng', 'Fnorm', 'hand'])}${[0, 1, 12, 17, frozen - 1, frozen, 500, 1000].map(k => { const r = L.rows[k], h = H.rows[k]; return row([r[0], `<span class="mono">${s17(r[1])}</span>`, `<span class="mono">${s17(r[2])}</span>`, h[1] === r[1] && h[2] === r[2] ? OK : BAD]); }).join('')}</table>`;
      $('.mn-bsum').innerHTML = `<p>All ${ok} of ${L.rows.length} rows identical; LAMMPS "${esc(L.stop)}", ${esc(L.it)}, hand ${H.niter} ${H.neval} ${L.it === `${H.niter} ${H.neval}` ? OK : BAD}. From iteration ${frozen} on, E and |F| do not change in a single bit. The last α₀ is ${e3(tr.alpha0)}; times the force (max ${fmax.toFixed(4)}) that is a move of ${e3(maxmove)}σ, below half the spacing of doubles near 5σ (4.4 × 10⁻¹⁶). By hand, x + α₀h equals x in ${12 - moved} of 12 coordinates. The atoms stay put, so f, β = 0, h = f and the next line search are the same: a fixed point, 7 force evaluations per iteration, forever.</p>
<p>The summary line is no help: <code>Final line search alpha, max atom move = ${esc(L.alpha)}</code>. That is the last <i>trial</i> α (${e3(tr.alpha)}) times the largest force component (finish.cpp:220). The step really taken, α₀, is not reported.</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.mn-stc'); b.disabled = true;
      starts = starts || randomStarts(NSTART, 1);
      const r = await batch(starts, {}, n => { st.textContent = `LAMMPS: ${n} of ${NSTART}…`; }); b.disabled = false;
      if (r.error) return err(st, r.error); st.textContent = '';
      let agree = 0; const cnt = {}; stalled = []; const stE = [];
      r.out.forEach((L, s) => {
        const H = minimizeCG(starts[s], {}), last = L.rows[L.rows.length - 1], hl = H.rows[H.rows.length - 1];
        if (L.stop === STOP_TEXT[H.stop] && L.it === `${H.niter} ${H.neval}` && last && last[1] === hl[1] && last[2] === hl[2]) agree++;
        cnt[L.stop] = (cnt[L.stop] || 0) + 1;
        if (L.stop === 'max iterations') { stalled.push(s); stE.push(last[2]); }
      });
      $('.mn-c').innerHTML = `<table>${head(['stopping criterion', 'starts'])}${Object.entries(cnt).map(([k, v]) => row([k, v])).join('')}</table>`;
      stE.sort((a, c) => a - c);
      $('.mn-csum').innerHTML = `<p>Hand code and LAMMPS agree on ${agree} of ${NSTART} starts: same stop, same iteration and evaluation counts, same final E and |F| to the last bit. ${stalled.length} starts (${(100 * stalled.length / NSTART).toFixed(1)} %) end on "max iterations" with |F| between ${stE[0].toFixed(3)} and ${stE[stE.length - 1].toFixed(2)}: every one frozen the same way, at a point where r = ½ (checked in the hand code). With etol = 10⁻⁴, the value in most example inputs, they would have reported "energy tolerance".</p>`;
      $('[data-k=d]').disabled = false;
    };

    $('[data-k=d]').onclick = async e => {
      const b = e.currentTarget, st = $('.mn-std'); b.disabled = true; const rows = [];
      const sub = stalled.map(s => starts[s]);
      for (const [style, line] of [['cg', 'backtrack'], ['cg', 'forcezero'], ['sd', ''], ['fire', ''], ['quickmin', '']]) {
        st.textContent = `${style} ${line}…`;
        const r = await batch(sub, { style, line, maxiter: 10000, maxeval: 100000 }); if (r.error) { b.disabled = false; return err(st, r.error); }
        const c = {}; let conv = 0, ev = 0;
        for (const L of r.out) { c[L.stop] = (c[L.stop] || 0) + 1; const last = L.rows[L.rows.length - 1]; if (last[2] < 1e-6) conv++; ev += +L.it.split(' ')[1]; }
        rows.push(row([`<code>${style}${line ? ' + ' + line : ''}</code>`, Object.entries(c).map(([k, v]) => `${k}: ${v}`).join('<br>'), `${conv} of ${sub.length}`, Math.round(ev / sub.length)]));
      }
      b.disabled = false; st.textContent = '';
      $('.mn-d').innerHTML = `<table>${head(['min_style (+ line)', 'stopping criterion', '|F| &lt; 10⁻⁶', 'evaluations (mean)'])}${rows.join('')}</table>`;
      $('.mn-dsum').innerHTML = '<p>Same starts, maxiter 10000. <code>forcezero</code> (a different line search) and the damped-dynamics minimizers <code>fire</code> and <code>quickmin</code> (page C.7) find E = −6 every time; plain backtracking stops early with |F| around 10⁻³ (page C.7 shows why); <code>sd</code> uses the same quadratic line search and freezes on some of them too. In this test, run over all 1000 starts, <code>forcezero</code> converged every one in about 80 evaluations; the default needed about 140 when it worked.</p>';
    };

    $('[data-k=e]').onclick = () => {
      starts = starts || randomStarts(NSTART, 1); const rows = [];
      for (const [name, v] of [['as in LAMMPS', false], ['strict: E(α₀) &lt; E₀', true], ['bracket guard', 'bracket']]) {
        const c = {}; let conv = 0, ev = 0;
        for (const P of starts) { const H = minimizeCG(P, { strict: v }); c[STOP_TEXT[H.stop]] = (c[STOP_TEXT[H.stop]] || 0) + 1; if (H.fnorm2 < 1e-6) conv++; ev += H.neval; }
        rows.push(row([name, Object.entries(c).map(([k, n]) => `${k}: ${n}`).join('<br>'), conv, Math.round(ev / starts.length)]));
      }
      $('.mn-e').innerHTML = `<table>${head(['line 458 variant', 'stopping criterion', '|F| &lt; 10⁻⁶', 'evaluations (mean)'])}${rows.join('')}</table>`;
      $('.mn-esum').innerHTML = `<p>The guard, <code>&amp;&amp; !(fh &lt; 0 &amp;&amp; fhprev &lt; 0)</code> added to the α₀ condition, ends all ${NSTART} starts on the force tolerance (|F| &lt; 10⁻¹⁰) and needs fewer evaluations on average, because the stalled starts no longer burn 7000 each. The strict test also escapes every stall, but then most starts stop with "linesearch alpha is zero" at |F| ≈ 10⁻⁷ instead of 10⁻¹⁰ (exercise 4). Four atoms are not a proof; in the hand code the guard also helped 7- and 13-atom clusters, which is the kind of evidence a bug report would start from.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
