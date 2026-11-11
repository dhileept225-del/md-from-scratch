// Research · B.3: styles, the plug-in design. Base class sockets, registration, the call order, and single() checked by hand through pair_write.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
// C's %.15g
const g15 = x => {
  if (x === 0) return '0';
  const ex = Math.floor(Math.log10(Math.abs(Number(x.toPrecision(15)))));
  if (ex < -4 || ex >= 15) { let [m, e] = x.toExponential(14).split('e'); m = m.replace(/\.?0+$/, ''); const n = Math.abs(+e); return `${m}e${+e < 0 ? '-' : '+'}${n < 10 ? '0' + n : n}`; }
  return String(Number(x.toPrecision(15)));
};

// the call order for pair_style lj/cut, with the real lines
const TRACE = [
  { when: 'pair_style lj/cut 2.5', fn: 'Input::pair_style() → Force::create_pair("lj/cut") → new PairLJCut → settings()', file: ['src/pair_lj_cut.cpp', 422],
    code: `422  void PairLJCut::settings(int narg, char **arg)
424    if (narg != 1) error->all(FLERR, "Illegal pair_style command");
426    cut_global = utils::numeric(FLERR, arg[0], false, lmp);   // 2.5` },
  { when: 'pair_coeff 1 1 1.0 1.0 2.5', fn: 'coeff(): store ε, σ, cutoff for the type pair', file: ['src/pair_lj_cut.cpp', 442],
    code: `451    double epsilon_one = utils::numeric(FLERR, arg[2], false, lmp);
452    double sigma_one = utils::numeric(FLERR, arg[3], false, lmp);
454    double cut_one = cut_global;
455    if (narg == 5) cut_one = utils::numeric(FLERR, arg[4], false, lmp);
       …
460        epsilon[i][j] = epsilon_one;
461        sigma[i][j] = sigma_one;
462        cut[i][j] = cut_one;
463        setflag[i][j] = 1;                  // "this pair is set" (page A.9, error 4)` },
  { when: 'run (setup, once)', fn: 'init_one(i, j): mix, precompute, offset; return the cutoff', file: ['src/pair_lj_cut.cpp', 501],
    code: `504    epsilon[i][j] = mix_energy(...)          // only if i j was not set (page A.7)
509    lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
       …                                       // lj2 = 24εσ⁶, lj3 = 4εσ¹², lj4 = 4εσ⁶
515    double ratio = sigma[i][j] / cut[i][j];
516    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));   // only with shift yes
556    return cut[i][j];                         // the neighbour list uses this` },
  { when: 'run (every step)', fn: 'compute(eflag, vflag): all forces, energy, virial', file: ['src/pair_lj_cut.cpp', 69],
    code: `the double loop over atoms and neighbours that pages A.7 and B.1 showed.
Verlet calls it without knowing which style it is (src/verlet.cpp:140):
140    if (pair_compute_flag) force->pair->compute(eflag,vflag);` },
  { when: 'pair_write …', fn: 'single(): energy and force of ONE pair at one distance', file: ['src/pair_lj_cut.cpp', 660],
    code: `665    r2inv = 1.0 / rsq;
666    r6inv = r2inv * r2inv * r2inv;
667    forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);
668    fforce = factor_lj * forcelj * r2inv;
670    philj = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];
671    return factor_lj * philj;` },
];

// three plugs for the same socket, each with single() and init_one() written out as in the C++
const PLUGS = {
  'lj/cut': { args: [['ε', 1.0], ['σ', 1.0]], rc: 2.5, file: ['src/pair_lj_cut.cpp', 660],
    formula: 'E = 4ε[(σ/r)¹² − (σ/r)⁶],  F = 24ε[2(σ/r)¹² − (σ/r)⁶]/r',
    coeff: (a, rc) => `${a[0]} ${a[1]} ${rc}`,
    cut: (a, rc) => rc,
    init: ([e, s], rc, sh) => ({ lj1: 48 * e * s ** 12, lj2: 24 * e * s ** 6, lj3: 4 * e * s ** 12, lj4: 4 * e * s ** 6, off: sh ? 4 * e * ((s / rc) ** 12 - (s / rc) ** 6) : 0 }),
    single: (P, rsq) => { const r2 = 1 / rsq, r6 = r2 * r2 * r2; return { f: r6 * (P.lj1 * r6 - P.lj2) * r2, e: r6 * (P.lj3 * r6 - P.lj4) - P.off }; } },
  'lj/expand': { args: [['ε', 1.0], ['σ', 1.0], ['Δ', 0.2]], rc: 2.5, file: ['src/pair_lj_expand.cpp', 396],
    formula: 'E = 4ε[(σ/(r−Δ))¹² − (σ/(r−Δ))⁶],  cutoff = rc + Δ',
    coeff: (a, rc) => `${a[0]} ${a[1]} ${a[2]} ${rc}`,
    cut: (a, rc) => rc + a[2],
    // the energy at the cutoff r = rc + Δ is 4ε[(σ/rc)¹² − (σ/rc)⁶]; LAMMPS line 239 uses σ/(rc + Δ) instead
    init: ([e, s, d], rc, sh, code) => { const q = code ? s / (rc + d) : s / rc; return { d, lj1: 48 * e * s ** 12, lj2: 24 * e * s ** 6, lj3: 4 * e * s ** 12, lj4: 4 * e * s ** 6, off: sh ? 4 * e * (q ** 12 - q ** 6) : 0 }; },
    single: (P, rsq) => { const r = Math.sqrt(rsq), rs = r - P.d, r2 = 1 / (rs * rs), r6 = r2 * r2 * r2; return { f: r6 * (P.lj1 * r6 - P.lj2) / rs / r, e: r6 * (P.lj3 * r6 - P.lj4) - P.off }; } },
  morse: { args: [['D₀', 1.0], ['α', 2.0], ['r₀', 1.1]], rc: 2.5, file: ['src/pair_morse.cpp', 338],
    formula: 'E = D₀[e^(−2α(r−r₀)) − 2e^(−α(r−r₀))],  F = 2αD₀[e^(−2α(r−r₀)) − e^(−α(r−r₀))]',
    coeff: (a, rc) => `${a[0]} ${a[1]} ${a[2]} ${rc}`,
    cut: (a, rc) => rc,
    init: ([D, al, r0], rc, sh) => { const adr = -al * (rc - r0); return { D, al, r0, m1: 2 * D * al, off: sh ? D * (Math.exp(2 * adr) - 2 * Math.exp(adr)) : 0 }; },
    single: (P, rsq) => { const r = Math.sqrt(rsq), dx = Math.exp(-P.al * (r - P.r0)); return { f: P.m1 * (dx * dx - dx) / r, e: P.D * (dx * dx - 2 * dx) - P.off }; } },
};

const script = s => `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10
create_box      1 box
mass            1 1.0
pair_style      ${s.style} ${s.rc}
pair_coeff      1 1 ${PLUGS[s.style].coeff(s.a, s.rc)}
${s.shift ? 'pair_modify     shift yes\n' : ''}pair_write      1 1 ${s.n} r ${s.lo} ${s.hi} pw.table PW
`;

export default {
  id: 'styles',
  title: 'Styles: the plug-in design',
  minutes: 20,
  summary: 'How LAMMPS can have 1400 styles without the main loop knowing any of them: every style is a C++ class that fills the same set of sockets (settings, coeff, init_one, compute, single). You follow one pair style through its sockets in the order LAMMPS calls them, then use pair_write to ask a style for its energy and force at chosen distances and match every number by hand, for three different styles.',
  learn: [
    'A base class declares the sockets; a style fills them (inheritance and virtual functions, in plain words)',
    'How a style registers its name, and how pair_style finds it',
    'The order the sockets are called: settings → coeff → init_one → compute every step',
    'Call single() with pair_write and match energy and force for lj/cut, lj/expand and morse by hand',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: sockets and plugs</h2>
    <p>Think of a wall socket. Any lamp, kettle or charger works in it, because the socket fixes the shape and every device is made to fit. In LAMMPS the <b>socket</b> is a <b>base class</b> such as <code>Pair</code>; each pair <b>style</b> is a <b>plug</b>: a class built from <code>Pair</code> that fills in its functions. The main loop only ever talks to the socket.</p>
    <p>The socket is a list of functions in ${gh('src/pair.h', 158)}. Those ending in <code>= 0</code> have no default: every style <i>must</i> write them.</p>
<pre class="lmp-mini">158  virtual void compute(int, int) = 0;                     // all forces, every step
164  virtual double single(int, int, int, int, double, ...)   // one pair (default: "not available")
187  virtual void settings(int, char **) = 0;                 // the pair_style arguments
188  virtual void coeff(int, char **) = 0;                    // the pair_coeff arguments
192  virtual double init_one(int, int) { return 0.0; }        // set up one type pair, return its cutoff</pre>
    <p>The plug says which socket it fits and which functions it fills, in ${gh('src/pair_lj_cut.h', 27)}:</p>
<pre class="lmp-mini">27  class PairLJCut : public Pair {          // "built from Pair"
31    void compute(int, int) override;        // "my own version of compute"
32    void settings(int, char **) override;
33    void coeff(int, char **) override;
35    double init_one(int, int) override;
42    double single(int, int, int, int, double, double, double, double &) override;</pre>
    <div class="tbl wide"><table>${head(['Base class (socket)', 'File', 'Styles (plugs) in 22 Jul 2025'])}
${row(['<code>Pair</code>', 'pair.h', 'pair_style: lj/cut, morse, eam, tersoff … (548)'])}
${row(['<code>Bond</code>, <code>Angle</code>, <code>Dihedral</code>, <code>Improper</code>', 'bond.h …', 'harmonic, fene, opls … (158)'])}
${row(['<code>KSpace</code>', 'kspace.h', 'ewald, pppm … (34)'])}
${row(['<code>Fix</code>', 'fix.h', 'nve, nvt, npt, langevin, rigid … (336)'])}
${row(['<code>Compute</code>', 'compute.h', 'temp, pe, rdf, msd … (176)'])}
${row(['<code>AtomVec</code>, <code>Region</code>, <code>Dump</code>, <code>Integrate</code>, <code>Min</code>, <code>Command</code>', '…', 'the rest of the 1412 (page B.1)'])}</table></div>
  </section>

  <section class="block prose">
    <h2>2 · Registration: how "lj/cut" finds PairLJCut</h2>
    <ol>
      <li>The header has one line: <code>PairStyle(lj/cut,PairLJCut);</code> (${gh('src/pair_lj_cut.h', 16)}).</li>
      <li>At build time CMake collects these lines from every header of the core and the switched-on packages into one generated file, <code>style_pair.h</code>. A package that is off contributes nothing: that is why <code>lj/smooth/linear</code> was missing on page B.1.</li>
      <li>When LAMMPS starts, ${gh('src/force.cpp', 92)} turns each line into an entry of a dictionary: name → "a function that does <code>new PairLJCut(lmp)</code>" (${gh('src/force.cpp', 41)}).</li>
      <li><code>pair_style lj/cut 2.5</code> looks the name up (${gh('src/force.cpp', 270)}). Found: the class is created and its <code>settings()</code> receives <code>2.5</code>. Not found: the "Unrecognized pair style" error.</li>
    </ol>
    <p>So adding a new style means: write one class, add one PairStyle line, rebuild. No other file of LAMMPS changes. Part E does exactly that.</p>
  </section>

  <section class="block">
    <h2>3 · The order LAMMPS calls the sockets</h2>
    <p class="muted">Step through the life of <code>pair_style lj/cut</code> in in.melt.</p>
    <div class="stepper st3"><div class="row"><button class="btn sm ghost" data-t="-1">◀</button><span class="st3-n mono"></span><button class="btn sm ghost" data-t="1">▶</button></div>
      <p class="st3-when"></p><pre class="lmp-mini st3-code"></pre></div>
  </section>

  <section class="block">
    <h2>4 · Modify: ask a style for E and F with pair_write</h2>
    <div class="prose"><p><code>pair_write</code> calls <code>single()</code> at evenly spaced distances and writes r, E, F to a file. It needs only a box, no atoms. Choose a plug, change its numbers, and predict the table by hand. The hand side below runs the same arithmetic as the C++ <code>init_one()</code> and <code>single()</code>.</p></div>
    <div class="row"><label class="lab" style="margin:0">pair style</label><select class="lc-u s3-style">${Object.keys(PLUGS).map(k => `<option>${k}</option>`).join('')}</select>
      <label class="rd-sp-lab" style="display:inline-block"><input type="checkbox" class="s3-shift">pair_modify shift yes</label>
      <label class="rd-sp-lab s3-codel" style="display:none"><input type="checkbox" class="s3-code">hand: use the offset exactly as written in pair_lj_expand.cpp line 239</label></div>
    <div class="lc-in s3-args"></div>
    <div class="lc-in">
      <label>cutoff r<sub>c</sub> <input type="number" class="lmp-in sm" data-q="rc" value="2.5" step="0.1" min="1" max="5"></label>
      <label>points <input type="number" class="lmp-in sm" data-q="n" value="6" min="2" max="10"></label>
      <label>from r <input type="number" class="lmp-in sm" data-q="lo" value="0.9" step="0.1" min="0.5"></label>
      <label>to r <input type="number" class="lmp-in sm" data-q="hi" value="2.6" step="0.1"></label>
    </div>
    <p class="muted s3-formula"></p>
    <pre class="lmp-mini s3-script"></pre>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run pair_write</button><span class="lmp-status muted s3-st"></span></div>
    <div class="tbl s3-cmp"></div>
  </section>

  <section class="block prose">
    <h2>A surprise: lj/expand with shift yes</h2>
    <p>Choose <b>lj/expand</b>, tick <b>shift yes</b> and run. The energies do not match. The hand side subtracts the energy at the cutoff, which is what <code>shift yes</code> is meant to do (the manual: the energy goes to 0 at the cutoff). For lj/expand the cutoff is r = r<sub>c</sub> + Δ, where r − Δ = r<sub>c</sub>, so the offset should be 4ε[(σ/r<sub>c</sub>)¹² − (σ/r<sub>c</sub>)⁶].</p>
    <p>The code computes something else (${gh('src/pair_lj_expand.cpp', 239)}):</p>
<pre class="lmp-mini">239    double ratio = sigma[i][j] / (cut[i][j] + shift[i][j]);     // σ/(rc + Δ), not σ/rc
240    offset[i][j] = 4.0 * epsilon[i][j] * (pow(ratio, 12.0) - pow(ratio, 6.0));</pre>
    <p>Tick "use the offset exactly as written" and the hand values match LAMMPS to 15 digits. With ε = σ = 1, Δ = 0.2, r<sub>c</sub> = 2.3, the shifted energy just inside the cutoff (r = 2.4999) is −0.0105 instead of 0, while lj/cut goes to 0 as it should. The forces are not affected, only the reported energy. The same line is in the current development version. This is what a hand calculation is for: it finds the one line where the code and the formula disagree.</p>
  </section>

  <section class="block prose">
    <h2>What the sockets show</h2>
    <ul>
      <li>The last row is always 0 when r reaches the cutoff: LAMMPS calls <code>single()</code> only if r² &lt; r<sub>c</sub>² (${gh('src/pair.cpp', 1943)}). At exactly r<sub>c</sub> the pair does not interact.</li>
      <li><code>shift yes</code> subtracts one constant, computed once in <code>init_one()</code>, from every energy. The force column does not change at all: a constant has no slope.</li>
      <li>The F column is the force itself. Inside <code>compute()</code> LAMMPS uses F/r (called <code>fpair</code>), so it can multiply by Δx, Δy, Δz without a square root; <code>pair_write</code> multiplies back by r.</li>
      <li>Three different plugs, one socket: <code>pair_write</code>, <code>compute()</code>'s callers and the rest of LAMMPS never needed to know which one was in.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> lj/cut with ε = 1, σ = 1, r<sub>c</sub> = 2.5, shift yes. Predict E at r = 1.0 before running.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Unshifted E(1) = 4(1 − 1) = 0. The offset is 4(2.5⁻¹² − 2.5⁻⁶) = −0.016316891136, so E = 0 − (−0.0163…) = <b>+0.016316891136</b>.</div></div>
    <div class="ex"><p><b>2.</b> For morse, which number in the table tells you r₀ without looking at the input?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The force changes sign at r = r₀: there e^(−α(r−r₀)) = 1 and F = 2αD₀(1 − 1) = 0. Put r₀ inside your range and look for where F goes from + (repulsive) to − (attractive).</div></div>
    <div class="ex"><p><b>3.</b> A new style forgets to write <code>init_one()</code>. What goes wrong?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The base version returns a cutoff of 0, so the neighbour list holds no pairs and the style never sees any neighbours. Writing <code>settings</code>, <code>coeff</code> and <code>compute</code> is forced by "= 0"; init_one is not forced but is needed in practice.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let t = 0;
    const showT = () => { const T = TRACE[t]; $('.st3-n').textContent = `${t + 1} / ${TRACE.length}`; $('.st3-when').innerHTML = `<b>${esc(T.when)}</b> → ${esc(T.fn)} · ${gh(...T.file)}`; $('.st3-code').textContent = T.code; };
    root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => { t = (t + +b.dataset.t + TRACE.length) % TRACE.length; showT(); }); showT();

    const argsUI = () => { const P = PLUGS[$('.s3-style').value]; $('.s3-args').innerHTML = P.args.map(([n, v], i) => `<label>${n} <input type="number" class="lmp-in sm" data-a="${i}" value="${v}" step="0.1"></label>`).join(''); root.querySelectorAll('[data-a]').forEach(i => { i.oninput = draw; }); };
    let S = null, H = null;
    const num = (sel, d) => { const v = parseFloat(root.querySelector(sel).value); return isFinite(v) ? v : d; };
    function draw() {
      const style = $('.s3-style').value, P = PLUGS[style];
      const a = P.args.map((_, i) => num(`[data-a="${i}"]`, P.args[i][1]));
      const n = Math.max(2, Math.min(10, Math.round(num('[data-q=n]', 6))));
      S = { style, a, rc: num('[data-q=rc]', 2.5), n, lo: num('[data-q=lo]', 0.9), hi: num('[data-q=hi]', 2.6), shift: $('.s3-shift').checked };
      $('.s3-formula').innerHTML = esc(P.formula + (S.shift ? '   − offset (E at the cutoff)' : '')) + ` · single() in ${gh(...P.file)}`;
      $('.s3-script').textContent = script(S);
      const C = P.init(a, S.rc, S.shift, $('.s3-code').checked), cut = P.cut(a, S.rc);
      H = []; for (let i = 0; i < n; i++) { const r = S.lo + (S.hi - S.lo) * i / (n - 1), rsq = r * r; if (rsq < cut * cut) { const o = P.single(C, rsq); H.push([g15(r), g15(o.e), g15(o.f * r)]); } else H.push([g15(r), '0', '0']); }
      H.off = C.off;
      $('.s3-cmp').innerHTML = `<table>${head(['r', 'E by hand', 'F by hand'])}${H.map(h => row(h.map(x => `<span class="mono">${x}</span>`))).join('')}</table>${S.shift ? `<p class="muted">offset from init_one = ${g15(C.off)}</p>` : ''}`;
      $('.s3-st').textContent = '';
    }
    const codeVis = () => { $('.s3-codel').style.display = $('.s3-style').value === 'lj/expand' && $('.s3-shift').checked ? 'block' : 'none'; };
    $('.s3-style').onchange = () => { argsUI(); codeVis(); draw(); }; $('.s3-shift').onchange = () => { codeVis(); draw(); }; $('.s3-code').onchange = draw;
    root.querySelectorAll('[data-q]').forEach(i => { i.oninput = draw; });
    argsUI(); draw();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.s3-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(script(S), { read: ['pw.table'] }); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      const txt = res.files && res.files['pw.table']; if (!txt) { st.textContent = 'No table was written.'; return; }
      st.textContent = 'pw.table read back:';
      const L = txt.split('\n').map(l => l.trim().split(/\s+/)).filter(w => w.length === 4 && /^\d+$/.test(w[0])).map(w => w.slice(1));
      const cell = (h, l) => { const ok = h === l; return `<span class="mono ${ok ? 'm-ok' : 'm-bad'}">${l ?? '—'}</span>`; };
      const all = H.every((h, i) => L[i] && h[1] === L[i][1] && h[2] === L[i][2]);
      $('.s3-cmp').innerHTML = `<table>${head(['r', 'E by hand', 'E LAMMPS', 'F by hand', 'F LAMMPS'])}${H.map((h, i) => row([`<span class="mono">${h[0]}</span>`, `<span class="mono">${h[1]}</span>`, cell(h[1], L[i] && L[i][1]), `<span class="mono">${h[2]}</span>`, cell(h[2], L[i] && L[i][2])])).join('')}</table>
        <p>${all ? '<span class="m-ok">✓</span> Every number agrees to all 15 printed digits.' : '<span class="m-bad">✗</span> Some digits differ: see the red cells.'}${S.shift ? ` Offset from init_one = ${g15(H.off)}.` : ''}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
