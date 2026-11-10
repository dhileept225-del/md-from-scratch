// Research · Part 0.5: reading C++ with school maths, by stepping through the real LAMMPS force loop.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/src/';

// the real lines of PairLJCut::compute that matter for atoms on a line (y = z = 0)
const CODE = [
  [93, 'for (ii = 0; ii < inum; ii++) {'],
  [94, '  i = ilist[ii];'],
  [95, '  xtmp = x[i][0];'],
  [102, '  for (jj = 0; jj < jnum; jj++) {'],
  [103, '    j = jlist[jj];'],
  [107, '    delx = xtmp - x[j][0];'],
  [110, '    rsq = delx * delx + dely * dely + delz * delz;'],
  [113, '    if (rsq < cutsq[itype][jtype]) {'],
  [114, '      r2inv = 1.0 / rsq;'],
  [115, '      r6inv = r2inv * r2inv * r2inv;'],
  [116, '      forcelj = r6inv * (lj1[itype][jtype] * r6inv - lj2[itype][jtype]);'],
  [117, '      fpair = factor_lj * forcelj * r2inv;'],
  [119, '      f[i][0] += delx * fpair;'],
  [123, '      f[j][0] -= delx * fpair;'],
  [129, '      evdwl = r6inv * (lj3[itype][jtype] * r6inv - lj4[itype][jtype]) - offset[itype][jtype];'],
  [133, '      if (evflag) ev_tally(i, j, nlocal, newton_pair, evdwl, …);'],
];
const SAYS = {
  93: 'loop over my atoms: ii = 0, 1, 2 (like i in Σ)', 94: 'which atom: the list says atom i', 95: 'remember its x position',
  102: 'loop over the neighbours of atom i', 103: 'which neighbour: atom j', 107: 'the gap along x: x_i − x_j', 110: 'the squared distance r² (dy = dz = 0 here)',
  113: 'piecewise: only if r² is below the cutoff² (2.5² = 6.25)', 114: '1/r²', 115: '1/r⁶', 116: '48r⁻¹² − 24r⁻⁶ = r·F(r)', 117: 'F(r)/r (factor_lj = 1)',
  119: 'add this pair\'s push to atom i', 123: 'Newton\'s third law: the opposite push on atom j', 129: 'the pair energy 4(r⁻¹² − r⁻⁶)', 133: 'add evdwl to the total energy',
};

function trace(xs) {
  const n = xs.length, f = xs.map(() => 0), steps = []; let E = 0, v = {};
  const snap = (line, extra = '') => steps.push({ line, v: { ...v }, f: f.slice(), E, extra });
  for (let ii = 0; ii < n; ii++) {
    v = { ii, i: ii }; snap(93); snap(94); v.xtmp = xs[ii]; snap(95);
    for (let jj = 0, j = ii + 1; j < n; jj++, j++) {
      v = { ii, i: ii, xtmp: xs[ii], jj }; snap(102); v.j = j; snap(103);
      v.delx = xs[ii] - xs[j]; snap(107); v.rsq = v.delx * v.delx; snap(110);
      if (!(v.rsq < 6.25)) { snap(113, 'false: too far apart, skip this pair'); continue; }
      snap(113, 'true');
      v.r2inv = 1 / v.rsq; snap(114); v.r6inv = v.r2inv ** 3; snap(115);
      v.forcelj = v.r6inv * (48 * v.r6inv - 24); snap(116); v.fpair = v.forcelj * v.r2inv; snap(117);
      f[ii] += v.delx * v.fpair; snap(119); f[j] -= v.delx * v.fpair; snap(123);
      v.evdwl = v.r6inv * (4 * v.r6inv - 4); snap(129); E += v.evdwl; snap(133);
    }
    if (ii === n - 1) snap(93, 'atom 3 has no neighbours left in the half list (each pair is stored once)');
  }
  steps.push({ line: 0, v: {}, f: f.slice(), E, extra: 'loop finished' });
  return steps;
}

const THREE = xs => `units           lj
atom_style      atomic
atom_modify     map yes
region          box block -5 5 -5 5 -5 5
create_box      1 box
${xs.map(x => `create_atoms    1 single ${x} 0.0 0.0`).join('\n')}
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
variable        f1 equal fx[1]
variable        f2 equal fx[2]
variable        f3 equal fx[3]
thermo_style    custom step pe v_f1 v_f2 v_f3
thermo_modify   norm no format float %.8g
run             0
`;

export default {
  id: 'read-cpp',
  title: 'Reading C++ with school maths',
  minutes: 25,
  summary: 'C++ looks foreign, but the parts of LAMMPS that do physics are ordinary maths: named numbers, tables, sums and cases. You learn the handful of ideas and symbols needed, then step line by line through the real LAMMPS force loop for three atoms, watching every number, and match the result to LAMMPS.',
  learn: [
    'The six C++ ideas behind every physics routine, each matched to school maths',
    'The symbols you will meet in LAMMPS and what they mean in words',
    'Step through the real force loop of pair lj/cut for three atoms and match LAMMPS',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: six ideas, all of them maths</h2>
    <div class="tbl wide"><table>${head(['C++', 'In words', 'In maths'])}
${row(['<code>double r2inv;</code>', 'make a box called r2inv that holds a decimal number', 'let r2inv be a real number'])}
${row(['<code>r2inv = 1.0 / rsq;</code>', 'put the value of 1/rsq into the box (a command, not an equation: <code>n = n + 1</code> is fine)', 'r2inv := 1/r²'])}
${row(['<code>x[i][0]</code>', 'row i, column 0 of the table x (counting starts at 0)', 'x-coordinate of atom i'])}
${row(['<code>for (j = 0; j &lt; n; j++) { … }</code>', 'repeat the lines in { } for j = 0, 1, …, n−1', 'Σ over j'])}
${row(['<code>if (rsq &lt; cutsq) { … }</code>', 'do the lines in { } only when the test is true', 'a case of a piecewise function'])}
${row(['<code>double init_one(int i, int j)</code>', 'a named recipe that takes i and j and returns a number', 'f(i, j)'])}</table></div>
    <p>A seventh idea, the <b>class</b>, is a form with blanks. LAMMPS has a form called <code>Pair</code> with a blank named <code>compute</code>; every force style fills in that blank its own way (section 4).</p>
  </section>

  <section class="block prose">
    <h2>2 · The symbols you will meet</h2>
    <div class="tbl wide"><table>${head(['Symbol', 'Read it as', 'Example from LAMMPS'])}
${row(['<code>;</code>', 'end of one instruction', '<code>r6inv = r2inv * r2inv * r2inv;</code>'])}
${row(['<code>{ }</code>', 'a group of lines that belong together', 'the body of a loop'])}
${row(['<code>//</code>', 'a comment for humans; the compiler ignores the rest of the line', '<code>// loop over neighbors of my atoms</code>'])}
${row(['<code>*</code>', 'multiply (in declarations such as <code>double **x</code> it marks a table: ignore it)', '<code>delx * fpair</code>'])}
${row(['<code>+=</code>, <code>-=</code>', 'add to, subtract from', '<code>f[i][0] += delx * fpair;</code>'])}
${row(['<code>++</code>', 'add one', '<code>ii++</code>'])}
${row(['<code>[ ]</code>', 'pick an entry from a table', '<code>lj1[itype][jtype]</code>'])}
${row(['<code>-&gt;</code>', 'the part called … inside …', '<code>atom-&gt;x</code>: the positions kept by the atom manager'])}
${row(['<code>::</code>', 'belongs to', '<code>PairLJCut::compute</code>: the compute of PairLJCut'])}
${row(['<code>&amp;</code>', 'here: keep only some bits of a number (a technical trick; skip it)', '<code>j &amp;= NEIGHMASK;</code>'])}
${row(['<code>||</code>, <code>&amp;&amp;</code>', 'or, and', '<code>if (newton_pair || j &lt; nlocal)</code>'])}
${row(['<code>#</code>', 'an instruction to the compiler, not to the computer', '<code>#include "atom.h"</code>'])}</table></div>
  </section>

  <section class="block">
    <h2>3 · Modify, observe, match: step through the real loop</h2>
    <p class="muted">These are the lines of <a href="${GH}pair_lj_cut.cpp#L93-L133" target="_blank" rel="noopener">pair_lj_cut.cpp, lines 93–133</a> that do the work, for three atoms on a line (ε = σ = 1, cutoff 2.5). Choose the positions, then press Step and watch every box change.</p>
    <div class="stepper" id="stepper">
      <div class="row pos-row">${[0, 1.1, 2.3].map((x, i) => `<label class="rlab">x<sub>${i + 1}</sub> <input type="number" class="r-in xin" step="0.05" value="${x}"></label>`).join('')}</div>
      <div class="row"><button class="btn sm primary" data-k="step">Step ▸</button><button class="btn sm ghost" data-k="pair">Next pair ▸▸</button><button class="btn sm ghost" data-k="end">To the end</button><button class="btn sm ghost" data-k="reset">Reset</button></div>
      <div class="step-grid">
        <pre class="code-lines"></pre>
        <div><div class="says"></div><div class="vars tbl"></div></div>
      </div>
      <div class="row"><button class="btn sm primary" data-k="lmp">▶ Run LAMMPS on the same three atoms</button><span class="lmp-status muted"></span></div>
      <div class="lmp-compare"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>What the loop is, in one line of maths</h2>
    <p>Line 119 adds each pair's push to atom i, and line 123 the opposite push to atom j. Over the whole loop that is</p>
    <div class="calc">F_i = Σ_j  F(r_ij) × (x_i − x_j) / r_ij        with   F(r) = 24(2r⁻¹³ − r⁻⁷)
E   = Σ_pairs  4(r⁻¹² − r⁻⁶)</div>
    <p>For x = 0, 1.1, 2.3 the three pairs are at r = 1.1, 2.3 and 1.2:</p>
    <div class="tbl"><table>${head(['pair', 'delx', 'fpair = F/r', 'evdwl'])}
${row(['1–2', '−1.1', '1.4437231', '−0.98337245'])}
${row(['1–3', '−2.3', '−0.030232992', '−0.026837948'])}
${row(['2–3', '−1.2', '−1.8430778', '−0.89096529'])}</table></div>
    <div class="calc">f₁ = (−1.1)(1.4437231) + (−2.3)(−0.030232992)        = <b>−1.5185595</b>
f₂ = −(−1.1)(1.4437231) + (−1.2)(−1.8430778)       = <b>+3.7997887</b>
f₃ = −(−2.3)(−0.030232992) − (−1.2)(−1.8430778)    = <b>−2.2812292</b>
E  = −0.98337245 − 0.026837948 − 0.89096529      = <b>−1.9011757</b>
check: f₁ + f₂ + f₃ = 0  (Newton's third law)</div>
  </section>

  <section class="block prose">
    <h2>4 · The class: a form with a blank</h2>
    <p>In <a href="${GH}pair.h#L158" target="_blank" rel="noopener">pair.h, line 158</a>, the general form for every force style says:</p>
    <div class="calc">virtual void compute(int, int) = 0;</div>
    <p>Read it as: "every pair style must have a recipe called compute; this form leaves it blank" (<code>= 0</code> means blank). In <a href="${GH}pair_lj_cut.h" target="_blank" rel="noopener">pair_lj_cut.h</a> the Lennard-Jones style fills it in:</p>
    <div class="calc">class PairLJCut : public Pair {          <b>PairLJCut is a filled-in Pair form</b>
 public:
  void compute(int, int) override;       <b>"I fill in the compute blank"</b>
  double init_one(int, int) override;    <b>"and the init_one blank"</b>
  …</div>
    <p>When your input says <code>pair_style lj/cut</code>, LAMMPS picks this filled-in form, and every timestep it calls its compute: the loop you just stepped through. Part B explains how the name lj/cut is connected to PairLJCut.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Translate into maths: <code>forcelj = r6inv * (lj1 * r6inv - lj2);</code> with lj1 = 48 and lj2 = 24.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>forcelj = r⁻⁶(48 r⁻⁶ − 24) = 48r⁻¹² − 24r⁻⁶ = r × F(r), where F(r) = 48r⁻¹³ − 24r⁻⁷ = −dE/dr.</div></div>
    <div class="ex"><p><b>2.</b> Move atom 3 to x = 3.0 in the stepper. Which pairs are skipped at line 113, and what is f₃?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Pair 1–3 (r = 3.0, r² = 9 &gt; 6.25) is skipped. Pair 2–3 has r = 1.9, so atom 3 still feels atom 2: fpair = F(1.9)/1.9 = −0.13530564, delx = 1.1 − 3.0 = −1.9, so f₃ = −(−1.9)(−0.13530564) = <b>−0.25708072</b> (pulled towards atom 2). Run LAMMPS in the stepper to check the digits.</div></div>
    <div class="ex"><p><b>3.</b> Why does LAMMPS write <code>f[j][0] -= delx * fpair</code> instead of looping over the pair a second time from atom j's side?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Newton's third law: the force of i on j is exactly minus the force of j on i. Computing it once and using it twice halves the work. That is why the neighbour list stores each pair once (a "half" list).</div></div>
  </section>`,

  mount(root) {
    const S = root.querySelector('#stepper'), codeEl = S.querySelector('.code-lines'), says = S.querySelector('.says'), vars = S.querySelector('.vars');
    let steps = [], k = 0;
    const xsNow = () => [...S.querySelectorAll('.xin')].map(i => parseFloat(i.value)).map(v => isFinite(v) ? v : 0);
    const fmt = v => v === undefined ? '' : typeof v === 'number' ? (Number.isInteger(v) && Math.abs(v) < 100 ? String(v) : g8(v)) : v;
    function draw() {
      const s = steps[k];
      codeEl.innerHTML = CODE.map(([n, c]) => `<span class="${n === s.line ? 'cur' : ''}"><i>${n}</i>${c.replace(/</g, '&lt;')}</span>`).join('');
      says.innerHTML = s.line ? `<b>line ${s.line}:</b> ${SAYS[s.line]}${s.extra ? ` <span class="m-ok">→ ${s.extra}</span>` : ''}` : `<b>${s.extra}.</b> Step ${k} of ${steps.length - 1}.`;
      const v = s.v, names = ['ii', 'i', 'xtmp', 'jj', 'j', 'delx', 'rsq', 'r2inv', 'r6inv', 'forcelj', 'fpair', 'evdwl'];
      vars.innerHTML = `<table>${names.map(n => `<tr><td class="mono">${n}</td><td class="mono">${fmt(v[n])}</td></tr>`).join('')}
        ${s.f.map((x, i) => `<tr class="fx"><td class="mono">f[${i}][0] (atom ${i + 1})</td><td class="mono">${g8(x)}</td></tr>`).join('')}
        <tr class="fx"><td class="mono">total energy</td><td class="mono">${g8(s.E)}</td></tr></table>`;
      const cur = codeEl.querySelector('.cur'); if (cur) codeEl.scrollTop = Math.max(0, cur.offsetTop - codeEl.clientHeight / 2);
    }
    const reset = () => { steps = trace(xsNow()); k = 0; draw(); S.querySelector('.lmp-compare').innerHTML = ''; };
    S.querySelectorAll('.xin').forEach(i => i.oninput = reset);
    S.querySelector('[data-k=step]').onclick = () => { if (k < steps.length - 1) k++; draw(); };
    S.querySelector('[data-k=pair]').onclick = () => { do { k++; } while (k < steps.length - 1 && steps[k].line !== 133 && !(steps[k].line === 113 && steps[k].extra.startsWith('false'))); k = Math.min(k, steps.length - 1); draw(); };
    S.querySelector('[data-k=end]').onclick = () => { k = steps.length - 1; draw(); };
    S.querySelector('[data-k=reset]').onclick = reset;
    S.querySelector('[data-k=lmp]').onclick = async e => {
      const b = e.currentTarget, st = S.querySelector('.lmp-status'), xs = xsNow(); b.disabled = true;
      st.textContent = 'Loading LAMMPS and running… (the first run downloads about 3 MB)';
      const res = await runLammps(THREE(xs)); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${res.error}</span>`; return; }
      st.textContent = `Done · LAMMPS ${(res.version.match(/\(([^)]+)\)/) || [])[1] || ''}`;
      const r0 = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols, last = trace(xs).at(-1);
      const line = (label, hv, lv) => { const m = matchMarkup(g8(hv), lv); return `<tr><td>${label}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; };
      S.querySelector('.lmp-compare').innerHTML = `<div class="tbl"><table><tr><th></th><th>the loop above</th><th>LAMMPS</th><th></th></tr>
        ${[1, 2, 3].map(i => line(`fx atom ${i}`, last.f[i - 1], r0[c.indexOf('v_f' + i)])).join('')}${line('energy', last.E, r0[c.indexOf('PotEng')])}</table></div>`;
    };
    reset();
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
