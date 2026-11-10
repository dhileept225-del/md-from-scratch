// Research · A.2: how LAMMPS reads one line (Input::file and Input::parse), with a line reader checked against LAMMPS.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/src/input.cpp';

// %.Ng the way C prints it
function gfmt(x, p) {
  if (!isFinite(x)) return String(x);
  if (x === 0) return '0';
  let s = x.toPrecision(p);
  const ex = Math.floor(Math.log10(Math.abs(Number(s))));
  if (ex < -4 || ex >= p) {
    let [m, e] = x.toExponential(p - 1).split('e');
    if (m.includes('.')) m = m.replace(/0+$/, '').replace(/\.$/, '');
    const n = Math.abs(+e); return `${m}e${+e < 0 ? '-' : '+'}${n < 10 ? '0' + n : n}`;
  }
  s = x.toFixed(Math.max(0, p - 1 - ex));
  return s.includes('.') ? s.replace(/0+$/, '').replace(/\.$/, '') : s;
}
function cfmt(x, f) {
  const m = f.match(/^%(\d*)\.(\d+)([fge])$/); if (!m) return gfmt(x, 20);
  const n = +m[2], s = m[3] === 'f' ? x.toFixed(n) : m[3] === 'e' ? x.toExponential(n).replace(/e([+-])(\d)$/, 'e$10$2') : gfmt(x, n || 1);
  return s.padStart(+m[1] || 0);
}

// a small evaluator for LAMMPS equal-style formulas: numbers, + - * / ^, ( ), v_name, PI and a few functions
function evaluate(expr, vars) {
  let i = 0; const s = expr.replace(/\s+/g, '');
  const fail = msg => { throw new Error(msg); };
  const num = () => { const m = s.slice(i).match(/^(\d+\.?\d*|\.\d+)([eE][+-]?\d+)?/); if (!m) fail(`cannot read the formula at "${s.slice(i)}"`); i += m[0].length; return parseFloat(m[0]); };
  const FN = { sqrt: Math.sqrt, exp: Math.exp, ln: Math.log, log: Math.log10, abs: Math.abs, sin: Math.sin, cos: Math.cos, tan: Math.tan };
  function atom() {
    if (s[i] === '(') { i++; const v = sum(); if (s[i] !== ')') fail('missing )'); i++; return v; }
    if (s[i] === '-') { i++; return -power(); }
    if (s[i] === '+') { i++; return power(); }
    const w = s.slice(i).match(/^[A-Za-z_]\w*/);
    if (w) {
      const name = w[0]; i += name.length;
      if (name.startsWith('v_')) return varValue(name.slice(2), vars);
      if (name === 'PI') return Math.PI;
      if (FN[name] && s[i] === '(') { i++; const v = sum(); if (s[i] !== ')') fail('missing )'); i++; return FN[name](v); }
      fail(`unknown name ${name} (this page knows v_name, PI, sqrt, exp, ln, log, abs, sin, cos, tan)`);
    }
    return num();
  }
  function power() { const b = atom(); if (s[i] === '^') { i++; return b ** power(); } return b; }
  function prod() { let v = power(); while (s[i] === '*' || s[i] === '/') { const o = s[i++]; const r = power(); v = o === '*' ? v * r : v / r; } return v; }
  function sum() { let v = prod(); while (s[i] === '+' || s[i] === '-') { const o = s[i++]; const r = prod(); v = o === '+' ? v + r : v - r; } return v; }
  const v = sum(); if (i < s.length) fail(`cannot read the formula at "${s.slice(i)}"`); return v;
}
// an equal-style variable is turned into text with 15 digits, and read back from that text
function varText(name, vars) {
  const d = vars[name]; if (!d) throw new Error(`variable ${name} is not defined`);
  return d.style === 'equal' ? gfmt(evaluate(d.value, vars), 15) : d.value;
}
function varValue(name, vars) { const t = parseFloat(varText(name, vars)); if (!isFinite(t)) throw new Error(`variable ${name} is not a number`); return t; }

function readDefs(text) {
  const vars = {};
  for (const line of text.split('\n')) {
    const m = line.trim().match(/^variable\s+(\w+)\s+(equal|string|index)\s+(.+)$/); if (!m) continue;
    let v = m[3].trim(); if (m[2] !== 'equal') v = v.split(/\s+/)[0].replace(/^["']|["']$/g, '');
    vars[m[1]] = { style: m[2], value: v.replace(/^"|"$/g, '') };
  }
  return vars;
}

// the six stages of Input::file / Input::parse
function stages(raw, vars) {
  const S = { raw };
  S.joined = raw.split('\n').reduce((acc, l) => { const t = acc.pending + l; if (/&\s*$/.test(t)) return { pending: t.replace(/&\s*$/, ' '), out: acc.out }; return { pending: '', out: [...acc.out, t] }; }, { pending: '', out: [] }).out[0] ?? raw;
  let q = null, cut = S.joined.length;
  for (let k = 0; k < S.joined.length; k++) { const c = S.joined[k]; if (q) { if (c === q) q = null; continue; } if (c === '"' || c === "'") q = c; else if (c === '#') { cut = k; break; } }
  S.stripped = S.joined.slice(0, cut);
  S.subErr = null;
  try {
    const subOne = txt => txt.replace(/\$\(([^()]*(?:\([^()]*\)[^()]*)*)\)|\$\{(\w+)\}|\$([A-Za-z0-9])/g, (all, ex, nm, ch) => {
      if (ex !== undefined) { const k = ex.lastIndexOf(':%'); const e = k >= 0 ? ex.slice(0, k) : ex; const v = evaluate(e, vars); return k >= 0 ? cfmt(v, ex.slice(k + 1)) : gfmt(v, 20); }
      return varText(nm || ch, vars);
    });
    // LAMMPS skips text inside quotes at this stage (input.cpp, parse)
    let outS = '', buf = '', qq0 = null;
    for (const c of S.stripped) {
      if (qq0) { buf += c; if (c === qq0) { outS += buf; buf = ''; qq0 = null; } continue; }
      if (c === '"' || c === "'") { outS += subOne(buf); buf = c; qq0 = c; continue; }
      buf += c;
    }
    S.subst = outS + (qq0 ? buf : subOne(buf));
  } catch (e) { S.subErr = e.message; S.subst = S.stripped; }
  const words = []; let cur = '', qq = null, has = false;
  for (const c of S.subst) { if (qq) { if (c === qq) qq = null; else cur += c; continue; } if (c === '"' || c === "'") { qq = c; has = true; continue; } if (/\s/.test(c)) { if (cur || has) { words.push(cur); cur = ''; has = false; } continue; } cur += c; }
  if (cur || has) words.push(cur);
  S.words = words;
  return S;
}

const PRESETS = [
  { name: 'variables in a pair_coeff', defs: 'variable eps string 1.0\nvariable rc equal 2.0', line: 'pair_coeff 1 1 ${eps} 1.0 $(v_rc+0.5)   # LJ for type 1' },
  { name: 'quotes protect # and $', defs: 'variable rc equal 2.0', line: 'print "cutoff ${rc} # inside quotes" # a real comment' },
  { name: 'a line continued with &', defs: 'variable n index 10', line: 'region box block 0 ${n} &\n                 0 ${n} 0 ${n}' },
  { name: 'the 15-digit trap', defs: 'variable third equal 1/3', line: 'variable show string $(1/3)/${third}/$(v_third)' },
  { name: 'choosing the format', defs: 'variable rc equal 2.0', line: 'variable show string $(v_rc/3:%.4f)_$(v_rc/3:%.3e)' },
];

const MELT_VARS = IN_MELT
  .replace('units           lj', 'variable        T equal 2*1.5\nvariable        rc equal 2.0\n\nunits           lj')
  .replace('create 3.0 87287', 'create ${T} 87287')
  .replace('pair_style      lj/cut 2.5', 'pair_style      lj/cut $(v_rc+0.5)')
  .replace('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 1 1.0 1.0 $(v_rc+0.5)   # same cutoff')
  .replace('run             250', 'run             0');

export default {
  id: 'reading-a-line',
  title: 'How LAMMPS reads one line',
  minutes: 20,
  summary: 'Before LAMMPS does anything with a line, it joins continued lines, cuts off comments, replaces variables with their values and splits the result into words. You follow these six steps in the real source, run any line through a line reader that predicts each step, and check the prediction against LAMMPS.',
  learn: [
    'The six steps between a line of text and a command being run, in input.cpp',
    'The rules for &, #, $x, ${name} and $(formula), including a 15-digit trap',
    'Use variables in the melt script and confirm by hand that nothing changes',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: six steps for every line</h2>
    <div class="tbl wide"><table>${head(['Step', 'What happens', 'Where in input.cpp'])}
${row(['1. read', 'take the next line of the file', '<code>Input::file()</code>, line 196'])}
${row(['2. join', 'a line ending in <code>&amp;</code> is joined to the next (the &amp; becomes a space)', 'line 205: "if line ends in continuation char \'&amp;\'"'])}
${row(['3. strip', 'everything from <code>#</code> on is a comment, unless the # is inside quotes', `<a href="${GH}#L444" target="_blank" rel="noopener">lines 444–472</a>`])}
${row(['4. substitute', '<code>$x</code>, <code>\${name}</code> and <code>$(formula)</code> are replaced by values, except inside quotes (a command such as print substitutes its quoted text itself, later)', `<a href="${GH}#L488" target="_blank" rel="noopener">line 488</a>: <code>substitute(copy,…)</code>`])}
${row(['5. split', 'cut into words at spaces; text in quotes stays one word (quotes removed)', '<code>nextword()</code>, line 525'])}
${row(['6. run', 'the first word is the command; LAMMPS looks it up and runs it with the other words', '<code>execute_command()</code>, line 764'])}</table></div>
    <p>The comment step in the real code, translated:</p>
    <div class="calc">  // strip a # comment by replacing it with 0
  // do not treat a # inside single/double/triple quotes as a comment
  while (*ptr) {                       <b>walk along the line, one character at a time</b>
    if (*ptr == '#') {                 <b>found a # outside quotes:</b>
      *ptr = '\\0';                     <b>end the line here</b>
      break;
    }
    if (*ptr == '\\'') { … }             <b>a quote: jump to the matching quote</b></div>
    <p>Step 6 is a long list of names (<code>if (mycmd == "clear") … else if (mycmd == "units") …</code>), followed by a lookup table for commands such as <code>run</code> and <code>velocity</code>. A name in neither gives the error "Unknown command".</p>
    <div class="tbl wide"><table>${head(['Write', 'Means', 'Example (rc = 2.0)'])}
${row(['<code>$x</code>', 'the variable with the one-letter name x', '<code>$n</code> → 10'])}
${row(['<code>\${name}</code>', 'the variable called name; an equal-style formula is printed with 15 digits', '<code>\${rc}</code> → 2'])}
${row(['<code>$(formula)</code>', 'work out the formula now, print with up to 20 digits', '<code>$(v_rc+0.5)</code> → 2.5'])}
${row(['<code>$(formula:%.4f)</code>', 'the same, with your own format', '<code>$(v_rc/3:%.4f)</code> → 0.6667'])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Modify and observe: the line reader</h2>
    <p class="muted">Pick an example or type your own line and variables. The page predicts every step; "Check with LAMMPS" asks the real LAMMPS to print the substituted line.</p>
    <div class="lmp" id="reader">
      <div class="row presets">${PRESETS.map((p, i) => `<button class="btn sm ghost" data-pre="${i}">${p.name}</button>`).join('')}</div>
      <label class="lab">Variables defined earlier in the script</label><textarea class="lmp-in small-ta defs" spellcheck="false" wrap="off"></textarea>
      <label class="lab">The line (or lines) LAMMPS reads</label><textarea class="lmp-in small-ta line" spellcheck="false" wrap="off"></textarea>
      <ol class="stages"></ol>
      <div class="row"><button class="btn sm primary" data-k="check">Check with LAMMPS</button><span class="lmp-status muted"></span></div>
      <div class="lmp-compare"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>3 · Match: variables change nothing when the values are the same</h2>
    <p>Here is melt with the temperature and the cutoff written as variables: <code>T = 2*1.5</code> and the cutoff <code>$(v_rc+0.5)</code> with rc = 2.0. Step 4 turns them into 3 and 2.5, so the hand prediction for step 0 is exactly the original line.</p>
    <div class="calc">variable        T equal 2*1.5
variable        rc equal 2.0
…
velocity        all create \${T} 87287 loop geom        <b>→ create 3 87287</b>
pair_style      lj/cut $(v_rc+0.5)                      <b>→ lj/cut 2.5</b>
pair_coeff      1 1 1.0 1.0 $(v_rc+0.5)   # same cutoff <b>→ 1 1 1.0 1.0 2.5</b>

predicted step 0:  E_pair −6.7733681   TotEng −2.2744931   Press −3.7033504</div>
    <div class="row"><button class="btn sm primary" data-k="melt">▶ Run this melt</button><span class="lmp-status muted melt-st"></span></div>
    <div class="lmp-compare melt-cmp"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> What does LAMMPS make of <code>pair_coeff 1 1 1.0 1.0 2.5#cutoff</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The # starts a comment even with no space before it, so the line becomes <code>pair_coeff 1 1 1.0 1.0 2.5</code>. Try it in the line reader.</div></div>
    <div class="ex"><p><b>2.</b> With <code>variable third equal 1/3</code>, why does <code>$(v_third)</code> print 0.33333333333333298176 but <code>$(1/3)</code> prints 0.33333333333333331483?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>An equal-style variable is handed over as text with 15 digits, 0.333333333333333, and read back as a number; the digits after the 15th are lost. $(1/3) computes 1/3 directly. For physics this is harmless (the difference is 3 × 10⁻¹⁶), but it explains odd digits in printed output.</div></div>
    <div class="ex"><p><b>3.</b> You write <code>variable rc equal 2.5</code> and then <code>pair_style lj/cut $rc</code>. What goes wrong?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>$rc</code> means the one-letter variable <code>r</code> followed by the letter c. If r is not defined LAMMPS stops with an error; if it is, you get a wrong number. Names longer than one letter need braces: <code>\${rc}</code>.</div></div>
  </section>`,

  mount(root) {
    const R = root.querySelector('#reader'), defs = R.querySelector('.defs'), line = R.querySelector('.line'), list = R.querySelector('.stages');
    let S = null;
    function draw() {
      const vars = readDefs(defs.value); S = stages(line.value, vars);
      const vis = t => esc(t).replace(/ /g, '<span class="sp">·</span>');
      list.innerHTML = `
        <li><b>read</b><div class="mono stage">${esc(S.raw).replace(/\n/g, '<span class="nl">↵</span>\n')}</div></li>
        <li><b>join &amp;</b><div class="mono stage">${vis(S.joined)}</div></li>
        <li><b>strip #</b><div class="mono stage">${vis(S.stripped)}${S.stripped.length < S.joined.length ? `<span class="m-bad">${vis(S.joined.slice(S.stripped.length))}</span>` : ''}</div></li>
        <li><b>substitute $</b><div class="mono stage">${S.subErr ? `<span class="lmp-err">${esc(S.subErr)}</span>` : vis(S.subst)}</div></li>
        <li><b>split into words</b><div class="stage words">${S.words.map((w, i) => `<span class="word${i ? '' : ' cmd'}">${esc(w) || '&nbsp;'}</span>`).join('')}</div></li>
        <li><b>run</b><div class="stage">${S.words.length ? `command <code>${esc(S.words[0])}</code> with ${S.words.length - 1} argument${S.words.length === 2 ? '' : 's'}` : 'an empty line: nothing to do'}</div></li>`;
      R.querySelector('.lmp-compare').innerHTML = '';
    }
    const load = p => { defs.value = p.defs; line.value = p.line; draw(); };
    R.querySelectorAll('[data-pre]').forEach(b => b.onclick = () => load(PRESETS[+b.dataset.pre]));
    defs.oninput = draw; line.oninput = draw;
    load(PRESETS[0]);
    R.querySelector('[data-k=check]').onclick = async e => {
      const b = e.currentTarget, st = R.querySelector('.lmp-status');
      b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(`${defs.value}\nxyzzy ${line.value}\n`); b.disabled = false; st.textContent = '';
      const pl = res.lines.find(l => l.startsWith('--> parsed line:'));
      if (!pl) { R.querySelector('.lmp-compare').innerHTML = `<p class="lmp-err">${esc(res.error || 'LAMMPS did not report a parsed line.')}</p>`; return; }
      const lw = []; let cur = '', q = null, has = false;
      for (const c of pl.slice(16)) { if (q) { if (c === q) q = null; else cur += c; continue; } if (c === '"' || c === "'") { q = c; has = true; continue; } if (/\s/.test(c)) { if (cur || has) { lw.push(cur); cur = ''; has = false; } continue; } cur += c; }
      if (cur || has) lw.push(cur);
      const got = lw.slice(1), ok = got.length === S.words.length && got.every((w, i) => w === S.words[i]);
      const chips = ws => ws.map((w, i) => `<span class="word${i ? '' : ' cmd'}">${esc(w) || '&nbsp;'}</span>`).join('');
      R.querySelector('.lmp-compare').innerHTML = `<div class="tbl wide"><table><tr><th>words predicted by this page</th><td><div class="words">${chips(S.words)}</div></td></tr><tr><th>words LAMMPS parsed</th><td><div class="words">${chips(got)}</div></td></tr></table></div>
        <p>${ok ? '<span class="m-ok">✓ identical, word by word.</span>' : '<span class="m-bad">They differ: this page\'s small formula reader may not know everything LAMMPS knows.</span>'}</p>
        <p class="muted small">How the check works: LAMMPS runs your variable lines, then your line with an unknown command <code>xyzzy</code> put in front. LAMMPS stops, and its error message shows the line exactly as it parsed it: <code>${esc(pl)}</code></p>`;
    };
    root.querySelector('[data-k=melt]').onclick = async e => {
      const b = e.currentTarget, st = root.querySelector('.melt-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(MELT_VARS); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const r = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols;
      const H = { E_pair: -6.7733681, TotEng: -2.2744931, Press: -3.7033504 };
      root.querySelector('.melt-cmp').innerHTML = `<div class="tbl"><table><tr><th>step 0</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${Object.entries(H).map(([k, h]) => { const lv = r[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
