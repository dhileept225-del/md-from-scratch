// Research · B.1: a tour of the source tree. Real folder sizes, the naming rule from a script line to a .cpp file, and a style swap checked by hand.
import { runLammps, parseThermo, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const TAG = 'stable_22Jul2025_update6', GH = `https://github.com/lammps/lammps/blob/${TAG}/`;
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;

const FOLDERS = {
  src: 'the C++ code of LAMMPS. The only folder that is compiled into <code>lmp</code>',
  lib: 'extra libraries some packages need (GPU, KOKKOS, ATC …), compiled separately',
  cmake: 'the build recipe used by <code>cmake</code> (page 0.4)',
  doc: 'the manual, written as text files that become docs.lammps.org',
  examples: 'input scripts for almost every feature, with sample logs. <code>examples/melt</code> is our in.melt',
  potentials: 'parameter files for many-body potentials (EAM, Tersoff, ReaxFF …)',
  tools: 'helper programs: converters, data file builders, the VS Code and Vim syntax files',
  python: 'the Python module that drives LAMMPS from a Python script',
  fortran: 'the same for Fortran',
  unittest: 'automatic tests: thousands of small checks run before every release',
  bench: 'benchmark inputs for measuring speed',
  third_party: 'notes on code from other projects',
};
const NOT_PKG = new Set(['MAKE', 'DEPEND', 'STUBS', 'fmt', 'nlohmann']);
const PREFIX = { pair_style: 'pair_', bond_style: 'bond_', angle_style: 'angle_', dihedral_style: 'dihedral_', improper_style: 'improper_', kspace_style: '', atom_style: 'atom_vec_', fix: 'fix_', compute: 'compute_', dump: 'dump_', region: 'region_', min_style: 'min_', run_style: '', command: '' };
const SUFFIX = { omp: ['_omp', 'OPENMP'], kk: ['_kokkos', 'KOKKOS'], 'kk/host': ['_kokkos', 'KOKKOS'], 'kk/device': ['_kokkos', 'KOKKOS'], gpu: ['_gpu', 'GPU'], intel: ['_intel', 'INTEL'], opt: ['_opt', 'OPT'] };
const predict = (kind, name) => PREFIX[kind] + name.replace(/[/-]/g, '_') + '.cpp';

let DATA = null;
const load = () => DATA ? Promise.resolve(DATA) : fetch('/data/lammps-src.json').then(r => r.json()).then(d => (DATA = d));

// Which source file does one input line lead to?
function find(line, D) {
  const w = line.replace(/#.*/, '').trim().split(/\s+/); if (!w[0]) return null;
  const at = { fix: 3, compute: 3, dump: 3, region: 2 };
  let kind = null, name = null;
  if (PREFIX[w[0]] !== undefined && w[0] !== 'command') { kind = w[0]; name = w[at[w[0]] || 1]; }
  const bi = D.builtin.find(b => b[0] === w[0]);
  const out = { cmd: w[0], kind, name, bi };
  if (!kind) { const st = D.styles.find(s => s[0] === 'command' && s[1] === w[0]); out.style = st; out.kind = st ? 'command' : null; out.name = st ? w[0] : null; return out; }
  out.style = name && D.styles.find(s => s[0] === kind && s[1] === name);
  out.pred = name && predict(kind, name);
  out.inBrowser = kind === 'command' || (D.browser[kind] || []).includes(name);
  return out;
}
function card(r) {
  if (!r) return '';
  if (!r.bi && !r.style) return `<p><span class="m-bad">✗</span> <code>${esc(r.cmd)}</code> is not a LAMMPS command. LAMMPS would stop with <code>ERROR: Unknown command</code> (page A.9).</p>`;
  const L = [];
  if (r.bi) L.push(`<b>1.</b> <code>${esc(r.cmd)}</code> is read in ${gh('src/input.cpp', r.bi[1])}, which calls ${r.bi[2] ? gh('src/input.cpp', r.bi[2]) : 'another file'}.`);
  if (r.kind === 'command' && r.style) L.push(`<b>${r.bi ? 2 : 1}.</b> <code>${esc(r.cmd)}</code> is a command class of its own: <code>${r.style[2]}</code> in ${gh(r.style[3])}.`);
  else if (r.kind && r.name) {
    if (!r.style) L.push(`<b>2.</b> No ${r.kind.replace('_style', '')} style called <code>${esc(r.name)}</code> exists in ${TAG}. LAMMPS would stop with "Unrecognized … style".`);
    else {
      const f = r.style[3].split('/').pop(), ok = f === r.pred;
      L.push(`<b>2.</b> The style name <code>${esc(r.name)}</code> is looked up in the list of registered styles and found in <code>${r.style[2]}</code>, file ${gh(r.style[3])} ${r.style[4] ? `in package <b>${r.style[4]}</b>` : 'in the <b>core</b>'}.`);
      L.push(`<span class="muted">Naming rule predicted <code>${r.pred}</code> ${ok ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span> (one of the exceptions)'} · in the LAMMPS running in this browser: ${r.inBrowser ? '<span class="m-ok">yes</span>' : '<span class="m-bad">no</span>'}</span>`);
    }
  }
  return L.map(x => `<p>${x}</p>`).join('');
}

// hand calculation for lj/expand on the melt crystal at step 0
function handExpand(D) {
  const rho = 0.8442, T = 3, N = 4000, a = Math.cbrt(4 / rho), rc = 2.5 + D, f = (N - 1) / N;
  const base = [[0, 0, 0], [.5, .5, 0], [.5, 0, .5], [0, .5, .5]], sh = new Map();
  for (let i = -5; i <= 5; i++) for (let j = -5; j <= 5; j++) for (let k = -5; k <= 5; k++) for (const b of base) {
    const r = a * Math.hypot(i + b[0], j + b[1], k + b[2]);
    if (r > 1e-9 && r < rc) { const key = r.toFixed(6); const s = sh.get(key) || { r, n: 0 }; s.n++; sh.set(key, s); }
  }
  let E = 0, W = 0; const rows = [];
  for (const { r, n } of [...sh.values()].sort((x, y) => x.r - y.r)) {
    const s = r - D, s6 = 1 / s ** 6, e = 0.5 * n * 4 * (s6 * s6 - s6), w = 0.5 * n * r * 24 * (2 * s6 * s6 - s6) / s;
    E += e; W += w; rows.push({ r, n, s, e });
  }
  return { rows, E, Tot: E + 1.5 * T * f, P: rho * T * f + rho * W / 3, kin: 1.5 * T * f, idp: rho * T * f };
}
const SCRIPT = D => IN_MELT.replace('pair_style      lj/cut 2.5', 'pair_style      lj/expand 2.5').replace('pair_coeff      1 1 1.0 1.0 2.5', `pair_coeff      1 1 1.0 1.0 ${D}`).replace('run             250', 'run             0');

export default {
  id: 'source-tree',
  title: 'A tour of the source tree',
  minutes: 20,
  summary: 'What is inside the LAMMPS download, folder by folder, with real file counts; which small part is actually compiled; and the naming rule that takes you from any line of an input script to the C++ file that carries it out. Then you swap one pair style for another, find its file, read its maths and check the new energy by hand.',
  learn: [
    'What each top-level folder holds and which one becomes the program',
    'Core code versus packages in src/',
    'The naming rule: from pair_style lj/cut to src/pair_lj_cut.cpp',
    'Swap lj/cut for lj/expand, read its 8 lines of maths, and match LAMMPS by hand',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: what you downloaded</h2>
    <p>The LAMMPS download (${TAG}) has twelve folders. Only one of them, <code>src</code>, becomes the program; the rest are the manual, examples, tests and helpers.</p>
    <div class="tbl wide st-top"><p class="muted">Loading the real counts…</p></div>
  </section>

  <section class="block prose">
    <h2>2 · Inside src: a core and packages</h2>
    <p>Files directly in <code>src/</code> are the <b>core</b>: always compiled. Each upper-case folder inside <code>src/</code> is a <b>package</b>, compiled only if you switch it on (page 0.4, and B.5). A build is the core plus the packages you chose.</p>
    <div class="st-pk"></div>
    <div class="st-pkinfo"></div>
  </section>

  <section class="block">
    <h2>3 · The naming rule: from a script line to a file</h2>
    <div class="prose">
      <p>Every command works the same way. The first word is looked up in <code>src/input.cpp</code>. If the command chooses a <b>style</b> (pair_style, fix, compute …), the style name is then looked up in a list of registered styles, and the matching class does the work. Its file name follows a rule:</p>
      <div class="tbl wide"><table>${head(['Script', 'Rule', 'File'])}
${row(['<code>pair_style lj/cut</code>', '<code>pair_</code> + name, with / → _', '<code>src/pair_lj_cut.cpp</code>'])}
${row(['<code>fix 1 all nve</code>', '<code>fix_</code> + the 4th word', '<code>src/fix_nve.cpp</code>'])}
${row(['<code>compute 1 all rdf 100</code>', '<code>compute_</code> + the 4th word', '<code>src/compute_rdf.cpp</code>'])}
${row(['<code>atom_style full</code>', '<code>atom_vec_</code> + name', '<code>src/MOLECULE/atom_vec_full.cpp</code>'])}
${row(['<code>create_atoms …</code>', 'the command name', '<code>src/create_atoms.cpp</code>'])}
${row(['<code>pair_style lj/cut/omp</code>', 'an accelerator suffix → the same name + <code>_omp</code>', '<code>src/OPENMP/pair_lj_cut_omp.cpp</code>'])}</table></div>
      <p class="st-rule muted"></p>
    </div>
    <label class="lab">Type any line of an input script</label>
    <div class="row"><input class="lmp-in st-q" value="pair_style lj/cut 2.5" spellcheck="false" autocapitalize="off" style="flex:1;min-width:0"><button class="btn sm primary" data-k="find">Find the file</button></div>
    <div class="st-res"></div>
    <div class="row"><button class="btn sm ghost" data-k="melt">Do it for every line of in.melt</button></div>
    <div class="tbl wide st-melt"></div>
  </section>

  <section class="block">
    <h2>4 · Modify: swap lj/cut for lj/expand</h2>
    <div class="prose">
      <p><code>lj/expand</code> is LJ with the whole curve moved outwards by a distance Δ: E(r) = 4ε[(σ/(r−Δ))¹² − (σ/(r−Δ))⁶], with the cutoff also moved to r<sub>c</sub> + Δ. It is used for particles of different sizes. Δ is the 5th word of pair_coeff. With Δ = 0 it must give exactly the lj/cut numbers.</p>
      <p>Its file (section 3 rule: <code>src/pair_lj_expand.cpp</code>) is 90 % the same as <code>pair_lj_cut.cpp</code>. The maths inside the double loop over atoms i and their neighbours j:</p>
<pre class="lmp-mini">106  if (rsq &lt; cutsq[itype][jtype]) {          // r² below (rc + Δ)²
107    r = sqrt(rsq);
108    rshift = r - shift[itype][jtype];       // r − Δ   ← the only new idea
109    rshiftsq = rshift * rshift;
110    r2inv = 1.0 / rshiftsq;
111    r6inv = r2inv * r2inv * r2inv;          // 1/(r−Δ)⁶
112    forcelj = r6inv * (lj1 * r6inv - lj2);  // lj1 = 48εσ¹², lj2 = 24εσ⁶
113    fpair = factor_lj * forcelj / rshift / r;
       …
125    evdwl = r6inv * (lj3 * r6inv - lj4) - offset;   // lj3 = 4εσ¹², lj4 = 4εσ⁶</pre>
      <p class="muted">${gh('src/pair_lj_expand.cpp', 106)}. Lines 112–113 are the force F = 24ε[2σ¹²/(r−Δ)¹³ − σ⁶/(r−Δ)⁷], divided by r so it can be multiplied by Δx, Δy, Δz. (The <code>[itype][jtype]</code> after lj1 … lj4 is dropped here for space.)</p>
    </div>
    <div class="row"><label class="lab" style="margin:0">Δ =</label><input class="lmp-in st-d" type="number" value="0.1" min="0" max="0.3" step="0.01" style="width:100px"><span class="muted">try 0, 0.1, 0.2</span></div>
    <pre class="lmp-mini st-script"></pre>
  </section>

  <section class="block prose">
    <h2>5 · Predict by hand, then run</h2>
    <p>At step 0 the atoms sit on the FCC crystal (page A.5), so every atom has the same neighbours: shells of 12, 6, 24, 12, 24 … atoms. Energy per atom = ½ Σ n·E(r) (the ½ because each pair is shared by two atoms).</p>
    <div class="tbl st-shells"></div>
    <div class="calc st-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run with LAMMPS</button><span class="lmp-status muted st-st"></span></div>
    <div class="lmp-compare st-cmp"></div>
  </section>

  <section class="block prose">
    <h2>6 · Observe: a style that is not in this build</h2>
    <p><code>lj/smooth/linear</code> is a real LAMMPS pair style. Find its file in section 3 first, then run it here.</p>
    <div class="row"><button class="btn sm primary" data-k="sl">▶ Run with lj/smooth/linear</button><span class="lmp-status muted sl-st"></span></div>
    <pre class="lmp-mini sl-out" hidden></pre>
    <p>The file is in the <b>EXTRA-PAIR</b> package, and the LAMMPS in this browser was built without it. Same command, same syntax: whether it works depends on how the program was built. On your own computer: <code>cmake -D PKG_EXTRA-PAIR=on</code> and build again (page 0.4).</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Which file carries out <code>fix 2 all npt temp 1 1 0.5 iso 0 0 5</code>? Predict, then check with the finder.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>src/fix_npt.cpp</code>. It is only a few lines long: it sets options and inherits everything from <code>fix_nh.cpp</code>, where the thermostat and barostat maths lives (inheritance: page B.3).</div></div>
    <div class="ex"><p><b>2.</b> With Δ = 0.2, E_pair rises from −7.47 to −0.62 and the pressure from 5 to 58. Which shell causes it? Use the shell table.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The first shell. Its 12 atoms now sit at r − Δ = 0.988, inside the repulsive wall (below 2<sup>1/6</sup> = 1.12), so their contribution turns from about −5.7 to about +2.0. The 5th shell (24 atoms at r = 2.6557) also enters the 2.7 cutoff, giving 78 neighbours instead of 54, but it adds only a small negative amount.</div></div>
    <div class="ex"><p><b>3.</b> The examples folder is the largest (about 217 MB). Does it make the program bigger?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. Only <code>src/</code> (the core plus switched-on packages) is compiled. examples, doc, potentials and tools are never part of <code>lmp</code>.</div></div>
  </section>`,

  mount(root) {
    load().then(D => setup(root, D)).catch(() => { root.querySelector('.st-top').innerHTML = '<p class="lmp-err">Could not load the source index.</p>'; });
    return () => {};
  },
};

function setup(root, D) {
    const $ = s => root.querySelector(s);
    $('.st-top').innerHTML = `<table>${head(['Folder', 'Files', 'MB', 'What it holds'])}${D.top.map(([d, n, mb]) => row([`<code>${d}/</code>`, n, mb, FOLDERS[d] || ''])).join('')}</table>`;
    const pk = D.pkgs.filter(p => !NOT_PKG.has(p[0]));
    const big = Math.max(...pk.map(p => p[1]));
    $('.st-pk').innerHTML = `<p><b>${D.core}</b> .cpp files in the core and <b>${pk.length}</b> packages. Tap a package:</p><div class="chips">${pk.map(([n, c]) => `<button class="chip" data-pk="${n}" title="${c} .cpp files" style="--w:${Math.max(4, 100 * c / big)}%">${n} <small>${c}</small></button>`).join('')}</div>`;
    root.querySelectorAll('[data-pk]').forEach(b => b.onclick = () => {
      const n = b.dataset.pk, st = D.styles.filter(s => s[4] === n), kinds = {};
      st.forEach(s => { (kinds[s[0]] = kinds[s[0]] || []).push(s[1]); });
      root.querySelectorAll('[data-pk]').forEach(x => x.classList.toggle('on', x === b));
      $('.st-pkinfo').innerHTML = `<p><b>${n}</b>: ${st.length} styles. ${Object.entries(kinds).map(([k, v]) => `<br><code>${k}</code>: ${v.slice(0, 30).map(esc).join(', ')}${v.length > 30 ? ` … (+${v.length - 30})` : ''}`).join('')}</p>`;
    });
    // how good is the rule?
    let exact = 0, suf = 0; const odd = [];
    for (const [k, n, , p] of D.styles) {
      const f = p.split('/').pop(); if (f === predict(k, n)) { exact++; continue; }
      const m = n.match(/^(.*)\/(omp|kk|kk\/host|kk\/device|gpu|intel|opt)$/);
      if (m && f === PREFIX[k] + m[1].replace(/[/-]/g, '_') + SUFFIX[m[2]][0] + '.cpp') suf++; else odd.push(n);
    }
    $('.st-rule').innerHTML = `Checked on all ${D.styles.length} styles registered in ${TAG}: the plain rule gives the right file for <b>${exact}</b>, the suffix rule for ${suf} more, and only <b>${odd.length}</b> are exceptions: old names kept so that old scripts still run, and a few files named differently (e.g. <code>lj/sdk</code> lives in <code>pair_lj_spica.cpp</code>).`;
    const go = () => { $('.st-res').innerHTML = card(find($('.st-q').value, D)); };
    $('[data-k=find]').onclick = go; $('.st-q').onkeydown = e => { if (e.key === 'Enter') go(); }; go();
    $('[data-k=melt]').onclick = () => {
      const lines = IN_MELT.split('\n').filter(l => l.trim() && !l.trim().startsWith('#'));
      $('.st-melt').innerHTML = `<table>${head(['Line', 'Read in', 'Done by'])}${lines.map(l => { const r = find(l, D); const by = r.style ? gh(r.style[3]) : r.bi && r.bi[2] ? gh('src/input.cpp', r.bi[2]) : ''; return row([`<code>${esc(l.trim().replace(/\s+/g, ' '))}</code>`, r.bi ? gh('src/input.cpp', r.bi[1]) : 'a command class', by]); }).join('')}</table>`;
    };
    // lj/expand
    let H = null;
    const draw = () => {
      let d = parseFloat($('.st-d').value); if (!isFinite(d) || d < 0) d = 0; if (d > 0.3) d = 0.3;
      $('.st-script').textContent = SCRIPT(d).split('\n').filter(l => /pair_style|pair_coeff|^run/.test(l)).join('\n');
      H = handExpand(d); H.d = d;
      $('.st-shells').innerHTML = `<table>${head(['shell r', 'atoms n', 'r − Δ', '½ n · 4[(1/(r−Δ))¹² − (1/(r−Δ))⁶]'])}${H.rows.map(x => row([g8(x.r), x.n, g8(x.s), g8(x.e)])).join('')}<tr><td colspan="3"><b>E_pair = sum</b></td><td><b>${g8(H.E)}</b></td></tr></table>`;
      $('.st-hand').innerHTML = esc(`cutoff rc + Δ = ${g8(2.5 + d)}: ${H.rows.length} shells, ${H.rows.reduce((s, x) => s + x.n, 0)} neighbours per atom
E_pair = `) + `<b>${g8(H.E)}</b>` + esc(`
TotEng = E_pair + 1.5 T (N−1)/N = ${g8(H.E)} + ${g8(H.kin)} = `) + `<b>${g8(H.Tot)}</b>` + esc(`
Press  = ρT(N−1)/N + (ρ/3) · ½ Σ n r F(r),  F = 24[2/(r−Δ)¹³ − 1/(r−Δ)⁷]
       = ${g8(H.idp)} + ${g8(H.P - H.idp)} = `) + `<b>${g8(H.P)}</b>`;
      $('.st-cmp').innerHTML = '';
    };
    $('.st-d').oninput = draw; draw();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.st-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(SCRIPT(H.d)); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const T = parseThermo(res.lines), r0 = T.rows.find(w => w[0] === '0'), c = T.cols;
      const want = [['E_pair', H.E], ['TotEng', H.Tot], ['Press', H.P]];
      $('.st-cmp').innerHTML = `<div class="tbl"><table><tr><th>Δ = ${H.d}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${want.map(([k, h]) => { const lv = r0[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
    };
    $('[data-k=sl]').onclick = async e => {
      const b = e.currentTarget, st = $('.sl-st'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(IN_MELT.replace('pair_style      lj/cut 2.5', 'pair_style      lj/smooth/linear 2.5').replace('pair_coeff      1 1 1.0 1.0 2.5', 'pair_coeff      1 1 1.0 1.0')); b.disabled = false; st.textContent = '';
      const i = res.lines.findIndex(x => /^ERROR/.test(x));
      const o = $('.sl-out'); o.textContent = i >= 0 ? res.lines.slice(i, i + 3).join('\n') : (res.error || 'It ran.'); o.hidden = false;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
}
