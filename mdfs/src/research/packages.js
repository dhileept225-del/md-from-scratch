// Research · B.5: packages. What is in this LAMMPS, counted by hand from the source tree; which packages a script needs, predicted and then run.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;

let DATA = null;
const load = () => DATA ? Promise.resolve(DATA) : fetch('/data/lammps-src.json').then(r => r.json()).then(d => (DATA = d));

const PKGS = [
  ['MOLECULE', 'bonds, angles, dihedrals, impropers; atom_style full, bond, angle', 'any molecule'],
  ['KSPACE', 'long-range Coulomb: Ewald, PPPM, and the pair styles that go with them (…/coul/long)', 'charged systems, water'],
  ['RIGID', 'rigid bodies, and fix shake / rattle for rigid bonds', 'rigid water models'],
  ['MANYBODY', 'many-body potentials: EAM, Tersoff, Stillinger–Weber …', 'metals, semiconductors'],
  ['EXTRA-PAIR, EXTRA-FIX, EXTRA-COMPUTE, EXTRA-MOLECULE', 'less common styles of each family', 'when the manual names one'],
  ['DRUDE', 'Drude-oscillator polarizable models', 'polarizable force fields'],
  ['ML-SNAP, ML-PACE, ML-IAP …', 'machine-learning potentials', 'MLIPs'],
  ['OPENMP, INTEL, OPT, GPU, KOKKOS', 'the same styles again, faster on threads or GPUs (suffix /omp, /intel, /opt, /gpu, /kk)', 'speed (page F.6)'],
];
const FAM = [['atom_style', 'atom styles'], ['bond_style', 'bond styles'], ['angle_style', 'angle styles'], ['dihedral_style', 'dihedral styles'], ['improper_style', 'improper styles'], ['pair_style', 'pair styles'], ['fix', 'fix styles'], ['compute', 'compute styles']];
const SEC = { Atom: 'atom_style', Pair: 'pair_style', Bond: 'bond_style', Angle: 'angle_style', Dihedral: 'dihedral_style', Improper: 'improper_style', Fix: 'fix', Compute: 'compute' };

const SCRIPTS = {
  water: `units           real
atom_style      full
region          box block 0 20 0 20 0 20
create_box      2 box bond/types 1 angle/types 1 extra/bond/per/atom 2 extra/angle/per/atom 1 extra/special/per/atom 2
pair_style      lj/cut/coul/long 10.0
kspace_style    pppm 1.0e-4
bond_style      harmonic
angle_style     harmonic
fix             1 all shake 1e-4 20 0 b 1 a 1
`,
  silicon: `units           metal
atom_style      atomic
lattice         diamond 5.431
region          box block 0 3 0 3 0 3
create_box      1 box
create_atoms    1 box
mass            1 28.0855
pair_style      tersoff
`,
  melt: `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 5 0 5 0 5
create_box      1 box
create_atoms    1 box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
compute         g all rdf 50
run             0
`,
};

// every style a script names, and the package it lives in
function needs(script, D) {
  const at = { fix: 3, compute: 3, dump: 3, region: 2 }, kinds = new Set(['pair_style', 'bond_style', 'angle_style', 'dihedral_style', 'improper_style', 'kspace_style', 'atom_style', 'fix', 'compute', 'dump', 'region', 'min_style', 'run_style']);
  const out = [];
  script.split('\n').forEach((l, i) => {
    const w = l.replace(/#.*/, '').trim().split(/\s+/); if (!kinds.has(w[0])) return;
    const name = w[at[w[0]] || 1]; if (!name || name === 'none') return;
    const st = D.styles.find(s => s[0] === w[0] && s[1] === name);
    out.push({ line: i + 1, kind: w[0], name, pkg: st ? (st[4] || 'core') : '?', file: st && st[3] });
  });
  return out;
}
const HERE = new Set(['core', 'MOLECULE']);

export default {
  id: 'packages',
  title: 'Packages',
  minutes: 15,
  summary: 'Most of LAMMPS is optional: 94 packages, each a folder of styles that is compiled only if you switch it on. You find out which packages the LAMMPS in this page was built with, count by hand how many styles that should give, and check the count against LAMMPS. Then you paste a script, predict which packages it needs and which line will fail, and run it to see.',
  learn: [
    'What a package is, and why LAMMPS is split into them',
    'Ask a LAMMPS binary which packages it has (info config, lmp -h)',
    'Count the styles a build should have, by hand, from the source tree',
    'Read a script and predict the packages it needs before you run it',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: a core and many packages</h2>
    <p>Page B.1 showed that <code>src/</code> holds the core (always compiled) and 94 upper-case folders. Each folder is a <b>package</b>: a group of styles that belong together. You choose them when you build (page 0.4):</p>
<pre class="lmp-mini">cmake -D PKG_MOLECULE=on -D PKG_KSPACE=on -D PKG_RIGID=on ../cmake
cmake --build . -j 4</pre>
    <p>Why not compile everything? Some packages need outside libraries (GPU drivers, Python, a machine-learning library), each one adds compile time, and every package switched on is code that must compile on your machine. Most users need five or six.</p>
    <div class="tbl wide"><table>${head(['Package', 'What it adds', 'You need it for'])}${PKGS.map(([p, w, n]) => row([`<code>${p}</code>`, w, n])).join('')}</table></div>
    <p>A style's package is written nowhere in the script. It is decided by which folder its file sits in. The registration step (page B.3) only collects styles from folders that were switched on; for all other styles the build writes a list of names and their packages (<code>packages_pair.h</code> and so on, read in ${gh('src/lammps.cpp', 1050)}), so that an unknown style can be answered with the package to add (${gh('src/utils.cpp', 444)}).</p>
  </section>

  <section class="block">
    <h2>2 · Observe: what is in this LAMMPS?</h2>
    <p class="muted">On your own computer, <code>lmp -h</code> prints the same list near the end, under "Installed packages".</p>
    <div class="row"><button class="btn sm primary" data-k="cfg">▶ Run info config</button><span class="lmp-status muted pk-st"></span></div>
    <pre class="lmp-mini pk-cfg" hidden></pre>
  </section>

  <section class="block">
    <h2>3 · Match: count the styles by hand</h2>
    <div class="prose"><p>If only the core and MOLECULE are built in, the styles available must be exactly those registered in <code>src/*.h</code> plus those in <code>src/MOLECULE/*.h</code>. Count them from the source tree of page B.1, then ask LAMMPS with <code>info styles</code>.</p></div>
    <div class="row"><button class="btn sm primary" data-k="count">▶ Count by hand and run info styles</button><span class="lmp-status muted pc-st"></span></div>
    <div class="tbl pc-cmp"></div>
    <div class="prose pc-extra"></div>
  </section>

  <section class="block">
    <h2>4 · Modify: which packages does a script need?</h2>
    <div class="prose"><p>Choose a script or paste your own. The page reads every style it names, finds its file in the source tree, and predicts the first line this LAMMPS (core + MOLECULE) will stop at. Then run it.</p></div>
    <div class="row"><button class="btn sm ghost" data-s="water">water (SPC/E-style setup)</button><button class="btn sm ghost" data-s="silicon">silicon (Tersoff)</button><button class="btn sm ghost" data-s="melt">LJ melt</button></div>
    <textarea class="lmp-in pn-in" spellcheck="false" wrap="off" style="height:230px"></textarea>
    <div class="tbl pn-need"></div>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run it</button><span class="lmp-status muted pn-st"></span></div>
    <pre class="lmp-mini pn-out" hidden></pre>
  </section>

  <section class="block prose">
    <h2>What the packages show</h2>
    <ul>
      <li>A rigid water model with long-range Coulomb needs three packages: MOLECULE for bonds, KSPACE for PPPM and the <code>coul/long</code> pair style, RIGID for <code>fix shake</code>. The water script stops at its first missing one; delete that line (or switch to <code>lj/cut/coul/cut</code>) and it stops at the next. LAMMPS reports one problem at a time; reading the script against the source tree finds them all at once.</li>
      <li>Some everyday styles are in the core even though they sound special: <code>fix nvt</code>, <code>fix npt</code>, <code>fix langevin</code>, <code>compute rdf</code>, <code>compute msd</code>.</li>
      <li>When you install LAMMPS from a package manager or use a cluster module, it was built by someone else. Run <code>lmp -h</code> once and keep the list.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Which packages does a script with <code>pair_style eam/alloy</code>, <code>fix nvt</code> and <code>compute msd</code> need?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only <b>MANYBODY</b> (for eam/alloy). fix nvt and compute msd are in the core. Check with the finder in section 4.</div></div>
    <div class="ex"><p><b>2.</b> The core has 24 pair styles and MOLECULE adds 7. KSPACE adds how many pair styles? Find out from the source tree of page B.1.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Tap KSPACE on page B.1: its pair styles are the <code>…/coul/long</code>, <code>…/coul/msm</code> and similar long-range versions. Building with KSPACE would raise the pair count by that number, which you can check with info styles on your own build.</div></div>
    <div class="ex"><p><b>3.</b> You switch on OPENMP. Does the number of physical models you can use go up?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. OPENMP adds about 200 styles, but each is a faster copy of an existing style (lj/cut/omp computes the same as lj/cut, using several threads). The models stay the same; the speed changes.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=cfg]').onclick = async e => {
      const b = e.currentTarget, st = $('.pk-st'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps('info config'); b.disabled = false; st.textContent = '';
      const L = res.lines, k = L.findIndex(l => /^LAMMPS version/.test(l)), j = L.findIndex(l => /^Installed packages/.test(l));
      const pick = [...L.slice(k, k + 1), ...L.filter(l => /^OS information|^Compiler/.test(l)), '', ...L.slice(j, j + 3).filter(l => l.trim())];
      const o = $('.pk-cfg'); o.textContent = pick.join('\n'); o.hidden = false;
    };
    $('[data-k=count]').onclick = async e => {
      const b = e.currentTarget, st = $('.pc-st'); b.disabled = true; st.textContent = 'Loading the source index and running LAMMPS…';
      const [D, res] = await Promise.all([load(), runLammps('info styles')]); b.disabled = false; st.textContent = '';
      const br = {}; let cur = null;
      for (const l of res.lines) { const m = l.match(/^(\w+) styles:/); if (m) { cur = SEC[m[1]] || null; if (cur) br[cur] = []; continue; } if (/^Info-Info/.test(l) || /^\w[\w ]* styles?:?$/.test(l.trim()) && !cur) { cur = null; continue; } if (cur && l.trim()) br[cur].push(...l.trim().split(/\s+/)); }
      const extra = {};
      const rows = FAM.map(([k, n]) => {
        const core = D.styles.filter(s => s[0] === k && s[4] === '').map(s => s[1]), mol = D.styles.filter(s => s[0] === k && s[4] === 'MOLECULE').map(s => s[1]);
        const want = new Set([...core, ...mol]), got = br[k] || [];
        extra[k] = got.filter(x => !want.has(x));
        const missing = [...want].filter(x => !got.includes(x));
        const ok = !missing.length && !extra[k].length;
        return row([n, core.length, mol.length, `<b>${want.size}</b>`, got.length, ok ? '<span class="m-ok">✓</span>' : `<span class="muted">+${extra[k].length}${missing.length ? ` −${missing.length}` : ''}</span>`]);
      });
      $('.pc-cmp').innerHTML = `<table>${head(['family', 'core', 'MOLECULE', 'by hand', 'LAMMPS', ''])}${rows.join('')}</table>`;
      const ex = Object.entries(extra).filter(([, v]) => v.length);
      $('.pc-extra').innerHTML = ex.length ? `<p>Where LAMMPS lists more, the extra names are styles that were <b>retired</b> and only print "no longer available", or were added after July 2025 (this browser runs the 30 Sep 2026 version):</p>
        <p>${ex.map(([k, v]) => `<code>${k}</code>: ${v.map(x => `<button class="btn sm ghost" data-ret="${k} ${x}">${esc(x)}</button>`).join(' ')}`).join('<br>')}</p><p class="muted">Tap a name to run it.</p><pre class="lmp-mini pc-ret" hidden></pre>` : '';
      root.querySelectorAll('[data-ret]').forEach(x => x.onclick = async () => {
        const [k, n] = x.dataset.ret.split(' '), s = `units real\natom_style full\nregion b block 0 10 0 10 0 10\ncreate_box 1 b\n${k === 'fix' ? `fix 1 all ${n}` : k === 'compute' ? `compute 1 all ${n}` : `${k} ${n}`}\n`;
        const r = await runLammps(s), i = r.lines.findIndex(l => /^ERROR/.test(l)), o = $('.pc-ret');
        o.textContent = s.split('\n').slice(-2, -1)[0] + '\n→ ' + (i >= 0 ? r.lines[i] : (r.error || 'accepted: this style works here')); o.hidden = false;
      });
    };
    const ta = $('.pn-in');
    const analyse = async () => {
      const D = await load(), N = needs(ta.value, D);
      const first = N.find(n => !HERE.has(n.pkg));
      const pk = [...new Set(N.map(n => n.pkg))];
      $('.pn-need').innerHTML = `<table>${head(['line', 'style', 'package', 'here?'])}${N.map(n => row([n.line, `<code>${n.kind} ${esc(n.name)}</code>`, n.pkg === '?' ? '<span class="m-bad">not a style in 22 Jul 2025</span>' : `<b>${n.pkg}</b>`, HERE.has(n.pkg) ? '<span class="m-ok">yes</span>' : '<span class="m-bad">no</span>'])).join('')}</table>
        <p>Packages needed: ${pk.map(p => `<b>${p}</b>`).join(', ') || 'core only'}. Prediction: ${first ? `LAMMPS stops at <b>line ${first.line}</b> (${first.kind} ${esc(first.name)}, package ${first.pkg}).` : 'it runs here.'}</p>`;
      $('.pn-out').hidden = true;
    };
    root.querySelectorAll('[data-s]').forEach(b => b.onclick = () => { ta.value = SCRIPTS[b.dataset.s]; analyse(); });
    ta.oninput = () => { clearTimeout(ta._t); ta._t = setTimeout(analyse, 300); };
    ta.value = SCRIPTS.water; analyse();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.pn-st'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(ta.value); b.disabled = false; st.textContent = '';
      const i = res.lines.findIndex(l => /^ERROR/.test(l)), o = $('.pn-out');
      const last = res.lines.filter(l => /^Last input line:/.test(l)).pop();
      o.textContent = i >= 0 ? [res.lines[i], last || ''].join('\n') : res.error ? res.error : 'No error: every style in the script is built into this LAMMPS.';
      o.hidden = false;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
