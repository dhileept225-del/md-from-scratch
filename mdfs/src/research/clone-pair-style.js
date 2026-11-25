// Research · E.2: clone a pair style. pair_lj_cut copied to pair_lj_mine with four renames, rebuilt natively;
// the 4000-atom melt with lj/mine (native) compared bit for bit with lj/cut run live in the browser.
import { runLammps } from './lammps-runner.js';
import { NATIVE } from './clone-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

const MELT = style => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 3.0 87287 loop geom
pair_style      ${style} 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nve
thermo_style    custom step temp epair emol etotal press
thermo_modify   format float %.17g
thermo          50
run             250`;
const COMMANDS = `cd lammps/src
cp pair_lj_cut.h   pair_lj_mine.h
cp pair_lj_cut.cpp pair_lj_mine.cpp
sed -i 's/PairLJCut/PairLJMine/g; s/PAIR_LJ_CUT_H/PAIR_LJ_MINE_H/g;
        s/lj\\/cut,/lj\\/mine,/; s/pair_lj_cut\\.h/pair_lj_mine.h/g' pair_lj_mine.h pair_lj_mine.cpp
cd ../build
make -j2`;
const changed = d => d.split('\n').filter(l => /^[-+][^-+]/.test(l)).join('\n');

export default {
  id: 'clone-pair-style',
  title: 'Clone a pair style',
  minutes: 15,
  summary: 'A new pair style is two files and one line that gives it a name. You copy lj/cut to lj/mine with four renames, watch CMake notice the new file and regenerate its list of styles, and check that the clone reproduces lj/cut exactly: a 4000-atom melt with lj/mine (native) against lj/cut run live in your browser, all thirty thermo numbers identical to the last bit. Then the two mistakes everyone makes, with the compiler errors they produce, and what the browser build says when asked for a style it was not compiled with.',
  learn: [
    'What makes a style: a class derived from Pair and the PairStyle(name, Class) line',
    'How CMake finds new files and generates style_pair.h',
    'Checking a clone: identical numbers, not just similar ones',
    'The header guard and the style name: two renames people forget',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: how LAMMPS finds a style</h2>
    <div class="prose"><p>Every pair style is a class derived from <code>Pair</code>, declared in a header that starts with one special line:</p>
<pre class="lmp-mini">#ifdef PAIR_CLASS
PairStyle(lj/cut,PairLJCut);       // the name in input scripts, the C++ class
#else
#ifndef LMP_PAIR_LJ_CUT_H          // header guard
… class PairLJCut : public Pair { … };</pre>
<p>At configure time CMake collects every <code>src/*.cpp</code> (<code>file(GLOB … CONFIGURE_DEPENDS …)</code>, cmake/CMakeLists.txt:243) and writes <code>style_pair.h</code>, an #include of every header that contains a PairStyle line. ${gh('src/force.cpp', 88)} includes that list twice: once for the class declarations, once with PAIR_CLASS defined, where each PairStyle line becomes</p>
<pre class="lmp-mini">(*pair_map)["lj/cut"] = &amp;style_creator&lt;Pair, PairLJCut&gt;;</pre>
<p><code>pair_style lj/cut</code> in a script looks the name up in this map and calls the creator. So a new style needs a new class name, a new style name, a new header guard, and a rebuild; nothing else in LAMMPS has to know about it.</p></div>
  </section>

  <section class="block">
    <h2>2 · Change: copy, rename, rebuild</h2>
    <pre class="lmp-mini">${esc(COMMANDS)}</pre>
    <div class="prose"><p>The header changes in four places:</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.hdiff.split('\n').slice(2).join('\n'))}</pre>
    <div class="prose"><p>The .cpp file only in the include and the class name, 21 lines:</p></div>
    <pre class="lmp-mini gd">${esc(changed(NATIVE.cppdiff).split('\n').slice(0, 10).join('\n'))}\n…</pre>
    <div class="prose"><p>The build (native, 13 s on two cores):</p></div>
    <pre class="lmp-mini">${esc([NATIVE.make[0], ...NATIVE.make.slice(3)].join('\n'))}</pre>
    <div class="prose"><p>"GLOB mismatch" is CMake seeing a .cpp file that was not there at the last configure; it re-runs itself, regenerates style_pair.h, and therefore recompiles force.cpp (which includes it) along with the new file. <code>lmp -h</code> now lists <code>lj/mine</code> among the pair styles.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: lj/mine against lj/cut</h2>
    <div class="prose"><p>The standard melt: 4000 atoms, T = 3.0 on an fcc lattice, 250 steps. Native LAMMPS ran it with <code>pair_style lj/mine</code>; your browser runs the original <code>lj/cut</code>.</p></div>
    <pre class="lmp-mini">${esc(MELT('lj/mine'))}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run lj/cut here and compare</button><span class="lmp-status muted cp-sta"></span></div>
    <div class="tbl cp-a"></div>
    <div class="prose cp-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: a style the program does not have</h2>
    <div class="prose"><p>The browser build was compiled from the original sources. Ask it for the clone:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run pair_style lj/mine here</button><span class="lmp-status muted cp-stb"></span></div>
    <pre class="lmp-mini cp-b" hidden></pre>
    <div class="prose cp-bsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You rename everything except the header guard (it stays LMP_PAIR_LJ_CUT_H). What happens?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>style_pair.h includes pair_lj_cut.h first; when it reaches pair_lj_mine.h the guard is already defined, so the class declaration is skipped, but the PairStyle line (outside the guard) is still read. force.cpp fails: <code>${esc(NATIVE.errGuard)}</code>. We tried it.</div></div>
    <div class="ex"><p><b>2.</b> You rename everything except the PairStyle line, which still says lj/cut. What happens?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>It compiles without a warning. Two lines now assign pair_map["lj/cut"]; the one from pair_lj_mine.h comes later in style_pair.h and wins, so <code>pair_style lj/cut</code> silently creates your class, and lj/mine does not exist (<code>lmp -h</code> does not list it; we checked). Harmless while the code is identical, very confusing once you change it.</div></div>
    <div class="ex"><p><b>3.</b> lj/cut has accelerated versions (lj/cut/omp, lj/cut/gpu, lj/cut/kk). Does <code>-sf omp</code> speed up lj/mine?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. With a suffix LAMMPS first tries "lj/mine/omp", does not find it, and falls back to plain lj/mine. Each accelerated variant is its own class in its package directory; you would clone those too.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.cp-sta'); b.disabled = true; st.textContent = 'Running 250 steps of 4000 atoms…';
      const res = await runLammps(MELT('lj/cut')); b.disabled = false;
      if (res.error) return err(st, res.error); st.textContent = '';
      const W = res.thermo.rows; let ok = 0, n = 0;
      const rows = NATIVE.mine.map((m, k) => { const w = W[k]; let same = 0; for (let c = 1; c < 6; c++) { n++; if (w && Number(w[c]) === Number(m[c])) { same++; ok++; } } return row([m[0], `<span class="mono">${esc(m[2])}</span>`, `<span class="mono">${esc(w ? w[2] : '')}</span>`, `<span class="mono">${esc(m[5])}</span>`, same === 5 ? `${OK} 5 of 5` : `${BAD} ${same} of 5`]); });
      $('.cp-a').innerHTML = `<table>${head(['step', 'E_vdwl, lj/mine (native)', 'E_vdwl, lj/cut (your browser)', 'Press, lj/mine', 'all 5 columns'])}${rows.join('')}</table>`;
      $('.cp-asum').innerHTML = `<p>${ok} of ${n} numbers identical, across two programs (native gcc and WebAssembly) and two class names. Native lj/cut gave the same thirty numbers as well, and <code>pair_write</code> tables of the two styles are identical. That is the test of a clone: not "close", identical, because the arithmetic is the same line for line. Timing is not a test: native ${esc(NATIVE.loop.cut.replace('Loop time of ', '').split(' on')[0])} s for lj/cut and ${esc(NATIVE.loop.mine.replace('Loop time of ', '').split(' on')[0])} s for lj/mine is run-to-run noise from the same code. Now lj/mine is yours to change (page E.3), and any difference from lj/cut will be your change and nothing else.</p>`;
    };
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.cp-stb'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(MELT('lj/mine')); b.disabled = false; st.textContent = '';
      const lines = res.lines.filter(l => /ERROR|Last input line/.test(l));
      $('.cp-b').hidden = false; $('.cp-b').textContent = (lines.length ? lines : [res.error || '(no error?)']).join('\n');
      $('.cp-bsum').innerHTML = '<p>The name is looked up at run time, so an unknown style is not a compile error, only an input error, raised by Force::create_pair when the map has no entry. You get the same error if you add the files but forget to rebuild, and a similar one if you read a restart file written with lj/mine into an executable that lacks it. When LAMMPS knows the name belongs to an uninstalled package it says which one; for your own style it cannot.</p>';
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
