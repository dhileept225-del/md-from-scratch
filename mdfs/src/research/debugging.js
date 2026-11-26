// Research · F.4: debugging. Three failures, three tools: an ERROR whose cause is elsewhere (gdb on Error::all,
// per-atom probes), a crash from a planted allocation bug (gdb, valgrind, AddressSanitizer), and a NaN that
// raises no error at all (floating-point traps). Each diagnosis is redone by hand.
import { NATIVE } from './debugging-native.js';
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const BASE = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 4 0 4 0 4
create_box 1 box
create_atoms 1 box
create_atoms 1 single 0.05 0 0
mass 1 1.0
velocity all create 1.0 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
neighbor 0.3 bin
neigh_modify every 1 delay 0 check yes
fix 1 all nve`;
const OV = `${BASE}
thermo_style custom step temp pe press
thermo 1
run 10`;
const PROBE = `${BASE}
compute pe all pe/atom
dump d all custom 1 a.dump id x y z fx c_pe
dump_modify d sort id format float %.17g
run 0`;
const SAME = OV.replace('create_atoms 1 single 0.05 0 0', 'create_atoms 1 single 0 0 0').replace('neighbor 0.3 bin\nneigh_modify every 1 delay 0 check yes\n', '').replace('run 10', 'run 3');
const GDBA = `gdb -q -batch -ex "set breakpoint pending on" \\
    -ex "break LAMMPS_NS::Error::all" \\
    -ex "run -in in.ov -log none -screen none" -ex "bt" lmp`;
const BUG = `--- a/src/pair_lj_cut.cpp        (PairLJCut::allocate, line 400)
-  int n = atom->ntypes + 1;
+  int n = atom->ntypes;`;
const CREATE = `// src/memory.h:144-160, memory->create(array, n1, n2, name), shortened
TYPE *data   = smalloc(sizeof(TYPE)   * n1 * n2);   // line 150: the numbers
TYPE **array = smalloc(sizeof(TYPE *) * n1);        // line 152: one pointer per row
for (i = 0; i < n1; i++) array[i] = &data[i * n2];`;
const COEFF = `// src/pair_lj_cut.cpp:458-466, pair_coeff * * eps sigma cut
for (int i = ilo; i <= ihi; i++)
  for (int j = MAX(jlo, i); j <= jhi; j++) {
    epsilon[i][j] = epsilon_one;      // line 460
    ...`;
const ASANB = `mkdir build-dbg && cd build-dbg
cmake ../cmake -D CMAKE_BUILD_TYPE=Debug -D ENABLE_SANITIZER=address \\
      -D BUILD_MPI=no -D PKG_MOLECULE=yes -D PKG_KSPACE=yes
make -j2 lmp                     # 172 s here, alongside the normal build`;
const FPEC = `// fpe.c: make the CPU raise SIGFPE at the first invalid, divide-by-zero or overflow
#define _GNU_SOURCE
#include <fenv.h>
__attribute__((constructor)) static void trap(void)
{ feenableexcept(FE_INVALID | FE_DIVBYZERO | FE_OVERFLOW); }

gcc -shared -fPIC -O2 -o libfpe.so fpe.c -lm
LD_PRELOAD=./libfpe.so lmp -in in.same          # Floating point exception, exit 136
gdb -ex "set environment LD_PRELOAD=./libfpe.so" -ex "run -in in.same" \\
    -ex bt -ex "print rsq" -ex "print atom->tag[i]" -ex "print atom->tag[j]" build-dbg/lmp`;
const CHEAT = `cmake -D CMAKE_BUILD_TYPE=RelWithDebInfo ...   # -O2 -g: fast, with file:line in gdb
gdb --args lmp -in in.x            # then: run, bt, frame 2, print x[i][0], up, list
break LAMMPS_NS::Error::all        # stop at any ERROR, see who raised it
break pair_lj_cut.cpp:114 if rsq < 0.5   # conditional breakpoint
watch atom->x[17][0]               # stop when this memory changes (slow, but exact)
valgrind lmp -in in.x              # no rebuild needed, about 25x slower
mpirun -np 4 xterm -e gdb --args lmp -in in.x   # one debugger per rank`;

export default {
  id: 'debugging',
  title: 'Debugging',
  minutes: 16,
  summary: 'Three ways a run goes wrong, each traced to its cause. An ERROR that LAMMPS raises far from the real problem; a segmentation fault from a one-character allocation bug, found by gdb, valgrind and AddressSanitizer; and a NaN that LAMMPS never reports, caught by making the CPU trap at the first invalid operation. Every diagnosis is then repeated by hand: the force that threw atoms out of the box, the exact byte offset the sanitizer reports, and the IEEE arithmetic that turned one zero into inf and nan.',
  learn: [
    'Where an ERROR was raised vs. where the problem is',
    'gdb, valgrind and -D ENABLE_SANITIZER=address on a real memory bug',
    'Finding the first NaN with floating-point traps',
    'Probing a bad state with compute pe/atom and a dump at step 0',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three symptoms</h2>
    <div class="tbl"><table>${head(['symptom', 'what it means', 'first tool'])}
    ${row(['ERROR: … (src/file.cpp:NNN)', 'LAMMPS noticed something; the file:line is where it noticed', 'look at earlier output; break on Error::all'])}
    ${row(['Segmentation fault', 'memory used outside what was allocated', 'gdb bt, then valgrind or AddressSanitizer'])}
    ${row(['nan or inf in thermo, no error', 'an arithmetic exception, silently propagated', 'thermo 1, then a floating-point trap'])}</table></div>
  </section>

  <section class="block">
    <h2>2 · An ERROR far from its cause</h2>
    <div class="prose"><p>A 256-atom crystal plus one more atom, placed 0.05 lattice units from atom 1 (a typical mistake when building a system by hand):</p></div>
    <pre class="lmp-mini">${esc(OV)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run it</button><span class="db-sta"></span></div>
    <div class="prose db-asum"></div>
    <div class="prose"><p>The message names ${gh('src/thermo.cpp', 526)}. That is only where LAMMPS counts the atoms. To see who called it, stop the program at the error:</p></div>
    <pre class="lmp-mini">${esc(GDBA)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.bt)}</pre>
    <div class="prose"><p>Thermo, called from the output at step 1, compared the atom count; the atoms were lost inside the step itself. gdb tells you where the error is raised, and here the cause is in the step-0 line: a potential energy of 10¹¹ per atom. The fix is in the input (<code>delete_atoms overlap 0.5 all all</code>, or a minimization first). It does not help to switch the check off: <code>thermo_modify lost ignore</code> lets the run continue with 255 atoms and no message.</p></div>
    <div class="prose"><p><b>Probe the bad state</b>: the same input with <code>run 0</code>, a per-atom energy and a dump at step 0. Then check by hand that this one pair explains everything.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a2">▶ Probe step 0 and check by hand</button><span class="db-sta2"></span></div>
    <div class="prose db-a2sum"></div>
  </section>

  <section class="block">
    <h2>3 · A crash: one character in allocate()</h2>
    <div class="prose"><p>Atom types count from 1, so every per-type array has ntypes + 1 rows and row 0 is not used. Forget the +1:</p></div>
    <pre class="lmp-mini">${esc(BUG)}</pre>
    <div class="prose"><p><b>Release build.</b> Any input with lj/cut stops at <code>pair_coeff</code>: <code>Segmentation fault</code>, exit 139, and the log ends at the LAMMPS version line. gdb on the Release binary names the function but not the line:</p></div>
    <pre class="lmp-mini">#0  LAMMPS_NS::PairLJCut::coeff(int, char**) () from liblammps.so.0
#1  LAMMPS_NS::Input::pair_coeff() () from liblammps.so.0</pre>
    <div class="prose"><p><b>valgrind</b>, same binary, no rebuild. It names the first bad access, not the crash:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.vg)}</pre>
    <div class="prose"><p><b>AddressSanitizer</b>: a second build in its own directory, with checks compiled into every memory access:</p></div>
    <pre class="lmp-mini">${esc(ASANB)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.a1)}</pre>
    <div class="prose"><p>The line, the array (allocated at pair_lj_cut.cpp:409, which is <code>epsilon</code>), and the size and position of the bad access. With two atom types and <code>pair_coeff * *</code> the report is different, a WRITE into a 32-byte block allocated at memory.h:150 instead of a READ at 152:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.a2)}</pre>
    <div class="prose"><p>Both reports follow from two short pieces of code:</p></div>
    <pre class="lmp-mini">${esc(CREATE)}</pre>
    <pre class="lmp-mini">${esc(COEFF)}</pre>
    <div class="row"><label class="lc-in">atom types <select class="lc-u db-nt"><option>1</option><option selected>2</option><option>3</option><option>4</option></select></label><button class="btn sm primary" data-k="b">▶ Predict the first bad access</button></div>
    <div class="prose db-bsum"></div>
    <div class="prose"><p>Cost on examples/melt here: Release 0.52 s, Debug + ASan 1.69 s, valgrind 13.2 s. The sanitizer build is cheap enough to keep next to the normal one while you develop a style.</p></div>
  </section>

  <section class="block">
    <h2>4 · A NaN that nobody reports</h2>
    <div class="prose"><p>Place the extra atom exactly on top of atom 1:</p></div>
    <pre class="lmp-mini">${esc(SAME)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.th)}</pre>
    <div class="prose"><p>No error, no warning, and from step 1 a reasonable-looking potential energy. Only Temp and Press give it away. To find the first bad operation, make the CPU stop at it. Normally an invalid operation just returns nan and execution continues. With traps enabled it raises a signal at the instruction that did it:</p></div>
    <pre class="lmp-mini">${esc(FPEC)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.fpe)}</pre>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run it here, and follow the arithmetic</button><span class="db-stc"></span></div>
    <div class="prose db-csum"></div>
    <div class="prose"><p>A caution with traps: some legitimate code produces inf on purpose and tests for it afterwards, and some libraries do the same. A trap that stops in unfamiliar code is not automatically the bug; look at the values first.</p></div>
  </section>

  <section class="block">
    <h2>5 · Commands worth remembering</h2>
    <pre class="lmp-mini">${esc(CHEAT)}</pre>
    <div class="prose"><p>When a run worked last month and fails now, there are two more tools before any debugger. To bisect the input, delete half of the commands and see whether it still fails. To bisect the code, <code>git bisect</code> between the two versions: LAMMPS builds in a few minutes, and about log₂(number of commits) builds find the commit.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> In section 2, why are exactly 2 atoms lost, not 1 or 257?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only the overlapping pair feels the huge force, equal and opposite on atoms 1 and 257. Every other atom moves normally in the first step. Both atoms of the pair travel about 10¹⁰ box lengths; Domain::pbc adds one box length (domain.cpp:807), so they are still outside, and in serial CommBrick::exchange drops atoms that leave the box (comm_brick.cpp:721–725).</div></div>
    <div class="ex"><p><b>2.</b> The Release binary with the allocation bug crashed in coeff() for both 1 and 2 types. Would it always crash?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. Reading or writing just past a heap block is undefined behaviour: it hits whatever the allocator put there. A crash is the lucky outcome. It can also overwrite another array silently and give wrong forces, or corrupt the allocator and crash much later in unrelated code. That is why valgrind or ASan, which report the first bad access, beat reading a backtrace of the eventual crash.</div></div>
    <div class="ex"><p><b>3.</b> After step 1 in section 4, atoms 1 and 257 have nan coordinates, yet PotEng stays finite. Why do they not make it nan?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Every comparison with nan is false, so <code>rsq &lt; cutsq[itype][jtype]</code> fails for every pair that contains them. They simply drop out of the energy sum. The kinetic energy has no such test, so Temp becomes nan. A finite PotEng proves nothing; check Temp and Press, or every atom's coordinates.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    const a = Math.cbrt(4 / 0.8442), dt = 0.005;

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.db-sta'); b.disabled = true; st.textContent = ' Running…';
      const R = await runLammps(OV); if (!alive) return; b.disabled = false; st.textContent = '';
      const t = R.thermo.rows.map(r => r.join('   ')).join('\n');
      $('.db-asum').innerHTML = `<pre class="lmp-mini">${esc(R.thermo.cols.join('   '))}\n${esc(t)}\n${esc(R.error || '(no error)')}</pre><p>${/Lost atoms: original 257 current 255/.test(R.error || '') ? `${OK} The same error as the native build, after one line of thermo.` : `${BAD} unexpected result`}</p>`;
    };

    $('[data-k=a2]').onclick = async e => {
      const b = e.currentTarget, st = $('.db-sta2'); b.disabled = true; st.textContent = ' Running…';
      const R = await runLammps(PROBE, { read: ['a.dump'] }); if (!alive) return; b.disabled = false;
      if (R.error || !R.files['a.dump']) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no dump')}</span>`; return; } st.textContent = '';
      const L = R.files['a.dump'].split('\n'), n = +L[3], A = [];
      for (let k = 0; k < n; k++) { const w = L[9 + k].trim().split(/\s+/).map(Number); A.push({ id: w[0], x: w[1], y: w[2], z: w[3], fx: w[4], pe: w[5] }); }
      const top = [...A].sort((p, q) => q.pe - p.pe).slice(0, 3), p1 = A.find(p => p.id === 1), p2 = A.find(p => p.id === 257);
      const dx = p2.x - p1.x, r = Math.hypot(dx, p2.y - p1.y, p2.z - p1.z), r6 = r ** -6, E = 4 * (r6 * r6 - r6), F = 24 * (2 * r6 * r6 - r6) / r, mv = 0.5 * dt * dt * F, Lb = 4 * a;
      $('.db-a2sum').innerHTML = `<div class="tbl"><table>${head(['atom', 'pe/atom', 'fx'])}${top.map(p => row([p.id, p.pe.toExponential(6), p.fx.toExponential(6)])).join('')}</table></div>
      <p>Atoms 1 and 257 stand out by twelve orders of magnitude, then normal atoms near −7. By hand: r = ${s17(r)} σ (0.05 × the lattice constant ${a.toFixed(6)}), E<sub>pair</sub> = 4(r⁻¹² − r⁻⁶) = ${E.toExponential(9)}. pe/atom gives each atom half of every pair it is in, so pe₁ + pe₂₅₇ = E<sub>pair</sub> + (their normal neighbours): ${(p1.pe + p2.pe).toExponential(9)}, ratio ${((p1.pe + p2.pe) / E).toFixed(12)}. Force F = 24(2r⁻¹³ − r⁻⁷) = ${F.toExponential(9)}; LAMMPS's fx on atom 257: ${p2.fx.toExponential(9)} (ratio ${(p2.fx / F).toFixed(12)}).</p>
      <p>The first half step moves each atom by ½·dt²·F/m = ${mv.toExponential(3)} σ, ${(mv / Lb).toExponential(2)} box lengths. Domain::pbc wraps a coordinate by one box length (${gh('src/domain.cpp', 807)}), the atom is still outside, and in serial CommBrick::exchange discards atoms that leave the box (${gh('src/comm_brick.cpp', 721)}): two atoms gone, "original 257 current 255" ${OK}.</p>`;
    };

    $('[data-k=b]').onclick = () => {
      const nt = +$('.db-nt').value, n = nt;          // the bug: n = ntypes instead of ntypes + 1
      let out = null;
      for (let i = 1; i <= nt && !out; i++) {
        if (i >= n) { out = ['READ', 8 * i - 8 * n, 8 * n, 'memory.h:152 (the row pointers)', `reading the pointer epsilon[${i}]`]; break; }
        for (let j = i; j <= nt; j++) { const k = i * n + j; if (k >= n * n) { out = ['WRITE', 8 * k - 8 * n * n, 8 * n * n, 'memory.h:150 (the numbers)', `writing epsilon[${i}][${j}] = data[${i}·${n}+${j}] = data[${k}]`]; break; } }
      }
      const rec = nt === 1 ? NATIVE.a1 : nt === 2 ? NATIVE.a2 : null;
      let m = '';
      if (rec) {
        const sz = +rec.match(/(\d+)-byte region/)[1], off = +rec.match(/located (\d+) bytes after/)[1], op = rec.match(/(READ|WRITE) of size/)[1], at = rec.match(/memory\.h:(\d+)/)[1];
        const ok = sz === out[2] && off === out[1] && op === out[0] && out[3].startsWith('memory.h:' + at);
        m = ` Recorded natively: ${op}, ${off} bytes after ${sz === 8 ? 'an' : 'a'} ${sz}-byte region, memory.h:${at} ${ok ? OK : BAD}`;
      } else m = ' (not recorded natively; try it with the sanitizer build)';
      $('.db-bsum').innerHTML = `<p>${nt} type${nt > 1 ? 's' : ''}: the arrays get ${n} × ${n} double${n > 1 ? 's' : ''} (${8 * n * n} bytes) and ${n} row pointer${n > 1 ? 's' : ''} (${8 * n} bytes). Walking the coeff loop for <code>pair_coeff * *</code>, the first access outside is ${out[4]}: a ${out[0]} of size 8, ${out[1]} bytes after ${out[2] === 8 ? 'an' : 'a'} ${out[2]}-byte block from ${out[3]}.${m}</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.db-stc'); b.disabled = true; st.textContent = ' Running…';
      const R = await runLammps(SAME); if (!alive) return; b.disabled = false; st.textContent = '';
      const nat = NATIVE.th.split('\n').slice(1).map(l => l.trim().split(/\s+/).map(v => Number(v.replace('-nan', 'NaN'))));
      const br = R.thermo.rows.map(r => r.map(v => Number(String(v).replace(/-?nan/i, 'NaN'))));
      let same = 0, tot = 0; nat.forEach((r, i) => r.forEach((v, c) => { tot++; const w = br[i] && br[i][c]; if (Object.is(v, w) || (Number.isNaN(v) && Number.isNaN(w)) || v === w) same++; }));
      // the IEEE chain in PairLJCut::compute for rsq = 0, delx = 0 (lj1 = 48, lj2 = 24, lj3 = 4, lj4 = 4)
      const rsq = 0, delx = 0, r2inv = 1.0 / rsq, r6inv = r2inv * r2inv * r2inv, forcelj = r6inv * (48 * r6inv - 24), fpair = forcelj * r2inv, fx = delx * fpair, ev = r6inv * (4 * r6inv - 4), vir = delx * delx * fpair;
      $('.db-csum').innerHTML = `<pre class="lmp-mini">${esc(R.thermo.cols.join('   '))}\n${esc(R.thermo.rows.map(r => r.join('   ')).join('\n'))}</pre>
      <p>Your browser: ${same} of ${tot} thermo values the same as native (nan counted equal to nan) ${same === tot ? OK : BAD}, and again no error. Now the same arithmetic as pair_lj_cut.cpp:113–131 with rsq = 0, in JavaScript's IEEE doubles: r2inv = 1/0 = ${r2inv}; r6inv = ${r6inv}; forcelj = ∞·(48·∞ − 24) = ${forcelj}; fpair = ${fpair}; f<sub>x</sub> += delx·fpair = 0·∞ = ${fx}; evdwl = ∞·(4·∞ − 4) = ${ev}; virial delx·delx·fpair = ${vir}.</p>
      <p>So PotEng = inf and Press = nan at step 0, exactly the native line. The forces on atoms 1 and 257 are nan, so at step 1 their velocities and positions are nan, and Temp is nan. Their pairs then fail <code>rsq &lt; cutsq</code> (nan &lt; 6.25 is ${NaN < 6.25}) and leave the energy sum. The trap stopped at the first of these operations, 1.0/rsq on line 114 with rsq = 0 and tags 1 and 257: the cause, not a symptom.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
