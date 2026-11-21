// Research · B.8: finding your way in the code. Four routes from what you see to the line that made it; Temp traced to compute_temp.cpp and recomputed from 4000 velocities.
import { runLammps } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const g15 = x => String(Number(x.toPrecision(15)));

const TRACE = [
  { what: 'You see the column Temp in the thermo output', where: ['src/thermo.cpp', 1430], code: `1430  } else if (word == "temp") {
1431    check_temp(word);
1432    compute_temp();` },
  { what: 'compute_temp() just copies a number from a compute', where: ['src/thermo.cpp', 1933], code: `1933  void Thermo::compute_temp()
1935    dvalue = temperature->scalar;` },
  { what: 'Which compute? One LAMMPS creates for you at start-up, called thermo_temp', where: ['src/output.cpp', 63], code: `63  modify->add_compute("thermo_temp all temp");
64  modify->add_compute("thermo_press all pressure thermo_temp");` },
  { what: '"temp" is a compute style, so by the naming rule (page B.1) its file is compute_temp.cpp. The sum:', where: ['src/compute_temp.cpp', 72], code: `72  double ComputeTemp::compute_scalar()
      …
92      t += (v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]) * mass[type[i]];
      …
99    scalar *= tfactor;` },
  { what: 'And tfactor: the degrees of freedom', where: ['src/compute_temp.cpp', 58], code: `58  void ComputeTemp::dof_compute()
62    dof = domain->dimension * natoms_temp;      // 3N
63    dof -= extra_dof + fix_dof;                 // minus 3, minus constraints
65      tfactor = force->mvv2e / (dof * force->boltz);` },
  { what: 'extra_dof is 3 unless you change it (compute.cpp)', where: ['src/compute.cpp', 85], code: `85  extra_dof = domain->dimension;     // the 3 momenta of the centre of mass` },
];

const SCRIPT = extra => IN_MELT.replace('thermo          50', 'thermo          50\nthermo_modify   format float %.15g')
  .replace('run             250', `${extra ? 'compute_modify  thermo_temp extra/dof 0\n' : ''}dump            1 all custom 100 v.dump id vx vy vz
dump_modify     1 format float %.17g
run             100`);

export default {
  id: 'finding-your-way',
  title: 'Finding your way in the code',
  minutes: 20,
  summary: 'LAMMPS is over a million lines of C++, but you never read it from the top. Four routes take you from what you see (an error message, a command, a thermo column) straight to the line that made it. You follow one of them from the Temp column to the eight lines that compute it, then recompute Temp yourself from the 4000 velocities LAMMPS dumps, and change the one setting that makes the naive formula right.',
  learn: [
    'Four routes into the code: the error location, the naming rule, the thermo keyword, grep',
    'Trace the Temp column to compute_temp.cpp',
    'Recompute T = Σ m v² / (dof · k_B) from 4000 velocities and match LAMMPS to 15 digits',
    'Why dof is 3N − 3, and what changes when you set it to 3N',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: four routes in</h2>
    <div class="tbl wide"><table>${head(['You have', 'Route', 'Example'])}
${row(['an error or a warning', 'read the end of the message: every ERROR and WARNING carries the file and line that raised it (' + gh('src/error.cpp', 139) + ', ' + gh('src/error.cpp', 228) + ')', '<code>ERROR: Illegal fix nvt command: missing argument(s) (src/fix_nh.cpp:61)</code>'])}
${row(['a command in the script', 'the naming rule of page B.1', '<code>fix 1 all nve</code> → <code>src/fix_nve.cpp</code>'])}
${row(['a column in the thermo output', 'find the keyword in <code>src/thermo.cpp</code>, follow it to the compute it reads', 'Temp → <code>compute_temp.cpp</code> (section 2)'])}
${row(['a word and nothing else', '<code>grep</code> (page 0.2) or <code>git grep</code>', '<code>grep -rn "Lost atoms" src/</code> → <code>src/thermo.cpp:526</code>'])}</table></div>
    <p>An editor with "go to definition" (VS Code with the C/C++ extension, opened on the <code>src</code> folder) then jumps from any name to where it is defined. You only ever read the few lines on the path.</p>
  </section>

  <section class="block">
    <h2>2 · Observe: from the Temp column to the formula</h2>
    <div class="stepper st3"><div class="row"><button class="btn sm ghost" data-t="-1">◀</button><span class="st3-n mono"></span><button class="btn sm ghost" data-t="1">▶</button></div>
      <p class="st3-when"></p><pre class="lmp-mini st3-code"></pre></div>
    <div class="prose"><p>Six stops, about 20 lines read. In LJ units <code>mvv2e</code> = 1 and <code>boltz</code> = 1 (page A.4), so the code says:</p>
<pre class="lmp-mini">T = Σ m (vx² + vy² + vz²) / dof,     dof = 3N − 3</pre>
    <p>The −3: <code>velocity … mom yes</code> (page A.6) and <code>fix nve</code> keep the total momentum at zero, so 3 of the 3N velocity components are not free. Most textbooks write 3N.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: Temp from 4000 velocities</h2>
    <div class="prose"><p>Run melt for 100 steps, dump every velocity with 17 digits, and apply the formula you just found. Thermo is printed with 15 digits so the two can be compared.</p></div>
    <label class="rd-sp-lab"><input type="checkbox" class="fw-x"> add <code>compute_modify thermo_temp extra/dof 0</code> (section 4)</label>
    <div class="row"><button class="btn sm primary" data-k="run">▶ Run melt and read the dump</button><span class="lmp-status muted fw-st"></span></div>
    <div class="calc fw-hand" hidden></div>
    <div class="tbl fw-cmp"></div>
  </section>

  <section class="block prose">
    <h2>4 · Modify: change one number in the code's path</h2>
    <p><code>extra_dof</code> is set in ${gh('src/compute.cpp', 85)}, but you do not need to recompile to change it: <code>compute_modify thermo_temp extra/dof 0</code> sets it from the script (${gh('src/compute.cpp', 149)}). Tick the box above and run again. Predict first:</p>
    <ul>
      <li>Temp becomes Σ m v² / 12000 instead of / 11997: smaller by the factor 11997/12000.</li>
      <li>Press does <b>not</b> change. ${gh('src/compute_pressure.cpp', 256)} uses <code>dof · boltz · T</code>, which is Σ m v² whatever dof is.</li>
      <li>TotEng does not change: the kinetic energy is ½ Σ m v², with no dof in it.</li>
    </ul>
    <p>Reading the code told you which numbers depend on the setting and which do not, before running anything.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Where does LAMMPS decide to stop with "Lost atoms"? Use route 4.</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>grep -rn "Lost atoms" src/</code> gives ${gh('src/thermo.cpp', 526)}, inside <code>Thermo::lost_check()</code> (line 495): it adds up the atoms on all processors every time thermo prints and compares with the number at the start. That is also why "Lost atoms" appears only at a thermo step.</div></div>
    <div class="ex"><p><b>2.</b> Follow the Press column the same way. Which file, and what is added to the kinetic part?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>thermo.cpp → <code>thermo_press</code> (output.cpp:64) → ${gh('src/compute_pressure.cpp', 256)}: P = (dof·k_B·T + virial_xx + virial_yy + virial_zz) / (3V) · nktv2p. The virial is Σ r·F over pairs, summed by the pair style (page C.5).</div></div>
    <div class="ex"><p><b>3.</b> For 2 atoms, what Temp does LAMMPS print if each moves with speed 1 in opposite directions (m = 1, LJ units)?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Σ m v² = 2, dof = 3·2 − 3 = 3, so T = <b>0.667</b>. With 3N it would be 0.333. For small systems the −3 matters a lot.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let t = 0;
    const showT = () => { const T = TRACE[t]; $('.st3-n').textContent = `${t + 1} / ${TRACE.length}`; $('.st3-when').innerHTML = `${esc(T.what)} · ${gh(...T.where)}`; $('.st3-code').textContent = T.code; };
    root.querySelectorAll('[data-t]').forEach(b => b.onclick = () => { t = (t + +b.dataset.t + TRACE.length) % TRACE.length; showT(); }); showT();
    $('[data-k=run]').onclick = async e => {
      const b = e.currentTarget, st = $('.fw-st'), extra = $('.fw-x').checked; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(SCRIPT(extra), { read: ['v.dump'] }); b.disabled = false;
      const txt = res.files && res.files['v.dump'];
      if (res.error || !txt) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump')}</span>`; return; }
      st.textContent = '';
      const fr = txt.split('ITEM: TIMESTEP').slice(1).map(f => { const ls = f.trim().split('\n'); return { step: ls[0].trim(), rows: ls.slice(8).map(l => l.trim().split(/\s+/).map(Number)) }; });
      const f = fr.find(x => x.step === '100'); let s = 0; for (const r of f.rows) s += (r[1] * r[1] + r[2] * r[2] + r[3] * r[3]) * 1.0;
      const N = f.rows.length, dof = extra ? 3 * N : 3 * N - 3, Th = s / dof, Tnaive = s / (3 * N), KE = 0.5 * s / N;
      const tl = res.lines.findIndex(l => /^\s+Step\s+Temp/.test(l)), cols = res.lines[tl].trim().split(/\s+/);
      const r100 = res.lines.slice(tl + 1).map(l => l.trim().split(/\s+/)).find(w => w[0] === '100'), L = k => r100[cols.indexOf(k)];
      const keL = +L('TotEng') - +L('E_pair');
      $('.fw-hand').hidden = false;
      $('.fw-hand').innerHTML = esc(`Σ m v² over ${N} atoms (m = 1)      = ${g15(s)}
dof = ${extra ? '3N' : '3N − 3'} = ${dof}
T = Σ m v² / dof                   = `) + `<b>${g15(Th)}</b>` + esc(extra ? '' : `
(the textbook 3N would give ${g15(Tnaive)})`) + esc(`
KE per atom = ½ Σ m v² / N         = ${g15(KE)}`);
      const ok = (h, l) => Math.abs(h - l) <= 2e-14 * Math.abs(l) ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>';
      $('.fw-cmp').innerHTML = `<table>${head(['step 100', 'by hand', 'LAMMPS', ''])}
        ${row(['Temp', `<span class="mono">${g15(Th)}</span>`, `<span class="mono">${L('Temp')}</span>`, ok(Th, +L('Temp'))])}
        ${row(['TotEng − E_pair (KE per atom)', `<span class="mono">${g15(KE)}</span>`, `<span class="mono">${g15(keL)}</span>`, ok(KE, keL)])}
        ${row(['Press', 'unchanged by dof', `<span class="mono">${L('Press')}</span>`, ''])}</table>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
