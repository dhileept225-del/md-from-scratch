// Research · A.1: what an input script is. Six rearranged versions of melt, predicted, run live and explained.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const PAIR = 'pair_style      lj/cut 2.5\npair_coeff      1 1 1.0 1.0 2.5\n';
const VEL = 'velocity        all create 3.0 87287 loop geom\n';

const VARIANTS = [
  { title: 'The force field before the box', change: 'the two pair lines moved above <code>lattice</code>',
    edit: s => s.replace(PAIR, '').replace('lattice         fcc 0.8442', PAIR + 'lattice         fcc 0.8442'), answer: 'error',
    why: '<code>pair_coeff 1 1</code> sets the force between atom types 1 and 1, but atom types do not exist until <code>create_box 1 box</code> says there is one type. <code>pair_style</code> alone is allowed early; <code>pair_coeff</code> is not.' },
  { title: 'No mass', change: 'the line <code>mass 1 1.0</code> deleted',
    edit: s => s.replace('mass            1 1.0\n', ''), answer: 'error',
    why: 'The error comes from <code>velocity</code>: to scale velocities to T = 3 it needs the kinetic energy ½mv², and m is unknown. LAMMPS names the file and line where it stopped (velocity.cpp).' },
  { title: 'Velocities before atoms', change: '<code>velocity</code> moved above <code>create_atoms</code>',
    edit: s => s.replace(VEL, '').replace('create_atoms    1 box', VEL + 'create_atoms    1 box'), answer: 'error',
    why: '<code>velocity</code> is an action: it gives velocities to the atoms that exist at that moment, and there are none yet.' },
  { title: 'Printing set after the run', change: '<code>thermo 50</code> moved below <code>run 250</code>',
    edit: s => s.replace('thermo          50\nrun             250', 'run             250\nthermo          50'), answer: 'different',
    why: '<code>thermo 50</code> is a setting for runs that come after it. During the run the default applies (print only the first and last step), so only steps 0 and 250 appear. The numbers on those two lines are identical to the original: the simulation itself did not change.',
    hand: { 0: { E_pair: -6.7733681, TotEng: -2.2744931, Press: -3.7033504 } } },
  { title: 'Moving atoms switched on too late', change: '<code>fix 1 all nve</code> moved below <code>run 250</code>',
    edit: s => s.replace('fix             1 all nve\n', '').replace('run             250', 'run             250\nfix             1 all nve'), answer: 'different',
    why: 'Nothing integrates the equations of motion during the run, so the atoms never move: every line repeats step 0. LAMMPS warns but does not stop, because a run without motion can be useful (for example to compute energies of a fixed structure).',
    hand: { 0: { E_pair: -6.7733681, TotEng: -2.2744931, Press: -3.7033504 }, 250: { E_pair: -6.7733681, TotEng: -2.2744931, Press: -3.7033504 } } },
  { title: 'No forces at all', change: 'both pair lines deleted',
    edit: s => s.replace(PAIR, ''), answer: 'different',
    why: 'LAMMPS happily simulates atoms that do not interact: an ideal gas. With no forces, E_pair = 0, velocities never change (Temp stays 3), and only the kinetic parts remain: TotEng = KE/N = 1.5 × 3 × 3999/4000 = 4.498875 and Press = ρT(N−1)/N = 0.8442 × 3 × 0.99975 = 2.5319668. Two warnings, but no error: LAMMPS cannot know you forgot the force field.',
    hand: { 0: { E_pair: 0, TotEng: 4.498875, Press: 2.5319668 }, 250: { E_pair: 0, TotEng: 4.498875, Press: 2.5319668 } } },
];

export default {
  id: 'input-script',
  title: 'What an input script is',
  minutes: 20,
  summary: 'An input script is a list of commands that LAMMPS carries out one at a time, top to bottom, the moment it reads each one. You take the melt example apart six ways, predict for each whether LAMMPS stops, runs normally or runs differently, then run it live and check the numbers by hand.',
  learn: [
    'The three rules of an input script: one command at a time, top to bottom, immediately',
    'Settings (remembered for later) versus actions (done now), and why order matters',
    'Which mistakes LAMMPS catches with an error, and which it cannot catch at all',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: three rules</h2>
    <ol>
      <li><b>One command per line.</b> The first word names the command; the rest are its settings (page A.2 shows how a line is read).</li>
      <li><b>Top to bottom.</b> LAMMPS never looks ahead. A command can only use what earlier lines created.</li>
      <li><b>Immediately.</b> Each command is carried out as soon as it is read, before the next line is even looked at.</li>
    </ol>
    <p>Commands come in two kinds. A <b>setting</b> is remembered and used later; an <b>action</b> does something now, with whatever exists at that moment.</p>
    <div class="tbl wide"><table>${head(['melt line', 'kind', 'what it needs before it'])}
${row(['<code>units lj</code>, <code>atom_style atomic</code>', 'setting', 'must come first: everything else depends on them'])}
${row(['<code>lattice fcc 0.8442</code>', 'setting', ''])}
${row(['<code>region box block …</code>', 'setting (a named shape)', 'the lattice (its numbers are in lattice cells)'])}
${row(['<code>create_box 1 box</code>', '<b>action</b>: makes the box and 1 atom type', 'the region'])}
${row(['<code>create_atoms 1 box</code>', '<b>action</b>: puts 4000 atoms in', 'the box and the lattice'])}
${row(['<code>mass 1 1.0</code>', 'setting', 'the atom type'])}
${row(['<code>velocity all create 3.0 …</code>', '<b>action</b>: gives velocities now', 'atoms and their mass'])}
${row(['<code>pair_style</code>, <code>pair_coeff</code>', 'setting', 'pair_coeff needs the atom types'])}
${row(['<code>neighbor</code>, <code>neigh_modify</code>, <code>fix</code>, <code>thermo</code>', 'setting', ''])}
${row(['<code>run 250</code>', '<b>action</b>: uses every setting made so far', 'all of the above'])}</table></div>
  </section>

  <section class="block">
    <h2>2 · Modify, predict, observe: six rearranged scripts</h2>
    <p class="muted">For each change, choose what you think LAMMPS will do, then run the real LAMMPS to see.</p>
    <div class="variants">${VARIANTS.map((v, i) => `
      <div class="exp var" data-v="${i}">
        <h3>${String.fromCharCode(65 + i)}. ${v.title}</h3>
        <p>Change: ${v.change}.</p>
        <div class="row choose">${[['ok', 'runs exactly as before'], ['different', 'runs, but the output differs'], ['error', 'stops with an error']].map(([k, t]) => `<button class="btn sm ghost" data-p="${k}">${t}</button>`).join('')}</div>
        <div class="row"><button class="btn sm primary" data-run>▶ Run it</button><span class="lmp-status muted"></span></div>
        <div class="var-out" hidden><pre class="lmp-mini"></pre><div class="verdict"></div><p class="why-txt"></p><div class="lmp-compare"></div></div>
      </div>`).join('')}
    </div>
  </section>

  <section class="block prose">
    <h2>3 · What the six results teach</h2>
    <div class="tbl wide"><table>${head(['Change', 'LAMMPS', 'Lesson'])}
${row(['A, B, C', 'stops with an error', 'an action or a setting used something that did not exist yet'])}
${row(['D', 'runs; fewer lines printed', 'a setting placed after an action does not affect it'])}
${row(['E', 'runs; a warning; nothing moves', 'a missing step is not always an error'])}
${row(['F', 'runs; two warnings; an ideal gas', 'LAMMPS checks the order, never your intention: only a hand calculation catches this'])}</table></div>
    <p>The last row is the reason for the hand calculations on every page of this section. An input that runs is not the same as an input that is right.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You add <code>timestep 0.001</code> after <code>run 250</code>. What happens to the run?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Nothing: the run uses the time step that was set when it started (the default 0.005 in LJ units). The new value applies only to later runs. Same rule as thermo in change D.</div></div>
    <div class="ex"><p><b>2.</b> In change F, why is the pressure exactly ρT(N−1)/N?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Pressure is a kinetic part plus a force (virial) part (Start page). With no forces the virial part is zero, leaving P = N k<sub>B</sub>T/V with the 3 removed degrees of freedom: (N−1)/N × ρT = 0.99975 × 0.8442 × 3 = 2.5319668, the ideal gas law.</div></div>
    <div class="ex"><p><b>3.</b> Can you put <code>pair_style lj/cut 2.5</code> as the very first line?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>In this version, yes: pair_style only chooses the kind of force, which needs no atoms. pair_coeff cannot go before create_box (change A). When in doubt, follow the order of the four blocks on page A.3.</div></div>
  </section>`,

  mount(root) {
    root.querySelectorAll('.var').forEach(box => {
      const v = VARIANTS[+box.dataset.v]; let pick = null;
      box.querySelectorAll('[data-p]').forEach(b => b.onclick = () => { pick = b.dataset.p; box.querySelectorAll('[data-p]').forEach(x => x.classList.toggle('on', x === b)); });
      box.querySelector('[data-run]').onclick = async e => {
        const btn = e.currentTarget, st = box.querySelector('.lmp-status'); btn.disabled = true;
        st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
        const res = await runLammps(v.edit(IN_MELT)); btn.disabled = false; st.textContent = '';
        const out = box.querySelector('.var-out'); out.hidden = false;
        const shown = res.lines.filter(l => /^(ERROR|WARNING)|^\s+Step|^\s+\d+\s+-?[\d.]/.test(l));
        box.querySelector('.lmp-mini').textContent = shown.join('\n') || res.lines.slice(-6).join('\n');
        const actual = res.error ? 'error' : v.answer === 'different' ? 'different' : 'ok';
        const names = { ok: 'runs exactly as before', different: 'runs, but the output differs', error: 'stops with an error' };
        box.querySelector('.verdict').innerHTML = `<p><b>LAMMPS ${names[actual]}.</b> ${pick ? (pick === actual ? '<span class="m-ok">Your prediction was right.</span>' : `<span class="m-bad">You predicted: ${names[pick]}.</span>`) : ''}</p>`;
        box.querySelector('.why-txt').innerHTML = v.why;
        if (v.hand && !res.error) {
          const { cols, rows } = res.thermo, rowsHtml = [];
          for (const [step, h] of Object.entries(v.hand)) {
            const r = rows.find(w => w[0] === step); if (!r) continue;
            for (const k of ['E_pair', 'TotEng', 'Press']) { const lv = r[cols.indexOf(k)], m = matchMarkup(g8(h[k]), lv); rowsHtml.push(`<tr><td>step ${step} ${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`); }
          }
          box.querySelector('.lmp-compare').innerHTML = `<div class="tbl"><table><tr><th></th><th>by hand</th><th>LAMMPS</th><th></th></tr>${rowsHtml.join('')}</table></div>`;
        }
      };
    });
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
