// Research · A.3: the four blocks of every input script; sort the melt lines, then build a new system from the same skeleton.
import { runLammps, g8, matchMarkup } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const BLOCKS = [['init', '1 · Initialization'], ['system', '2 · System'], ['settings', '3 · Settings'], ['run', '4 · Run']];
const LINES = IN_MELT.split('\n').filter(l => l.trim() && !l.startsWith('#'));
const ANSWER = l => /^(units|atom_style|dimension|boundary)\b/.test(l) ? 'init' : /^(lattice|region|create_box|create_atoms|read_data)\b/.test(l) ? 'system' : /^run\b/.test(l) ? 'run' : 'settings';

const TWO_D = `# 1 · Initialization
units           lj
dimension       2
atom_style      atomic

# 2 · System: a 3 x 3 square crystal, spacing 1.2
lattice         sq $(1/1.2^2)
region          box block 0 3 0 3 -0.5 0.5
create_box      1 box
create_atoms    1 box

# 3 · Settings
mass            1 1.0
velocity        all create 0.72 87287 loop geom
pair_style      lj/cut 1.8
pair_coeff      1 1 1.0 1.0 1.8
fix             1 all nve
fix             2 all enforce2d

# 4 · Run
run             0
`;

export default {
  id: 'four-blocks',
  title: 'The four blocks of every script',
  minutes: 20,
  summary: 'Almost every LAMMPS input has the same skeleton: initialization, system, settings, run. You sort the melt lines into the four blocks, then reuse the skeleton for a different system, a two-dimensional crystal of 9 atoms, predict its first line by hand and check it.',
  learn: [
    'The four blocks and which commands belong in each',
    'Sort every line of the melt example into its block',
    'Build a new 2D system from the same skeleton and match its first line by hand',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: one skeleton for every script</h2>
    <div class="tbl wide"><table>${head(['Block', 'Answers the question', 'Typical commands'])}
${row(['1 · Initialization', 'what kind of world? units, 2D or 3D, what an atom carries', '<code>units</code>, <code>dimension</code>, <code>boundary</code>, <code>atom_style</code>'])}
${row(['2 · System', 'where are the atoms?', '<code>lattice</code>, <code>region</code>, <code>create_box</code>, <code>create_atoms</code>, or <code>read_data</code> (page A.8)'])}
${row(['3 · Settings', 'how do they interact and move, and what is printed?', '<code>mass</code>, <code>velocity</code>, <code>pair_style</code>, <code>pair_coeff</code>, <code>neighbor</code>, <code>fix</code>, <code>compute</code>, <code>thermo</code>'])}
${row(['4 · Run', 'go!', '<code>run</code>, <code>minimize</code>'])}</table></div>
    <p>The order follows the rules of page A.1: block 1 changes the meaning of everything after it, so it comes first; block 3 needs the atom types created in block 2; block 4 uses everything. <code>velocity</code> is strictly an action, but it is usually written with the settings, just before the run.</p>
  </section>

  <section class="block">
    <h2>2 · Sort the melt lines</h2>
    <p class="muted">Choose a block for each line, then check.</p>
    <div class="sorter">${LINES.map((l, i) => `<div class="srow" data-i="${i}"><code>${esc(l.replace(/\s+/g, " "))}</code><div class="seg mini">${BLOCKS.map(([k, t]) => `<button data-b="${k}" title="${t}">${t.split(' · ')[0]}</button>`).join('')}</div></div>`).join('')}</div>
    <div class="row"><button class="btn sm primary" data-k="sort">Check</button><button class="btn sm ghost" data-k="show">Show the answer</button><span class="sort-res muted"></span></div>
    <p class="muted small">1 = Initialization, 2 = System, 3 = Settings, 4 = Run.</p>
  </section>

  <section class="block prose">
    <h2>3 · Modify: a new system from the same skeleton</h2>
    <p>Keep the four blocks and change their contents: a flat (2D) crystal of 3 × 3 atoms on a square grid with spacing 1.2, at temperature 0.72, cutoff 1.8. These are the 9 atoms of Lab 3.9.</p>
    <div class="calc">${esc(TWO_D)}</div>
    <p>New commands: <code>dimension 2</code> (atoms move in a plane), <code>lattice sq</code> (a square grid; in LJ units the number is the density 1/1.2² = 0.69444, so the spacing is 1.2) and <code>fix enforce2d</code> (keeps the motion exactly flat).</p>
    <h3>Predict by hand</h3>
    <p>Every atom has 4 neighbours at 1.2 and 4 on the diagonals at 1.2√2 = 1.6971; the next ones, at 2.4, are beyond the cutoff 1.8.</p>
    <div class="calc">u(r) = 4(r⁻¹² − r⁻⁶)
u(1.2)    = 4(0.11215665 − 0.33489798) = −0.89096529
u(1.6971) = 4(0.00175245 − 0.04186225) = −0.16043920
E_pair = ½ × (4 × u(1.2) + 4 × u(1.6971))   (each pair shared by 2)
       = 2 × (−0.89096529 − 0.16043920)     = <b>−2.102809</b>

in 2D an atom has 2 degrees of freedom; 2 are removed for the whole crystal:
KE/N   = (2N − 2)/2 × T / N = 8 × 0.72 / 9 = 0.64
TotEng = −2.102809 + 0.64                   = <b>−1.462809</b>

area A = 3.6 × 3.6 = 12.96;  pressure in 2D = (kinetic + virial)/(2A)
kinetic: (2N − 2) × T /(2A) = 16 × 0.72 / 25.92          = 0.44444444
virial : N × ½ × 4 × (rF(1.2) + rF(1.6971)) /(2A)
       = 9 × 2 × (−2.65403201 − 0.92057644) / 25.92     = −2.48236700
Press  = 0.44444444 − 2.48236700            = <b>−2.0379225</b>
       (rF = 24(2r⁻¹² − r⁻⁶) = r times the force)</div>
    <div class="row"><button class="btn sm primary" data-k="run2d">▶ Run the 2D script</button><span class="lmp-status muted st2d"></span></div>
    <div class="lmp-compare cmp2d"></div>
  </section>

  <section class="block prose">
    <h2>What the new system shows</h2>
    <ul>
      <li>The skeleton did not change; only the contents of each block did. Every script in this course, and most in research, has this shape.</li>
      <li>A change in block 1 (<code>dimension 2</code>) changed rules everywhere: 2 degrees of freedom per atom instead of 3, an area instead of a volume.</li>
      <li>The hand calculation is the same method as for melt: shells of neighbours, half of each pair per atom.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> In which block do <code>compute</code> and <code>dump</code> belong?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Settings: they say what to calculate and save during the run. Like thermo, they must come before the run they are meant for.</div></div>
    <div class="ex"><p><b>2.</b> Change the 2D cutoff to 1.5 (both lines). Which neighbours remain, and what is E_pair?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Only the 4 at 1.2 (the diagonals at 1.6971 are now outside). E_pair = ½ × 4 × (−0.89096529) = <b>−1.7819306</b>. Edit the 2D script in the Start page's runner or on your own computer to check.</div></div>
    <div class="ex"><p><b>3.</b> Why does the 2D kinetic energy per atom use (2N − 2)/2 and not 3N/2?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>In 2D each atom moves in x and y only: 2 degrees of freedom, each with ½T on average. LAMMPS removes the motion of the whole crystal (2 more), leaving 2N − 2 = 16 for 9 atoms.</div></div>
  </section>`,

  mount(root) {
    const rows = [...root.querySelectorAll('.srow')], picks = rows.map(() => null);
    rows.forEach((r, i) => r.querySelectorAll('[data-b]').forEach(b => b.onclick = () => { picks[i] = b.dataset.b; r.querySelectorAll('[data-b]').forEach(x => x.classList.toggle('on', x === b)); r.classList.remove('ok', 'bad'); }));
    root.querySelector('[data-k=sort]').onclick = () => {
      let right = 0; rows.forEach((r, i) => { const ok = picks[i] === ANSWER(LINES[i]); if (ok) right++; r.classList.toggle('ok', ok); r.classList.toggle('bad', !!picks[i] && !ok); });
      root.querySelector('.sort-res').textContent = `${right} of ${rows.length} right${right === rows.length ? ': every line in its block.' : '.'}`;
    };
    root.querySelector('[data-k=show]').onclick = () => { rows.forEach((r, i) => { picks[i] = ANSWER(LINES[i]); r.querySelectorAll('[data-b]').forEach(x => x.classList.toggle('on', x.dataset.b === picks[i])); }); root.querySelector('[data-k=sort]').click(); };
    root.querySelector('[data-k=run2d]').onclick = async e => {
      const b = e.currentTarget, st = root.querySelector('.st2d'); b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(TWO_D); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      const created = res.lines.filter(l => /Created|Lattice spacing/.test(l)).map(l => l.trim()).join(' · ');
      st.textContent = created;
      const r = res.thermo.rows.find(w => w[0] === '0'), c = res.thermo.cols;
      const H = { E_pair: -2.102809, TotEng: -1.462809, Press: -2.0379225 };
      root.querySelector('.cmp2d').innerHTML = `<div class="tbl"><table><tr><th>step 0, 2D</th><th>by hand</th><th>LAMMPS</th><th></th></tr>${Object.entries(H).map(([k, h]) => { const lv = r[c.indexOf(k)], m = matchMarkup(g8(h), lv); return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${lv}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('')}</table></div>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
