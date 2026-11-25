// Research · C.3: hooks and fix masks. setmask() bits, the per-hook lists built in fix order, and two live experiments where the order of fixes changes the numbers.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;

// bit, name, where in the loop (page C.2 stage)
const HOOKS = [
  [0, 'INITIAL_INTEGRATE', '2 · first half'],
  [1, 'POST_INTEGRATE', '2 · right after'],
  [2, 'PRE_EXCHANGE', '3 · rebuild steps only'],
  [3, 'PRE_NEIGHBOR', '3 · rebuild steps only'],
  [4, 'POST_NEIGHBOR', '3 · rebuild steps only'],
  [5, 'PRE_FORCE', '4 · before forces'],
  [6, 'PRE_REVERSE', '4 · before reverse_comm'],
  [7, 'POST_FORCE', '5 · after forces'],
  [8, 'FINAL_INTEGRATE', '5 · second half'],
  [9, 'END_OF_STEP', '5 · end, every N steps'],
  [10, 'POST_RUN', 'after the whole run'],
];
const H = Object.fromEntries(HOOKS.map(([b, n]) => [n, 1 << b]));
const R = { INITIAL_INTEGRATE_RESPA: 1 << 11, PRE_FORCE_RESPA: 1 << 13, POST_FORCE_RESPA: 1 << 14, FINAL_INTEGRATE_RESPA: 1 << 15, MIN_PRE_REVERSE: 1 << 20, MIN_POST_FORCE: 1 << 21 };
const FIXES = {
  'nve': { f: 'src/fix_nve.cpp', l: 44, h: ['INITIAL_INTEGRATE', 'FINAL_INTEGRATE'], x: ['INITIAL_INTEGRATE_RESPA', 'FINAL_INTEGRATE_RESPA'] },
  'nvt': { f: 'src/fix_nh.cpp', l: 614, h: ['INITIAL_INTEGRATE', 'FINAL_INTEGRATE'], x: ['INITIAL_INTEGRATE_RESPA', 'FINAL_INTEGRATE_RESPA', 'PRE_FORCE_RESPA'] },
  'langevin': { f: 'src/fix_langevin.cpp', l: 217, h: ['POST_FORCE'], x: ['POST_FORCE_RESPA'] },
  'temp/rescale': { f: 'src/fix_temp_rescale.cpp', l: 97, h: ['END_OF_STEP'], x: [] },
  'setforce': { f: 'src/fix_setforce.cpp', l: 115, h: ['POST_FORCE'], x: ['POST_FORCE_RESPA', 'MIN_POST_FORCE'] },
  'addforce': { f: 'src/fix_addforce.cpp', l: 129, h: ['POST_FORCE'], x: ['POST_FORCE_RESPA', 'MIN_POST_FORCE'] },
  'spring/self': { f: 'src/fix_spring_self.cpp', l: 124, h: ['POST_FORCE'], x: ['POST_FORCE_RESPA', 'MIN_POST_FORCE'] },
  'shake': { f: 'src/RIGID/fix_shake.cpp', l: 353, h: ['PRE_NEIGHBOR', 'POST_FORCE', 'POST_RUN'], x: ['POST_FORCE_RESPA', 'MIN_PRE_REVERSE', 'MIN_POST_FORCE'] },
  'recenter': { f: 'src/fix_recenter.cpp', l: 109, h: ['INITIAL_INTEGRATE'], x: ['INITIAL_INTEGRATE_RESPA'] },
  'momentum': { f: 'src/fix_momentum.cpp', l: 74, h: ['END_OF_STEP'], x: [] },
  'ave/time': { f: 'src/fix_ave_time.cpp', l: 467, h: ['END_OF_STEP'], x: [] },
  'print': { f: 'src/fix_print.cpp', l: 110, h: ['END_OF_STEP'], x: [] },
};
const mask = k => FIXES[k].h.reduce((m, n) => m | H[n], 0) | FIXES[k].x.reduce((m, n) => m | R[n], 0);
const rowsOf = res => res.thermo.rows.map(r => Object.fromEntries(res.thermo.cols.map((c, i) => [c, r[i]])));
const DEF = ['nve', 'langevin', 'setforce', 'ave/time'];

const BASE = `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 3.0 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no
fix             1 all nve
`;
const ADD = 'fix             a all addforce 1.0 0.0 0.0', SET = 'fix             s all setforce 0.0 NULL NULL';
const ORDERS = { as: [ADD, SET], sa: [SET, ADD], a: [ADD], s: [SET] };
const PRED = { as: [0, 0, 0], sa: [1, 1, 4000], a: ['?', '?', 4000], s: [0, 0, 0] };
const scriptA = o => `${BASE}${ORDERS[o].join('\n')}
compute         mn all reduce min fx
compute         mx all reduce max fx
compute         sm all reduce sum fx
thermo_style    custom step c_mn c_mx c_sm
thermo_modify   norm no
thermo          1
run             2`;
const PRINT = 'fix             p all print 5 "$(step) $(temp:%.17g)" screen no file p.txt', RESC = 'fix             r all temp/rescale 10 1.0 1.0 0.0 1.0';
const scriptB = pf => `${BASE}${(pf ? [PRINT, RESC] : [RESC, PRINT]).join('\n')}
thermo_style    custom step temp f_r
thermo_modify   norm no format float %.17g
thermo          10
run             30`;

export default {
  id: 'hooks',
  title: 'Hooks and fix masks',
  minutes: 20,
  summary: 'A fix does not choose when it runs; it fills in a bit mask, and LAMMPS builds one list per hook, in the order you defined the fixes. You build those lists for any set of fixes, then change only the order of two lines and watch the forces go from 0 to 1, and a thermostat\'s printed temperature jump from 2.103 to 1.000. The energy fix temp/rescale removes is matched by hand to 17 digits.',
  learn: [
    'setmask(): 11 hook bits, and which fixes set which',
    'Modify::init builds one list per hook, in input-script order',
    'post_force order: addforce then setforce is not setforce then addforce',
    'end_of_step every N steps, and the energy fix temp/rescale takes out, by hand',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a mask and a list per hook</h2>
    <div class="prose"><p>Every fix has a <code>setmask()</code> that returns a sum of bits (${gh('src/fix.h', 308)}). For fix nve (${gh('src/fix_nve.cpp', 44)}):</p>
<pre class="lmp-mini">int FixNVE::setmask()
{
  int mask = 0;
  mask |= INITIAL_INTEGRATE;        // 1 &lt;&lt; 0 = 1
  mask |= FINAL_INTEGRATE;          // 1 &lt;&lt; 8 = 256
  mask |= INITIAL_INTEGRATE_RESPA;  // 1 &lt;&lt; 11 = 2048
  mask |= FINAL_INTEGRATE_RESPA;    // 1 &lt;&lt; 15 = 32768
  return mask;                      // 35073
}</pre>
<p>At <code>run</code>, ${gh('src/modify.cpp', 1692)} turns the masks into one list per hook, going through the fixes in the order they were defined:</p>
<pre class="lmp-mini">for (int i = 0; i &lt; nfix; i++)
  if (fmask[i] &amp; mask) list[n++] = i;</pre>
<p>and each stage of the loop (page C.2) just walks its list, ${gh('src/modify.cpp', 465)}:</p>
<pre class="lmp-mini">for (int i = 0; i &lt; n_post_force; i++) fix[list_post_force[i]]-&gt;post_force(vflag);</pre>
<p>Pick four fixes in your script's order:</p></div>
    <div class="row hk-sel">${DEF.map((d, i) => `<label class="lc-in">fix ${i + 1}<select class="lc-u" data-f="${i}">${Object.keys(FIXES).map(k => `<option${k === d ? ' selected' : ''}>${k}</option>`).join('')}</select></label>`).join('')}</div>
    <div class="tbl hk-masks"></div>
    <div class="tbl hk-lists"></div>
    <div class="prose"><p>fix nvt is not "fix nve plus a thermostat fix": it takes the same two hooks and does its own kicks, with the thermostat half-steps around them (${gh('src/fix_nh.cpp', 834)} then nve_v at 859, nve_x at 865). Use one integrator per atom, or the atoms are moved twice.</p></div>
  </section>

  <section class="block">
    <h2>2 · Modify: swap two post_force fixes</h2>
    <div class="prose"><p>Both fixes sit in the post_force list. setforce sets fx to 0, addforce adds 1 to fx. Predict min, max and sum of fx over the 4000 atoms, then run. (The pair forces alone sum to zero: Newton's third law.)</p></div>
    <div class="row"><label class="lc-in">order<select class="lc-u hk-o"><option value="as">addforce, then setforce</option><option value="sa">setforce, then addforce</option><option value="a">addforce only</option><option value="s">setforce only</option></select></label><button class="btn sm primary" data-k="a">▶ Run 2 steps</button><span class="lmp-status muted hk-sta"></span></div>
    <pre class="lmp-mini hk-scra"></pre>
    <div class="tbl hk-outa"></div>
  </section>

  <section class="block">
    <h2>3 · Match: end_of_step and the thermostat's energy</h2>
    <div class="prose"><p><code>fix temp/rescale 10 1.0 1.0 0.0 1.0</code> runs at end_of_step on steps 10, 20, 30 (${gh('src/modify.cpp', 493)}: <code>ntimestep % nevery == 0</code>). It scales every v by √(T<sub>target</sub>/T) and keeps a running total of the kinetic energy it removed (${gh('src/fix_temp_rescale.cpp', 165)}):</p>
<pre class="lmp-mini">factor  = sqrt(t_target / t_current);
energy += (t_current - t_target) * 0.5 * boltz * dof;   // dof = 3·4000 − 3 = 11997</pre>
<p>fix print sits in the same end_of_step list. Put it before or after the thermostat:</p></div>
    <div class="row"><label class="lc-in">order<select class="lc-u hk-p"><option value="1">print, then temp/rescale</option><option value="0">temp/rescale, then print</option></select></label><button class="btn sm primary" data-k="b">▶ Run 30 steps</button><span class="lmp-status muted hk-stb"></span></div>
    <div class="tbl hk-outb"></div>
    <div class="calc hk-hand" hidden></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe</h2>
    <ul>
      <li><b>Step 0 already has the fixes' forces.</b> In section 2, fx at step 0 is 0 or 1 too: Verlet::setup calls <code>modify->setup()</code> (${gh('src/verlet.cpp', 159)}), and fix setforce's setup just calls its own post_force (${gh('src/fix_setforce.cpp', 200)}).</li>
      <li><b>The thermostat's T is not exactly 1.</b> After rescaling, thermo prints 1.0000000000000024: the sum over 4000 atoms of the scaled v² rounds differently from T · factor². The 15th digit is round-off, not physics.</li>
      <li><b>Redefining a fix keeps its place.</b> <code>fix s all setforce …</code> again with the same ID replaces it in the same slot (${gh('src/modify.cpp', 849)}). To move it later in the list, <code>unfix</code> it and define it again.</li>
      <li><b>Without <code>thermo_modify norm no</code></b>, LJ units divide extensive values by N: the sum of fx would print as 1, not 4000, and f_r as 1.654, not 6617.6.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> fix shake sets PRE_NEIGHBOR, POST_FORCE, POST_FORCE_RESPA, MIN_PRE_REVERSE, MIN_POST_FORCE and POST_RUN. What is its mask?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>8 + 128 + 16384 + 1048576 + 2097152 + 1024 = <b>3163272</b>.</div></div>
    <div class="ex"><p><b>2.</b> You have <code>fix 1 all nve</code>, <code>fix 2 all langevin …</code>, then <code>fix 3 all setforce 0 0 0</code> on a frozen wall group. Does the wall feel the Langevin random force?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. Both are post_force; langevin (fix 2) adds first, then setforce (fix 3) sets the wall's forces to 0. Swap them and the wall atoms keep the random and friction force.</div></div>
    <div class="ex"><p><b>3.</b> With <code>fix temp/rescale 10 1.0 1.0 0.5 1.0</code> (window 0.5), would anything happen at step 20 of section 3?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. T before rescaling at step 20 is 1.0082; |1.0082 − 1| &lt; 0.5, so the if at line 165 is false and v is not touched. At step 10 (T = 2.10) it still rescales.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const lists = () => {
      const sel = [...root.querySelectorAll('[data-f]')].map(s => s.value);
      $('.hk-masks').innerHTML = `<table>${head(['#', 'fix', 'setmask()', 'mask'])}${sel.map((k, i) => row([i + 1, k, gh(FIXES[k].f, FIXES[k].l), `<span class="mono">${mask(k)}</span>`])).join('')}</table>`;
      const rows = HOOKS.map(([b, n, w]) => { const on = sel.map((k, i) => (mask(k) & (1 << b)) ? `${i + 1} ${k}` : null).filter(Boolean); return on.length ? row([n, w, on.join(' → ')]) : ''; }).join('');
      const dup = sel.filter(k => k === 'nve' || k === 'nvt').length;
      $('.hk-lists').innerHTML = `<table>${head(['hook', 'stage', 'called in this order'])}${rows}</table>${dup > 1 ? '<p class="m-bad">Two integrators on the same atoms: each step they would both move them.</p>' : ''}`;
    };
    root.querySelectorAll('[data-f]').forEach(s => s.onchange = lists); lists();

    const showA = () => $('.hk-scra').textContent = ORDERS[$('.hk-o').value].join('\n');
    $('.hk-o').onchange = () => { showA(); $('.hk-outa').innerHTML = ''; }; showA();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.hk-sta'), o = $('.hk-o').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptA(o)); b.disabled = false;
      if (res.error || !res.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; return; }
      st.textContent = '';
      const P = PRED[o], near = (p, v) => p === '?' || Math.abs(v - p) < 1e-9;
      $('.hk-outa').innerHTML = `<table>${head(['step', 'min fx', 'max fx', 'sum fx', 'your prediction'])}${rowsOf(res).map(t => { const v = [t.c_mn, t.c_mx, t.c_sm]; const ok = v.every((x, i) => near(P[i], x)); return row([t.Step, ...v.map(x => `<span class="mono">${x}</span>`), `${P.join(', ')} ${ok ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}`]); }).join('')}</table>
      <p class="muted">${o === 'a' ? 'min and max depend on the pair forces; the sum is 4000 to round-off.' : o === 'as' ? 'setforce runs last, so it wins.' : o === 'sa' ? 'addforce runs last: every atom feels exactly 1, and nothing else.' : ''}</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.hk-stb'), pf = $('.hk-p').value === '1'; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptB(pf), { read: ['p.txt'] }); b.disabled = false;
      const txt = res.files?.['p.txt'];
      if (res.error || !txt) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; return; }
      st.textContent = '';
      const pr = {}; for (const l of txt.trim().split('\n')) { const w = l.trim().split(/\s+/); if (w.length === 2) pr[+w[0]] = w[1]; }
      const th = Object.fromEntries(rowsOf(res).map(t => [t.Step, t]));
      $('.hk-outb').innerHTML = `<table>${head(['step', 'fix print says T', 'thermo says T', 'f_r'])}${[5, 10, 15, 20, 25, 30].map(s => row([s, `<span class="mono">${pr[s] ?? ''}</span>`, th[s] ? `<span class="mono">${th[s].Temp}</span>` : '', th[s] ? `<span class="mono">${th[s].f_r}</span>` : ''])).join('')}</table>
      <p class="muted">${pf ? 'print runs first: on steps 10, 20, 30 it sees T before rescaling.' : 'temp/rescale runs first: on steps 10, 20, 30 print already sees ≈ 1.'}</p>`;
      const hand = $('.hk-hand'); hand.hidden = false;
      if (!pf) { hand.textContent = 'Switch to "print, then temp/rescale" to get T before rescaling, which the hand sum needs.'; return; }
      let E = 0; const L = [];
      for (const s of [10, 20, 30]) { const T = parseFloat(pr[s]); const dE = (T - 1) * 0.5 * 1 * 11997; E += dE; const lm = +th[s].f_r;
        L.push(`step ${s}:  (${pr[s]} − 1) × ½ × 11997 = ${dE.toPrecision(17)}`);
        L.push(`         f_r = ${E.toPrecision(17)}   LAMMPS ${lm}  ${Math.abs(E - lm) <= 1e-12 * Math.abs(lm) ? '✓' : '✗'}`);
        L.push(`         factor = √(1 / ${T.toPrecision(6)}) = ${Math.sqrt(1 / T).toPrecision(10)}`); }
      hand.textContent = L.join('\n');
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
