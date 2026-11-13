// Research · B.6: how atoms are stored. Per-atom arrays, owned and ghost atoms counted by hand, sorting, image flags checked by hand.
import { runLammps, g8 } from './lammps-runner.js';
import { IN_MELT } from './terminal.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const A = Math.cbrt(4 / 0.8442), H = A / 2, L = 10 * A, NP = 20;   // melt: 10 fcc cells, planes every a/2

// ghosts by hand: planes within the communication cutoff of each face
function ghosts(rc, skin) {
  const cut = rc + skin;
  const lo = [], hi = [];
  for (let p = 0; p < NP; p++) { if (p * H <= cut) lo.push(p); if (p * H >= L - cut) hi.push(p); }
  const idx = []; for (let i = -hi.length; i < NP + lo.length; i++) idx.push(i);
  let n = 0; for (const i of idx) for (const j of idx) for (const k of idx) if (((i + j + k) % 2 + 2) % 2 === 0) n++;
  return { cut, lo: lo.length, hi: hi.length, side: idx.length, sites: n, ghost: n - 4000 };
}
const GSCRIPT = (rc, skin) => IN_MELT.replace('pair_style      lj/cut 2.5', `pair_style      lj/cut ${rc}`).replace('pair_coeff      1 1 1.0 1.0 2.5', `pair_coeff      1 1 1.0 1.0 ${rc}`).replace('neighbor        0.3 bin', `neighbor        ${skin} bin`).replace('run             250', 'run             0');
const ISCRIPT = sortEvery => IN_MELT.replace('run             250', `atom_modify     sort ${sortEvery} 0.0
dump            1 all custom 100 a.dump id x xu ix
dump_modify     1 format float %.15g
run             250`);

export default {
  id: 'atom-storage',
  title: 'How atoms are stored',
  minutes: 20,
  summary: 'Inside LAMMPS an atom is not an object: it is one row number in a set of long arrays (x, v, f, type, tag, image …). You see which arrays an atom_style creates, why every processor also stores copies of other atoms (ghosts) and count them by hand, watch LAMMPS reorder its atoms in memory, and check by hand how image flags rebuild the true path of an atom that crossed the periodic box.',
  learn: [
    'Per-atom arrays, the local index i and the atom ID (tag)',
    'Owned and ghost atoms: count the ghosts of the melt crystal by hand',
    'Why LAMMPS sorts atoms in memory, and what that does to the order you see',
    'Image flags: unwrapped x = x + ix·L, checked against LAMMPS',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: an atom is a row number</h2>
    <p>The class <code>Atom</code> (${gh('src/atom.h', 72)}) holds one array per property. Atom number <code>i</code> is simply row <code>i</code> of every array:</p>
<pre class="lmp-mini">46  int nlocal, nghost;      // atoms this processor owns, and copies of other atoms
47  int nmax;                // length the arrays are allocated for
72  tagint *tag;             // the atom ID you see in data files and dumps
73  int *type, *mask;        // type 1, 2 …; mask = which groups it belongs to (one bit per group)
74  imageint *image;         // how many times it crossed the periodic box (section 5)
75  double **x, **v, **f;    // position, velocity, force: x[i][0], x[i][1], x[i][2]</pre>
    <p>Rows <code>0 … nlocal−1</code> are the atoms this processor <b>owns</b> and moves. Rows <code>nlocal … nlocal+nghost−1</code> are <b>ghosts</b>: copies of atoms that sit just outside its region, kept so every owned atom sees all its neighbours (section 3). The row number <code>i</code> is not the atom ID: <code>tag[i]</code> is. They start equal and drift apart (section 4).</p>
    <p>Which arrays exist is decided by <code>atom_style</code>. Every style has the defaults (${gh('src/atom_vec.cpp', 33)}); each adds its own list (<code>fields_grow</code>):</p>
    <div class="tbl wide"><table>${head(['atom_style', 'arrays per atom', 'bytes per atom (approx.)'])}
${row(['<code>atomic</code>', 'tag, type, mask, image, x, v, f', '4 × 4 + 9 × 8 = <b>88</b>'])}
${row(['<code>charge</code>', '+ q (charge)', '88 + 8 = 96'])}
${row(['<code>full</code>', '+ q, molecule, and lists of the bonds, angles, dihedrals, impropers and special neighbours of each atom', 'a few hundred, depending on the per-atom limits'])}</table></div>
    <p class="muted">${gh('src/atom_vec_charge.cpp', 33)} · ${gh('src/MOLECULE/atom_vec_full.cpp', 35)}. Integers are 4 bytes, doubles 8, in the default build.</p>
  </section>

  <section class="block">
    <h2>2 · Predict: how many ghost atoms?</h2>
    <div class="prose">
      <p>In serial, the box is the processor's region, and its neighbours across a periodic face are copies of its own atoms. LAMMPS makes the copies one direction at a time (${gh('src/comm_brick.cpp', 775)}). For each face it copies every atom within the <b>communication cutoff</b> of the opposite face (page B.2: pair cutoff + skin):</p>
<pre class="lmp-mini">429  slabhi[iswap] = sublo[dim] + cutghost[dim];   // atoms near the low face → ghosts beyond the high face
456  slablo[iswap] = subhi[dim] - cutghost[dim];   // atoms near the high face → ghosts beyond the low face
842  if (x[i][dim] >= lo && x[i][dim] <= hi)       // inside the slab: send it</pre>
      <p>Because x, then y, then z are done in turn, ghosts of ghosts fill the edges and corners: the result is the whole box grown by the cutoff on every side. At step 0 the melt atoms sit on FCC planes every a/2 = ${g8(H)} apart (page A.5), 20 planes per direction, so the count can be done by hand.</p>
    </div>
    <div class="lc-in">
      <label>pair cutoff <input type="number" class="lmp-in sm" data-g="rc" value="2.5" step="0.1" min="1.2" max="4"></label>
      <label>skin <input type="number" class="lmp-in sm" data-g="skin" value="0.3" step="0.1" min="0" max="1.5"></label>
    </div>
    <div class="calc as-hand"></div>
    <div class="row"><button class="btn sm primary" data-k="ghost">▶ Run with LAMMPS</button><span class="lmp-status muted as-st"></span></div>
    <div class="lmp-compare as-cmp"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe: the count moves in steps</h2>
    <p>Try skin 0.3, 0.5 and 0.9 with cutoff 2.5. The first two give the same 5841 ghosts; the third jumps. Ghosts come in whole planes, so nothing changes until the cutoff passes the next plane at a multiple of ${g8(H)} (3 × ${g8(H)} = ${g8(3 * H)}, 4 × ${g8(H)} = ${g8(4 * H)}). After melting there are no planes, and the count changes smoothly: 5506 at the end of melt (page A.9).</p>
    <p>Ghosts cost memory and time: here there are more ghosts (5841) than owned atoms (4000). A small box, or many processors each holding a small region, makes this worse (page B.7).</p>
  </section>

  <section class="block">
    <h2>4 · Observe: LAMMPS reorders its atoms</h2>
    <div class="prose"><p>Every 1000 steps by default (${gh('src/atom.cpp', 103)}), <code>Atom::sort()</code> (${gh('src/atom.cpp', 2379)}) reorders all arrays so atoms close in space sit close in memory, which makes the force loop faster. A dump without <code>dump_modify sort id</code> writes atoms in memory order, so the IDs come out shuffled. Here sorting is set to every 100 steps:</p></div>
    <div class="row"><button class="btn sm primary" data-k="img">▶ Run melt and read the dump</button><span class="lmp-status muted ai-st"></span></div>
    <div class="tbl ai-order"></div>
  </section>

  <section class="block">
    <h2>5 · Match: image flags by hand</h2>
    <div class="prose">
      <p>When an atom leaves the box through a face, LAMMPS puts it back on the other side, so <code>x</code> always stays between 0 and L = ${g8(L)}, and it adds or subtracts 1 from the atom's image flag <code>ix</code>. The true, unwrapped position is then (${gh('src/dump_custom.cpp', 2665)}):</p>
<pre class="lmp-mini">2665  buf[n] = x[j][0] + ((image[j] &amp; IMGMASK) - IMGMAX) * xprd;    // xu = x + ix · L</pre>
      <p>Diffusion (the mean squared displacement lesson) must use the unwrapped xu: an atom that steps across a face would otherwise appear to jump by a whole box length. Below, the atoms in the step-200 dump with ix ≠ 0: predict xu from x and ix.</p>
    </div>
    <div class="tbl ai-img"><p class="muted">Run section 4 first.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why are there 4 ghost planes beyond the high face but only 3 beyond the low face, with cutoff 2.8?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Near the low face the planes at 0, 0.84, 1.68 and 2.52 are all within 2.8, including the plane at x = 0 itself: 4 planes, which become ghosts at L, L + 0.84 …. Near the high face the planes are at L − 0.84, L − 1.68, L − 2.52: only 3, because there is no plane at x = L (that is the plane at 0 again).</div></div>
    <div class="ex"><p><b>2.</b> About how much memory do the per-atom arrays of melt need at step 0?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>(4000 owned + 5841 ghosts) × 88 bytes ≈ <b>0.87 MB</b>. LAMMPS prints a larger "Per MPI rank memory allocation" because it adds the neighbour lists, communication buffers and spare room (nmax &gt; nlocal + nghost).</div></div>
    <div class="ex"><p><b>3.</b> An atom has x = 0.4 and ix = −2 in a box of length 20. Where has it really gone?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>xu = 0.4 + (−2) × 20 = <b>−39.6</b>: it left through the low face twice.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let G = null;
    const num = (k, d, lo, hi) => { const v = parseFloat(root.querySelector(`[data-g=${k}]`).value); return isFinite(v) ? Math.min(hi, Math.max(lo, v)) : d; };
    const draw = () => {
      const rc = num('rc', 2.5, 1.2, 4), skin = num('skin', 0.3, 0, 1.5); G = { ...ghosts(rc, skin), rc, skin };
      $('.as-hand').innerHTML = esc(`cutoff = ${rc} + ${skin} = ${g8(G.cut)}
planes within it:
  of the low face  ${G.lo}
  of the high face ${G.hi}
per direction: ${G.hi} + 20 + ${G.lo} = ${G.side}
FCC sites in ${G.side}³ = ${G.sites}
  (grid points with i + j + k even)
ghosts = ${G.sites} − 4000 = `) + `<b>${G.ghost}</b>`;
      $('.as-cmp').innerHTML = ''; $('.as-st').textContent = '';
    };
    root.querySelectorAll('[data-g]').forEach(i => { i.oninput = draw; }); draw();
    $('[data-k=ghost]').onclick = async e => {
      const b = e.currentTarget, st = $('.as-st'), g = G; b.disabled = true; st.textContent = 'Running the real LAMMPS… (the first run downloads about 3 MB)';
      const res = await runLammps(GSCRIPT(g.rc, g.skin)); b.disabled = false;
      if (res.error) { st.innerHTML = `<span class="lmp-err">${esc(res.error)}</span>`; return; }
      st.textContent = '';
      const get = re => { const l = res.lines.find(x => re.test(x)); return l ? l.match(re)[1] : '—'; };
      const nl = get(/^Nlocal:\s+(\d+)/), ng = get(/^Nghost:\s+(\d+)/);
      $('.as-cmp').innerHTML = `<div class="tbl"><table><tr><th>cutoff ${g8(g.cut)}</th><th>by hand</th><th>LAMMPS</th><th></th></tr>
        <tr><td>Nlocal (owned)</td><td class="mono">4000</td><td class="mono">${nl}</td><td>${nl === '4000' ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>
        <tr><td>Nghost</td><td class="mono">${g.ghost}</td><td class="mono">${ng}</td><td>${ng === String(g.ghost) ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr></table></div>`;
    };
    $('[data-k=img]').onclick = async e => {
      const b = e.currentTarget, st = $('.ai-st'); b.disabled = true; st.textContent = 'Running melt with a dump every 100 steps…';
      const res = await runLammps(ISCRIPT(100), { read: ['a.dump'] }); b.disabled = false;
      const txt = res.files && res.files['a.dump'];
      if (res.error || !txt) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no dump written')}</span>`; return; }
      st.textContent = '';
      const frames = txt.split('ITEM: TIMESTEP').slice(1).map(f => { const ls = f.trim().split('\n'); return { step: ls[0].trim(), rows: ls.slice(8).map(l => l.trim().split(/\s+/)) }; });
      $('.ai-order').innerHTML = `<table>${head(['step', 'IDs in the first 8 rows of memory (rows i = 0 … 7)'])}${frames.map(f => row([f.step, `<span class="mono">${f.rows.slice(0, 8).map(r => r[0]).join(', ')}</span>`])).join('')}</table>
        <p class="muted">At step 0 row i holds ID i + 1. After each sort the rows hold the atoms of the first spatial bin, whatever their IDs.</p>`;
      const f2 = frames.find(f => f.step === '200') || frames.at(-1), cross = f2.rows.filter(r => +r[3] !== 0);
      const show = cross.slice(0, 8);
      $('.ai-img').innerHTML = `<p>Step ${f2.step}: <b>${cross.length}</b> of 4000 atoms have ix ≠ 0. L = ${g8(L)}.</p><table>${head(['ID', 'x', 'ix', 'x + ix·L by hand', 'xu from LAMMPS', ''])}${show.map(r => { const h = +r[1] + (+r[3]) * L, lv = +r[2], ok = Math.abs(h - lv) < 1e-9 * Math.max(1, Math.abs(lv)); return row([r[0], `<span class="mono">${g8(+r[1])}</span>`, r[3], `<span class="mono">${g8(h)}</span>`, `<span class="mono">${g8(lv)}</span>`, ok ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>']); }).join('')}</table>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
