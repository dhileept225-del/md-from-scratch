// Research · E.1: hello from inside the engine. One logmesg line added to PairLJCut::compute, rebuilt natively,
// and every printed rsq, fpair and evdwl recomputed by hand from the dumped positions. Two side effects observed:
// a dump every step makes LAMMPS compute energies every step, and with MPI the screen shows every rank but the log only rank 0.
import { NATIVE } from './hello-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const clean = s => s.replace(/\/home\/[^ ]*?\/src\//g, '…/src/');

const INPUT = `units           lj
atom_style      atomic
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 single 4.0 5.0 5.0
create_atoms    1 single 5.2 5.1 5.0
create_atoms    1 single 4.5 6.0 5.4
mass            1 1.0
velocity        all set 0.0 0.0 0.0
set             atom 1 vx 0.3
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
fix             1 all nve
dump            d all custom 1 x.dump id x y z
dump_modify     d format float %.17g sort id
thermo_style    custom step pe ke etotal
thermo_modify   norm no format float %.17g
thermo          2
run             3`;

// the HELLO lines and thermo lines of an output
const parse = lines => {
  const hello = [], thermo = {};
  for (const l of lines) {
    const m = l.match(/^HELLO rank (\d+) step (\d+) pair (\d+)-(\d+) rsq (\S+) fpair (\S+) evdwl (.*)$/);
    if (m) hello.push({ rank: +m[1], step: +m[2], i: +m[3], j: +m[4], rsq: +m[5], fpair: +m[6], evdwl: m[7] === 'not computed' ? null : +m[7] });
    else { const w = l.trim().split(/\s+/); if (w.length === 4 && /^\d+$/.test(w[0])) thermo[+w[0]] = +w[1]; }
  }
  return { hello, thermo };
};
const dumpFrames = text => {
  const L = text.split('\n'), F = {};
  for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const s = +L[i + 1], n = +L[i + 3], X = {}; for (let k = 0; k < n; k++) { const w = L[i + 9 + k].trim().split(/\s+/).map(Number); X[w[0]] = w.slice(1); } F[s] = X; i += 8 + n; }
  return F;
};
// the lines of PairLJCut::compute for one pair, ε = σ = 1 (lj1 = 48, lj2 = 24, lj3 = lj4 = 4, offset 0, factor_lj 1)
function pairByHand(xi, xj) {
  const delx = xi[0] - xj[0], dely = xi[1] - xj[1], delz = xi[2] - xj[2];
  const rsq = delx * delx + dely * dely + delz * delz, r2inv = 1.0 / rsq, r6inv = r2inv * r2inv * r2inv;
  const forcelj = r6inv * (48.0 * r6inv - 24.0), fpair = 1.0 * forcelj * r2inv;
  let evdwl = r6inv * (4.0 * r6inv - 4.0) - 0.0; evdwl *= 1.0;
  return { rsq, fpair, evdwl };
}

export default {
  id: 'hello-engine',
  title: 'Hello from inside the engine',
  minutes: 15,
  summary: 'The first change to LAMMPS itself: one line in PairLJCut::compute that prints every pair the force loop handles. You read the diff, see which files the rebuild touches, and check all twelve printed lines (rsq, force, energy) against a hand calculation from the dumped positions, to the last bit. The print also reveals two things you could not see from outside: a dump every step makes LAMMPS compute energies every step, and with MPI the screen shows every processor but log.lammps only processor 0. Then the change is undone and the original numbers come back.',
  learn: [
    'Where in the pair loop a print goes, and why {:.17g}',
    'utils::logmesg: screen and log file, and what MPI ranks do with it',
    'eflag: energies are computed only on steps that need them',
    'Rebuild only what changed; undo with git checkout',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the loop and where to listen</h2>
    <div class="prose"><p>${gh('src/pair_lj_cut.cpp', 69)} visits every pair in the neighbour list. Inside <code>if (rsq &lt; cutsq)</code> it has everything about one pair in local variables: <code>rsq</code>, <code>fpair</code>, and <code>evdwl</code> when <code>eflag</code> is set. Just before <code>ev_tally</code> is the place to print them. Three choices:</p>
<ul><li><b>utils::logmesg(lmp, …)</b>, LAMMPS's own print: it writes to the screen and to log.lammps, with fmt-style <code>{}</code> placeholders (${gh('src/utils.cpp', 255)}).</li>
<li><b>{:.17g}</b>: 17 significant digits, enough to recover the exact double, so the hand check can be bit for bit.</li>
<li><b>eflag</b> is 0 on steps where nobody asked for an energy; evdwl then still holds a stale value from an earlier pair, so the line prints "not computed" instead.</li></ul>
<p>The test: three atoms in a big box, atom 1 moving at 0.3σ/τ, 3 steps, thermo every 2 steps, positions dumped every step with 17 digits.</p></div>
    <pre class="lmp-mini">${esc(INPUT)}</pre>
  </section>

  <section class="block">
    <h2>2 · Change and rebuild</h2>
    <pre class="lmp-mini gd">${esc(NATIVE.diff)}</pre>
    <div class="prose"><p>Then <code>cd build; make -j2</code>. Make compares file times, so only the changed file is compiled again, then the library and the executable are linked:</p></div>
    <pre class="lmp-mini">${esc(clean(NATIVE.make.join('\n')))}</pre>
    <div class="prose"><p>About 10–25 s on two cores, against 266 s for the full build (page 0.4). "Generating lmpgitversion.h" is CMake noticing that the git tree now has local changes; the version banner of this binary says "${esc(NATIVE.version)}", the build still carries its git hash and a "modified" flag.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: twelve lines from inside, by hand</h2>
    <div class="prose"><p>The run (native, serial) prints three pairs per step. Your browser takes the dumped positions of the two atoms and repeats the lines of compute() for each pair.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Check every HELLO line</button></div>
    <div class="tbl he-a"></div>
    <div class="prose he-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: what the print reveals</h2>
    <div class="prose"><p><b>Energies every step?</b> Thermo asks for the energy at steps 0 and 2 (and the last step, 3), yet the run above printed evdwl at step 1 too. Same input without the dump:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.nodump.filter(l => /step 1 /.test(l)).join('\n'))}</pre>
    <div class="prose"><p>Now step 1 has no energy. A dump may contain computes, so Output::setup and Output::write register every dump step with all computes that track time (modify->addstep_compute, output.cpp:278), compute pe among them; Integrate::ev_set then sees that pe "matches" the step and turns eflag on. A dump every step therefore costs an energy evaluation every step, even when the dump writes only positions.</p>
<p><b>Two processors.</b> Same input with atom 1 moved to z = 4.6 and <code>mpirun -np 2</code>: the box is split at z = 5, atom 1 on rank 0, atoms 2 and 3 on rank 1. Step 0 on the screen and in log.lammps:</p></div>
    <div class="tbl he-b"></div>
    <div class="prose he-bsum"></div>
  </section>

  <section class="block">
    <h2>5 · Undo</h2>
    <div class="prose"><p><code>git checkout src/pair_lj_cut.cpp</code>, <code>make -j2</code> (same two files again), and the original binary is back: no HELLO lines, and the step-0 line of the input above is again <code>0 -2.6718824134574581 0.044999999999999998 -2.6268824134574582</code>, identical to the run before the change. The print changed nothing in the physics: it only reads variables.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> A 32 000-atom LJ liquid at ρ = 0.8442: how many HELLO lines per step?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Pairs inside 2.5σ only, each once (newton on): about 32 000 × 4π/3 × 2.5³ × 0.8442 / 2 ≈ 880 000 lines per step, roughly 100 MB. Guard the print, e.g. <code>if (atom->tag[i] == 1)</code>, or print only on one step.</div></div>
    <div class="ex"><p><b>2.</b> Why does rank 0 compute both pairs that cross the boundary, 1-2 and 1-3, and rank 1 none?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>With newton on, each pair is computed once, from a half stencil that looks "up" (higher z, then y, then x; page D.3). Rank 0 is below the boundary, so its ghost copies of atoms 2 and 3 are in the bins above atom 1. Rank 1's ghost of atom 1 lies below its atoms, outside their half stencil.</div></div>
    <div class="ex"><p><b>3.</b> Where would you print the total force on each atom after all pair forces are added?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Not in compute(): with newton on, ghost atoms still hold forces that reverse communication adds to their owners later. Print in Verlet::run after <code>comm->reverse_comm()</code> and before <code>modify->post_force()</code> (stage 4 on page C.2), or from a fix's post_force.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = () => {
      const { hello, thermo } = parse(NATIVE.serial), F = dumpFrames(NATIVE.dump);
      let ok = 0, n = 0; const sums = {};
      const rows = hello.map(h => {
        const H = pairByHand(F[h.step][h.i], F[h.step][h.j]);
        const same = H.rsq === h.rsq && H.fpair === h.fpair && (h.evdwl === null || H.evdwl === h.evdwl);
        n++; ok += same;
        if (h.evdwl !== null) sums[h.step] = (sums[h.step] ?? 0.0) + h.evdwl;
        return row([h.step, `${h.i}-${h.j}`, `<span class="mono">${s17(h.rsq)}</span>`, `<span class="mono">${s17(h.fpair)}</span>`, h.evdwl === null ? '–' : `<span class="mono">${s17(h.evdwl)}</span>`, same ? OK : BAD]);
      });
      $('.he-a').innerHTML = `<table>${head(['step', 'pair', 'rsq', 'fpair', 'evdwl', 'hand'])}${rows.join('')}</table>`;
      const chk = Object.keys(thermo).map(s => `step ${s}: ${s17(0.0 + sums[s])} ${sums[s] === thermo[s] ? OK : BAD}`);
      $('.he-asum').innerHTML = `<p>${ok} of ${n} lines identical to the last bit in rsq, fpair and evdwl. Added up in the printed order (the order ev_tally adds them to eng_vdwl), the three evdwl of a thermo step give PotEng: ${chk.join(', ')}. fpair is the force divided by r: the force on atom i is fpair·(x_i − x_j), and negative means attraction; all three pairs sit beyond the minimum at 2<sup>1/6</sup> (rsq &gt; 1.26), so the cluster is contracting.</p>`;
    };
    // MPI: screen vs log
    const S = parse(NATIVE.screen).hello.filter(h => h.step === 0), L = parse(NATIVE.log).hello.filter(h => h.step === 0);
    const X0 = { 1: [4.0, 5.0, 4.6], 2: [5.2, 5.1, 5.0], 3: [4.5, 6.0, 5.4] };
    $('.he-b').innerHTML = `<table>${head(['pair', 'computed on', 'on the screen', 'in log.lammps', 'rsq, fpair, evdwl by hand'])}${S.map(h => { const H = pairByHand(X0[h.i], X0[h.j]); return row([`${h.i}-${h.j}`, `rank ${h.rank}`, OK, L.some(l => l.i === h.i && l.j === h.j) ? OK : '–', H.rsq === h.rsq && H.fpair === h.fpair && H.evdwl === h.evdwl ? OK : BAD]); }).join('')}</table>`;
    const tot = parse(NATIVE.screen).hello.length, totL = parse(NATIVE.log).hello.length;
    $('.he-bsum').innerHTML = `<p>${tot} HELLO lines on the screen, ${totL} in the log. LAMMPS sets <code>screen = stdout</code> on every rank (lammps.cpp:543) but opens the log file only on rank 0, so logmesg from rank 1 reaches the terminal (in any order relative to rank 0) and never the file. LAMMPS's own messages are printed by rank 0 only (<code>if (comm->me == 0)</code>), which is why you normally never notice. For debugging, label each line with <code>comm->me</code> as here, and read the screen, not the log. The last column recomputes each pair from the create_atoms coordinates: identical, ghost atoms included. Splitting the work does not change the arithmetic of a pair, only who does it and in which order the sums are added.</p>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
