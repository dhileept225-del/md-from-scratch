// Research · D.18: dumps, restarts and output. Dump columns (xu, xs) reproduced from x and image flags; what the default %g
// dump loses; four ways to continue a run compared bit for bit (two runs, binary restart, data file, uninterrupted);
// pair coefficients written with %g in a data file.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const SYS = eps => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 4 0 4 0 4
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 ${eps} 1.0 2.5
fix             1 all nve`;
const THERMO = `thermo_style    custom step pe ke press
thermo_modify   format float %.17g
thermo          50`;

const DUMP = `${SYS('1.0')}
${THERMO}
dump            d1 all custom 100 full.dump id x y z ix iy iz xu yu zu xs ys zs
dump_modify     d1 format float %.17g sort id
dump            d2 all custom 100 short.dump id x y z
dump_modify     d2 sort id
run             300`;

const ROUTES = eps => ({
  A: { label: 'one run of 200 steps', s: `${SYS(eps)}\n${THERMO}\nrun             200` },
  B: { label: 'run 100, run 100', s: `${SYS(eps)}\n${THERMO}\nrun             100\nrun             100` },
  C: { label: 'run 100, write_restart, clear, read_restart, run 100', s: `${SYS(eps)}\n${THERMO}\nrun             100\nwrite_restart   r.bin\nclear\nread_restart    r.bin\nfix             1 all nve\n${THERMO}\nrun             100` },
  D: { label: 'run 100, write_data, clear, read_data, run 100', s: `${SYS(eps)}\n${THERMO}\nrun             100\nwrite_data      d.data\nclear\nunits           lj\natom_style      atomic\npair_style      lj/cut 2.5\nread_data       d.data\nfix             1 all nve\n${THERMO}\nrun             100` },
});

function frames(text) {
  const L = text.split('\n'), out = [];
  for (let i = 0; i < L.length; i++) {
    if (!L[i].startsWith('ITEM: TIMESTEP')) continue;
    const step = +L[i + 1], n = +L[i + 3], box = [0, 1, 2].map(k => L[i + 5 + k].trim().split(/\s+/).map(Number)), A = [];
    for (let k = 0; k < n; k++) A.push(L[i + 9 + k].trim().split(/\s+/).map(Number));
    out.push({ step, box, A }); i += 8 + n;
  }
  return out;
}
// plain LJ energy, all pairs, nearest image (box > 2 rc)
function ljEnergy(A, L, eps = 1.0) {
  let e = 0.0;
  for (let i = 0; i < A.length; i++) for (let j = i + 1; j < A.length; j++) {
    let dx = A[i][1] - A[j][1], dy = A[i][2] - A[j][2], dz = A[i][3] - A[j][3];
    dx -= L * Math.round(dx / L); dy -= L * Math.round(dy / L); dz -= L * Math.round(dz / L);
    const rsq = dx * dx + dy * dy + dz * dz; if (rsq < 6.25) { const r6 = 1 / (rsq * rsq * rsq); e += 4 * eps * r6 * (r6 - 1); }
  }
  return e / A.length;
}
// thermo rows of every run in the output (step, pe, ke, press as numbers)
function allRows(lines) {
  const out = [];
  for (let i = 0; i < lines.length; i++) {
    if (!/^\s*Step\s+PotEng\s+KinEng\s+Press/.test(lines[i])) continue;
    for (let k = i + 1; k < lines.length; k++) { const w = lines[k].trim().split(/\s+/); if (w.length !== 4 || !/^\d+$/.test(w[0])) break; out.push(w.map(Number)); }
  }
  return out;
}

export default {
  id: 'output',
  title: 'Dumps, restarts and output',
  minutes: 18,
  summary: 'What leaves LAMMPS, and how much of it comes back. You reproduce the unwrapped and scaled coordinates of a dump from x and the image flags, measure what the default 6-digit dump format costs, and continue a run four ways: two run commands, a binary restart and a data file all give the same trajectory to the last bit; one uninterrupted run does not. Then a data file quietly rounds a pair coefficient to 6 digits, and you predict the new energy by hand.',
  learn: [
    'dump custom: when it writes, xu = x + ix·Lx, xs = (x − xlo)/Lx, and the default %g',
    'write_restart vs write_data: what each stores, exactly or not',
    'Why "run 100; run 100" is not "run 200", and a restart is the former',
    'Pair Coeffs in a data file are written with %g',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three ways out</h2>
    <div class="prose"><p><b>dump</b> writes per-atom columns every N steps, starting at the setup of each run (step 0 here). Unwrapped and scaled coordinates are computed when written (${gh('src/dump_custom.cpp', 2661)}, ${gh('src/dump_custom.cpp', 2565)}):</p>
<pre class="lmp-mini">xu = x + ((image &amp; IMGMASK) − IMGMAX) · xprd          // x + ix·Lx
xs = (x − boxxlo) · (1.0/xprd)</pre>
<p>The default format for every floating-point column is <code>%g</code>: 6 significant digits (line 145). Atoms come in memory order unless <code>dump_modify sort id</code>.</p>
<p><b>write_data</b> is text. Positions and velocities are written with fmt's shortest exact form (<code>{}</code>), so they read back bit for bit; masses and scalar per-atom values (charge, diameter) with 16 digits (<code>{:.16}</code>, ${gh('src/atom_vec.cpp', 1810)}), which is not always enough; pair coefficients with the pair style's own <code>%g</code> (${gh('src/pair_lj_cut.cpp', 644)}). It does not store the step number, fixes or computes.</p>
<p><b>write_restart</b> is binary: the doubles themselves, plus the step number, units, pair style and coefficients, neighbour settings and the state of some fixes. Fixes still have to be defined again after read_restart.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: what a dump contains</h2>
    <div class="prose"><p>256 LJ atoms, 300 steps. One dump with 17 digits and every coordinate flavour, one with the defaults.</p></div>
    <pre class="lmp-mini">${esc(DUMP)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run and check by hand</button><span class="lmp-status muted ou-sta"></span></div>
    <div class="tbl ou-a"></div>
    <div class="prose ou-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: four ways to continue a run</h2>
    <div class="prose"><p>Same system, 200 steps in total, split four ways. Predict which agree to the last bit.</p></div>
    <pre class="lmp-mini">${esc(ROUTES('1.0').C.s.split('\n').slice(-9).join('\n'))}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run all four</button><span class="lmp-status muted ou-stb"></span></div>
    <div class="tbl ou-b"></div>
    <div class="prose ou-bsum"></div>
    <div class="prose"><p>Now the same with ε = 0.123456789 instead of 1.0, and the data-file route with no <code>pair_coeff</code> after read_data, so the coefficient comes from the file's Pair Coeffs section.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Run with ε = 0.123456789</button><span class="lmp-status muted ou-stc"></span></div>
    <pre class="lmp-mini ou-file" hidden></pre>
    <div class="tbl ou-c"></div>
    <div class="prose ou-csum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> An atom has x = 0.2 and ix = −3 in a box of length 6.718: xu?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>0.2 − 3 × 6.718 = −19.954. It has crossed the lower x face three times more often than the upper one. Image flags are what make MSD (page D.15) possible from wrapped coordinates.</div></div>
    <div class="ex"><p><b>2.</b> Which file would you use to continue a run on a different machine, or to edit the system by hand?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The data file: text, portable, editable, exact for positions and velocities. Binary restarts can depend on the LAMMPS version and platform; convert them with <code>lmp -restart2data</code>. But check that the data file's Pair Coeffs (6 digits) are what you want, or set pair_coeff again after read_data.</div></div>
    <div class="ex"><p><b>3.</b> Why does read_data start at step 0 and read_restart at step 100?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The restart file stores update->ntimestep; the data file only mentions it in the comment line. Use <code>reset_timestep 100</code> after read_data if output files must continue the numbering.</div></div>
    <div class="ex"><p><b>4.</b> Postprocessing a default-format dump: how large is the position error, and does it matter for g(r) or MSD?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Half a unit in the 6th digit: up to 5 × 10⁻⁶ × x, about 3 × 10⁻⁵σ near the box edge. Invisible in g(r) or MSD. It matters for anything recomputed from positions and compared closely: energies and forces (rerun), or bit-for-bit checks like the ones on this site.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.ou-sta'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(DUMP, { read: ['full.dump', 'short.dump'] }); b.disabled = false;
      if (res.error || !res.files['full.dump']) return err(st, res.error); st.textContent = '';
      const F = frames(res.files['full.dump']), S = frames(res.files['short.dump']), rows = allRows(res.lines);
      let n = 0, ok = 0, maxim = 0;
      for (const f of F) {
        const lo = f.box.map(r => r[0]), prd = f.box.map(r => r[1] - r[0]);
        for (const a of f.A) for (let k = 0; k < 3; k++) {
          const x = a[1 + k], im = a[4 + k], xu = x + im * prd[k], xs = (x - lo[k]) * (1.0 / prd[k]);
          n += 2; ok += (xu === a[7 + k]) + (xs === a[10 + k]); maxim = Math.max(maxim, Math.abs(im));
        }
      }
      const L = F[0].box[0][1] - F[0].box[0][0], tab = [];
      let maxdx = 0;
      F.forEach((f, k) => {
        const s = S[k]; for (let i = 0; i < f.A.length; i++) for (let c = 1; c <= 3; c++) maxdx = Math.max(maxdx, Math.abs(f.A[i][c] - s.A[i][c]));
        const eFull = ljEnergy(f.A, L), eShort = ljEnergy(s.A, L), lmp = rows.find(r => r[0] === f.step)[1];
        tab.push(row([f.step, `<span class="mono">${s17(lmp)}</span>`, (eFull - lmp).toExponential(1), (eShort - lmp).toExponential(1)]));
      });
      const a0 = F[F.length - 1].A[0], s0 = S[S.length - 1].A[0];
      $('.ou-a').innerHTML = `<table>${head(['step', 'PotEng (LAMMPS)', 'hand from %.17g dump − LAMMPS', 'hand from %g dump − LAMMPS'])}${tab.join('')}</table>`;
      $('.ou-asum').innerHTML = `<p>xu and xs by hand from x, the image flags and the box in the dump header: ${ok} of ${n} values identical to the last bit (image flags up to ±${maxim} after 300 steps). Atom 1 at step 300 in the two dumps: <code>${esc(a0.slice(1, 4).map(s17).join(' '))}</code> and <code>${esc(s0.slice(1, 4).join(' '))}</code>; the largest position difference over all atoms and frames is ${maxdx.toExponential(2)}σ.</p><p>The table recomputes the energy per atom from each dump (plain double loop, all pairs, so the pair order differs from LAMMPS's neighbour list). From the 17-digit dump it agrees to round-off; from the default dump it is off by as much as 10⁻⁵ per atom (a few parts per million), except at step 0, where the lattice positions happen to need few digits. Fine for pictures, g(r) and MSD; not fine for <code>rerun</code> checks or hand calculations: always set <code>dump_modify format float %.17g</code> (or <code>format line</code>) when the numbers will be used again.</p>`;
    };

    const runRoutes = async (eps, st) => {
      const R = ROUTES(eps), out = {};
      if (eps !== '1.0') R.D.label = 'run 100, write_data, clear, read_data (no pair_coeff), run 100';
      for (const k of Object.keys(R)) {
        st.textContent = `Running ${k}…`;
        const res = await runLammps(R[k].s, { read: k === 'D' ? ['d.data'] : [] });
        if (res.error) return { error: res.error };
        out[k] = { rows: allRows(res.lines), file: res.files['d.data'], label: R[k].label };
      }
      st.textContent = ''; return out;
    };
    const last = o => o.rows[o.rows.length - 1];
    const table = O => {
      const ref = last(O.B);
      return `<table>${head(['route', 'steps printed', 'PotEng at the end', 'same as B'])}${Object.entries(O).map(([k, o]) => { const l = last(o), same = l[1] === ref[1] && l[2] === ref[2] && l[3] === ref[3]; return row([`<b>${k}</b> ${esc(o.label)}`, `${o.rows[o.rows.length - 3][0]}…${l[0]}`, `<span class="mono">${s17(l[1])}</span>`, same ? OK : BAD]); }).join('')}</table>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.ou-stb'); b.disabled = true;
      const O = await runRoutes('1.0', st); b.disabled = false; if (O.error) return err(st, O.error);
      $('.ou-b').innerHTML = table(O);
      const a = last(O.A), bb = last(O.B), d100 = O.B.rows.filter(r => r[0] === 100);
      $('.ou-bsum').innerHTML = `<p>B, C and D end on the same numbers to the last bit: the restart and the data file both bring back exactly the positions and velocities B had after its first run, and all three then start a new run from them. A, the uninterrupted run, differs from step 150 on (PotEng ${(a[1] - bb[1]).toExponential(1)} apart at step 200). The difference is the second setup, not the files: B prints step 100 twice, ${s17(d100[0][1])} at the end of the first run and ${s17(d100[1][1])} at the start of the second, because setup sorts the atoms and rebuilds the neighbour list, which changes the order of the pair sums (page C.6). Note also D's steps: read_data starts at 0, read_restart at 100.</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.ou-stc'); b.disabled = true;
      const O = await runRoutes('0.123456789', st); b.disabled = false; if (O.error) return err(st, O.error);
      const pc = O.D.file.split('\n'), i = pc.findIndex(l => l.startsWith('Pair Coeffs'));
      $('.ou-file').hidden = false; $('.ou-file').textContent = pc.slice(i, i + 3).join('\n') + '\n…';
      $('.ou-c').innerHTML = table(O);
      const e100 = O.B.rows.filter(r => r[0] === 100)[1][1], d0 = O.D.rows.filter(r => r[0] === 0)[1][1], pred = e100 * (0.123457 / 0.123456789);
      $('.ou-csum').innerHTML = `<p>The data file holds ε = 0.123457: written with %g. D's first energy is therefore B's energy at step 100 times 0.123457/0.123456789 (every pair term is proportional to ε): by hand ${s17(pred)}, LAMMPS ${s17(d0)} (${Math.abs(pred - d0) < 1e-14 ? 'equal to round-off' : 'differ'}, the multiplication order is not LAMMPS's). The run then follows a slightly different potential, ${((0.123457 / 0.123456789 - 1) * 1e6).toFixed(2)} parts per million stronger. The default <code>write_data</code> writes only the diagonal coefficients and lets the mixing rule rebuild the rest; <code>write_data … pair ij</code> writes all of them, also with %g. For a force field whose numbers have more than 6 significant digits (unit conversions, fitted parameters), give pair_coeff again after read_data, or use the binary restart.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
