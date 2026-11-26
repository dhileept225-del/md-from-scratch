// Research · F.6: restarts and system surgery. A binary restart written natively is continued in the browser,
// with and without the thermostat's stored state (decoded from the file by hand); then four edits to a system
// (delete_atoms overlap, replicate, change_box remap) are predicted atom by atom from LAMMPS's own source.
import { NATIVE } from './restarts-native.js';
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const bin = () => Uint8Array.from(atob(NATIVE.bin), c => c.charCodeAt(0));
const readDump = txt => {   // first frame: {box:[[lo,hi]x3], at: Map(id -> [cols...])}
  const L = txt.split('\n'), n = +L[3], box = [5, 6, 7].map(k => L[k].trim().split(/\s+/).map(Number)), at = new Map();
  for (let k = 0; k < n; k++) { const w = L[9 + k].trim().split(/\s+/).map(Number); at.set(w[0], w.slice(1)); }
  return { box, at, n };
};

const STORED = [
  ['units, atom style, box, timestep number and dt', 'yes', 'read_restart starts at step 100'],
  ['positions, velocities, image flags, types, masses', 'yes, as raw doubles', 'exact'],
  ['groups', 'yes', ''],
  ['pair style and coefficients, special_bonds, neighbor settings', 'yes (most pair styles)', 'some styles ask for pair_coeff again'],
  ['fixes', 'only their internal state', 'you must define the fix again, with the same ID and style'],
  ['computes, variables, thermo_style, dumps, outputs', 'no', 'set them again'],
];
const SURGERY = [
  ['delete_atoms overlap / region / group', 'remove atoms; overlap keeps one of each close pair'],
  ['replicate nx ny nz', 'tile the box; new IDs are offsets of the old ones'],
  ['change_box … remap', 'resize or tilt the box; remap scales positions with it'],
  ['displace_atoms move / ramp / rotate', 'move a group rigidly'],
  ['set atom/type/group … type, charge, vx …', 'edit per-atom values, also a random fraction of types'],
  ['velocity … create / scale / zero linear', 'new or rescaled velocities, remove drift'],
  ['read_data … add append shift', 'merge a second system into the current box'],
];
const DEL = `units lj
atom_style atomic
region box block 0 8 0 8 0 8 units box
create_box 1 box
create_atoms 1 random 400 4321 NULL
mass 1 1.0
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
write_dump all custom before.dump id x y z modify sort id format float %.17g
delete_atoms overlap 0.8 all all compress no
write_dump all custom after.dump id modify sort id`;
const DEL2 = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 4 0 4 0 4
create_box 1 box
create_atoms 1 box
create_atoms 1 single 0.05 0 0      # page F.4's mistake, atom 257
mass 1 1.0
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
delete_atoms overlap 0.5 all all compress no
write_dump all custom after.dump id modify sort id`;
const MELT = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 4 0 4 0 4
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.44 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
fix 1 all nve
thermo_style custom step pe
thermo_modify format float %.17g
run 300`;
const REP = `${MELT}
write_dump all custom a.dump id x y z ix iy iz modify sort id format float %.17g
replicate 2 2 2
write_dump all custom b.dump id x y z modify sort id format float %.17g
run 0`;
const CB = `${MELT}
write_dump all custom a.dump id x y z modify sort id format float %.17g
change_box all x scale 1.1 y final 0 7 remap units box
write_dump all custom b.dump id x y z modify sort id format float %.17g`;

export default {
  id: 'restarts',
  title: 'Restarts and system surgery',
  minutes: 15,
  summary: 'Page D.18 showed that a restart brings positions and velocities back to the last bit. This page is about the rest. A binary restart written by native LAMMPS is continued in your browser, with and without the thermostat state it carries; that state is decoded from the file by hand. Then the system is edited between runs: delete_atoms overlap, replicate and change_box. For each one you predict, from LAMMPS\'s source, which atoms survive, what their new IDs are and where they end up, and check every atom.',
  learn: [
    'What a restart stores, and the "same fix ID and style" rule',
    'Decoding fix nvt\'s chain variables from the binary file',
    'delete_atoms overlap: which atom of a close pair goes',
    'replicate and change_box remap, atom by atom',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: what survives a restart</h2>
    <div class="tbl rs-st"></div>
    <div class="prose"><p>A fix that has state (a thermostat's chain variables, a running average, fix spring/self's anchor points) writes it into the restart file under its ID and style. After read_restart that state waits. It is given to the first fix defined with the <i>same ID and the same style</i>, and the log says "Resetting global fix info from restart file". Otherwise the run ends with "Unused restart file global fix info", and the fix starts from zero.</p></div>
  </section>

  <section class="block">
    <h2>2 · A native restart, continued here</h2>
    <div class="prose"><p>Native LAMMPS ran 256 LJ atoms for 100 steps under fix nvt and wrote <code>r.bin</code> (${(bin().length / 1024).toFixed(1)} kB):</p></div>
    <pre class="lmp-mini">${esc(NATIVE.inw)}</pre>
    <div class="prose"><p>Your browser receives that binary file unchanged and continues it twice: with <code>fix 1 all nvt</code> (same ID) and with <code>fix 2 all nvt</code> (new ID).</p></div>
    <pre class="lmp-mini">${esc(NATIVE.inr)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Continue r.bin both ways</button><span class="rs-sta"></span></div>
    <div class="prose rs-asum"></div>
    <div class="prose"><p><b>What fix 1 got back.</b> ${gh('src/fix_nh.cpp', 1278)} packs its state as [tstat_flag, mtchain, η₁…η₃, η̇₁…η̇₃, pstat_flag] after the fix ID and style. Find those nine doubles in r.bin and compute the thermostat energy the way compute_scalar does: L·kT·η₁ + ½Q₁η̇₁² + Σ<sub>k≥2</sub>(kT·η<sub>k</sub> + ½Q<sub>k</sub>η̇<sub>k</sub>²), with L = 3N − 3 = 765, Q₁ = L·kT/ω², Q<sub>k</sub> = kT/ω², ω = 1/T<sub>damp</sub> = 2.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Decode the thermostat from r.bin</button></div>
    <div class="prose rs-bsum"></div>
    <div class="prose"><p>For a file you need to read or edit, convert it: <code>lmp -restart2data r.bin r.data</code> (beginning):</p></div>
    <pre class="lmp-mini">${esc(NATIVE.data)}</pre>
    <div class="prose"><p>The data file keeps positions and velocities exactly (page D.18) but has no fix state at all, so a thermostatted run continued from it restarts its chain at η = 0, as with the new ID above.</p></div>
  </section>

  <section class="block">
    <h2>3 · Surgery: the commands</h2>
    <div class="tbl rs-su"></div>
    <div class="prose"><p>Each one rewrites the per-atom arrays between runs. The three below are the ones whose results most often surprise people; for each, the rule is read from the source and then applied by hand to every atom.</p></div>
  </section>

  <section class="block">
    <h2>4 · delete_atoms overlap: which one goes?</h2>
    <div class="prose"><p>${gh('src/delete_atoms.cpp', 350)} loops over the owned atoms in their storage order with a full neighbor list. Atom I is deleted if some neighbour J within the cutoff qualifies. An owned J qualifies if it has not been deleted already. A ghost J qualifies only if I has the smaller tag, so that two processors agree. The first atom of a close pair in storage order goes, the second stays. Before any run, storage order is creation order. One more default to know: <code>compress yes</code> renumbers the survivors 1…N afterwards, so to see which IDs went, add <code>compress no</code>.</p></div>
    <pre class="lmp-mini">${esc(DEL)}</pre>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Delete, and predict every deletion</button><span class="rs-stc"></span></div>
    <div class="prose rs-csum"></div>
    <div class="prose"><p>The same rule applied to page F.4's mistake:</p></div>
    <pre class="lmp-mini">${esc(DEL2)}</pre>
    <div class="row"><button class="btn sm primary" data-k="c2">▶ Which atom is deleted?</button></div>
    <div class="prose rs-c2sum"></div>
  </section>

  <section class="block">
    <h2>5 · replicate: new IDs, new positions</h2>
    <div class="prose"><p>${gh('src/replicate.cpp', 179)} first unwraps every atom with its image flags, then for each copy (i<sub>x</sub>, i<sub>y</sub>, i<sub>z</sub>) adds i<sub>x</sub>·L<sub>x</sub> (etc.), wraps the result into the new box (${gh('src/domain.cpp')} remap: add or subtract the new period until inside) and offsets the tag by (i<sub>z</sub>·n<sub>y</sub>n<sub>x</sub> + i<sub>y</sub>·n<sub>x</sub> + i<sub>x</sub>)·maxtag.</p></div>
    <pre class="lmp-mini">${esc(REP)}</pre>
    <div class="row"><button class="btn sm primary" data-k="d">▶ Replicate and check 2048 atoms</button><span class="rs-std"></span></div>
    <div class="prose rs-dsum"></div>
  </section>

  <section class="block">
    <h2>6 · change_box remap: positions follow the box</h2>
    <div class="prose"><p>${gh('src/change_box.cpp', 243)}: <code>scale</code> keeps the midpoint, lo′ = mid + s·(lo − mid), hi′ = mid + s·(hi − mid). <code>remap</code> converts each atom to fractional coordinates with the old box, λ = (x − lo)·(1/L), and back with the new one, x′ = L′·λ + lo′ (${gh('src/domain.cpp', 2314)}). Two roundings, both reproducible:</p></div>
    <pre class="lmp-mini">${esc(CB)}</pre>
    <div class="row"><button class="btn sm primary" data-k="e">▶ Resize and predict every coordinate</button><span class="rs-ste"></span></div>
    <div class="prose rs-esum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> You restart an NPT run, but name the fix <code>npt1</code> instead of <code>1</code>. What happens to the box at the first steps?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The barostat starts with ω̇ = 0 and η = 0 instead of its stored values: the box velocity is suddenly zero and the thermostat chain is reset. The log ends with "Unused restart file global fix info". Usually you see a small jolt in pressure and volume, and Econserve jumps by the stored thermostat and barostat energy. Keep fix IDs fixed across restarts.</div></div>
    <div class="ex"><p><b>2.</b> After delete_atoms overlap, the system in section 4 has gaps in its IDs. Does that matter?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Not for LAMMPS; IDs need not be contiguous. It matters for tools that assume 1…N, and for replicate, which offsets by the largest tag rather than the count. <code>reset_atoms id</code> renumbers 1…N (with <code>sort yes</code> it also orders them spatially).</div></div>
    <div class="ex"><p><b>3.</b> Section 4 kept atom 257, the one you added by mistake, and deleted lattice atom 1. How would you make it delete atom 257 instead?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Use the two groups: <code>group new id 257</code>, then <code>delete_atoms overlap 0.5 new all</code>. Only atoms in the first group are candidates for deletion (the check on group1bit for I), so the lattice atom stays. In general, put the atoms you are inserting in group 1.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true;
    $('.rs-st').innerHTML = `<table>${head(['', 'in a binary restart', 'note'])}${STORED.map(row).join('')}</table>`;
    $('.rs-su').innerHTML = `<table>${head(['command', 'what it does'])}${SURGERY.map(row).join('')}</table>`;
    const busy = (b, st, t) => { b.disabled = true; st.textContent = t; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.rs-sta'); busy(b, st, ' Same ID…');
      const R1 = await runLammps(NATIVE.inr, { files: { 'r.bin': bin() } }); if (!alive) return; st.textContent = ' New ID…';
      const R2 = await runLammps(NATIVE.inr.replace('fix 1 all nvt', 'fix 2 all nvt').replace('f_1', 'f_2'), { files: { 'r.bin': bin() } }); if (!alive) return; b.disabled = false; st.textContent = '';
      if (R1.error || R2.error) { st.innerHTML = `<span class="lmp-err">${esc(R1.error || R2.error)}</span>`; return; }
      const cmp = (R, N) => { let s = 0, n = 0; R.thermo.rows.forEach((r, i) => r.forEach((v, c) => { n++; if (N[i] && String(v) === N[i][c]) s++; })); return [s, n]; };
      const [s1, n1] = cmp(R1, NATIVE.r), [s2, n2] = cmp(R2, NATIVE.r2), m1 = R1.lines.some(l => /Resetting global fix info/.test(l)), m2 = R2.lines.some(l => /Unused restart file global fix info/.test(l));
      const tab = (R, lab) => `<tr><th colspan="6">${lab}</th></tr>${R.thermo.rows.map(r => row(r)).join('')}`;
      const jump = Number(NATIVE.r2[0][5]) - Number(NATIVE.r[0][5]);
      $('.rs-asum').innerHTML = `<div class="tbl"><table>${head(R1.thermo.cols)}${tab(R1, 'fix 1 (same ID)')}${tab(R2, 'fix 2 (new ID)')}</table></div>
      <p>Same ID: ${s1}/${n1} strings identical to the native continuation ${s1 === n1 ? OK : BAD}, and the browser log says "Resetting global fix info from restart file" ${m1 ? OK : BAD}. New ID: ${s2}/${n2} identical to native ${s2 === n2 ? OK : BAD}, "Unused restart file global fix info" ${m2 ? OK : BAD}. A binary file written on x86-64 Linux, read by WebAssembly: same format, same doubles.</p>
      <p>With the new ID, f_2 starts at 0 and the trajectory changes from the first step (Temp at step 150: ${R2.thermo.rows[1][1]} instead of ${R1.thermo.rows[1][1]}), because the chain starts from η̇ = 0 instead of the stored state. Econserve at step 100 jumps by ${s17(jump)}, which is minus the stored f_1, ${NATIVE.r[0][4]}, up to the last bit of the subtraction. The "conserved" quantity is only conserved while the thermostat keeps its memory.</p>`;
    };

    $('[data-k=b]').onclick = () => {
      const u = bin(), dv = new DataView(u.buffer), txt = new TextDecoder('latin1').decode(u);
      const p = txt.indexOf('nvt\u0000'), size = dv.getInt32(p + 4, true), v = []; for (let k = 0; k < size / 8; k++) v.push(dv.getFloat64(p + 8 + 8 * k, true));
      const idlen = dv.getInt32(p - 4 - 2 - 4, true);
      const kT = 1.0, L = 765, w = 2.0, Q0 = L * 1.0 * 1.0 / (w * w), Qk = 1.0 * 1.0 / (w * w), ke = L * 1.0 * 1.0;
      let E = ke * v[2] + 0.5 * Q0 * v[5] * v[5]; for (let k = 1; k < 3; k++) E += kT * v[2 + k] + 0.5 * Qk * v[5 + k] * v[5 + k];
      const f1 = E / 256, nat = Number(NATIVE.w[NATIVE.w.length - 1][4]);
      $('.rs-bsum').innerHTML = `<p>At byte ${p} of r.bin: the style string "nvt" (after the ID "1", ${idlen} bytes with its terminator), an int ${size} (= ${size / 8} doubles), then ${v.map(x => s17(x)).join(', ')}. That is tstat_flag = ${v[0]}, mtchain = ${v[1]}, η = (${v.slice(2, 5).map(x => x.toFixed(6)).join(', ')}), η̇ = (${v.slice(5, 8).map(x => x.toFixed(6)).join(', ')}), pstat_flag = ${v[8]}.</p>
      <p>Thermostat energy / N = ${s17(f1)}; the f_1 native LAMMPS printed at step 100 before writing the file: ${s17(nat)} ${f1 === nat ? OK + ' identical' : BAD}. And it is the value the same-ID continuation printed after reading it. Nothing else of fix 1 is in the file: T<sub>start</sub>, T<sub>damp</sub> and the group come from the new fix command, so a restart can also change them.</p>`;
    };

    $('[data-k=c]').onclick = async e => {
      const b = e.currentTarget, st = $('.rs-stc'); busy(b, st, ' Running…');
      const R = await runLammps(DEL, { read: ['before.dump', 'after.dump'] }); if (!alive) return; b.disabled = false;
      if (R.error || !R.files['after.dump']) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no dump')}</span>`; return; } st.textContent = '';
      const B = readDump(R.files['before.dump']), A = readDump(R.files['after.dump']), Lb = 8, cut2 = 0.8 * 0.8, ids = [...B.at.keys()].sort((p, q) => p - q);
      // the loop of delete_atoms.cpp:350-412, storage order = creation order = tag order
      const del = new Set(); let ghostPairs = 0;
      for (const i of ids) {
        const xi = B.at.get(i);
        for (const j of ids) {
          if (j === i) continue;
          const xj = B.at.get(j); let ghost = false; const d = [0, 1, 2].map(c => { let dd = xi[c] - xj[c]; if (dd > Lb / 2) { dd -= Lb; ghost = true; } else if (dd < -Lb / 2) { dd += Lb; ghost = true; } return dd; });
          if (d[0] * d[0] + d[1] * d[1] + d[2] * d[2] >= cut2) continue;
          if (!ghost) { if (del.has(j)) continue; } else { ghostPairs++; if (i > j) continue; }
          del.add(i); break;
        }
      }
      const surv = ids.filter(i => !del.has(i)), got = [...A.at.keys()];
      const same = surv.length === got.length && surv.every((x, k) => x === got[k]);
      const msg = R.lines.find(l => /Deleted \d+ atoms/.test(l)) || '';
      $('.rs-csum').innerHTML = `<p>LAMMPS: "${esc(msg.trim())}". By hand, walking the atoms in creation order: ${del.size} deleted, ${surv.length} remain, and the surviving IDs are ${same ? 'identical, ID by ID' : 'different'} ${same ? OK : BAD}. Every remaining pair is now at least 0.8 apart. Note that the outcome depends on the storage order. After a run has sorted the atoms spatially, the same command on the same positions can keep the other atom of a pair, and the count can change too (three atoms in a row, each close only to the middle one: if the middle atom comes first in storage, only it is deleted; if an end atom comes first, two are).</p>`;
    };

    $('[data-k=c2]').onclick = async () => {
      const R = await runLammps(DEL2, { read: ['after.dump'] }); if (!alive) return;
      const A = readDump(R.files['after.dump'] || '\n\n\n0\n\n\n\n\n\n'), has1 = A.at.has(1), has257 = A.at.has(257);
      $('.rs-c2sum').innerHTML = `<p>${A.n} atoms remain; atom 1 ${has1 ? 'kept' : 'deleted'}, atom 257 ${has257 ? 'kept' : 'deleted'} ${!has1 && has257 ? OK : BAD}. Atom 1 is first in storage order, its neighbour 257 has not been deleted, so atom 1 goes. The lattice site is now occupied by the atom you misplaced, 0.084 σ away from it. Usually harmless, but not what most people expect; exercise 3 shows how to choose.</p>`;
    };

    $('[data-k=d]').onclick = async e => {
      const b = e.currentTarget, st = $('.rs-std'); busy(b, st, ' Running…');
      const R = await runLammps(REP, { read: ['a.dump', 'b.dump'] }); if (!alive) return; b.disabled = false;
      if (R.error || !R.files['b.dump']) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no dump')}</span>`; return; } st.textContent = '';
      const A = readDump(R.files['a.dump']), B = readDump(R.files['b.dump']), Lo = A.box.map(r => r[1] - r[0]), Pn = B.box.map(r => r[1] - r[0]), lo = B.box.map(r => r[0]), hi = B.box.map(r => r[1]);
      let ok = 0, n = 0, img = 0; const maxtag = 256;
      for (let ix = 0; ix < 2; ix++) for (let iy = 0; iy < 2; iy++) for (let iz = 0; iz < 2; iz++) {
        const off = iz * 2 * 2 * maxtag + iy * 2 * maxtag + ix * maxtag, sh = [ix, iy, iz];
        for (const [t, w] of A.at) {
          if (w[3] || w[4] || w[5]) img++;
          const x = [0, 1, 2].map(c => { let v = w[c] + w[3 + c] * Lo[c]; v = v + sh[c] * Lo[c]; while (v < lo[c]) v += Pn[c]; while (v >= hi[c]) v -= Pn[c]; return Math.max(v, lo[c]); });
          const got = B.at.get(t + off); n++;
          if (got && got[0] === x[0] && got[1] === x[1] && got[2] === x[2]) ok++;
        }
      }
      const pe = R.thermo.rows, pe0 = Number(pe[pe.length - 1][1]);
      const i2 = R.lines.findIndex(l => /^\s*Step\s/.test(l) && false);
      const all = []; let on = false; for (const l of R.lines) { if (/^\s*Step\s/.test(l)) { on = true; continue; } if (on && /^\s*\d+\s/.test(l)) all.push(l.trim().split(/\s+/)); else if (on && /^Loop/.test(l)) on = false; }
      const peA = Number(all[all.length - 2][1]), peB = Number(all[all.length - 1][1]);
      $('.rs-dsum').innerHTML = `<p>${ok} of ${n} atoms found under the predicted ID at the predicted position, to the last bit ${ok === n ? OK : BAD}. ${img / 8} of the 256 originals had non-zero image flags after 300 steps; unwrapping them first is what keeps each copy contiguous.</p>
      <p>Potential energy per atom, last step before replicate: ${s17(peA)}; after replicate (run 0): ${s17(peB)} (difference ${(peB - peA).toExponential(1)}). Each copy has the same neighbours at the same distances, so the energy per atom is the same; any difference is the summation order, since setup re-sorts 2048 atoms and rebuilds the neighbour list. A replicated box is the same state, not just a similar one. Its trajectory will still diverge from the small box's (page F.3), and it can now hold wavelengths twice as long.</p>`;
      void pe0; void i2;
    };

    $('[data-k=e]').onclick = async e => {
      const b = e.currentTarget, st = $('.rs-ste'); busy(b, st, ' Running…');
      const R = await runLammps(CB, { read: ['a.dump', 'b.dump'] }); if (!alive) return; b.disabled = false;
      if (R.error || !R.files['b.dump']) { st.innerHTML = `<span class="lmp-err">${esc(R.error || 'no dump')}</span>`; return; } st.textContent = '';
      const A = readDump(R.files['a.dump']), B = readDump(R.files['b.dump']);
      // box: x scale 1.1 about the midpoint, then y final 0 7
      const [xlo, xhi] = A.box[0], mid = 0.5 * (xlo + xhi), nxlo = mid + 1.1 * (xlo - mid), nxhi = mid + 1.1 * (xhi - mid);
      const nb = [[nxlo, nxhi], [0, 7], A.box[2]], boxOK = nb.every((r, c) => r[0] === B.box[c][0] && r[1] === B.box[c][1]);
      const hinv = A.box.map(r => 1.0 / (r[1] - r[0])), h = nb.map(r => r[1] - r[0]);
      let ok = 0, n = 0, wrapped = 0;
      for (const [t, w] of A.at) {
        // change_box.cpp:211 first wraps atoms into the old box (Domain::pbc), then remaps
      const lo0 = A.box.map(r => r[0]), hi0 = A.box.map(r => r[1]), p0 = A.box.map(r => r[1] - r[0]);
      const v = [0, 1, 2].map(c => { let q = w[c]; if (q < lo0[c]) q += p0[c]; if (q >= hi0[c]) { q -= p0[c]; q = Math.max(q, lo0[c]); } if (q !== w[c]) wrapped++; return q; });
      const lam = [0, 1, 2].map(c => hinv[c] * (v[c] - A.box[c][0]));
        const x = [0, 1, 2].map(c => h[c] * lam[c] + nb[c][0]);
        const g = B.at.get(t); n += 3; ok += (g[0] === x[0]) + (g[1] === x[1]) + (g[2] === x[2]);
      }
      const t1 = A.at.get(1), g1 = B.at.get(1);
      $('.rs-esum').innerHTML = `<p>New box: x [${s17(nxlo)}, ${s17(nxhi)}], y [0, 7], z unchanged; the dump header agrees ${boxOK ? OK : BAD}. Atom 1: x ${s17(t1[0])} → ${s17(g1[0])}, y ${s17(t1[1])} → ${s17(g1[1])}. All ${ok} of ${n} coordinates predicted to the last bit ${ok === n ? OK : BAD}, including z, which goes through λ and back with the same box. One detail was needed for that: ${wrapped} coordinates were slightly outside the box in the dump (atoms are only wrapped at reneighboring), and change_box first puts them back inside (${gh('src/change_box.cpp', 211)}, Domain::pbc), then remaps; remapping first and wrapping afterwards gives different last bits. Without <code>remap</code> the box would change and the atoms would stay where they were: a gap in x, and in y (box 6.72 → 7) a slab of empty space, or overlaps if you shrink.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
