// Research · F.7: rerun and post-processing. A trajectory recomputed from its dump: what comes back to the last
// bit (an RDF), what to round-off (energies), what the default dump format costs (predicted from the forces), a
// different potential evaluated on old frames (checked against an integral over g(r)), and the image-flag trap.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const N = 500, T0 = 0.72;

const MD = `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 5 0 5 0 5
create_box 1 box
create_atoms 1 box
mass 1 1.0
velocity all create 1.44 87287 loop geom
pair_style lj/cut 2.5
pair_coeff 1 1 1.0 1.0 2.5
fix 1 all nve
compute r all rdf 50
compute m all msd
fix ar all ave/time 100 1 100 c_r[*] file rdf.md mode vector format ' %.17g'
thermo_style custom step temp pe press c_m[4]
thermo_modify format float %.17g
thermo 100
dump d all custom 100 t17.dump id type x y z ix iy iz vx vy vz
dump_modify d sort id format float %.17g     # all digits
dump g all custom 100 tg.dump id type x y z ix iy iz vx vy vz
dump_modify g sort id                        # default format
run 1000`;
const RE = (file, extra = '', pair = 'lj/cut 2.5', fields = 'x y z ix iy iz vx vy vz', tag = 're') => `units lj
atom_style atomic
lattice fcc 0.8442
region box block 0 5 0 5 0 5
create_box 1 box
mass 1 1.0
pair_style ${pair}
pair_coeff 1 1 1.0 1.0 ${pair.split(' ')[1]}
compute r all rdf 50
compute m all msd
fix ar all ave/time 100 1 100 c_r[*] file rdf.${tag} mode vector format ' %.17g'
thermo_style custom step temp pe press c_m[4]
thermo_modify format float %.17g
thermo 100${extra}
rerun ${file} dump ${fields} add keep`;
const thermoAll = lines => { const out = []; let on = false; for (const l of lines) { if (/^\s*Step\s/.test(l)) { on = true; continue; } if (on && /^\s*-?\d/.test(l)) out.push(l.trim().split(/\s+/)); else if (on && /^Loop/.test(l)) on = false; } return out; };
const frames = txt => { const L = txt.split('\n'), F = []; let i = 0; while (i < L.length) { if (L[i].startsWith('ITEM: TIMESTEP')) { const st = +L[i + 1], n = +L[i + 3], box = [5, 6, 7].map(k => L[i + k].trim().split(/\s+/).map(Number)), at = []; for (let a = 0; a < n; a++) at.push(L[i + 9 + a].trim().split(/\s+/).map(Number)); F.push({ st, box, at }); i += 9 + n; } else i++; } return F; };

export default {
  id: 'rerun',
  title: 'Rerun and post-processing',
  minutes: 14,
  summary: 'A dump file is a trajectory you can run again. rerun reads it frame by frame, recomputes forces and energies, possibly with a different potential, and feeds any compute you like, without moving an atom. You check what comes back exactly (an RDF, bit for bit), what only to round-off (the energy), and what the default dump format costs (predicted to two digits from the forces). Then you evaluate a longer cutoff on old frames, check the result against an integral over g(r), and fall into the image-flag trap with an MSD.',
  learn: [
    'What rerun does per frame (setup_minimal: neighbours, forces, output)',
    'Which results are exact, which round-off, which format-limited',
    'ΔU of a changed potential on stored frames, and a check via g(r)',
    'Why unwrapped analysis needs ix iy iz in the dump',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a run without time steps</h2>
    <div class="prose"><p>${gh('src/rerun.cpp', 159)}: for each snapshot, read_dump replaces the positions (and whatever else is listed), the timestep is set to the snapshot's, and <code>integrate->setup_minimal(1)</code> (line 170) rebuilds the neighbor list, computes forces and energies, and runs output: thermo, dumps, fix ave/time. No fix moves any atom. What it cannot recover is what the dump did not store: velocities (Temp is 0 without vx vy vz), image flags, and digits.</p></div>
    <pre class="lmp-mini">${esc(MD)}</pre>
    <div class="prose"><p>Then rerun it, with nothing but the dump:</p></div>
    <pre class="lmp-mini">${esc(RE('t17.dump'))}</pre>
    <div class="prose"><p><code>add keep</code> creates the atoms from the first frame with their IDs; without it, read_dump only updates atoms that already exist (none here), and every thermo line is 0.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run, then rerun the 17-digit dump</button><span class="rr-sta"></span></div>
    <div class="prose rr-asum"></div>
  </section>

  <section class="block">
    <h2>2 · Observe: what the default dump format costs</h2>
    <div class="prose"><p>The second dump used the default format, %g: 6 significant digits, so a position near 5 σ is rounded to about ±5 × 10⁻⁶. To first order the energy changes by ΔU ≈ −Σ<sub>i</sub> F<sub>i</sub>·δx<sub>i</sub>. Both factors are known exactly: δx from the two dumps, F from a 17-digit rerun that also dumps the forces.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Rerun the %g dump and predict its error</button><span class="rr-stb"></span></div>
    <div class="prose rr-bsum"></div>
    <div class="tbl rr-btab"></div>
  </section>

  <section class="block">
    <h2>3 · Modify the potential, keep the frames</h2>
    <div class="prose"><p>The classic use of rerun: how would the energy differ with another potential, on the same configurations? Rerun the 17-digit dump with <code>pair_style lj/cut 3.0</code>. The difference per frame, ΔU = U<sub>3.0</sub> − U<sub>2.5</sub>, is the energy of all pairs between 2.5 and 3.0 σ. On average it must equal 2πρ∫<sub>2.5</sub><sup>3.0</sup> g(r) u(r) r² dr per atom, with g(r) measured on the same frames.</p></div>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Rerun with cutoff 3.0 and integrate g(r)</button><span class="rr-stc"></span></div>
    <div class="prose rr-csum"></div>
    <div class="tbl rr-ctab"></div>
  </section>

  <section class="block">
    <h2>4 · The image-flag trap</h2>
    <div class="prose"><p>Positions in a dump are wrapped into the box; the image flags say how many box lengths to add back. Rerun the same dump reading only <code>x y z</code>: every atom starts with image flags 0, and an atom that crossed the boundary appears to have jumped by a box length. What exactly happens is a little stranger, as you will see.</p></div>
    <div class="row"><button class="btn sm primary" data-k="d">▶ MSD with and without ix iy iz</button><span class="rr-std"></span></div>
    <div class="prose rr-dsum"></div>
    <div class="tbl rr-dtab"></div>
  </section>

  <section class="block prose">
    <h2>5 · Rerun or a script?</h2>
    <ul><li><b>rerun</b> when the analysis needs LAMMPS's own physics: energies or forces with another potential, per-atom stress, anything a compute already does correctly (RDF with the right normalization, MSD with image flags, order parameters). Options <code>first/last/every/skip/start/stop</code> select frames.</li>
    <li><b>A script</b> (Python with numpy, MDAnalysis, OVITO) when the analysis is geometric and you want to iterate quickly. The dump is a text table; page D.18 showed the columns.</li>
    <li>In both cases <b>write the dump with the digits you will need</b> (<code>dump_modify format float %.17g</code> or a binary/compressed format), and <b>with ix iy iz</b> (or xu yu zu) if anything depends on unwrapped positions.</li></ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why is the RDF from rerun identical to the last bit, while the energy is not?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The RDF is a histogram: integer counts per bin, then the same normalization. Adding counts is exact in any order, and with identical positions every pair lands in the same bin. The energy is a sum of ~19 000 floating-point terms, and rerun visits the atoms in the dump's (ID) order while the MD run had them sorted spatially: same terms, different order, different last bits (page C.6).</div></div>
    <div class="ex"><p><b>2.</b> In section 3, can you use ΔU to say what the free energy change of the cutoff is?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Yes, by free-energy perturbation: ΔA = −kT ln⟨exp(−ΔU/kT)⟩ over frames of the old potential. It works when ΔU fluctuates by much less than kT, so that the old frames are representative of the new potential. Here ΔU (total, for 500 atoms) fluctuates by about the σ in the table, and the estimate is close to ⟨ΔU⟩ − σ²/2kT. For a larger change (a new charge model, say) the old frames stop being representative and the estimate is dominated by a few frames; then you use intermediate potentials, or run the new potential.</div></div>
    <div class="ex"><p><b>3.</b> Rerun a dump that has only <code>id type x y z ix iy iz</code>. Temp is 0 on every line, but Press is not. What is it then?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Pressure is P = (N k T + W)/V with W the virial, ⅓Σ r·F. Without velocities the kinetic term is zero and Press is the virial part alone, the "configurational" pressure. Add N k T/V with the known temperature by hand, or dump vx vy vz, if you need the full pressure from a rerun.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    let alive = true, md = null;
    const need = async st => { if (md) return md; st.textContent = ' Running the MD first…'; const R = await runLammps(MD, { read: ['t17.dump', 'tg.dump', 'rdf.md'] }); if (R.error) throw new Error(R.error); md = { th: thermoAll(R.lines), t17: R.files['t17.dump'], tg: R.files['tg.dump'], rdf: R.files['rdf.md'] }; return md; };
    const go = (k, fn) => { $(`[data-k=${k}]`).onclick = async e => { const b = e.currentTarget, st = $(`.rr-st${k}`); b.disabled = true; try { await fn(st); st.textContent = ''; } catch (err) { st.innerHTML = `<span class="lmp-err">${esc(err.message)}</span>`; } b.disabled = false; }; };

    go('a', async st => {
      const M = await need(st); if (!alive) return; st.textContent = ' Rerun…';
      const R = await runLammps(RE('t17.dump'), { files: { 't17.dump': M.t17 }, read: ['rdf.re'] }); if (!alive) return; if (R.error) throw new Error(R.error);
      const t = thermoAll(R.lines), cols = ['Temp', 'PotEng', 'Press', 'MSD'];
      const res = cols.map((c, k) => { let same = 0, mx = 0; M.th.forEach((r, i) => { const a = Number(r[k + 1]), b = Number(t[i][k + 1]); if (r[k + 1] === t[i][k + 1]) same++; if (a) mx = Math.max(mx, Math.abs(b / a - 1)); }); return [c, `${same}/${M.th.length}`, mx ? mx.toExponential(1) : '0']; });
      const rdfSame = M.rdf.split('\n').filter(l => !l.startsWith('#')).join('\n') === (R.files['rdf.re'] || '').split('\n').filter(l => !l.startsWith('#')).join('\n');
      $('.rr-asum').innerHTML = `<div class="tbl"><table>${head(['quantity', 'identical strings', 'largest relative difference'])}${res.map(row).join('')}${row(['RDF, 11 frames × 50 bins', rdfSame ? 'all' : 'no', rdfSame ? '0' : '?'])}</table></div><p>The RDF comes back identical, every digit of every bin ${rdfSame ? OK : BAD}. Energy, pressure, temperature and MSD agree to round-off (largest relative difference above, at the level of 10⁻¹⁵), with only a few lines identical. Rerun recomputes them from the same doubles in a different order: atoms in ID order as read, not sorted spatially as in the MD run. For an analysis the difference is irrelevant; for a regression test (page F.3) it is a different number.</p>`;
    });

    go('b', async st => {
      const M = await need(st); if (!alive) return; st.textContent = ' Rerun %g dump…';
      const Rg = await runLammps(RE('tg.dump'), { files: { 'tg.dump': M.tg } }); if (!alive) return; st.textContent = ' Rerun 17-digit dump with forces…';
      const Rf = await runLammps(RE('t17.dump', '\ndump f all custom 100 f.dump id fx fy fz\ndump_modify f sort id format float %.17g'), { files: { 't17.dump': M.t17 }, read: ['f.dump'] }); if (!alive) return;
      if (Rg.error || Rf.error) throw new Error(Rg.error || Rf.error);
      const tg = thermoAll(Rg.lines), t17 = thermoAll(Rf.lines), F17 = frames(M.t17), Fg = frames(M.tg), FF = frames(Rf.files['f.dump']);
      const rows = []; let maxdx = 0, worst = 0, jumps = 0; const Lb = F17[0].box[0][1] - F17[0].box[0][0], uc = 4 * (2.5 ** -12 - 2.5 ** -6);
      const r2 = (P, a, b) => { let s = 0; for (let c = 0; c < 3; c++) { let d = P[a][2 + c] - P[b][2 + c]; d -= Lb * Math.round(d / Lb); s += d * d; } return s; };
      F17.forEach((fr, k) => {
        let s = 0; for (let a = 0; a < N; a++) { const x = fr.at[a], g = Fg[k].at[a], f = FF[k].at[a]; for (let c = 0; c < 3; c++) { const dx = g[2 + c] - x[2 + c]; maxdx = Math.max(maxdx, Math.abs(dx)); s -= f[1 + c] * dx; } }
        // pairs that cross r_c = 2.5 because of the rounding: the energy jumps by -u(r_c) (in) or +u(r_c) (out)
        let nin = 0, nout = 0; for (let a = 0; a < N; a++) for (let b = a + 1; b < N; b++) { const p = r2(fr.at, a, b) < 6.25, q = r2(Fg[k].at, a, b) < 6.25; if (q && !p) nin++; if (p && !q) nout++; }
        jumps += nin + nout;
        const pred = (s + (nin - nout) * uc) / N, obs = Number(tg[k][2]) - Number(t17[k][2]);
        if (k > 0) worst = Math.max(worst, Math.abs(pred - obs) / Math.abs(obs));
        rows.push([fr.st, obs.toExponential(3), (s / N).toExponential(3), nin || nout ? `${nin} in, ${nout} out` : '', pred.toExponential(3), k ? (pred / obs).toFixed(4) : '(lattice, F = 0)']);
      });
      const dT = Math.max(...tg.map((r, k) => Math.abs(Number(r[1]) / Number(t17[k][1]) - 1)));
      $('.rr-btab').innerHTML = `<table>${head(['step', 'ΔPE per atom, observed', '−ΣF·δx / N', 'pairs crossing r<sub>c</sub>', 'predicted', 'ratio'])}${rows.map(row).join('')}</table>`;
      $('.rr-bsum').innerHTML = `<p>The %g positions differ from the exact ones by up to ${maxdx.toExponential(1)} σ. The resulting energy error per atom, ~10⁻⁶ to 10⁻⁵, is ten orders of magnitude above round-off, and it is predicted frame by frame to within ${(100 * worst).toFixed(1)} % ${worst < 0.01 ? OK : BAD} by two terms: the first-order change −ΣF·δx, and ${jumps === 1 ? 'one pair' : jumps + ' pairs'} that the rounding moved across the cutoff, ${jumps === 1 ? 'which changes' : 'each changing'} the energy by the full u(r<sub>c</sub>) = ${uc.toFixed(6)} because lj/cut is not shifted. At step 800 that one pair outweighs everything else almost forty times. The rest is second order, ½δx·H·δx, which is also why the lattice frame (zero forces) has a tiny error left. It is not random noise to average away: it is a definite, computable error that depends on the forces. Temperature from %g velocities: relative error up to ${dT.toExponential(1)}. Whether that matters depends on the analysis: not for an RDF, yes for energy differences between potentials, which is what section 3 does.</p>`;
    });

    go('c', async st => {
      const M = await need(st); if (!alive) return; st.textContent = ' Rerun lj/cut 3.0…';
      const R3 = await runLammps(RE('t17.dump', '', 'lj/cut 3.0', 'x y z ix iy iz vx vy vz', 'c3').replace('compute r all rdf 50', 'compute r all rdf 120 cutoff 3.0'), { files: { 't17.dump': M.t17 }, read: ['rdf.c3'] }); if (!alive) return; st.textContent = ' Rerun lj/cut 2.5…';
      const R2 = await runLammps(RE('t17.dump'), { files: { 't17.dump': M.t17 } }); if (!alive) return;
      if (R3.error || R2.error) throw new Error(R3.error || R2.error);
      const t3 = thermoAll(R3.lines), t2 = thermoAll(R2.lines), dU = t3.map((r, k) => Number(r[2]) - Number(t2[k][2])).slice(1);   // skip the lattice frame
      const mean = dU.reduce((a, b) => a + b, 0) / dU.length, sd = Math.sqrt(dU.reduce((a, b) => a + (b - mean) ** 2, 0) / (dU.length - 1));
      // g(r) from fix ave/time: blocks of 120 rows: bin r g(r) coord; average frames 100..1000
      const L = R3.files['rdf.c3'].split('\n').filter(l => l && !l.startsWith('#')), g = new Float64Array(120), r = new Float64Array(120); let nb = 0;
      for (const l of L) { const w = l.trim().split(/\s+/).map(Number); if (w.length === 2) { nb++; continue; } if (nb >= 2) { g[w[0] - 1] += w[2] / 10; r[w[0] - 1] = w[1]; } }
      const rho = N / (5 * Math.cbrt(4 / 0.8442)) ** 3, dr = 3.0 / 120, u = x => 4 * (x ** -12 - x ** -6);
      let I = 0; for (let b = 0; b < 120; b++) if (r[b] > 2.5) I += g[b] * u(r[b]) * r[b] * r[b] * dr;
      const pred = 2 * Math.PI * rho * I, uni = 2 * Math.PI * rho * (() => { let s = 0; for (let b = 0; b < 120; b++) if (r[b] > 2.5) s += u(r[b]) * r[b] * r[b] * dr; return s; })();
      const lse = (() => { const x = dU.map(v => -v * N / T0), m = Math.max(...x); return -(T0 / N) * (m + Math.log(x.reduce((a, b) => a + Math.exp(b - m), 0) / x.length)); })();
      $('.rr-ctab').innerHTML = `<table>${head(['step', 'U(2.5) per atom', 'U(3.0) per atom', 'ΔU per atom'])}${t3.slice(1).map((rr, k) => row([rr[0], Number(t2[k + 1][2]).toFixed(9), Number(rr[2]).toFixed(9), dU[k].toFixed(6)])).join('')}</table>`;
      $('.rr-csum').innerHTML = `<p>Over the ten liquid frames, ⟨ΔU⟩ = ${mean.toFixed(5)} per atom, frame-to-frame σ = ${sd.toFixed(5)}. From g(r) on the same frames (120 bins of 0.025 σ, midpoint rule), 2πρ∫ g u r² dr = ${pred.toFixed(5)} per atom: ${(100 * Math.abs(pred / mean - 1)).toFixed(1)} % apart ${Math.abs(pred / mean - 1) < 0.02 ? OK : BAD}, the size of the midpoint-rule error with this bin width. With g = 1 (the usual tail-correction assumption) the integral would be ${uni.toFixed(5)}: ${(100 * Math.abs(uni / pred - 1)).toFixed(0)} % off, because the liquid's structure between 2.5 and 3.0 σ (g below 1 near its second minimum) is not negligible. As a free-energy perturbation (T ≈ ${T0}, 500 atoms), −kT ln⟨e<sup>−ΔU/kT</sup>⟩ = ${lse.toFixed(5)} per atom, close to ⟨ΔU⟩ because the fluctuations are small compared with kT.</p>`;
    });

    go('d', async st => {
      const M = await need(st); if (!alive) return; st.textContent = ' Rerun without image flags…';
      const Rn = await runLammps(RE('t17.dump', '', 'lj/cut 2.5', 'x y z vx vy vz'), { files: { 't17.dump': M.t17 } }); if (!alive) return;
      if (Rn.error) throw new Error(Rn.error);
      const tn = thermoAll(Rn.lines), F = frames(M.t17), lo = F[0].box[0][0], hi = F[0].box[0][1], Lb = hi - lo, x0 = F[0].at;
      // read_dump without ix iy iz keeps each atom's current image flags (read_dump.cpp:843) and remaps atoms that lie
      // outside the box into it, counting the crossing (line 1144); the flags persist from frame to frame
      const img = Array.from({ length: N }, () => [0, 0, 0]); let ref = null;
      const rows = F.map((fr, k) => {
        let sw = 0, su = 0, crossed = 0; const cur = [];
        for (let a = 0; a < N; a++) { const p = fr.at[a]; const v = [0, 1, 2].map(c => { let q = p[2 + c]; while (q < lo) { q += Lb; img[a][c]--; } while (q >= hi) { q -= Lb; img[a][c]++; } return Math.max(q, lo) + img[a][c] * Lb; }); cur.push(v); }
        if (!ref) ref = cur;
        for (let a = 0; a < N; a++) { const p = fr.at[a], q = x0[a]; let any = false; for (let c = 0; c < 3; c++) { const dw = cur[a][c] - ref[a][c], du = (p[2 + c] + p[5 + c] * Lb) - (q[2 + c] + q[5 + c] * Lb); sw += dw * dw; su += du * du; if (p[5 + c] !== q[5 + c]) any = true; } crossed += any; }
        return [fr.st, Number(md.th[k][4]).toFixed(5), (su / N).toFixed(5), Number(tn[k][4]).toFixed(5), (sw / N).toFixed(5), crossed];
      });
      const okU = rows.every((r, k) => Math.abs(Number(r[2]) - Number(md.th[k][4])) < 1e-5), okW = rows.every(r => Math.abs(Number(r[4]) - Number(r[3])) < 1e-5);
      $('.rr-dtab').innerHTML = `<table>${head(['step', 'MSD, MD run', 'by hand, unwrapped', 'MSD, rerun x y z only', 'by hand, as read_dump does it', 'atoms that crossed'])}${rows.map(row).join('')}</table>`;
      $('.rr-dsum').innerHTML = `<p>With image flags the rerun MSD is the MD run's, and both equal ⟨|x + n·L − x₀ − n₀·L|²⟩ computed here ${okU ? OK : BAD}. Without them the rerun MSD is neither the true one nor that of the wrapped coordinates. read_dump keeps each atom's previous image flags (${gh('src/read_dump.cpp', 843)}) and moves atoms that lie outside the box back in, counting that as a crossing (${gh('src/read_dump.cpp', 1144)}). So only atoms caught just outside the box at a dump step get credit, and the credit accumulates. Reproduced here frame by frame ${okW ? OK : BAD}. Every atom that crossed a face without being caught contributes about L² = ${(Lb * Lb).toFixed(1)} σ², and the "MSD" is a hundred times too large from the first frame on. Nothing warns you; the only defence is to dump ix iy iz (or xu yu zu) whenever a later analysis may need unwrapped positions.</p>`;
    });

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => { alive = false; };
  },
};
