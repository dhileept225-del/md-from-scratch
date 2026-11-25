// Research · C.5: energy and virial bookkeeping. Where pe, ke, etotal, ecouple, econserve and press come from; tail and shift matched by hand; the thermostat's energy account.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
const rowsOf = res => res.thermo.rows.map(r => Object.fromEntries(res.thermo.cols.map((c, i) => [c, r[i]])));

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
`;
const scriptA = m => `${BASE}${m ? `pair_modify     ${m}\n` : ''}fix             1 all nve
compute         pk all pressure thermo_temp ke
compute         pv all pressure thermo_temp virial
thermo_style    custom step temp pe ke etotal press c_pk c_pv vol
thermo_modify   norm no format float %.17g
run             0`;
const THERMO = {
  nve: 'fix             1 all nve',
  rescale: 'fix             1 all nve\nfix             2 all temp/rescale 20 1.0 1.0 0.0 1.0',
  nvt: 'fix             1 all nvt temp 1.0 1.0 0.5',
  langevin: 'fix             1 all nve\nfix             2 all langevin 1.0 1.0 0.5 48279 tally yes',
};
const scriptB = k => `${BASE}${THERMO[k]}
thermo_style    custom step temp pe ke etotal ecouple econserve
thermo_modify   norm no format float %.17g
thermo          20
run             200`;

const SRC = [
  ['pe', 'Σ pair (eng_vdwl + eng_coul) + bond + angle + dihedral + improper + kspace + etail/V + fixes with fix_modify energy yes', gh('src/compute_pe.cpp', 91)],
  ['ke', '½ · dof · kB · T (from compute temp, not from ½mv²)', gh('src/thermo.cpp', 1958)],
  ['etotal', 'pe + ke', gh('src/thermo.cpp', 1964)],
  ['ecouple', 'Σ compute_scalar() of thermostat and barostat fixes: the energy they have taken out', gh('src/modify.cpp', 505)],
  ['econserve', 'etotal + ecouple: what should stay constant', gh('src/thermo.cpp', 2129)],
  ['press', '(dof · kB · T + Wxx + Wyy + Wzz) / 3V; W = Σ x·f of pair + bonded + kspace + fixes with fix_modify virial yes, + ptail/V', gh('src/compute_pressure.cpp', 256)],
];

export default {
  id: 'energy-bookkeeping',
  title: 'Energy and virial bookkeeping',
  minutes: 20,
  summary: 'Each thermo number is a sum someone chose. You see what goes into pe, ke, etotal, ecouple, econserve and press, then turn on the tail correction and the energy shift and predict both changes by hand: one from a formula, one by counting 108 000 neighbour pairs in the FCC lattice. Finally four thermostats: TotEng falls by thousands, but Econserve, which adds back what the thermostat took, stays flat.',
  learn: [
    'Which C++ sum is behind each thermo keyword',
    'pair_modify tail: etail/V and ptail/V², by hand to the last bit',
    'pair_modify shift: ΔPE = number of pairs × offset, with the pairs counted from FCC shells',
    'ecouple and econserve for temp/rescale, nvt and langevin',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: who adds what</h2>
    <div class="tbl"><table>${head(['keyword', 'sum', 'source'])}${SRC.map(r => row([`<code>${r[0]}</code>`, esc(r[1]), r[2]])).join('')}</table></div>
    <div class="prose"><p>In LJ units every extensive one of these is divided by N unless <code>thermo_modify norm no</code> (page C.3). This page sets <code>norm no</code> so the sums stay sums.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: tail and shift at step 0</h2>
    <div class="prose"><p>The same melt box, <code>run 0</code>, with three versions of the pair: plain, <code>pair_modify tail yes</code> (adds the missing energy and pressure of pairs beyond 2.5, assuming g(r) = 1 there) and <code>pair_modify shift yes</code> (moves each pair's E(r) up so it is 0 at the cutoff). Tail formula, ${gh('src/pair_lj_cut.cpp', 551)}:</p>
<pre class="lmp-mini">prefactor = 8π N_i N_j ε σ⁶ / (9 rc⁹)
etail = prefactor · (σ⁶ − 3 rc⁶)          // added to pe as etail / V
ptail = 2 · prefactor · (2σ⁶ − 3 rc⁶)     // added to each W_kk as ptail / V</pre></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run all three</button><span class="lmp-status muted ec-sta"></span></div>
    <div class="tbl ec-outa"></div>
    <div class="calc ec-hand" hidden></div>
  </section>

  <section class="block">
    <h2>3 · Modify: pick a thermostat</h2>
    <div class="prose"><p>Melt starts at T = 3. Each thermostat pulls it to 1 and keeps the books in <code>ecouple</code>. Watch TotEng and Econserve.</p></div>
    <div class="row"><label class="lc-in">integrator<select class="lc-u ec-t"><option value="nve">fix nve (no thermostat)</option><option value="rescale" selected>nve + temp/rescale every 20</option><option value="nvt">fix nvt, Tdamp 0.5</option><option value="langevin">nve + langevin, tally yes</option></select></label><button class="btn sm primary" data-k="b">▶ Run 200 steps</button><span class="lmp-status muted ec-stb"></span></div>
    <pre class="lmp-mini ec-scrb"></pre>
    <div class="tbl ec-outb"></div>
    <div class="prose ec-sumb"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe</h2>
    <ul>
      <li><b>Shift changes pe, never forces or pressure.</b> Press is identical with and without shift; the offset is a constant, so its derivative is 0.</li>
      <li><b>Tail changes pe and press, never forces.</b> Atoms move exactly the same with tail on. Only what you report changes, not the trajectory.</li>
      <li><b>Press is not exactly c_pk + c_pv.</b> The two parts add to …067 in the last digit, LAMMPS's own total is …062: (a + b)/3V and a/3V + b/3V round differently.</li>
      <li><b>Even fix nve's Econserve falls at the start.</b> From −9097.97 to −9142.89 in the first 20 steps: melting a perfect lattice from T = 3 with dt = 0.005 is a hard start. After that, Econserve changes by at most about 20 in every run, while with a thermostat TotEng falls by thousands.</li>
      <li><b>temp/rescale at step 20 gives the nve number back.</b> Up to step 20 the two runs are the same atoms; rescaling moves energy from etotal into ecouple, so Econserve at 20 matches the nve TotEng at 20 to the round-off.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Double the cutoff to 5.0 (and the box stays the same). By what factor does etail shrink, roughly?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>etail ≈ −(8π/3) N² ε σ⁶ / rc³ for rc ≫ σ, so about <b>8 times</b> smaller (2³). The exact ratio, with the σ⁶ term: 7.99.</div></div>
    <div class="ex"><p><b>2.</b> How many pairs would <code>shift yes</code> shift with cutoff 2.7 instead of 2.5, in the perfect lattice?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The 5th FCC shell, a·√(5/2) = 2.656, comes in: 12 + 6 + 24 + 12 + 24 = 78 neighbours per atom, 4000 × 78 / 2 = <b>156 000 pairs</b>.</div></div>
    <div class="ex"><p><b>3.</b> In a run with <code>fix nvt</code>, Econserve drifts up by 2 units over a long run. Is the thermostat wrong?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Not necessarily: Econserve carries the integrator's error, exactly like TotEng under nve. Check with a smaller dt; if the drift falls as dt², it is the integration, not the thermostat.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.ec-sta'); b.disabled = true; const R = {};
      for (const [k, m] of [['plain', ''], ['tail', 'tail yes'], ['shift', 'shift yes']]) {
        st.textContent = `Running the real LAMMPS (${k})…`;
        const res = await runLammps(scriptA(m));
        if (res.error || !res.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; b.disabled = false; return; }
        R[k] = rowsOf(res)[0];
      }
      b.disabled = false; st.textContent = '';
      const f = k => Object.fromEntries(Object.entries(R[k]).map(([a, v]) => [a, parseFloat(v)]));
      const p = f('plain'), t = f('tail'), s = f('shift');
      $('.ec-outa').innerHTML = `<table>${head(['', 'PotEng', 'Press', 'c_pk (ke part)', 'c_pv (virial part)'])}${[['plain', R.plain], ['tail yes', R.tail], ['shift yes', R.shift]].map(([n, r]) => row([n, ...['PotEng', 'Press', 'c_pk', 'c_pv'].map(c => `<span class="mono">${r[c]}</span>`)])).join('')}</table>`;
      // hand
      const N = 4000, V = p.Volume, rc3 = 2.5 * 2.5 * 2.5, rc6 = rc3 * rc3, rc9 = rc3 * rc6;
      const pre = 8.0 * Math.PI * N * N * 1.0 * 1.0 / (9.0 * rc9), et = pre * (1.0 - 3.0 * rc6), pt = 2.0 * pre * (2.0 - 3.0 * rc6);
      const peT = p.PotEng + et / V, dP = pt / V / V;
      const a = Math.cbrt(4 / 0.8442), shells = [[0.5, 12], [1, 6], [1.5, 24], [2, 12], [2.5, 24]].map(([k, n]) => [a * Math.sqrt(k), n]);
      const nn = shells.filter(([d]) => d < 2.5).reduce((m, [, n]) => m + n, 0), pairs = N * nn / 2;
      const off = 4.0 * (Math.pow(1 / 2.5, 12) - Math.pow(1 / 2.5, 6)), peS = p.PotEng - pairs * off;
      const pk = 11997 * 1.0 * p.Temp / 3.0 * (1.0 / V) * 1.0;
      const ok = (h, l, tol) => (tol ? Math.abs(h - l) <= tol * Math.abs(l) : h === l) ? '✓' : '✗';
      const L = [
        `TAIL   V = ${s17(V)}, N = 4000, rc = 2.5`,
        `  prefactor = 8π · 4000² / (9 · 2.5⁹) = ${s17(pre)}`,
        `  etail = prefactor · (1 − 3 · 2.5⁶) = ${s17(et)}`,
        `  pe + etail/V = ${s17(p.PotEng)} + ${s17(et / V)} = ${s17(peT)}   LAMMPS ${s17(t.PotEng)}  ${ok(peT, t.PotEng)}`,
        `  ptail = ${s17(pt)};  ΔPress = 3 · (ptail/V) / 3V = ptail/V² = ${s17(dP)}`,
        `  LAMMPS ΔPress = ${s17(t.Press - p.Press)}  ${ok(dP, t.Press - p.Press, 1e-12)} (to 12 digits; the virial sum rounds in between)`, '',
        `SHIFT  offset = 4(2.5⁻¹² − 2.5⁻⁶) = ${s17(off)}`,
        `  FCC shells, a = (4/0.8442)^(1/3) = ${a.toFixed(4)}:  ${shells.map(([d, n]) => `${d.toFixed(3)}×${n}`).join('  ')}`,
        `  inside 2.5: ${shells.filter(([d]) => d < 2.5).map(([, n]) => n).join(' + ')} = ${nn} neighbours → 4000 × ${nn} / 2 = ${pairs} pairs`,
        `  pe − pairs · offset = ${s17(peS)}   LAMMPS ${s17(s.PotEng)}  ${ok(peS, s.PotEng, 1e-13)}`,
        `  (LAMMPS adds −offset once per pair, 108 000 times; you multiply once: same to 13 digits)`, '',
        `KE PART  c_pk = dof · T / 3V = 11997 × ${s17(p.Temp)} / 3 / V = ${s17(pk)}   LAMMPS ${s17(p.c_pk)}  ${ok(pk, p.c_pk)}`,
        `  c_pk + c_pv = ${s17(p.c_pk + p.c_pv)}   Press = ${s17(p.Press)}`];
      $('.ec-hand').hidden = false; $('.ec-hand').innerHTML = L.map(esc).join('\n').replace(/✓/g, '<span class="m-ok">✓</span>');
    };

    const showB = () => { $('.ec-scrb').textContent = THERMO[$('.ec-t').value]; $('.ec-outb').innerHTML = ''; $('.ec-sumb').innerHTML = ''; };
    $('.ec-t').onchange = showB; showB();
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.ec-stb'), k = $('.ec-t').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptB(k)); b.disabled = false;
      if (res.error || !res.thermo.rows.length) { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; return; }
      st.textContent = '';
      const rs = rowsOf(res); let same = 0;
      const short = x => (+x).toFixed(4);
      $('.ec-outb').innerHTML = `<table>${head(['step', 'Temp', 'TotEng', 'Ecouple', 'Econserve', 'TotEng + Ecouple'])}${rs.map(r => { const h = parseFloat(r.Ecouple) + parseFloat(r.TotEng), ok = h === parseFloat(r.Econserve); if (ok) same++; return row([r.Step, (+r.Temp).toFixed(4), short(r.TotEng), short(r.Ecouple), short(r.Econserve), ok ? '<span class="m-ok">✓ bit for bit</span>' : `✗ ${s17(h)}`]); }).join('')}</table>`;
      const after = rs.filter(r => +r.Step >= 20), span = c => { const v = after.map(r => +r[c]); return Math.max(...v) - Math.min(...v); };
      $('.ec-sumb').innerHTML = `<p>From step 20 to 200: TotEng spans <b>${span('TotEng').toFixed(2)}</b>, Econserve spans <b>${span('Econserve').toFixed(2)}</b>. ${same} of ${rs.length} Econserve values are TotEng + Ecouple to the last bit.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
