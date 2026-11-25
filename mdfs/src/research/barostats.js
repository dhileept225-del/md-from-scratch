// Research · D.11: barostats. fix press/berendsen reproduced by hand on an ideal gas (pressure from the kinetic term alone), its closed-form end point, and a live comparison with the Martyna–Tobias–Klein barostat of fix npt on the LJ melt.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const M = 39.948, MVV2E = 48.88821291 * 48.88821291, KB = 0.0019872067, NKTV2P = 68568.415, N = 10, DOF = 3 * N - 3;
const V0 = Array.from({ length: N }, (_, i) => [0, 1, 2].map(k => +(0.004 * Math.sin(1.7 * i + 2.3 * k + 0.5) + 0.0011 * Math.cos(3.1 * i * k + i)).toFixed(6)));
const DATA = `ideal gas

${N} atoms
1 atom types

0 60 xlo xhi
0 60 ylo yhi
0 60 zlo zhi

Masses

1 ${M}

Atoms # atomic

${V0.map((_, i) => `${i + 1} 1 ${5 + 10 * (i % 5)} ${10 + 20 * Math.floor(i / 5)} 30`).join('\n')}

Velocities

${V0.map((v, i) => `${i + 1} ${v.join(' ')}`).join('\n')}
`;
const scriptA = (P, tp, K) => `units           real
atom_style      atomic
read_data       gas.data
timestep        1.0
comm_modify     cutoff 2.0
fix             1 all nve
fix             2 all press/berendsen iso ${P} ${P} ${tp} modulus ${K}
thermo_style    custom step temp press vol lx
thermo_modify   format float %.17g
thermo          100
run             1000`;
const temp = v => { let t = 0; for (const q of v) t += (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) * M; return t * (MVV2E / (DOF * KB)); };
function handA(P0, tp, K) {
  const T = temp(V0); let lo = 0, hi = 60; const out = [];
  const press = () => { const L = hi - lo, inv = 1.0 / (L * L * L); return (DOF * KB * T + 0 + 0 + 0) / 3.0 * inv * NKTV2P; };
  out.push({ s: 0, P: press(), V: (hi - lo) ** 3, L: hi - lo });
  for (let s = 1; s <= 1000; s++) {
    const pc = press(), dil = Math.pow(1.0 - 1.0 / tp * (P0 - pc) / K, 1.0 / 3.0), ctr = 0.5 * (lo + hi);
    lo = (lo - ctr) * dil + ctr; hi = (hi - ctr) * dil + ctr;
    if (s % 100 === 0) { const L = hi - lo; out.push({ s, P: press(), V: L * L * L, L }); }
  }
  return { out, T, Vinf: DOF * KB * T * NKTV2P / (3 * P0) };
}

const scriptB = k => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 6 0 6 0 6
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.5 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        0.3 bin
neigh_modify    every 1 delay 0 check yes
${k === 'ber' ? 'fix             1 all nvt temp 1.5 1.5 0.5\nfix             2 all press/berendsen iso 1.0 1.0 0.5 modulus 10' : 'fix             1 all npt temp 1.5 1.5 0.5 iso 1.0 1.0 0.5'}
thermo_style    custom step temp press vol
thermo          50
run             3000`;

export default {
  id: 'barostats',
  title: 'Barostats',
  minutes: 20,
  summary: 'A barostat changes the box. fix press/berendsen does it the simple way: once per step, scale every length by a cube root that depends on how far the pressure is from the target, divided by a bulk modulus you guess. On an ideal gas the pressure is pure kinetic energy, so you can follow it by hand for 1000 steps and predict where the volume ends up. Then the LJ melt with Berendsen and with fix npt: one relaxes, the other rings, and both reach the same density.',
  learn: [
    'Pressure = kinetic part + virial; for an ideal gas only the first',
    'press/berendsen: dilation = [1 − Δt/τ_P (P₀ − P)/K]^(1/3), box scaled about its centre',
    'The bulk modulus K is only a rate setting, not physics',
    'MTK (fix npt): the box as a particle with mass W, and why it oscillates',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: scaling the box</h2>
    <div class="prose"><p>${gh('src/fix_press_berendsen.cpp', 314)}, at end_of_step:</p>
<pre class="lmp-mini">temperature-&gt;compute_scalar();  pressure-&gt;compute_scalar();     // P = (N_f k T + W)/3V
dilation = pow(1.0 - dt/p_period * (p_target - p_current)/bulkmodulus, 1.0/3.0);
remap();     // 385: x → λ, boxlo/hi = (old − ctr)·dilation + ctr, λ → x</pre>
<p>Velocities are not touched. The bulk modulus (default 10, in pressure units) only sets how fast the box moves: if it is much smaller than the real one, the box overshoots; much larger, and nothing happens. <code>fix npt</code> instead (${gh('src/fix_nh.cpp')}, Martyna–Tobias–Klein) gives the box a velocity ω̇ and a mass W = (N+1) k T₀ τ_P², and scales velocities too: a proper dynamics with the right fluctuations, and oscillations of period about 2π τ_P.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: an ideal gas, 1000 steps by hand</h2>
    <div class="prose"><p>The ten argon atoms of page D.10, no forces: T never changes, the virial is 0, and P = N_f k T / (3V) × nktv2p. Predict P, V and the box length every 100 steps; the end point is V∞ = N_f k T nktv2p / (3P₀).</p></div>
    <div class="row">
      <label class="lc-in">P₀ (atm)<select class="lc-u bs-p"><option>1</option><option selected>5</option><option>20</option></select></label>
      <label class="lc-in">τ_P (fs)<select class="lc-u bs-t"><option>50</option><option selected>100</option><option>1000</option></select></label>
      <label class="lc-in">modulus K (atm)<select class="lc-u bs-k"><option>1</option><option selected>10</option><option>100</option></select></label>
      <button class="btn sm primary" data-k="a">▶ Run 1000 steps</button><span class="lmp-status muted bs-sta"></span></div>
    <div class="calc bs-pred"></div>
    <div class="tbl bs-outa"></div>
    <div class="prose bs-suma"></div>
  </section>

  <section class="block">
    <h2>3 · Observe: Berendsen versus fix npt on a liquid</h2>
    <div class="prose"><p>864 LJ atoms from the fcc lattice at density 0.844, T = 1.5, target P = 1, τ_P = 0.5 (100 steps). The lattice starts at P ≈ −5, so the box must grow by almost half. Left: nvt + press/berendsen (modulus 10). Right: npt (MTK). 3000 steps each.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run both</button><span class="lmp-status muted bs-stb"></span></div>
    <div class="tbl bs-outb"></div>
    <div class="prose bs-sumb"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> P = 2.63 atm, P₀ = 5, τ_P = 100 fs, K = 10 atm, Δt = 1 fs: the first dilation?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>[1 − 0.01 × (5 − 2.63)/10]^(1/3) = 0.99763^(1/3) = <b>0.999210</b>: the box shrinks by 0.079 % in length.</div></div>
    <div class="ex"><p><b>2.</b> Why does the ideal-gas volume approach V∞ only slowly once P is close to P₀?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Linearize: ΔV/V ≈ −(Δt/τ_P)(P₀ − P)/K, and P ∝ 1/V, so the deviation decays with rate (Δt/τ_P)(P/K) per step: at P₀ = 5, K = 10 that is half of 1/τ_P. Raising K slows it, lowering K speeds it until it overshoots.</div></div>
    <div class="ex"><p><b>3.</b> When is Berendsen pressure coupling acceptable?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>For equilibration: it relaxes the density without ringing. It suppresses volume fluctuations, so do not compute compressibility or other fluctuation properties from it; switch to npt (MTK) or a stochastic barostat for production.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const pA = () => ({ P: +$('.bs-p').value, tp: +$('.bs-t').value, K: +$('.bs-k').value });
    const predA = () => { const p = pA(), H = handA(p.P, p.tp, p.K), P0 = H.out[0].P; $('.bs-pred').textContent = `T = ${H.T.toFixed(4)} K (constant), N_f = ${DOF}\nstart: V = 216000 Å³, P = ${P0.toFixed(5)} atm\nend point: V∞ = ${DOF} × ${KB} × ${H.T.toFixed(3)} × ${NKTV2P} / (3 × ${p.P}) = ${H.Vinf.toFixed(1)} Å³ (L = ${Math.cbrt(H.Vinf).toFixed(3)} Å)`; $('.bs-outa').innerHTML = ''; $('.bs-suma').innerHTML = ''; return H; };
    root.querySelectorAll('.bs-p,.bs-t,.bs-k').forEach(e => e.onchange = predA); predA();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.bs-sta'), p = pA(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptA(p.P, p.tp, p.K), { files: { 'gas.data': DATA } }); b.disabled = false;
      if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = predA(), c = res.thermo.cols; let ex = 0, nr = 0, n = 0;
      const cmp = (h, l) => { n++; if (h === l) { ex++; return OK; } if (Math.abs(h - l) <= 1e-13 * Math.abs(l)) { nr++; return '◐'; } return BAD; };
      $('.bs-outa').innerHTML = `<table>${head(['step', 'P by hand (atm)', '', 'V by hand (Å³)', '', 'Lx', ''])}${res.thermo.rows.map(r => { const h = H.out.find(o => o.s === +r[0]); const g = k => parseFloat(r[c.indexOf(k)]); return row([r[0], h.P.toPrecision(12), cmp(h.P, g('Press')), h.V.toPrecision(12), cmp(h.V, g('Volume')), h.L.toPrecision(12), cmp(h.L, g('Lx'))]); }).join('')}</table>`;
      const last = H.out[H.out.length - 1];
      $('.bs-suma').innerHTML = `<p>${ex} of ${n} identical to the last bit${nr ? `, ${nr} (◐) to 13 digits (the cube root is computed with pow, which the browser's C library may round differently from JavaScript)` : ''}. After 1000 steps V = ${last.V.toFixed(1)} Å³, ${(100 * (last.V / H.Vinf - 1)).toFixed(2)} % from V∞${p.K <= 1 ? ': with K this small the box moves fast (fine for an ideal gas; with real forces it can overshoot)' : p.K >= 100 ? ': with K this large the box creeps, and 1000 steps are not enough' : ''}.</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.bs-stb'); b.disabled = true; const R = {};
      for (const k of ['ber', 'npt']) { st.textContent = `Running ${k === 'ber' ? 'Berendsen' : 'npt'}…`; const res = await runLammps(scriptB(k)); if (res.error || !res.thermo.rows.length) { b.disabled = false; return err(st, res); } const c = res.thermo.cols; R[k] = res.thermo.rows.map(r => ({ s: +r[0], V: +r[c.indexOf('Volume')], P: +r[c.indexOf('Press')] })); }
      b.disabled = false; st.textContent = '';
      const all = [...R.ber, ...R.npt].map(o => o.V), lo = Math.min(...all), hi = Math.max(...all), bar = V => '█'.repeat(Math.max(1, Math.round(1 + 24 * (V - lo) / (hi - lo))));
      $('.bs-outb').innerHTML = `<table>${head(['step', 'V Berendsen', '', 'V npt', ''])}${R.ber.filter((_, i) => i % 3 === 0).map((o, i) => { const q = R.npt[3 * i]; return row([o.s, o.V.toFixed(1), `<span class="mono">${bar(o.V)}</span>`, q.V.toFixed(1), `<span class="mono">${bar(q.V)}</span>`]); }).join('')}</table>`;
      const avg = (a, k) => { const t = a.filter(o => o.s >= 1500); return t.reduce((s, o) => s + o[k], 0) / t.length; }, sd = a => { const t = a.filter(o => o.s >= 1500), m = avg(a, 'V'); return Math.sqrt(t.reduce((s, o) => s + (o.V - m) ** 2, 0) / t.length); };
      $('.bs-sumb').innerHTML = `<p>Steps 1500–3000: Berendsen ⟨V⟩ = ${avg(R.ber, 'V').toFixed(1)} (σ ${sd(R.ber).toFixed(1)}), npt ⟨V⟩ = ${avg(R.npt, 'V').toFixed(1)} (σ ${sd(R.npt).toFixed(1)}); ⟨P⟩ ${avg(R.ber, 'P').toFixed(3)} and ${avg(R.npt, 'P').toFixed(3)}. Both find the same density; Berendsen gets there without the overshoot, and with smaller volume fluctuations than the ensemble really has.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
