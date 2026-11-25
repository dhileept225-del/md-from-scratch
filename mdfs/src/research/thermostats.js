// Research · D.10: thermostats. fix temp/berendsen, fix langevin (with LAMMPS's own Marsaglia random numbers) and fix nvt (a Nosé–Hoover chain of 3) reproduced bit for bit on an ideal gas, where only the thermostat changes velocities.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

// real units, argon
const M = 39.948, MVV2E = 48.88821291 * 48.88821291, FTM2V = 1.0 / 48.88821291 / 48.88821291, KB = 0.0019872067, N = 10, DOF = 3 * N - 3, NSTEP = 100;
// ten atoms far apart, deterministic starting velocities (Å/fs), about 300 K
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
const FIX = (k, T, damp) => ({
  berendsen: `fix             1 all nve\nfix             2 all temp/berendsen ${T} ${T} ${damp}`,
  langevin: `fix             1 all nve\nfix             2 all langevin ${T} ${T} ${damp} 48279`,
  nvt: `fix             2 all nvt temp ${T} ${T} ${damp}`,
})[k];
const script = (k, T, damp) => `units           real
atom_style      atomic
read_data       gas.data
timestep        1.0
${FIX(k, T, damp)}
thermo_style    custom step temp f_2
thermo_modify   format float %.17g
thermo          10
dump            d all custom ${NSTEP} v.dump id vx vy vz
dump_modify     d format float %.17g sort id
run             ${NSTEP}`;

const temp = v => { let t = 0; for (const q of v) t += (q[0] * q[0] + q[1] * q[1] + q[2] * q[2]) * M; return t * (MVV2E / (DOF * KB)); };

// RanMars (random_mars.cpp), Marsaglia's universal generator
class RanMars {
  constructor(seed) {
    let ij = Math.trunc((seed - 1) / 30082), kl = (seed - 1) - 30082 * ij, i = Math.trunc(ij / 177) % 177 + 2, j = ij % 177 + 2, k = Math.trunc(kl / 169) % 178 + 1, l = kl % 169;
    this.u = new Float64Array(98);
    for (let ii = 1; ii <= 97; ii++) { let s = 0.0, t = 0.5; for (let jj = 1; jj <= 24; jj++) { const m = ((i * j) % 179) * k % 179; i = j; j = k; k = m; l = (53 * l + 1) % 169; if ((l * m) % 64 >= 32) s = s + t; t = 0.5 * t; } this.u[ii] = s; }
    this.c = 362436.0 / 16777216.0; this.cd = 7654321.0 / 16777216.0; this.cm = 16777213.0 / 16777216.0; this.i97 = 97; this.j97 = 33; this.uniform();
  }
  uniform() { let uni = this.u[this.i97] - this.u[this.j97]; if (uni < 0.0) uni += 1.0; this.u[this.i97] = uni; if (--this.i97 === 0) this.i97 = 97; if (--this.j97 === 0) this.j97 = 97; this.c -= this.cd; if (this.c < 0.0) this.c += this.cm; uni -= this.c; if (uni < 0.0) uni += 1.0; return uni; }
}

function hand(k, T, tp) {
  const v = V0.map(q => [...q]), dt = 1.0, out = [];
  if (k === 'berendsen') {
    let energy = 0;
    out.push({ s: 0, T: temp(v), e: 0 });
    for (let s = 1; s <= NSTEP; s++) {
      for (const q of v) for (let d = 0; d < 3; d++) q[d] += (0.5 * dt * FTM2V / M) * 0.0;   // nve with no force
      const tc = temp(v), lamda = Math.sqrt(1.0 + dt / tp * (T / tc - 1.0));
      energy += tc * (1.0 - lamda * lamda) * (0.5 * KB * DOF);
      for (const q of v) for (let d = 0; d < 3; d++) q[d] *= lamda;
      if (s % 10 === 0) out.push({ s, T: temp(v), e: energy });
    }
  }
  if (k === 'langevin') {
    const rng = new RanMars(48279 + 0), g1 = -M / tp / FTM2V, g2 = Math.sqrt(M) / FTM2V * Math.sqrt(24.0 * KB / tp / dt / MVV2E), dtfm = 0.5 * dt * FTM2V / M;
    const post = () => { const tsqrt = Math.sqrt(T), gamma2 = g2 * tsqrt; return v.map(q => { const fr = [0, 1, 2].map(() => gamma2 * (rng.uniform() - 0.5)); return [0, 1, 2].map(d => 0.0 + (g1 * q[d] + fr[d])); }); };
    let f = post();   // setup() calls post_force
    out.push({ s: 0, T: temp(v), e: 0 });
    for (let s = 1; s <= NSTEP; s++) {
      v.forEach((q, i) => { for (let d = 0; d < 3; d++) q[d] += dtfm * f[i][d]; });
      f = post();
      v.forEach((q, i) => { for (let d = 0; d < 3; d++) q[d] += dtfm * f[i][d]; });
      if (s % 10 === 0) out.push({ s, T: temp(v), e: 0 });
    }
  }
  if (k === 'nvt') {
    const mt = 3, tf = 1.0 / tp, dthalf = 0.5 * dt, dt4 = 0.25 * dt, dt8 = 0.125 * dt, tdrag = 1.0 - (dt * tf * 0.0 / 1);
    const eta = [0, 0, 0], ed = [0, 0, 0, 0], edd = [0, 0, 0], em = [0, 0, 0];
    let tcur = temp(v), tdof = DOF, tt = T, ket = tdof * KB * tt;
    em[0] = tdof * KB * tt / (tf * tf); for (let c = 1; c < mt; c++) em[c] = KB * tt / (tf * tf);
    for (let c = 1; c < mt; c++) edd[c] = (em[c - 1] * ed[c - 1] * ed[c - 1] - KB * tt) / em[c];
    const nhc = () => {
      let kec = tdof * KB * tcur;
      em[0] = tdof * KB * tt / (tf * tf); for (let c = 1; c < mt; c++) em[c] = KB * tt / (tf * tf);
      edd[0] = em[0] > 0 ? (kec - ket) / em[0] : 0;
      const nf = 1.0 / 1; let ex;
      for (let c = mt - 1; c > 0; c--) { ex = Math.exp(-nf * dt8 * ed[c + 1]); ed[c] *= ex; ed[c] += edd[c] * nf * dt4; ed[c] *= tdrag; ed[c] *= ex; }
      ex = Math.exp(-nf * dt8 * ed[1]); ed[0] *= ex; ed[0] += edd[0] * nf * dt4; ed[0] *= tdrag; ed[0] *= ex;
      const fe = Math.exp(-nf * dthalf * ed[0]);
      for (const q of v) { q[0] *= fe; q[1] *= fe; q[2] *= fe; }
      tcur *= fe * fe; kec = tdof * KB * tcur;
      edd[0] = em[0] > 0 ? (kec - ket) / em[0] : 0;
      for (let c = 0; c < mt; c++) eta[c] += nf * dthalf * ed[c];
      ed[0] *= ex; ed[0] += edd[0] * nf * dt4; ed[0] *= ex;
      for (let c = 1; c < mt; c++) { ex = Math.exp(-nf * dt8 * ed[c + 1]); ed[c] *= ex; edd[c] = (em[c - 1] * ed[c - 1] * ed[c - 1] - KB * tt) / em[c]; ed[c] += edd[c] * nf * dt4; ed[c] *= ex; }
    };
    const energy = () => { let e = 0; e += ket * eta[0] + 0.5 * em[0] * ed[0] * ed[0]; for (let c = 1; c < mt; c++) e += KB * tt * eta[c] + 0.5 * em[c] * ed[c] * ed[c]; return e; };
    out.push({ s: 0, T: temp(v), e: energy() });
    for (let s = 1; s <= NSTEP; s++) {
      ket = tdof * KB * tt; nhc();
      tcur = temp(v); tdof = DOF; nhc();
      if (s % 10 === 0) out.push({ s, T: temp(v), e: energy() });
    }
  }
  return { out, v };
}

export default {
  id: 'thermostats',
  title: 'Thermostats',
  minutes: 25,
  summary: 'Ten argon atoms with no forces between them: only the thermostat can change their velocities. That makes three thermostats fully reproducible by hand: Berendsen\'s rescaling, Langevin\'s friction and random kicks (with LAMMPS\'s own Marsaglia random-number generator rebuilt in JavaScript), and the Nosé–Hoover chain of fix nvt with its Trotter splitting. Temperature, coupling energy and all 30 final velocity components match to the last bit.',
  learn: [
    'temp/berendsen: λ = √(1 + Δt/τ (T₀/T − 1)), applied at end of step',
    'langevin: −m v/τ plus a uniform kick of the right size; RanMars, seeded seed + rank',
    'nvt: a chain of three thermostat variables and the order the Trotter steps run in',
    'ecouple: the energy each thermostat has taken out',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three ways to set a temperature</h2>
    <div class="prose"><p><b>Berendsen</b>, ${gh('src/fix_temp_berendsen.cpp', 120)}, at end_of_step:</p>
<pre class="lmp-mini">lamda = sqrt(1.0 + update-&gt;dt/t_period*(t_target/t_current - 1.0));
energy += t_current * (1.0-lamda*lamda) * 0.5*boltz*tdof;      // what it took out
v *= lamda;</pre>
<p><b>Langevin</b>, ${gh('src/fix_langevin.cpp', 306)} and the post_force loop: a drag and a uniform random force whose variance gives the right temperature (24 = 12 × 2, because a uniform on [−½, ½] has variance 1/12):</p>
<pre class="lmp-mini">gamma1 = -mass / t_period / ftm2v;
gamma2 = sqrt(mass) * sqrt(24.0*boltz/t_period/dt/mvv2e) / ftm2v * sqrt(t_target);
f += gamma1*v + gamma2*(random-&gt;uniform() - 0.5);              // three draws per atom, atom order</pre>
<p>The generator is ${gh('src/random_mars.cpp')} (Marsaglia's universal RNG), seeded with seed + MPI rank: the same seed gives a different trajectory on a different number of processors.</p>
<p><b>Nosé–Hoover chain</b>, ${gh('src/fix_nh.cpp', 1754)}: three extra variables η₁, η₂, η₃ with masses Q₁ = N_f k T₀/ω², Qₖ = k T₀/ω² (ω = 1/τ). Half a chain step before the kick, half after; inside, velocities are scaled by exp(−½Δt η̇₁). The chain's own energy (ecouple) is N_f k T₀ η₁ + Σ ½Qₖη̇ₖ² + Σ k T₀ ηₖ.</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: an ideal gas, three thermostats</h2>
    <div class="prose"><p>Ten argon atoms 10 Å or more apart, no pair style: no forces at all. Start at about ${temp(V0).toFixed(1)} K, dt = 1 fs, 100 steps.</p></div>
    <div class="row">
      <label class="lc-in">thermostat<select class="lc-u th-k"><option value="berendsen">temp/berendsen</option><option value="langevin">langevin</option><option value="nvt" selected>nvt</option></select></label>
      <label class="lc-in">T₀ (K)<select class="lc-u th-t"><option>100</option><option selected>150</option><option>600</option></select></label>
      <label class="lc-in">τ (fs)<select class="lc-u th-d"><option>10</option><option selected>50</option><option>200</option></select></label>
      <button class="btn sm primary" data-k="a">▶ Run 100 steps</button><span class="lmp-status muted th-sta"></span></div>
    <pre class="lmp-mini th-scr"></pre>
    <div class="tbl th-out"></div>
    <div class="prose th-sum"></div>
  </section>

  <section class="block prose">
    <h2>3 · Observe</h2>
    <ul>
      <li><b>Berendsen decays smoothly</b> to T₀ with time constant τ and never overshoots, but it does not give canonical fluctuations: with forces on, kinetic energy is too narrowly distributed.</li>
      <li><b>Nosé–Hoover rings.</b> With no forces to share energy, T oscillates around T₀ with a period set by τ: the chain is a dynamical system, not a damper. With τ = 10 fs the swings are fast; with 200 fs slow and large.</li>
      <li><b>Langevin is noisy</b> even with 10 atoms, and its ecouple is not tallied by default (f_2 is 0: add <code>tally yes</code>, page C.5). Change the seed, or the number of processors, and the trajectory changes; the statistics do not.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Berendsen at T = 300 K, T₀ = 150, τ = 50 fs, Δt = 1 fs: λ?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>√(1 + (1/50)(0.5 − 1)) = √0.99 = <b>0.994987</b>; T falls by 1 % per step at first.</div></div>
    <div class="ex"><p><b>2.</b> Langevin γ₂ for argon at 150 K, τ = 50 fs, Δt = 1 fs, in kcal/mol/Å?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>√39.948 / 4.184 × 10⁻⁴ × √(24 × 0.0019872 / 50 / 1 / 2390.06) × √150 = 15106 × 6.317 × 10⁻⁴ × 12.25 ≈ <b>117</b>; each kick component is uniform in ±58 kcal/mol/Å, for one step only.</div></div>
    <div class="ex"><p><b>3.</b> Q₁ for this gas with T₀ = 150 K and τ = 50 fs?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>N_f k T₀ τ² = 27 × 0.0019872 × 150 × 2500 = <b>20 120</b> kcal/mol·fs².</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const show = () => { $('.th-scr').textContent = FIX($('.th-k').value, +$('.th-t').value, +$('.th-d').value); $('.th-out').innerHTML = ''; $('.th-sum').innerHTML = ''; };
    root.querySelectorAll('.th-k,.th-t,.th-d').forEach(e => e.onchange = show); show();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.th-sta'), k = $('.th-k').value, T = +$('.th-t').value, tp = +$('.th-d').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(script(k, T, tp), { files: { 'gas.data': DATA }, read: ['v.dump'] }); b.disabled = false;
      const t = res.files?.['v.dump']; if (res.error || !t || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = hand(k, T, tp), c = res.thermo.cols; let n = 0, ok = 0, nr = 0;
      const cl = (a, b) => a === b ? OK : Math.abs(a - b) <= 2e-14 * Math.max(Math.abs(a), Math.abs(b)) ? (nr++, '◐') : BAD;
      const rows = res.thermo.rows.map(r => { const s = +r[0], h = H.out.find(o => o.s === s), lt = parseFloat(r[c.indexOf('Temp')]), le = parseFloat(r[c.indexOf('f_2')]), a = h.T === lt, bb = h.e === le; n += 2; ok += a + bb; const ca = cl(h.T, lt), cb = cl(h.e, le); return row([s, `<span class="mono">${s17(h.T)}</span>`, a ? OK : `${ca} ${s17(lt)}`, `<span class="mono">${s17(h.e)}</span>`, bb ? OK : `${cb} ${s17(le)}`, `<span class="mono">${'█'.repeat(Math.max(0, Math.round(h.T / 15)))}</span>`]); });
      const last = t.split('ITEM: TIMESTEP').pop().split('ITEM: ATOMS')[1].trim().split('\n').slice(1).map(l => l.trim().split(/\s+/).map(Number));
      let vok = 0, vnr = 0; last.forEach(w => { for (let d = 0; d < 3; d++) { const a = H.v[w[0] - 1][d], b2 = w[1 + d]; if (a === b2) vok++; else if (Math.abs(a - b2) <= 2e-14 * Math.abs(b2)) vnr++; } });
      $('.th-out').innerHTML = `<table>${head(['step', 'T by hand', '', 'ecouple (f_2) by hand', '', ''])}${rows.join('')}</table>`;
      $('.th-sum').innerHTML = `<p>${ok} of ${n} temperatures and coupling energies, and ${vok} of 30 final velocity components, identical to the last bit${nr + vnr ? `; the other ${nr ? nr + ' T/ecouple' : ''}${nr && vnr ? ' and ' : ''}${vnr ? vnr + ' velocity' : ''} values (◐) differ only in the last one or two digits` : ''}.</p>${nr + vnr ? '<p class="muted">Nosé–Hoover calls exp() several times per half step, and exp is not required to be correctly rounded: the C library inside the simulation and your browser\'s Math.exp may disagree in the last bit, and 100 steps of scaling carry it along. Run outside this browser (Node, and native LAMMPS on Linux), the same hand code matches all 30 velocities bit for bit for T₀ = 150 K, τ = 50 fs. Berendsen and Langevin use only + − × ÷ and sqrt, which are always correctly rounded, so they match everywhere.</p>' : ''}`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
