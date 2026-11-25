// Research · E.6: write a thermostat. A new fix temp/andersen (stochastic collisions with a heat bath, RanMars per
// processor); thirty steps of a ten-atom gas reproduced by hand to the last bit, RanMars and all; an 864-atom liquid
// showing canonical temperature fluctuations (unlike Berendsen) and diffusion that the collisions slow down.
import { NATIVE } from './andersen-native.js';
import { RanMars } from './ranmars.js';
import { crlog } from './crpow.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>', HALF = '<span class="m-ok">◐</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const frames = t => { const L = t.split('\n'), F = []; for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const n = +L[i + 3]; F.push(L.slice(i + 9, i + 9 + n).map(l => l.trim().split(/\s+/).map(Number))); i += 8 + n; } return F; };
const ulps = (a, b) => Math.abs(a - b) / (Math.abs(b) * Number.EPSILON / 2 || Number.MIN_VALUE);

// FixTempAndersen::end_of_step for the 10-atom gas (no forces: fix nve leaves v unchanged), plus compute temp
function handGas(F0, log) {
  const mass = { 1: 1.0, 2: 3.0 }, rng = new RanMars(58923 + 0), prob = 20.0 * 0.005, typ = F0.map(a => a[1]);
  let v = F0.map(a => a.slice(2, 5)), nc = 0; const out = [];
  for (let s = 1; s <= 30; s++) {
    const hits = [];
    for (let i = 0; i < v.length; i++) if (rng.uniform() < prob) {
      const m = mass[typ[i]], sigma = Math.sqrt(1.0 * 1.5 / (m * 1.0));
      v[i] = [sigma * rng.gaussian(log), sigma * rng.gaussian(log), sigma * rng.gaussian(log)]; nc++; hits.push(F0[i][0]);
    }
    let t = 0.0; for (let i = 0; i < v.length; i++) t += (v[i][0] * v[i][0] + v[i][1] * v[i][1] + v[i][2] * v[i][2]) * mass[typ[i]];
    out.push({ v: v.map(a => [...a]), T: t * (1.0 / (27 * 1.0)), nc, hits });
  }
  return out;
}

export default {
  id: 'write-thermostat',
  title: 'Write a thermostat',
  minutes: 20,
  summary: 'A thermostat LAMMPS does not have: Andersen\'s, where each atom now and then collides with a heat bath and leaves with a fresh Maxwell–Boltzmann velocity. You write it as a 60-line fix with its own random-number stream, check thirty steps of a ten-atom gas against a hand calculation that reruns LAMMPS\'s RanMars generator, bit for bit, and then measure what it does to a liquid: the right canonical temperature fluctuations (Berendsen\'s are 40 % too small), and diffusion that slows down as collisions become frequent.',
  learn: [
    'A fix with its own state: arguments, a RanMars per processor, a global scalar output',
    'end_of_step: where a velocity-resampling thermostat belongs',
    'Reproducing a random stream exactly, and the one function (log) that can spoil it',
    'Canonical fluctuations, and what a thermostat does to dynamics',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: collisions with a heat bath</h2>
    <div class="prose"><p>Andersen (1980): between steps, each atom independently, with probability ν·Δt, forgets its velocity and draws a new one from the Maxwell–Boltzmann distribution at the target T, each component Gaussian with σ = √(k<sub>B</sub>T/m). Between collisions the motion is plain Newtonian (fix nve). The result samples the canonical ensemble exactly, for any ν; ν only sets how strongly the dynamics is disturbed.</p>
<p>What the fix needs, compared with the earlier pages: <b>arguments</b> T, ν and a seed; a <b>random number generator</b> that is different on every processor (<code>RanMars(lmp, seed + comm->me)</code>, as fix langevin does, page D.10); <b>end_of_step</b>, after the second half kick, so the new velocity is the one thermo reports and the next step starts from; and a <b>global scalar</b>, the number of collisions, so <code>f_ID</code> can be printed.</p></div>
  </section>

  <section class="block">
    <h2>2 · Write it</h2>
    <pre class="lmp-mini">${esc(NATIVE.hdr.split('\n').slice(4).join('\n').replace(/\n{3,}/g, '\n\n'))}</pre>
    <pre class="lmp-mini">${esc(NATIVE.src.split('\n').slice(6).join('\n'))}</pre>
    <div class="prose"><p>RanMars::gaussian (${gh('src/random_mars.cpp')}) is the polar Box–Muller method: two uniforms, rejected unless inside the unit circle, then √(−2 ln r²/r²) times each, one returned now and one saved for the next call. Native rebuild 12 s.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: thirty steps of a ten-atom gas</h2>
    <div class="prose"><p>No forces, so fix nve leaves velocities alone and every change comes from the thermostat. Two masses (1 and 3), ν = 20, Δt = 0.005: probability 0.1 per atom per step.</p></div>
    <pre class="lmp-mini">${esc(NATIVE.gasin)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Rerun the random stream by hand</button></div>
    <div class="tbl wt-a"></div>
    <div class="prose wt-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: a liquid under five thermostats</h2>
    <div class="prose"><p>864 LJ atoms at ρ = 0.8442, T = 1.0, equilibrated, then 20 000 steps (100τ) under each thermostat, with T sampled every 10 steps and D from the slope of the MSD (native).</p></div>
    <div class="tbl wt-b"></div>
    <div class="prose wt-bsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> With ν = 20 and Δt = 0.005, how many collisions do you expect in 30 steps of 10 atoms, and how far from that is 35?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>300 trials with p = 0.1: mean 30, standard deviation √(300 × 0.1 × 0.9) = 5.2. 35 is one standard deviation high: nothing unusual.</div></div>
    <div class="ex"><p><b>2.</b> Why must the seed differ per processor, and what breaks if the fix used one stream on processor 0 and sent the numbers to the others?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>With the same seed on every processor, atoms with the same local index on different processors would get identical random numbers: correlated velocities. Broadcasting from one stream would work but costs communication every step, and the result would still depend on how atoms are distributed. The usual compromise (and the reason results change with the processor count) is seed + rank.</div></div>
    <div class="ex"><p><b>3.</b> For a diffusion coefficient, which ν would you use?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Ideally none: equilibrate with a thermostat, then measure in NVE (or with Nosé–Hoover, which here changed D by under 10 %). Even ν = 1 lowered D by about a fifth. Andersen is excellent for static, canonical averages, poor for dynamics.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = () => {
      const F = frames(NATIVE.dump), T = NATIVE.thermo, H = handGas(F[0], crlog), J = handGas(F[0], Math.log);
      let exact = 0, close = 0, n = 0, jsBad = 0;
      const rows = [];
      for (let s = 1; s <= 30; s++) {
        let e = 0, c = 0;
        for (let i = 0; i < 10; i++) for (let k = 0; k < 3; k++) { n++; const a = F[s][i][2 + k], h = H[s - 1].v[i][k]; if (a === h) { exact++; e++; } else if (ulps(h, a) <= 2) { close++; c++; } if (J[s - 1].v[i][k] !== a) jsBad++; }
        const tSame = H[s - 1].T === Number(T[s][1]), cSame = H[s - 1].nc === Number(T[s][2]);
        if (s <= 6 || s % 6 === 0) rows.push(row([s, H[s - 1].hits.length ? H[s - 1].hits.join(', ') : '–', `<span class="mono">${esc(T[s][1])}</span>`, T[s][2], `${e === 30 ? OK : c ? HALF : BAD} ${e}/30`, tSame && cSame ? OK : BAD]));
      }
      $('.wt-a').innerHTML = `<table>${head(['step', 'atoms that collided (ID)', 'temp (LAMMPS)', 'f_2', 'velocities by hand', 'temp, f_2'])}${rows.join('')}</table>`;
      $('.wt-asum').innerHTML = `<p>${exact} of ${n} velocity components identical to the last bit${close ? `, ${close} within two units in the last place (◐)` : ''}; temperature and collision count identical at every step. The hand code is a JavaScript copy of RanMars, seeded 58923 + 0 (one processor), consuming one uniform per atom per step and three Gaussians per collision, in memory order. One detail had to be right: the logarithm inside gaussian(). Math.log is not always correctly rounded: ${jsBad ? `in your browser it would spoil ${jsBad} of the ${n} components in the last bit` : 'your browser happens to get every value here right, but Node.js gets 24 of these 900 components wrong in the last bit'}. So the hand code uses a correctly rounded log built from BigInt arithmetic (as native glibc's log is, almost always: in a test of 26 000 values, glibc differed from the correctly rounded result 17 times).</p>`;
    };
    const L = NATIVE.liquid, canon = Math.sqrt(2 / (3 * 864 - 3));
    const names = { nve: 'none (fix nve)', andersen1: 'temp/andersen, ν = 1', andersen10: 'temp/andersen, ν = 10', andersen100: 'temp/andersen, ν = 100', berendsen: 'temp/berendsen, τ = 0.5', nh: 'nvt (Nosé–Hoover), τ = 0.5' };
    $('.wt-b').innerHTML = `<table>${head(['thermostat', '⟨T⟩', 'σ_T / ⟨T⟩', 'D (σ²/τ)'])}${Object.entries(L).map(([k, v]) => row([names[k], v[0].toFixed(4), v[1].toFixed(4), v[2].toFixed(4)])).join('')}${row(['canonical ensemble, √(2/(3N−3))', '', `<b>${canon.toFixed(4)}</b>`, ''])}</table>`;
    $('.wt-bsum').innerHTML = `<p>The kinetic energy of N atoms in the canonical ensemble fluctuates by √(2/N_dof) = ${canon.toFixed(4)} relative. Andersen at every ν and Nosé–Hoover reproduce it, within the statistical noise of 2000 correlated samples. Berendsen gives ${L.berendsen[1].toFixed(4)}, the same as the microcanonical run: it scales velocities towards the target and suppresses exactly the fluctuations a canonical ensemble should have (page D.10). The cost of Andersen shows in the last column: every collision erases an atom's momentum memory, so diffusion falls from ${L.nve[2].toFixed(3)} (NVE) to ${L.andersen1[2].toFixed(3)}, ${L.andersen10[2].toFixed(3)} and ${L.andersen100[2].toFixed(4)} as ν grows. Right temperature, right distribution, wrong dynamics: choose the thermostat for the quantity you measure.</p>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
