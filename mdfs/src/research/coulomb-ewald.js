// Research · D.7: Coulomb and Ewald. Why a cut-off Coulomb sum fails, Wolf and DSF live, and kspace_style ewald read line by line: g_ewald, kmax, k-vectors, real and reciprocal energies by hand against native runs, all against the Madelung energy of rock salt.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const grab = (lines, re) => { for (const l of lines) { const m = l.match(re); if (m) return m; } return null; };
const f8 = x => String(+x.toPrecision(8));

const QQ = 332.06371, A = 5.64, H = A / 2, MAD = 1.74756459463, MY_PIS = Math.sqrt(Math.PI);
const EP = 0.3275911, A1 = 0.254829592, A2 = -0.284496736, A3 = 1.421413741, A4 = -1.453152027, A5 = 1.061405429;   // ewald_const.h
const erfcPoly = x => { const ex = Math.exp(-x * x), t = 1.0 / (1.0 + EP * x); return t * (A1 + t * (A2 + t * (A3 + t * (A4 + t * A5)))) * ex; };
// accurate erfc (what the C library gives coul/wolf): series below 2.5, continued fraction above
function erfc(x) {
  if (x < 2.5) { let s = 0, term = x, n = 0; const x2 = x * x; while (Math.abs(term) > 1e-17 * Math.abs(s) || n < 3) { s += term / (2 * n + 1); n++; term *= -x2 / n; if (n > 200) break; } return 1 - 2 / MY_PIS * s; }
  let f = 0; for (let k = 60; k >= 1; k--) f = (k / 2) / (x + f); return Math.exp(-x * x) / MY_PIS / (x + f);
}
const madelung = n => -MAD * QQ / H * (8 * n ** 3) / 2;
// all neighbours of one ion in rock salt, within rc, as half-lattice offsets: sign +1 same species, −1 other
function shell(rc) { const R = Math.floor(rc / H) + 1, o = []; for (let i = -R; i <= R; i++) for (let j = -R; j <= R; j++) for (let k = -R; k <= R; k++) { if (!i && !j && !k) continue; const r = H * Math.sqrt(i * i + j * j + k * k); if (r * r < rc * rc) o.push([r, (i + j + k) % 2 ? -1 : 1]); } return o; }

// ---------- live: plain cut, Wolf, DSF ----------
const scriptL = (st, rc, al) => `units           real
atom_style      charge
lattice         fcc 5.64
region          box block 0 3 0 3 0 3
create_box      2 box
create_atoms    1 box
lattice         fcc 5.64 origin 0.5 0.0 0.0
create_atoms    2 box
mass            1 22.99
mass            2 35.45
set             type 1 charge 1.0
set             type 2 charge -1.0
pair_style      ${st === 'cut' ? `coul/cut ${rc}` : `coul/${st} ${al} ${rc}`}
pair_coeff      * *
thermo_style    custom step pe
thermo_modify   format float %.17g
run             0`;
function handL(st, rc, al) {
  const N = 216, sh = shell(rc); let per = 0, self = 0;
  if (st === 'cut') for (const [r, s] of sh) per += QQ * s / r;
  if (st === 'wolf') { const es = erfc(al * rc) / rc; self = -(es / 2.0 + al / MY_PIS) * QQ; for (const [r, s] of sh) per += (erfc(al * r) - es * r) * (QQ * s / r); }
  if (st === 'dsf') { const ec = erfc(al * rc), ed = Math.exp(-al * al * rc * rc), fs = -(ec / (rc * rc) + 2.0 / MY_PIS * al * ed / rc), es = ec / rc - fs * rc; self = -(es / 2.0 + al / MY_PIS) * QQ; for (const [r, s] of sh) per += (QQ * s / r) * (erfcPoly(al * r) - r * es - r * r * fs); }
  return N * self + 0.5 * N * per;
}

// ---------- native Ewald: by hand ----------
function ewald(n, acc, rc = 10.0) {
  const L = n * A, N = 8 * n ** 3, q2 = N * 1.0 * QQ, accuracy = acc * QQ, V = L * L * L;
  let g = accuracy * Math.sqrt(N * rc * L * L * L) / (2.0 * q2);
  g = g >= 1.0 ? (1.35 - 0.15 * Math.log(accuracy)) / rc : Math.sqrt(-Math.log(g)) / rc;
  const rms = km => 2.0 * q2 * g / L * Math.sqrt(1.0 / (Math.PI * km * N)) * Math.exp(-Math.PI * Math.PI * km * km / (g * g * L * L));
  let km = 1; while (rms(km) > accuracy) km++;
  const uk = 2.0 * Math.PI / L, gsq = uk * uk * km * km * 1.00001, ks = [];
  for (let kx = 0; kx <= km; kx++) for (let ky = -km; ky <= km; ky++) for (let kz = -km; kz <= km; kz++) {
    if (!kx && !ky && !kz) continue; if (kx === 0 && (ky < 0 || (ky === 0 && kz < 0))) continue;
    const s = (uk * kx) ** 2 + (uk * ky) ** 2 + (uk * kz) ** 2; if (s <= gsq) ks.push([kx, ky, kz, s]);
  }
  const lpr = rms(km), spr = 2.0 * (q2 / Math.sqrt(N * rc * V)) * Math.exp(-g * g * rc * rc), est = Math.sqrt(lpr * lpr + spr * spr);
  let per = 0; for (const [r, s] of shell(rc)) per += QQ * s / r * erfcPoly(g * r);
  const ecoul = 0.5 * N * per;
  // structure factor of rock salt: ions on all half-lattice points, Na where i+j+k even
  const m = 2 * n; let E = 0;
  for (const [kx, ky, kz, s] of ks) {
    let re = 0, im = 0; const ph = uk * H;
    for (let i = 0; i < m; i++) for (let j = 0; j < m; j++) for (let k = 0; k < m; k++) { const q = (i + j + k) % 2 ? -1 : 1, a = ph * (kx * i + ky * j + kz * k); re += q * Math.cos(a); im += q * Math.sin(a); }
    E += (4.0 * Math.PI / V) * Math.exp(-0.25 * s / (g * g)) / s * (re * re + im * im);
  }
  const selfE = g * N / MY_PIS * QQ, elong = (E - g * N / MY_PIS) * QQ;
  return { g, km, kcount: ks.length, kmax3d: 4 * km ** 3 + 6 * km ** 2 + 3 * km, est, ecoul, elong, selfE, pe: ecoul + elong, N };
}
// measured natively: LAMMPS 22 Jul 2025, pair coul/long 10.0, pair_modify table 0, kspace_style ewald acc
const NATIVE = [
  [2, 1e-4, 0.26827392, 0.033302057, [61, 3, 171], -6584.9537162344586, -3368.2955532190085, -3216.6581630154506],
  [2, 1e-5, 0.30821542, 0.003321771, [128, 4, 364], -6584.961598412985, -2890.1073901617674, -3694.8542082512172],
  [2, 1e-6, 0.34354417, 0.00033207714, [257, 5, 665], -6584.9973981586409, -2470.6403520627669, -4114.3570460958736],
  [3, 1e-4, 0.26827392, 0.039019911, [128, 4, 364], -22224.218792289656, -11367.99749211251, -10856.221300177147],
  [3, 1e-5, 0.30821542, 0.003321771, [462, 6, 1098], -22224.245394642428, -9754.1124417945721, -12470.132952847858],
  [3, 1e-6, 0.34354417, 0.00033305655, [709, 7, 1687], -22224.366218783565, -8338.4111882099933, -13885.955030573574],
  [4, 1e-4, 0.26827392, 0.033302057, [462, 6, 1098], -52679.629729872133, -26946.364425748532, -25733.265304123604],
  [4, 1e-5, 0.30821542, 0.0038674915, [709, 7, 1687], -52679.692787294232, -23120.859121284495, -29558.833666009737],
  [4, 1e-6, 0.34354417, 0.00033965882, [1535, 9, 3429], -52679.979185259814, -19765.122816492829, -32914.856368766988],
];

export default {
  id: 'coulomb-ewald',
  title: 'Coulomb and Ewald',
  minutes: 30,
  summary: 'Rock salt has an exact answer: the Madelung energy. A plain Coulomb cutoff misses it by 80 %, and you see why by summing shells by hand. Wolf and DSF damp and shift the sum and land within 0.2 % at rc = 10 Å, checked live. Then kspace_style ewald, read line by line: you compute g_ewald, kmax, the number of k-vectors, the accuracy estimate, the real-space and the reciprocal energies by hand, and match nine native runs to 12–15 digits; the total converges to Madelung as the accuracy tightens.',
  learn: [
    'Why 1/r cannot be cut: charged shells and conditional convergence',
    'Wolf and DSF: damped, shifted, and a self term',
    'Ewald in LAMMPS: g_ewald, kmax, k-vectors, self and neutralizing terms',
    'The polynomial erfc of coul/long, and what it costs',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the sum that will not converge</h2>
    <div class="prose"><p>Rock salt, a = 5.64 Å, ions ±1, nearest distance r₀ = 2.82 Å. The exact energy per ion pair is −M q²/(4πε₀ r₀) with the Madelung constant M = 1.747565; in real units −1.747565 × 332.06371 / 2.82 = <b>−205.781 kcal/mol</b>. Summing 1/r up to a cutoff adds whole shells of alternately charged ions, and the running sum swings by more than the answer. The order of summation changes the result: the sum is conditionally convergent.</p></div>
    <div class="tbl ce-shells"></div>
  </section>

  <section class="block">
    <h2>2 · Match: cut, Wolf and DSF, live</h2>
    <div class="prose"><p>216 ions (3 × 3 × 3 cells). Wolf (${gh('src/pair_coul_wolf.cpp', 55)}) damps each pair with erfc(αr), shifts it to zero at rc, and adds a self term per ion; DSF (${gh('src/pair_coul_dsf.cpp', 56)}) also shifts the force, and uses the polynomial erfc of section 3:</p>
<pre class="lmp-mini">Wolf:  E = Σ_pairs qᵢqⱼ (erfc(αr)/r − erfc(αrc)/rc) − Σᵢ qᵢ² (erfc(αrc)/(2rc) + α/√π)
DSF:   E = Σ_pairs qᵢqⱼ (erfc(αr)/r − e_shift − r f_shift) − Σᵢ qᵢ² (e_shift/2 + α/√π)</pre></div>
    <div class="row">
      <label class="lc-in">pair style<select class="lc-u ce-st"><option value="cut">coul/cut</option><option value="wolf" selected>coul/wolf</option><option value="dsf">coul/dsf</option></select></label>
      <label class="lc-in">rc (Å)<select class="lc-u ce-rc"><option>6</option><option>8</option><option selected>10</option><option>12</option></select></label>
      <label class="lc-in">α (1/Å)<select class="lc-u ce-al"><option>0.15</option><option selected>0.2</option><option>0.25</option><option>0.3</option></select></label>
      <button class="btn sm primary" data-k="l">▶ run 0</button><span class="lmp-status muted ce-stl"></span></div>
    <div class="calc ce-predl"></div>
    <div class="prose ce-outl"></div>
  </section>

  <section class="block">
    <h2>3 · Understand: Ewald in LAMMPS</h2>
    <div class="prose"><p>Split 1/r = erfc(gr)/r + erf(gr)/r. The first part is short-ranged: <code>pair_style coul/long</code> sums it to rc (${gh('src/KSPACE/pair_coul_long.cpp', 115)}), with erfc from a 5-term polynomial (Abramowitz–Stegun 7.1.26, error about 10⁻⁷; ${gh('src/ewald_const.h', 19)}) when <code>pair_modify table 0</code>. The second part is smooth: <code>kspace_style ewald</code> sums it over k-vectors. ${gh('src/KSPACE/ewald.cpp', 96)}:</p>
<pre class="lmp-mini">accuracy = acc_relative * qqr2e                     // force units, one e at 1 Å
q2 = Σqᵢ² · qqrd2e
g_ewald = accuracy*sqrt(N*rc*V) / (2*q2);           // 169
g_ewald = sqrt(-log(g_ewald)) / rc;                 // 171 (if &lt; 1)
rms(k) = 2 q2 g/L · sqrt(1/(π k N)) · exp(−π² k² / (g² L²));   // 341
kxmax = smallest k with rms(k) ≤ accuracy;          // 240
keep k-vectors in half space with |k|² ≤ (2π kmax/L)² · 1.00001
E_long = qqrd2e [ Σ_k (4π/V) exp(−k²/4g²)/k² |S(k)|² − g/√π Σqᵢ² − π/2 (Σqᵢ)²/(g² V) ]   // 457</pre>
<p>The browser build has no KSPACE package, so the runs below are native (mpirun, 1 process), embedded on the page. Everything else is computed now, in your browser.</p></div>
  </section>

  <section class="block">
    <h2>4 · Match: nine native Ewald runs</h2>
    <div class="row"><button class="btn sm primary" data-k="e">▶ Compute all nine by hand</button><span class="lmp-status muted ce-ste"></span></div>
    <div class="tbl ce-oute"></div>
    <div class="tbl ce-oute2"></div>
    <div class="prose ce-sume"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> g_ewald for 1e-5, rc = 10, 216 ions, L = 16.92 Å?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>accuracy = 1e-5 × 332.06 = 3.3206e-3; q2 = 216 × 332.06 = 71726; first estimate 3.3206e-3 × √(216 × 10 × 16.92³) / (2 × 71726) = 7.49 × 10⁻⁵; g = √(−ln 0.07479)/10 = <b>0.30822 Å⁻¹</b>, the same for every box here because N/V is fixed.</div></div>
    <div class="ex"><p><b>2.</b> Why does E_coul fall and E_long rise in magnitude as the accuracy tightens, while their sum stays put?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A tighter accuracy gives a larger g: erfc(gr) decays faster, so less is left in real space and more moves to k-space (more k-vectors). The split is a choice; the total is the physics.</div></div>
    <div class="ex"><p><b>3.</b> The 1e-6 runs land closer to Madelung than the 1e-4 runs, but not closer than about 10⁻⁶. What limits them?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The polynomial erfc in coul/long has an absolute error of about 1.5 × 10⁻⁷, and the k-space cutoff error estimate is an RMS force, not an energy bound. Tables (<code>pair_modify table 12</code>, the default) or a smaller error setting change the last digits again.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    // shell table
    { const sh = shell(13.0).sort((a, b) => a[0] - b[0]), byR = []; for (const [r, s] of sh) { const k = r.toFixed(4); const last = byR[byR.length - 1]; if (last && last.k === k) last.n++; else byR.push({ k, r, s, n: 1 }); }
      let run = 0; $('.ce-shells').innerHTML = `<table>${head(['shell r (Å)', 'ions', 'charge', 'adds (kcal/mol)', 'running sum', 'running sum ÷ exact (−205.781)'])}${byR.slice(0, 10).map(b => { const add = b.n * b.s * QQ / b.r; run += add; return row([b.r.toFixed(3), b.n, b.s < 0 ? 'opposite' : 'same', add.toFixed(1), run.toFixed(1), `${(run / -205.7811).toFixed(2)} ×`]); }).join('')}</table><p class="muted">Σⱼ qᵢqⱼ/rᵢⱼ around one ion = the energy per ion pair (each pair is shared by two ions). Opposite charges lower it, equal charges raise it.</p>`; }

    const pL = () => ({ st: $('.ce-st').value, rc: +$('.ce-rc').value, al: +$('.ce-al').value });
    const predL = () => { const p = pL(), E = handL(p.st, p.rc, p.al), M = madelung(3); $('.ce-al').disabled = p.st === 'cut'; $('.ce-predl').textContent = `hand: E = ${E.toPrecision(15)} kcal/mol   (Madelung: ${M.toPrecision(10)}, off by ${(100 * (E / M - 1)).toFixed(3)} %)`; $('.ce-outl').innerHTML = ''; return E; };
    root.querySelectorAll('.ce-st,.ce-rc,.ce-al').forEach(e => e.onchange = predL); predL();
    $('[data-k=l]').onclick = async e => {
      const b = e.currentTarget, st = $('.ce-stl'), p = pL(); b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptL(p.st, p.rc, p.al)); b.disabled = false; if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const E = predL(), l = parseFloat(res.thermo.rows[0][1]), rel = Math.abs(E - l) / Math.abs(l);
      $('.ce-outl').innerHTML = `<p>LAMMPS: PotEng = <b>${res.thermo.rows[0][1]}</b> ${rel < 1e-12 ? `${OK} (to ${rel === 0 ? 'the last bit' : rel.toExponential(0)})` : BAD}.</p>`;
    };

    $('[data-k=e]').onclick = () => {
      const st = $('.ce-ste'); st.textContent = 'Summing k-vectors…';
      setTimeout(() => {
        let ok = 0, tot = 0; const r1 = [], r2 = [];
        for (const [n, acc, gL, estL, [kc, k1, k3], pe, ec, el] of NATIVE) {
          const h = ewald(n, acc), c = [f8(h.g) === f8(gL), f8(h.est) === f8(estL), h.kcount === kc && h.km === k1 && h.kmax3d === k3];
          const e1 = Math.abs(h.ecoul - ec) / Math.abs(ec), e2 = Math.abs(h.elong - el) / Math.abs(el), M = madelung(n);
          c.push(e1 < 1e-12, e2 < 1e-12); c.forEach(x => { tot++; if (x) ok++; });
          r1.push(row([`${n}³ cells, ${8 * n ** 3} ions`, acc, `${f8(h.g)} ${c[0] ? OK : BAD}`, `${f8(h.est)} ${c[1] ? OK : BAD}`, `${h.kcount} ${h.km} ${h.kmax3d} ${c[2] ? OK : BAD}`]));
          r2.push(row([`${n}³`, acc, `${h.ecoul.toFixed(6)} ${c[3] ? OK : BAD}`, `${h.elong.toFixed(6)} ${c[4] ? OK : BAD}`, pe.toFixed(4), M.toFixed(4), `${((pe - M) / M).toExponential(1)}`]));
        }
        $('.ce-oute').innerHTML = `<table>${head(['system', 'accuracy', 'g_ewald (log: G vector)', 'estimated RMS force accuracy', 'k-vectors, kmax, max3d'])}${r1.join('')}</table>`;
        $('.ce-oute2').innerHTML = `<table>${head(['system', 'accuracy', 'E_coul (real space)', 'E_long (k-space + self)', 'LAMMPS total', 'Madelung', 'relative error'])}${r2.join('')}</table>`;
        $('.ce-sume').innerHTML = `<p>${ok} of ${tot} checks match the native log (g and accuracy to the 8 printed digits, counts exactly, energies to 10⁻¹² or better). The self term alone is −g N/√π × 332.06: for 216 ions at 1e-5, ${(-ewald(3, 1e-5).selfE).toFixed(3)} kcal/mol, inside E_long.</p>`;
        st.textContent = '';
      }, 30);
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
