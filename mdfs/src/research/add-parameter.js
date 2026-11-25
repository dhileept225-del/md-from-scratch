// Research · E.4: add a parameter. A fresh clone of lj/cut (lj/mine/soft) gets a per-pair lambda and a soft core;
// every place a per-pair parameter must be added; 18 two-atom values by hand; melts and a fix adapt lambda ramp
// bit for bit against LAMMPS's lj/cut/soft; restart and data files carry lambda.
import { NATIVE } from './softcore-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

// the hunks of the diff that add lambda, without the deleted rRESPA / born_matrix bodies
function keyHunks(diff) {
  return diff.split(/\n(?=@@)/).slice(1)
    .map(h => { const k = h.search(/\n-void PairLJMineSoft::(compute_inner|born_matrix)/); return k > 0 ? h.slice(0, k) : h; })
    .filter(h => /^\+/m.test(h) && !/^@@[^\n]*\n[^\n]*compute_(inner|middle|outer)/.test(h))
    .map(h => h.split('\n').filter(l => !/^ \s*$/.test(l)).join('\n')).join('\n');
}
// PairLJMineSoft::compute for one pair along x (eps = sigma = 1, n = 2, alpha = 0.5)
function softPair(lam, r) {
  const lj1 = lam * lam, lj2 = 1.0, lj3 = 0.5 * (1.0 - lam) * (1.0 - lam), eps = 1.0;   // pow(lam, 2.0) is correctly rounded, so = lam·lam
  const x2 = 5.0 + Number(r), delx = 5.0 - x2, rsq = delx * delx + 0.0 * 0.0 + 0.0 * 0.0;
  const r4sig6 = rsq * rsq / lj2, denlj = lj3 + rsq * r4sig6;
  const fpair = 1.0 * (lj1 * eps * (48.0 * r4sig6 / (denlj * denlj * denlj) - 24.0 * r4sig6 / (denlj * denlj)));
  let ev = lj1 * 4.0 * eps * (1.0 / (denlj * denlj) - 1.0 / denlj) - 0.0; ev *= 1.0;
  return { pe: 0.0 + ev, fx: 0.0 + delx * fpair };
}

export default {
  id: 'add-parameter',
  title: 'Add a parameter (soft core)',
  minutes: 20,
  summary: 'A per-pair parameter touches more of a pair style than the formula: storage, input, mixing, restart files, data files, and the extract() hook that lets fix adapt change it during a run. You add λ to a fresh clone of lj/cut and turn it into the soft-core potential used in free-energy calculations, finite even at r = 0. Eighteen two-atom energies and forces match a hand calculation, and melts plus a fix adapt run that switches the atoms off match LAMMPS\'s own lj/cut/soft, all to the last bit; a restart continues bit for bit and the data file carries λ.',
  learn: [
    'The soft core: E = λⁿ 4ε [1/D² − 1/D], D = α(1−λ)² + (r/σ)⁶',
    'The seven places a per-pair parameter lives in a pair style',
    'extract(): how fix adapt and FEP computes reach inside a style',
    'Restart and data files: what you must write, and how to check it',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: a parameter that switches atoms off</h2>
    <div class="prose"><p>Free-energy calculations make molecules appear or disappear gradually. Scaling LJ by λ does not work: at small λ atoms overlap freely and the r⁻¹² wall makes energies explode. The soft core (Beutler et al. 1994) moves the singularity away:</p>
<pre class="lmp-mini">E = λⁿ · 4ε [ 1/D² − 1/D ],    D = α(1 − λ)² + (r/σ)⁶          (n = 2, α = 0.5 here)
λ = 1: D = (r/σ)⁶, plain Lennard-Jones      λ = 0: E = 0, no interaction
r = 0: D = α(1−λ)², a finite energy for every λ &lt; 1</pre>
<p>LAMMPS has it as <code>lj/cut/soft</code> in the FEP package (${gh('src/FEP/pair_lj_cut_soft.cpp')}), with n and α on the pair_style line. Here n and α are fixed constants, so the only new input is λ, one per type pair, after σ: <code>pair_coeff 1 1 ε σ λ [cutoff]</code>. Where it has to go:</p>
<ol><li><b>Header</b>: a new array <code>double **lambda</code>.</li>
<li><b>allocate() / destructor</b>: create and destroy it.</li>
<li><b>coeff()</b>: read it (now 5 or 6 arguments), check 0 ≤ λ ≤ 1.</li>
<li><b>init_one()</b>: a mixing rule for i ≠ j (here: only equal λ may mix, else an error) and the derived constants lj1 = λⁿ, lj2 = σ⁶, lj3 = α(1−λ)², plus the shift and tail.</li>
<li><b>write_restart / read_restart</b>: one more double per type pair, in the same position in both.</li>
<li><b>write_data / write_data_all</b>: one more column.</li>
<li><b>extract("lambda")</b>: return the array, so fix adapt and compute fep can change it.</li></ol>
<p>Plus the formula in compute() and single(), and, as in E.3, the rRESPA functions and born_matrix removed.</p></div>
  </section>

  <section class="block">
    <h2>2 · Change</h2>
    <div class="prose"><p>A new clone of lj/cut (recipe of page E.2) named <code>lj/mine/soft</code>; then the header:</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.hdiff.split('\n').slice(2).join('\n'))}</pre>
    <div class="prose"><p>and the .cpp (48 lines added; 255 removed, mostly the deleted functions):</p></div>
    <pre class="lmp-mini gd">${esc(keyHunks(NATIVE.cppdiff))}</pre>
  </section>

  <section class="block">
    <h2>3 · Match: two atoms, three λ, six distances</h2>
    <div class="prose"><p>Atoms at x = 5 and 5 + r, ε = σ = 1, <code>run 0</code>, PotEng and the force on atom 1 with 17 digits. No pow() in the loop this time, so the hand check needs nothing special.</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Compute by hand</button></div>
    <div class="tbl ap-a"></div>
    <div class="prose ap-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Match: against lj/cut/soft, including fix adapt</h2>
    <div class="prose"><p>lj/cut/soft cloned as <code>lj/soft/ref</code> (the FEP package is not in either build), set to the same n = 2, α = 0.5: <code>pair_style lj/soft/ref 2.0 0.5 2.5</code>. The 4000-atom melt at λ = 0.6 with shift, with tail, and a third run that starts at λ = 1 and switches every pair off over 250 steps:</p>
<pre class="lmp-mini">variable        lam equal ramp(1.0,0.0)
fix             ad all adapt 1 pair lj/mine/soft lambda 1 1 v_lam</pre></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Compare</button></div>
    <div class="tbl ap-b"></div>
    <div class="prose ap-bsum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: restart and data files</h2>
    <div class="prose ap-c"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Energy of two atoms on top of each other (r = 0) at λ = 0.7?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>D = 0.5 × 0.3² = 0.045, E = 0.49 × 4 × (1/0.045² − 1/0.045) = 0.49 × 4 × (493.83 − 22.22) = 924.3. The run at r = 0.1 printed 924.30: the core is almost flat there (force −2.5 against −4.8 × 10¹⁴ for plain LJ).</div></div>
    <div class="ex"><p><b>2.</b> Types 1 and 2 have λ = 1 and λ = 0.5 and no <code>pair_coeff 1 2</code>. What happens?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>init_one stops with "different lambda values in mix": there is no sensible mixing rule for λ, so the style asks you to set the 1-2 pair explicitly. lj/cut/soft does the same. ε and σ still mix by the usual rule.</div></div>
    <div class="ex"><p><b>3.</b> You forget the extract() line. Which command fails, and with what?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>fix adapt (and compute fep) ask the style for a pointer to the named array; without it they stop with the message we got by asking for "alpha", which this style does not export: <code>${esc(NATIVE.adaptErr)}</code>.</div></div>
    <div class="ex"><p><b>4.</b> And if you write λ in write_restart but forget to read it in read_restart?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Every value after it is read one position too early: cut gets the λ bytes, the next type pair's setflag reads the first half of the real cutoff, and so on. Usually garbage or a crash, sometimes worse: plausible but wrong numbers. Write and read must mirror each other exactly; the restart check in section 5 is how you catch it.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = () => {
      let ok = 0;
      const rows = NATIVE.two.map(w => { const h = softPair(Number(w[1]), w[3]), s = h.pe === Number(w[5]) && h.fx === Number(w[7]); ok += s; return row([w[1], w[3], `<span class="mono">${esc(w[5])}</span>`, `<span class="mono">${esc(w[7])}</span>`, s ? OK : BAD]); });
      $('.ap-a').innerHTML = `<table>${head(['λ', 'r', 'PotEng (LAMMPS)', 'fx on atom 1 (LAMMPS)', 'both by hand'])}${rows.join('')}</table>`;
      $('.ap-asum').innerHTML = `<p>${ok} of ${NATIVE.two.length} identical to the last bit. At λ = 1 the soft core is plain LJ (r = 0.1 gives 4 × 10¹², E(σ) = 0, E(2σ) = −0.0615…), at λ = 0.7 the same overlap costs 924, at λ = 0.3 only 4.5, and the force near r = 0 goes to zero: atoms can pass through each other, which is the point.</p>`;
    };
    $('[data-k=b]').onclick = () => {
      const rows = []; let all = 0, same = 0;
      for (const [lab, a, b] of [['λ = 0.6, shift yes', 'mS', 'rS'], ['λ = 0.6, tail yes', 'mT', 'rT'], ['fix adapt, λ 1 → 0', 'mA', 'rA']]) NATIVE[a].forEach((r, k) => { const q = NATIVE[b][k], s = r.every((v, c) => v === q[c]); all++; same += s; if (k === 0 || k === NATIVE[a].length - 1 || a === 'mA') rows.push(row([k ? '' : lab, r[0], a === 'mA' ? r[1] : '', `<span class="mono">${esc(r[a === 'mA' ? 3 : 2])}</span>`, s ? OK : BAD])); });
      $('.ap-b').innerHTML = `<table>${head(['run', 'step', 'λ', 'E_vdwl, lj/mine/soft', 'same line in lj/soft/ref'])}${rows.join('')}</table>`;
      const a0 = NATIVE.mA[0], aN = NATIVE.mA[NATIVE.mA.length - 1];
      $('.ap-bsum').innerHTML = `<p>${same} of ${all} thermo lines identical in every column. The fix adapt run proves the extract() path: each step fix adapt evaluates v_lam, writes it into the lambda array it got from extract("lambda") and calls init_one again, in both styles alike. As λ falls the potential energy rises from ${Number(a0[3]).toFixed(3)} to ${Number(aN[3])}, and the temperature ends at ${Number(aN[2]).toFixed(2)}, against 1.66 for the same melt with plain LJ (page E.2): weakening the attraction while the atoms sit in their wells lifts potential energy that is paid out of kinetic energy (fix adapt changes the Hamiltonian, so total energy is not conserved). At λ = 1 the first line equals lj/cut's in E but not in the pressure's last digit: the soft-core arithmetic is different arithmetic for the same function.</p>`;
    };
    const R = NATIVE.restart;
    $('.ap-c').innerHTML = `<p>256 atoms at λ = 0.6: run 50, <code>write_restart</code>, run 50; then <code>clear</code>, <code>read_restart</code>, run 50. PotEng at step 50 when the second run starts: ${R.before} both times; at step 100: ${R.end1} and ${R.end2} ${R.end1 === R.end2 ? OK : BAD}. λ survived the binary file bit for bit. <code>write_data d.data pair ij</code> writes the new column (ε σ λ cutoff):</p><pre class="lmp-mini">PairIJ Coeffs # lj/mine/soft\n\n${esc(NATIVE.dataLine)}</pre><p>with %g, so λ values with more than 6 significant digits would be rounded there (page D.18).</p>`;
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
