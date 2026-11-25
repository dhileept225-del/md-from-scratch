// Research · E.3: change the formula. lj/mine (the clone from E.2) turned into a Mie n-m potential with the exponents
// on the pair_style line; two-atom energies and forces checked by hand with a correctly rounded pow; a melt checked
// bit for bit against LAMMPS's own mie/cut; the price of pow().
import { NATIVE } from './mie-native.js';
import { crpow } from './crpow.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

// show the parts of the .cpp diff that matter (compute, init_one, settings), skip the deleted rRESPA/born_matrix bodies
function keyHunks(diff) {
  const hunks = diff.split(/\n(?=@@)/).slice(1);
  return hunks.map(h => { const k = h.search(/\n.*(\/\* -{10,}[^\n]*\n[^\n]*\n?[^\n]*)?born_matrix|\n-void PairLJMine::compute_inner/); return k > 0 ? h.slice(0, k) : h; })
    .filter(h => /pow|n_rep|m_att|cmie|evdwl = \(lj3|philj = \(lj3/.test(h) && !/^@@[^\n]*\n[^\n]*compute_(inner|middle|outer)/.test(h))
    .map(h => h.split('\n').filter(l => !/^ \s*$/.test(l)).join('\n')).join('\n');
}
const cmieOf = (n, m) => n / (n - m) * Math.pow(n / m, m / (n - m));

export default {
  id: 'change-formula',
  title: 'Change the formula (Mie potential)',
  minutes: 20,
  summary: 'The clone from E.2 becomes a different potential: Mie n–m, E = Cε[(σ/r)ⁿ − (σ/r)ᵐ], with the exponents given on the pair_style line. You predict its constants by hand, read every place the formula lives (compute, single, the coefficients, the shift, the tail correction, the restart file) and the two features that had to go, then check: eight two-atom energies and forces against a correctly rounded hand calculation, and a 4000-atom melt against LAMMPS\'s own mie/cut, both to the last bit. The price shows in the timing table: pow() makes the pair loop 3.5 times slower than lj/cut\'s r⁻⁶·r⁻⁶.',
  learn: [
    'Mie n–m: the constant C, the minimum at σ(n/m)^(1/(n−m)), depth −ε',
    'Every copy of a formula in a pair style: compute, single, init_one, tail, restart',
    'Checking against an independent implementation, bit for bit',
    'Why lj/cut never calls pow() in its loop',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: the formula and its constants</h2>
    <div class="prose"><p>The Mie potential generalises 12–6:</p>
<pre class="lmp-mini">E(r) = C ε [ (σ/r)ⁿ − (σ/r)ᵐ ],   C = n/(n−m) · (n/m)^(m/(n−m))
F(r)·r = C ε [ n (σ/r)ⁿ − m (σ/r)ᵐ ]          (what compute() needs: fpair = F/r)</pre>
<p>C is chosen so that the well depth is −ε for any n, m, at r_min = σ (n/m)^(1/(n−m)). For n = 12, m = 6: C = 2 · 2¹ = 4, r_min = 2<sup>1/6</sup>σ, plain Lennard-Jones. Softer repulsion (smaller n) is common in coarse-grained models (SAFT-γ Mie uses n from 8 to 30 with m = 6).</p></div>
    <div class="tbl mf-c"></div>
  </section>

  <section class="block">
    <h2>2 · Change: every place the formula lives</h2>
    <div class="prose"><p>Usage becomes <code>pair_style lj/mine cutoff [n m]</code>, defaults 12 and 6. The header gets the two exponents and loses four functions:</p></div>
    <pre class="lmp-mini gd">${esc(NATIVE.hdiff.split('\n').slice(2).join('\n'))}</pre>
    <div class="prose"><p>In the .cpp (39 lines added, 252 removed, most of them the deleted functions), the parts that carry the formula:</p></div>
    <pre class="lmp-mini gd mf-diff">${esc(keyHunks(NATIVE.cppdiff))}</pre>
    <div class="prose"><p>The checklist a new formula needs, and what was done here:</p>
<ul><li><b>compute()</b> and <b>single()</b>: the same lines twice (single is used by pair_write, compute pair/local and others).</li>
<li><b>init_one()</b>: the coefficients lj1–lj4, the <b>shift</b> (offset at the cutoff) and the <b>tail correction</b>, which for Mie is ∫ r²E(r) dr from r_c to ∞ (and the matching pressure integral), written in the same form as mie/cut.</li>
<li><b>settings()</b> reads n and m; <b>write/read_restart_settings</b> must store them, or a restarted run silently falls back to 12–6.</li>
<li><b>Removed</b>: the three rRESPA functions (respa_enable = 0) and born_matrix (born_matrix_enable = 0), which still contained 12–6. Leaving them would give wrong numbers without an error; removing them gives a clear error if anyone asks for them.</li></ul>
<p>Native rebuild: about 10 s.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: two atoms, by hand</h2>
    <div class="prose"><p>Atoms at x = 5 and 5 + r, <code>pair_style lj/mine 2.5 14 6</code>, ε = σ = 1, <code>run 0</code> and print PotEng and the force on atom 1 with 17 digits, for eight distances. The hand code repeats init_one and compute; its pow is a correctly rounded one built from BigInt arithmetic (glibc's pow, used by native LAMMPS, is correctly rounded in practice; browsers' Math.pow is not always).</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Compute by hand</button></div>
    <div class="tbl mf-a"></div>
    <div class="prose mf-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Match: a melt against LAMMPS's own mie/cut</h2>
    <div class="prose"><p>LAMMPS has a Mie style, <code>mie/cut</code> in the EXTRA-PAIR package, with per-pair exponents: <code>pair_coeff 1 1 ε σ n m</code>. It is an independent implementation to test against. (The browser build lacks EXTRA-PAIR. Copying its two files into src/ is how the old make-based build installs a package, and CMake refuses that: "Found package(s) installed by the make-based build system". So it was cloned under the name mie/ref, as on page E.2.) The 4000-atom melt, 250 steps, with <code>pair_modify shift yes</code> and, separately, <code>tail yes</code>:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Compare</button></div>
    <div class="tbl mf-b"></div>
    <div class="prose mf-bsum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: what pow() costs</h2>
    <div class="tbl mf-t"></div>
    <div class="prose mf-tsum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Where is the minimum of the 14–6 potential, and what is the energy there?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>r_min = (14/6)^(1/8) = 1.11178σ, E = −ε. The run at r = 1.1118 printed E = −0.9999998 and a force of 0.005 (r is a little beyond the minimum, so the force is attractive: positive on atom 1, which sits at the lower x).</div></div>
    <div class="ex"><p><b>2.</b> You forgot to add n and m to write_restart_settings. What happens after read_restart?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The constructor sets 12 and 6, nothing overwrites them, and the run continues with Lennard-Jones: no error, different physics. With the two lines added, a restart at step 100 continues bit for bit like a second <code>run</code> command (we checked: PotEng −4.0841443568763918 at step 100 either way). Note that neigh_modify settings are not in a restart file and must be given again (documented in read_restart).</div></div>
    <div class="ex"><p><b>3.</b> How could lj/mine keep pow() out of the loop for integer exponents?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Precompute in init_one whether n/2 and m/2 are integers and use repeated multiplication (r2inv·r2inv·…), as lj/cut does for 6, or a small switch over common exponents. That is what makes lj/cut 3.5 times faster than lj/mine 12 6, and also why the two differ in the last bits.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    // constants table, by hand
    $('.mf-c').innerHTML = `<table>${head(['n–m', 'C', 'r_min / σ', 'E(r_min) / ε'])}${[[12, 6], [14, 6], [16, 6], [24, 12], [9, 6]].map(([n, m]) => { const C = cmieOf(n, m), r = Math.pow(n / m, 1 / (n - m)); return row([`${n}–${m}`, C.toFixed(5), r.toFixed(5), (C * (Math.pow(r, -n) - Math.pow(r, -m))).toFixed(6)]); }).join('')}</table>`;

    $('[data-k=a]').onclick = () => {
      const n = 14, m = 6, cmie = (n / (n - m) * crpow(n / m, 3, 4));       // pow(14/6, 6/8)
      const lj1 = cmie * n * 1.0 * 1.0, lj2 = cmie * m * 1.0 * 1.0, lj3 = cmie * 1.0 * 1.0, lj4 = cmie * 1.0 * 1.0;
      let ok = 0, jsDiff = 0;
      const rows = NATIVE.two.map(w => {
        const r = w[1], x2 = 5.0 + Number(r), delx = 5.0 - x2, rsq = delx * delx + 0.0 * 0.0 + 0.0 * 0.0, r2inv = 1.0 / rsq;
        const rm = crpow(r2inv, 3), rn = crpow(r2inv, 7), fpair = 1.0 * (lj1 * rn - lj2 * rm) * r2inv, fx = 0.0 + delx * fpair;
        let ev = (lj3 * rn - lj4 * rm) - 0.0; ev *= 1.0; const pe = 0.0 + ev;
        if (Math.pow(r2inv, 3) !== rm || Math.pow(r2inv, 7) !== rn) jsDiff++;
        const same = pe === Number(w[5]) && fx === Number(w[7]); ok += same;
        return row([r, `<span class="mono">${esc(w[5])}</span>`, `<span class="mono">${esc(w[7])}</span>`, same ? OK : BAD]);
      });
      $('.mf-a').innerHTML = `<table>${head(['r', 'PotEng (LAMMPS)', 'fx on atom 1 (LAMMPS)', 'both by hand'])}${rows.join('')}</table>`;
      $('.mf-asum').innerHTML = `<p>${ok} of ${NATIVE.two.length} distances identical to the last bit, with C = ${s17(cmie)} (1.75 × (7/3)^0.75). At r = σ the energy is exactly 0, as for any n–m. In this browser, Math.pow differs from the correctly rounded value for ${jsDiff} of the ${NATIVE.two.length} distances${jsDiff ? ': with it, the hand check would miss in the last bit there' : ''}.</p>`;
    };

    $('[data-k=b]').onclick = () => {
      const rows = []; let all = 0, same = 0;
      for (const [lab, a, b] of [['shift yes', 'mineS', 'refS'], ['tail yes', 'mineT', 'refT']]) NATIVE[a].forEach((r, k) => { const q = NATIVE[b][k], s = r.every((v, c) => v === q[c]); all++; same += s; rows.push(row([k ? '' : lab, r[0], `<span class="mono">${esc(r[2])}</span>`, `<span class="mono">${esc(q[2])}</span>`, s ? OK : BAD])); });
      $('.mf-b').innerHTML = `<table>${head(['', 'step', 'E_vdwl, lj/mine 2.5 14 6', 'E_vdwl, mie/ref, pair_coeff … 14 6', 'all 5 columns'])}${rows.join('')}</table>`;
      const a = NATIVE.mine126, c = NATIVE.cut, la = a[a.length - 1], lc = c[c.length - 1];
      $('.mf-bsum').innerHTML = `<p>${same} of ${all} thermo lines identical in all five columns, including the shifted energies and the tail-corrected energy and pressure: lj/mine does the same arithmetic as mie/cut in the same order.</p><p>Now <code>lj/mine 2.5</code> (12–6) against <code>lj/cut</code>: the same potential, but step 0 already differs in the pressure's last digit, and after 250 steps E_vdwl is ${s17(Number(la[2]))} against ${s17(Number(lc[2]))} (${(Number(la[2]) - Number(lc[2])).toExponential(1)} apart). pow(r2inv, 6.0) and r6inv·r6inv round differently, and the trajectory carries the difference forward. Equal physics, not equal bits: which is why the checks above compare with an implementation that also uses pow().</p>`;
    };

    const T = NATIVE.pairTime;
    $('.mf-t').innerHTML = `<table>${head(['pair style', 'Pair time, 1000 steps, 4000 atoms (s)', 'relative'])}${Object.entries(T).map(([k, v]) => row([`<code>${k}</code>`, v.toFixed(3), (v / T['lj/cut']).toFixed(2)])).join('')}</table>`;
    $('.mf-tsum').innerHTML = '<p>Native, one core, best of three. Two pow() calls per pair cost far more than the rest of the loop: lj/mine is 3.5 times slower than lj/cut even for 12–6, where the answer is the same. mie/cut pays the same price. A production style would special-case common exponents (exercise 3).</p>';

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
