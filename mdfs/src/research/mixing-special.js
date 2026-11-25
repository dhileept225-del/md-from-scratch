// Research · D.6: mixing rules and special bonds. mix_energy / mix_distance (including cutoffs), special_lj/special_coul weights, which pairs leave the neighbour list, and the special_bonds line that resets everything.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const grab = (lines, re) => { for (const l of lines) { const m = l.match(re); if (m) return m; } return null; };
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));
function g15(v) {
  if (v === 0) return '0';
  const ex = v.toExponential(14), e = +ex.split('e')[1];
  if (e < -4 || e >= 15) { let [m, x] = ex.split('e'); m = m.replace(/\.?0+$/, ''); const xs = Math.abs(+x); return `${m}e${+x < 0 ? '-' : '+'}${xs < 10 ? '0' : ''}${xs}`; }
  return v.toFixed(14 - e).replace(/(\.\d*?)0+$/, '$1').replace(/\.$/, '');
}
const near15 = (h, l) => g15(h) === l || Math.abs(h - +l) <= 6e-15 * Math.abs(h);

// ---------- A: mixing ----------
const T1 = { eps: 0.2381, sig: 3.405, cut: 8.5 }, T2 = { eps: 0.0914, sig: 3.95, cut: 12.0 };
const powint = (x, n) => { let yy = 1.0, ww = x, nn = n; for (; nn !== 0; nn >>= 1, ww *= ww) if (nn & 1) yy *= ww; return yy; };   // math_special.h
function mix(rule) {
  const e = (a, b, s1, s2) => rule === 'sixthpower' ? 2.0 * Math.sqrt(a * b) * powint(s1, 3) * powint(s2, 3) / (powint(s1, 6) + powint(s2, 6)) : Math.sqrt(a * b);
  const d = (a, b) => rule === 'geometric' ? Math.sqrt(a * b) : rule === 'arithmetic' ? 0.5 * (a + b) : Math.pow(0.5 * (powint(a, 6) + powint(b, 6)), 1.0 / 6.0);
  return { eps: e(T1.eps, T2.eps, T1.sig, T2.sig), sig: d(T1.sig, T2.sig), cut: d(T1.cut, T2.cut) };
}
function table(p, n, inner, outer) {
  const lj1 = 48.0 * p.eps * Math.pow(p.sig, 12.0), lj2 = 24.0 * p.eps * Math.pow(p.sig, 6.0), lj3 = 4.0 * p.eps * Math.pow(p.sig, 12.0), lj4 = 4.0 * p.eps * Math.pow(p.sig, 6.0), cutsq = p.cut * p.cut, out = [];
  for (let i = 0; i < n; i++) { const r = inner + (outer - inner) * i / (n - 1), rsq = r * r; let e = 0, f = 0; if (rsq < cutsq) { const r2 = 1.0 / rsq, r6 = r2 * r2 * r2; e = r6 * (lj3 * r6 - lj4); f = r6 * (lj1 * r6 - lj2) * r2 * r; } out.push([i + 1, r, e, f]); }
  return out;
}
const parseTable = (txt, key) => { const ls = txt.split('\n'), i = ls.findIndex(l => l.trim() === key); const o = []; for (let k = i + 3; k < ls.length; k++) { const w = ls[k].trim().split(/\s+/); if (w.length !== 4 || !/^\d+$/.test(w[0])) break; o.push(w); } return o; };
const scriptA = rule => `units           real
atom_style      atomic
region          box block 0 30 0 30 0 30 units box
create_box      2 box
mass            * 39.948
pair_style      lj/cut 10.0
pair_coeff      1 1 ${T1.eps} ${T1.sig} ${T1.cut}
pair_coeff      2 2 ${T2.eps} ${T2.sig} ${T2.cut}
pair_modify     mix ${rule}
pair_write      1 2 12 r 3.0 13.0 w.table MIX`;

// ---------- B: special bonds on a 4-atom chain ----------
const DATA = `four united atoms in a chain

4 atoms
3 bonds
1 atom types
1 bond types

0 20 xlo xhi
0 20 ylo yhi
0 20 zlo zhi

Masses

1 14.027

Atoms # full

1 1 1  0.2  9.0 10.2 10.0
2 1 1 -0.2 10.3  9.4 10.0
3 1 1 -0.2 11.6 10.2 10.0
4 1 1  0.2 12.9  9.4 10.0

Bonds

1 1 1 2
2 1 2 3
3 1 3 4
`;
const XYZ = [[9.0, 10.2, 10.0], [10.3, 9.4, 10.0], [11.6, 10.2, 10.0], [12.9, 9.4, 10.0]], Q = [0.2, -0.2, -0.2, 0.2];
const PAIRS = [[0, 1, 1], [1, 2, 1], [2, 3, 1], [0, 2, 2], [1, 3, 2], [0, 3, 3]];   // i, j, 1-2 / 1-3 / 1-4
const SB = {
  charmm: { line: 'special_bonds   charmm', lj: [0, 0, 0], co: [0, 0, 0] },
  amber: { line: 'special_bonds   amber', lj: [0, 0, 0.5], co: [0, 0, 5.0 / 6.0] },
  dreiding: { line: 'special_bonds   dreiding', lj: [0, 0, 1], co: [0, 0, 1] },
  fene: { line: 'special_bonds   fene', lj: [0, 1, 1], co: [0, 1, 1] },
  all: { line: 'special_bonds   lj/coul 1 1 1', lj: [1, 1, 1], co: [1, 1, 1] },
  half13: { line: 'special_bonds   lj/coul 0 0.5 1', lj: [0, 0.5, 1], co: [0, 0.5, 1] },
  oneline: { line: 'special_bonds   lj 0 0 0.5 coul 0 0 0.8333', lj: [0, 0, 0.5], co: [0, 0, 0.8333] },
  twolines: { line: 'special_bonds   lj 0 0 0.5\nspecial_bonds   coul 0 0 0.8333', lj: [0, 0, 0], co: [0, 0, 0.8333] },
};
const scriptB = k => `units           real
atom_style      full
read_data       chain.data
pair_style      lj/cut/coul/cut 10.0
pair_coeff      1 1 0.0914 3.95
bond_style      harmonic
bond_coeff      1 260.0 1.526
${SB[k].line}
thermo_style    custom step evdwl ecoul
thermo_modify   format float %.17g
run             0`;
function handB(k) {
  const w = SB[k], lj3 = 4.0 * 0.0914 * Math.pow(3.95, 12.0), lj4 = 4.0 * 0.0914 * Math.pow(3.95, 6.0), qq = 332.06371;
  let ev = 0, ec = 0, nn = 0; const rows = [];
  for (const [i, j, s] of PAIRS) {
    const fl = w.lj[s - 1], fc = w.co[s - 1], d = [0, 1, 2].map(c => XYZ[i][c] - XYZ[j][c]), rsq = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
    const kept = !(fl === 0 && fc === 0);
    const r2 = 1.0 / rsq, r6 = r2 * r2 * r2, elj = r6 * (lj3 * r6 - lj4), eco = qq * Q[i] * Q[j] * Math.sqrt(r2);
    if (kept) { nn++; ev += fl * elj; ec += fc * eco; }
    rows.push({ i, j, s, r: Math.sqrt(rsq), elj, eco, fl, fc, kept });
  }
  return { ev, ec, nn, rows };
}

export default {
  id: 'mixing-special',
  title: 'Mixing rules and special bonds',
  minutes: 20,
  summary: 'Two small functions decide most of a force field\'s cross terms: mix_energy and mix_distance fill every pair_coeff you did not type, cutoffs included; special_lj and special_coul scale or delete bonded neighbours. You reproduce mixed pair_write tables for three rules, predict LJ and Coulomb energies and the neighbour count of a 4-atom chain under eight special_bonds settings, and catch the one that silently undoes the line before it.',
  learn: [
    'geometric, arithmetic and sixthpower mixing, for ε, σ and the cutoff',
    'special_lj / special_coul: 1-2, 1-3, 1-4 weights and the bits that carry them',
    'Weight 0 for both LJ and Coulomb removes the pair from the neighbour list',
    'Each special_bonds command starts from all zeros',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: mixing</h2>
    <div class="prose"><p>When <code>pair_coeff 1 2</code> is missing, init_one (page D.5) asks ${gh('src/pair.cpp', 720)}:</p>
<pre class="lmp-mini">mix_energy:   geometric, arithmetic:  sqrt(eps1*eps2)
              sixthpower:  2 sqrt(eps1 eps2) σ1³ σ2³ / (σ1⁶ + σ2⁶)
mix_distance: geometric:   sqrt(sig1*sig2)
              arithmetic:  0.5 * (sig1+sig2)        // Lorentz–Berthelot = arithmetic
              sixthpower:  (0.5 (σ1⁶ + σ2⁶))^(1/6)
cut[i][j] = mix_distance(cut[i][i], cut[j][j]);     // the cutoff is mixed too</pre>
<p>Type 1: argon (ε 0.2381, σ 3.405, rc 8.5). Type 2: a CH₂ bead (ε 0.0914, σ 3.95, rc 12). No <code>pair_coeff 1 2</code>.</p></div>
    <div class="row"><label class="lc-in">pair_modify mix<select class="lc-u ms-r"><option>geometric</option><option>arithmetic</option><option>sixthpower</option></select></label><button class="btn sm primary" data-k="a">▶ pair_write 1 2</button><span class="lmp-status muted ms-sta"></span></div>
    <div class="calc ms-mix"></div>
    <div class="tbl ms-outa"></div>
    <div class="prose ms-suma"></div>
  </section>

  <section class="block">
    <h2>2 · Understand: special bonds</h2>
    <div class="prose"><p>Atoms 1 and 2 bonded are a <b>1-2</b> pair, separated by one more bond <b>1-3</b>, two more <b>1-4</b>. The neighbour list stores j with two extra bits saying which (${gh('src/pair.h', 306)}: <code>sbmask(j) = j &gt;&gt; 30 &amp; 3</code>), and every pair style multiplies:</p>
<pre class="lmp-mini">factor_lj   = special_lj[sbmask(j)];      // [0] = 1 for ordinary pairs
factor_coul = special_coul[sbmask(j)];
evdwl *= factor_lj;   ecoul = factor_coul * qqrd2e * qi * qj / r;</pre>
<p>${gh('src/neighbor.cpp', 538)}: if both weights of a kind are 0, the pair is <b>not stored at all</b>; if both are 1, it is stored as an ordinary pair; otherwise it is stored with its bits (and always stored when a KSpace style is active, which must subtract the full Coulomb). And ${gh('src/force.cpp', 704)} begins every <code>special_bonds</code> command with:</p>
<pre class="lmp-mini">special_lj[1] = special_lj[2] = special_lj[3] = 0.0;
special_coul[1] = special_coul[2] = special_coul[3] = 0.0;</pre></div>
  </section>

  <section class="block">
    <h2>3 · Match: a 4-atom chain</h2>
    <div class="prose"><p>Four united atoms in a zigzag (bonds 1.526 Å), charges +0.2, −0.2, −0.2, +0.2, <code>lj/cut/coul/cut 10</code>. Every pair here is bonded: three 1-2, two 1-3, one 1-4. Pick a setting and predict E_vdwl, E_coul and the number of neighbours.</p></div>
    <div class="row"><label class="lc-in">special_bonds<select class="lc-u ms-sb">${Object.entries(SB).map(([k, v]) => `<option value="${k}">${esc(v.line.replace(/special_bonds\s+/g, '').replace('\n', ' ; then '))}</option>`).join('')}</select></label><button class="btn sm primary" data-k="b">▶ run 0</button><span class="lmp-status muted ms-stb"></span></div>
    <div class="tbl ms-pairs"></div>
    <div class="prose ms-outb"></div>
  </section>

  <section class="block prose">
    <h2>4 · Observe</h2>
    <ul>
      <li><b>Two lines are not one.</b> <code>special_bonds lj 0 0 0.5</code> followed by <code>special_bonds coul 0 0 0.8333</code> leaves LJ 1-4 at <b>0</b>: the second command reset it. Put both keywords on one line, or use a preset (amber).</li>
      <li><b>Mixing mixes cutoffs.</b> With per-type cutoffs 8.5 and 12, geometric mixing gives the 1-2 pair a cutoff of 10.0995 Å, arithmetic 10.25 Å. If you meant one global cutoff, do not give per-type ones.</li>
      <li><b>sixthpower changes ε too.</b> For these two types ε drops from 0.1475 (geometric) to 0.1340 kcal/mol, because 2σ₁³σ₂³/(σ₁⁶ + σ₂⁶) &lt; 1 when the σ differ; σ becomes 3.7266 and the cutoff 10.9048 Å.</li>
      <li><b>All 1 1 1 is not "no special bonds".</b> At 1.53 Å apart, LJ between two bonded atoms is about 3 × 10⁴ kcal/mol: that is why force fields exclude 1-2 and 1-3 and let the bond term do the work.</li>
    </ul>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Geometric mixing of σ 3.405 and 3.95, and of ε 0.2381 and 0.0914?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>σ = √(3.405 × 3.95) = <b>3.6674 Å</b>; ε = √(0.2381 × 0.0914) = <b>0.14752 kcal/mol</b>.</div></div>
    <div class="ex"><p><b>2.</b> With <code>special_bonds amber</code>, how many pairs of this chain are in the neighbour list, and why?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><b>1</b>: the 1-4 pair (LJ 0.5, Coulomb 5/6). 1-2 and 1-3 have both weights 0, so they are never stored.</div></div>
    <div class="ex"><p><b>3.</b> You add <code>kspace_style pppm</code> to a CHARMM system. Do 1-2 pairs come back into the list?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Yes, with weight 0: PPPM computes Coulomb between all charges, bonded or not, so the pair style must subtract the excluded part in real space; it needs those pairs in the list (${gh('src/neighbor.cpp', 556)}).</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, res) => { st.innerHTML = `<span class="lmp-err">${esc(res.error || 'no output')}</span>`; };
    const showMix = () => { const r = $('.ms-r').value, m = mix(r); $('.ms-mix').textContent = `${r}:  ε₁₂ = ${s17(m.eps)}   σ₁₂ = ${s17(m.sig)}   rc₁₂ = ${s17(m.cut)}`; $('.ms-outa').innerHTML = ''; $('.ms-suma').innerHTML = ''; return m; };
    $('.ms-r').onchange = showMix; showMix();
    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.ms-sta'), r = $('.ms-r').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptA(r), { read: ['w.table'] }); b.disabled = false;
      const t = res.files?.['w.table']; if (res.error || !t) return err(st, res); st.textContent = '';
      const m = showMix(), H = table(m, 12, 3.0, 13.0), L = parseTable(t, 'MIX'), gen = grab(res.lines, /Generated (\d+) of (\d+) mixed pair_coeff terms from (\S+) mixing rule/); let ok = 0, n = 0;
      $('.ms-outa').innerHTML = `<table>${head(['r', 'E by hand', 'F by hand', 'LAMMPS'])}${H.map((h, i) => { const l = L[i] || [], s = [near15(h[2], l[2]), near15(h[3], l[3])]; s.forEach(x => { n++; if (x) ok++; }); return row([g15(h[1]), g15(h[2]), g15(h[3]), s.every(Boolean) ? OK : `${BAD} ${l.slice(2).join(' ')}`]); }).join('')}</table>`;
      $('.ms-suma').innerHTML = `<p>${ok} of ${n} values match to 15 digits (as on page D.5, the browser's pow may change the 15th). The table is zero beyond the mixed cutoff ${g15(m.cut)} Å. Log: <i>${gen ? esc(gen[0]) : ''}</i>.</p>`;
    };

    const showB = () => {
      const k = $('.ms-sb').value, H = handB(k), w = SB[k];
      $('.ms-pairs').innerHTML = `<table>${head(['pair', 'kind', 'r (Å)', 'E_LJ', 'E_coul', 'weights LJ / coul', 'in list?'])}${H.rows.map(p => row([`${p.i + 1}–${p.j + 1}`, `1-${p.s + 1}`, p.r.toFixed(4), p.elj.toPrecision(6), p.eco.toPrecision(6), `${+p.fl.toFixed(4)} / ${+p.fc.toFixed(4)}`, p.kept ? 'yes' : 'no'])).join('')}</table><p class="muted">Effective weights after the whole input: LJ ${w.lj.map(x => +x.toFixed(4)).join(' ')}, Coulomb ${w.co.map(x => +x.toFixed(4)).join(' ')}. Prediction: E_vdwl = ${s17(H.ev)}, E_coul = ${s17(H.ec)}, Total # of neighbors = ${H.nn}.</p>`;
      $('.ms-outb').innerHTML = ''; return H;
    };
    $('.ms-sb').onchange = showB; showB();
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.ms-stb'), k = $('.ms-sb').value; b.disabled = true; st.textContent = 'Running the real LAMMPS…';
      const res = await runLammps(scriptB(k), { files: { 'chain.data': DATA } }); b.disabled = false;
      if (res.error || !res.thermo.rows.length) return err(st, res); st.textContent = '';
      const H = showB(), c = res.thermo.cols, r0 = res.thermo.rows[0], ev = parseFloat(r0[c.indexOf('E_vdwl')]), ec = parseFloat(r0[c.indexOf('E_coul')]), nb = grab(res.lines, /Total # of neighbors = (\d+)/);
      const cl = (h, l) => h === l || Math.abs(h - l) <= 1e-14 * Math.max(Math.abs(h), 1e-300);
      $('.ms-outb').innerHTML = `<p>LAMMPS: E_vdwl = ${s17(ev)} ${cl(H.ev, ev) ? OK : BAD}, E_coul = ${s17(ec)} ${cl(H.ec, ec) ? OK : BAD}, Total # of neighbors = ${nb ? nb[1] : '?'} ${nb && +nb[1] === H.nn ? OK : BAD}.${k === 'twolines' ? ' The LJ 1-4 weight is 0: the second line reset it.' : ''}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
