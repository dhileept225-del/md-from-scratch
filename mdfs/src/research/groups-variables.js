// Research · D.17: groups, regions and variables. Region inside() tests and group counts reproduced from dumped positions,
// the boundary cases (shared faces, a sphere through a neighbour shell), and the variable parser: precedence, left-associative ^,
// unary minus above ^, and the 15-digit rounding of ${name} and of v_name inside another formula.
import { runLammps } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

// regions: [name, LAMMPS args, hand test (lattice constant a, position p after remap)]
const R2 = Math.SQRT1_2;   // 0.70710678118654757, the fcc neighbour distance in lattice units
const sph = (c, r) => (a, p) => { const d = [0, 1, 2].map(k => p[k] - a * c[k]); return Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]) <= a * r; };
const REG = [
  ['left', 'block 0 2.5 INF INF INF INF', (a, p) => p[0] >= 0 * a && p[0] <= 2.5 * a],
  ['right', 'block 2.5 5 INF INF INF INF', (a, p) => p[0] >= 2.5 * a && p[0] <= 5 * a],
  ['ball', 'sphere 2.5 2.5 2.5 1.0', sph([2.5, 2.5, 2.5], 1.0)],
  ['rod', 'cylinder z 2.5 2.5 1.0 INF INF', (a, p) => { const dx = p[0] - 2.5 * a, dy = p[1] - 2.5 * a; return Math.sqrt(dx * dx + dy * dy) <= 1.0 * a; }],
  ['outside', 'sphere 2.5 2.5 2.5 1.0 side out', (a, p) => !sph([2.5, 2.5, 2.5], 1.0)(a, p)],
  ['ball1', 'sphere 2 2 2 1.0', sph([2, 2, 2], 1.0)],
  ['ballnn', `sphere 2 2 2 ${s17(R2)}`, sph([2, 2, 2], R2)],
];
const GROUPS = [...REG.map(r => [r[0], `region ${r[0]}`]), ['both', 'intersect left right'], ['lr', 'union left right'], ['ring', 'subtract rod ball']];
const GSCRIPT = `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 5 0 5 0 5
create_box      1 box
create_atoms    1 box
mass            1 1.0
${REG.map(r => `region          ${r[0].padEnd(8)}${r[1]}`).join('\n')}
${GROUPS.map(g => `group           ${g[0].padEnd(8)}${g[1]}`).join('\n')}
write_dump      all custom x.dump id x y z modify format float %.17g
print           "XCM $(xcm(ball,x):%.17g) $(xcm(ballnn,y):%.17g) A $(xlat:%.17g) L $(lx:%.17g)"`;

// ---------- a hand copy of the variable parser (Variable::evaluate, operators only) ----------
const PREC = { '||': 1, '|^': 1, '&&': 2, '==': 3, '!=': 3, '<': 4, '<=': 4, '>': 4, '>=': 4, '+': 5, '-': 5, '*': 6, '/': 6, '%': 6, '^': 7, 'u-': 8, '!': 8 };
function evaluate(str, vars) {
  const ops = [], args = [];
  const apply = op => {
    if (op === 'u-') { args.push(-args.pop()); return; }
    if (op === '!') { args.push(args.pop() === 0.0 ? 1.0 : 0.0); return; }
    const b = args.pop(), a = args.pop();
    switch (op) {
      case '+': args.push(a + b); break; case '-': args.push(a - b); break; case '*': args.push(a * b); break;
      case '/': if (b === 0.0) throw new Error('Divide by 0'); args.push(a / b); break;
      case '%': if (b === 0.0) throw new Error('Modulo 0'); args.push(a % b); break;          // fmod
      case '^': if (b === 0.0) throw new Error('Power by 0'); if (a < 0 && Math.round(b) !== b) throw new Error('Invalid power'); args.push(Math.pow(a, b)); break;
      case '==': args.push(a === b ? 1 : 0); break; case '!=': args.push(a !== b ? 1 : 0); break;
      case '<': args.push(a < b ? 1 : 0); break; case '<=': args.push(a <= b ? 1 : 0); break;
      case '>': args.push(a > b ? 1 : 0); break; case '>=': args.push(a >= b ? 1 : 0); break;
      case '&&': args.push(a !== 0 && b !== 0 ? 1 : 0); break; case '||': args.push(a !== 0 || b !== 0 ? 1 : 0); break;
      case '|^': args.push((a === 0) !== (b === 0) ? 1 : 0); break;
    }
  };
  let i = 0, expectArg = true;
  while (i < str.length) {
    const c = str[i];
    if (c === ' ') { i++; continue; }
    if (c === '(') { const start = ++i; let depth = 1; while (depth) { if (str[i] === '(') depth++; if (str[i] === ')') depth--; i++; } args.push(evaluate(str.slice(start, i - 1), vars)); expectArg = false; continue; }
    if (/[0-9.]/.test(c)) { const m = str.slice(i).match(/^(\d*\.?\d*(?:[eE][-+]?\d+)?)/)[1]; args.push(Number(m)); i += m.length; expectArg = false; continue; }
    if (str.startsWith('v_', i)) { const m = str.slice(i + 2).match(/^\w+/)[0]; args.push(Number(vars[m].toPrecision(15))); i += 2 + m.length; expectArg = false; continue; }   // retrieve(): "%.15g", then stod
    let op = ['||', '|^', '&&', '==', '!=', '<=', '>=', '<', '>', '+', '-', '*', '/', '%', '^', '!'].find(o => str.startsWith(o, i)); i += op.length;
    if (expectArg) { if (op === '-') op = 'u-'; ops.push(op); continue; }      // unary minus and NOT are pushed, not applied
    while (ops.length && PREC[ops[ops.length - 1]] >= PREC[op]) apply(ops.pop());   // >= : every operator is left-associative
    ops.push(op); expectArg = true;
  }
  while (ops.length) apply(ops.pop());
  return args[0];
}
const EXPR = ['-2^2', '2^3^2', '2*3^2', '-2*-3', '7%3', '-7%3', '2^-1', '1-2-3', '8/4/2', '!0+1', '1<2==1', '3>2>1', '1+1==2&&0||1'];
const VSCRIPT = `variable        a equal 0.1+0.2
variable        b equal 1/3
variable        c equal \${b}
variable        d equal v_b*3
print           "A \${a} $(v_a) $(v_a==0.3) $(0.1+0.2==0.3)"
print           "B \${b} $(v_b) $(v_b*3==1) $((1/3)*3==1) $(v_d)"
print           "C \${c} $(v_c)"
region          box block 0 1 0 1 0 1
create_box      1 box
mass            1 1.0
thermo_style    custom step v_a v_b v_d
thermo_modify   format float %.20g
run             0`;
const ESCRIPT = `print           "E ${EXPR.map(e => `$(${e})`).join(' ')}"`;

export default {
  id: 'groups-variables',
  title: 'Groups, regions and variables',
  minutes: 18,
  summary: 'Regions are a yes/no test on a point; groups are a bit in a 32-bit mask, set once; variables are formulas parsed by LAMMPS\'s own little calculator. You reproduce every group count of an fcc crystal from the inside() tests (and find atoms counted twice on a shared face, and a sphere through a neighbour shell that catches 3 of 12), then a hand copy of the variable parser that agrees with LAMMPS on −2^2 = 4 and 2^3^2 = 64, and the places where a value is rounded to 15 digits: ${a}, and v_a inside any formula, $(v_a) included.',
  learn: [
    'Region inside() tests: inclusive faces, ≤ on a computed radius, lattice units',
    'Groups: a bitmask set at definition; union, intersect, subtract',
    'The variable parser: precedence, left-associative ^, unary minus above ^',
    'Where values pass through a 15-digit string, and where they do not',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: three small machines</h2>
    <div class="prose"><p><b>Regions</b> answer one question, inside(x, y, z). Region::match (${gh('src/region.cpp', 143)}) first remaps the point into the periodic box, then flips the answer for <code>side out</code>. A block is <code>x &gt;= xlo &amp;&amp; x &lt;= xhi …</code> (${gh('src/region_block.cpp', 288)}), both faces included; a sphere is <code>sqrt(dx² + dy² + dz²) &lt;= radius</code>. Numbers are in <b>lattice units</b> unless you write <code>units box</code>: 2.5 means 2.5 × xlat.</p>
<p><b>Groups</b> are one bit each in atom->mask, at most 32 including <code>all</code> (group.h:26). <code>group g region r</code> sets the bit for atoms inside r <i>now</i>; it is not re-checked later unless the group is <code>dynamic</code>.</p>
<p><b>Variables</b>: an equal-style variable stores the formula as text and evaluates it each time it is used (${gh('src/variable.cpp', 1081)}). How you use it decides what you get:</p>
<pre class="lmp-mini">\${a}          text substitution when the line is read: the value printed with %.15g
$(expr)       immediate, printed with %.20g (or $(expr:%.17g))
v_a in thermo, fix ave/time, dump …      the double itself (compute_equal)
v_a inside any formula, even $(v_a)      retrieve(): "%.15g" string, then stod  (variable.cpp:2173)</pre>
<p>The parser is a two-stack shunting-yard (${gh('src/variable.cpp', 2548)}): before an operator is pushed, every operator on the stack with <b>equal or higher</b> precedence is applied. Precedence: || |^ (1) &lt; &amp;&amp; (2) &lt; == != (3) &lt; &lt; &lt;= &gt; &gt;= (4) &lt; + − (5) &lt; * / % (6) &lt; ^ (7) &lt; unary −, ! (8).</p></div>
  </section>

  <section class="block">
    <h2>2 · Match: every group count, by hand</h2>
    <div class="prose"><p>A 5 × 5 × 5 fcc crystal (500 atoms), seven regions and ten groups. The hand code applies each inside() test to the dumped positions.</p></div>
    <pre class="lmp-mini">${esc(GSCRIPT)}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Run and count by hand</button><span class="lmp-status muted gv-sta"></span></div>
    <div class="tbl gv-a"></div>
    <div class="prose gv-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Match: the variable parser</h2>
    <div class="prose"><p>Predict each value before you run. The hand parser is 40 lines in this page's source and follows the rules above.</p></div>
    <pre class="lmp-mini">${esc(ESCRIPT)}\n${esc(VSCRIPT)}</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Run and compare</button><span class="lmp-status muted gv-stb"></span></div>
    <div class="tbl gv-b"></div>
    <div class="prose gv-bsum"></div>
    <div class="tbl gv-c"></div>
    <div class="prose gv-csum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Two regions <code>block 0 2.5 …</code> and <code>block 2.5 5 …</code> in this crystal: how many atoms are in neither?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>None: every atom is in at least one (lr = 500), and the 50 on the plane x = 2.5 are in both. To split a crystal into halves use <code>block 0 2.49 …</code> or give one block a slightly smaller upper bound.</div></div>
    <div class="ex"><p><b>2.</b> Why does <code>left</code> have 300 atoms and <code>right</code> only 250, if both include both faces?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The plane x = 5 is the same as x = 0 in a periodic box, and the atoms there are stored at x = 0. left gets x = 0 and x = 2.5 (two planes of 50 at the faces); right gets x = 2.5 but nothing at x = 5.</div></div>
    <div class="ex"><p><b>3.</b> How do you write −2² the mathematical way?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>-(2^2)</code> or <code>0-2^2</code>. And 2<sup>3²</sup> is <code>2^(3^2)</code>.</div></div>
    <div class="ex"><p><b>4.</b> A loop stops when <code>v_t &gt;= v_tmax</code>, both equal-style. Can the 15-digit rounding bite?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Yes, near equality: both are rounded to 15 digits before the comparison, so values that differ only past the 15th digit compare equal. For counters and times built from integers, compare integers (step numbers), not accumulated doubles.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    const err = (st, e) => { st.innerHTML = `<span class="lmp-err">${esc(e || 'no output')}</span>`; };

    $('[data-k=a]').onclick = async e => {
      const b = e.currentTarget, st = $('.gv-sta'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(GSCRIPT, { read: ['x.dump'] }); b.disabled = false;
      if (res.error || !res.files['x.dump']) return err(st, res.error); st.textContent = '';
      const lmp = {}; for (const l of res.lines) { const m = l.match(/^(\d+) atoms in group (\w+)/); if (m) lmp[m[2]] = +m[1]; }
      const xl = res.lines.find(l => l.startsWith('XCM')).split(/\s+/), a = +xl[4], L = +xl[6];
      const D = res.files['x.dump'].split('\n'), n = +D[3], P = [];
      for (let k = 0; k < n; k++) { const w = D[9 + k].trim().split(/\s+/).map(Number); P.push([w[1], w[2], w[3]].map(v => (v >= L ? v - L : v < 0 ? v + L : v))); }
      const mem = {}; for (const [name, , test] of REG) mem[name] = P.map(p => test(a, p));
      mem.both = P.map((_, i) => mem.left[i] && mem.right[i]); mem.lr = P.map((_, i) => mem.left[i] || mem.right[i]); mem.ring = P.map((_, i) => mem.rod[i] && !mem.ball[i]);
      let ok = 0;
      const rows = GROUPS.map(([g, def]) => { const h = mem[g].filter(Boolean).length; ok += h === lmp[g]; return row([`<code>${g}</code>`, esc(def.startsWith('region') ? REG.find(r => r[0] === g)[1] : def), lmp[g], h, h === lmp[g] ? OK : BAD]); });
      $('.gv-a').innerHTML = `<table>${head(['group', 'definition', 'LAMMPS', 'hand', ''])}${rows.join('')}</table>`;
      // xcm by hand (Group::xcm: sum of x·m in memory order, divided by the group mass)
      const xcm = (g, k) => { let s = 0.0, m = 0.0; P.forEach((p, i) => { if (mem[g][i]) { s += p[k] * 1.0; m += 1.0; } }); return s / m; };
      const hx = xcm('ball', 0), hy = xcm('ballnn', 1);
      // the neighbour shell around (2,2,2)
      const c = [2 * a, 2 * a, 2 * a], R = a * Math.SQRT1_2, shell = [];
      P.forEach(p => { const d = [0, 1, 2].map(k => p[k] - c[k]), r = Math.sqrt(d[0] * d[0] + d[1] * d[1] + d[2] * d[2]); if (Math.abs(r - R) < 1e-9) shell.push({ r, d }); });
      const inn = shell.filter(s => s.r <= R).length, ulp = R * Number.EPSILON / 2;
      const offs = [...new Set(shell.map(s => Math.round((s.r - R) / ulp)))].sort((x, y) => x - y);
      $('.gv-asum').innerHTML = `<p>${ok} of ${GROUPS.length} counts identical; xcm(ball, x) = ${s17(+xl[1])} by LAMMPS, ${s17(hx)} by hand ${hx === +xl[1] ? OK : BAD}; xcm(ballnn, y) ${hy === +xl[2] ? OK : BAD}.</p>
<p><b>Shared face.</b> left + right = ${lmp.left + lmp.right}, not 500: the 50 atoms on the plane x = 2.5 lattice are in both (<code>both</code>). Inclusive faces are the rule for block, so adjacent blocks always overlap on their common face.</p>
<p><b>A radius on a neighbour shell.</b> ball1 (radius 1.0 around a lattice site) has 1 + 12 + 6 = ${lmp.ball1}: the six second neighbours at exactly one lattice constant are all in. ballnn uses the nearest-neighbour distance √½ itself: the 12 neighbours are at r = R in exact arithmetic, but computed r − R comes out as ${offs.map(o => `${o > 0 ? '+' : ''}${o}`).join(', ')} half-ulps (1 half-ulp = ${ulp.toExponential(2)}σ), and only ${inn} of 12 pass <code>r &lt;= radius</code>. The group has ${lmp.ballnn} atoms where a reader expects 13 (or 1). Never put a region boundary exactly on a lattice distance; add or subtract a little (0.01σ).</p>`;
    };

    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.gv-stb'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps(`${ESCRIPT}\n${VSCRIPT}`); b.disabled = false;
      if (res.error) return err(st, res.error); st.textContent = '';
      const get = k => res.lines.find(l => l.startsWith(k + ' ')).split(/\s+/).slice(1);
      const E = get('E'); let ok = 0;
      const rows = EXPR.map((x, k) => { const h = evaluate(x, {}); ok += h === Number(E[k]); return row([`<code>${esc(x)}</code>`, E[k], h, h === Number(E[k]) ? OK : BAD]); });
      $('.gv-b').innerHTML = `<table>${head(['expression', 'LAMMPS', 'hand parser', ''])}${rows.join('')}</table>`;
      $('.gv-bsum').innerHTML = `<p>${ok} of ${EXPR.length} identical. The surprises all come from two rules. Unary minus has the highest precedence, above ^, so −2^2 = (−2)² = 4 (most languages and all of mathematics give −4). And ^ is popped by an incoming ^ like any other operator, so 2^3^2 = (2³)² = 64, not 2⁹ = 512. The rest is C: % is fmod (sign of the left operand), comparisons give 1 or 0 and chain left to right, so 3 &gt; 2 &gt; 1 is (1) &gt; 1 = 0.</p>`;
      const A = get('A'), B = get('B'), C = get('C'), T = res.thermo.rows[0];
      const va = 0.1 + 0.2, vb = 1 / 3, vars = { a: va, b: vb }, vc = Number(Number(vb.toPrecision(15)).toString()); vars.d = evaluate('v_b*3', vars);
      const p20 = x => x.toPrecision(20);
      const hb = [vb.toPrecision(15).replace(/0+$/, ''), evaluate('v_b*3==1', vars), (1 / 3) * 3 === 1 ? 1 : 0];
      const lines = [
        ['<code>${a}</code> (text, %.15g)', A[0], Number(va.toPrecision(15)).toString()],
        ['<code>$(v_a)</code> (immediate, %.20g)', A[1], p20(evaluate('v_a', vars))],
        ['<code>$(v_a==0.3)</code>', A[2], evaluate('v_a==0.3', vars)],
        ['<code>$(0.1+0.2==0.3)</code>', A[3], evaluate('0.1+0.2==0.3', vars)],
        ['<code>${b}</code>', B[0], hb[0]],
        ['<code>$(v_b)</code>', B[1], p20(evaluate('v_b', vars))],
        ['<code>$(v_b*3==1)</code>', B[2], hb[1]],
        ['<code>$((1/3)*3==1)</code>', B[3], hb[2]],
        ['<code>$(v_d)</code>, d = v_b*3', B[4], p20(evaluate('v_d', vars))],
        ['<code>$(v_c)</code>, c = ${b}', C[1], p20(vc)],
        ['thermo <code>v_a</code>', T[1], p20(va)],
        ['thermo <code>v_b</code>', T[2], p20(vb)],
        ['thermo <code>v_d</code>', T[3], p20(vars.d)],
      ];
      let ok2 = 0;
      $('.gv-c').innerHTML = `<table>${head(['', 'LAMMPS', 'hand', ''])}${lines.map(([l, x, h]) => { const same = Number(x) === Number(h); ok2 += same; return row([l, `<span class="mono">${esc(x)}</span>`, `<span class="mono">${esc(String(h))}</span>`, same ? OK : BAD]); }).join('')}</table>`;
      $('.gv-csum').innerHTML = `<p>${ok2} of ${lines.length} identical. 0.1 + 0.2 is 0.30000000000000004 as a double, and only the thermo column shows it. Everywhere else v_a is a reference <i>inside a formula</i> (even <code>$(v_a)</code> is the formula "v_a"), fetched through retrieve(), printed with %.15g ("0.3") and parsed back. So <code>v_a==0.3</code> is <b>true</b> while <code>0.1+0.2==0.3</code> is false, and <code>v_b*3==1</code> is false although (1/3)·3 is exactly 1 in doubles. v_d is rounded twice: thermo evaluates d's own formula with b already rounded (0.99999999999999889), and $(v_d) rounds that again (0.999999999999999).</p><p>None of this matters for physics (1 part in 10¹⁵), but it decides <code>if</code> tests and loop conditions that compare variables for equality, and it is why a hand check of a nested variable must round the same way. An internal-style variable (<code>variable x internal 0.1</code>) is the exception: it is read as a double (variable.cpp:2160).</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
