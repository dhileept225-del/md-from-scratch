// Research · Part 0.6: git, and reading a diff, with a real one-character change to LAMMPS compiled and run.
import { g8 } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const RHO = 0.8442, X = RHO * RHO / 2, S3 = 1991 / 144, S6 = 3018233 / 248832, F = 3999 / 4000;

// melt step 0 for coefficients lj1..lj4 (epsilon = sigma = 1): E = ½(c x² S6 − d x S3), P_vir = (ρ/6)(a x² S6 − b x S3)
const step0 = (a, b, c, d) => { const E = 0.5 * (c * X * X * S6 - d * X * S3); return { E, TotEng: E + 1.5 * 3 * F, Press: RHO * 3 * F + RHO / 6 * (a * X * X * S6 - b * X * S3) }; };
const LINES = n => [
  '',
  `  lj1[i][j] = ${n[0]} * epsilon[i][j] * pow(sigma[i][j], 12.0);`,
  `  lj2[i][j] = ${n[1]} * epsilon[i][j] * pow(sigma[i][j], 6.0);`,
  `  lj3[i][j] = ${n[2]} * epsilon[i][j] * pow(sigma[i][j], 12.0);`,
  `  lj4[i][j] = ${n[3]} * epsilon[i][j] * pow(sigma[i][j], 6.0);`,
  '',
  '  if (offset_flag && (cut[i][j] > 0.0)) {',
];
const ORIG = ['48.0', '24.0', '4.0', '4.0'];
const esc = s => s.replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

export default {
  id: 'git-diff',
  title: 'Git: seeing exactly what changed',
  minutes: 20,
  summary: 'Git keeps every version of LAMMPS and shows any change as a diff: the lines removed in red, the lines added in green. You learn to read one, then follow a real experiment: one character of the LAMMPS source changed, recompiled, run, matched by hand, and undone.',
  learn: [
    'The git commands that download a version, show changes and undo them',
    'How to read a diff, line by line',
    'A real one-character change to LAMMPS: predict its effect, compare with the run, undo it',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: git remembers every version</h2>
    <p>LAMMPS is developed with <b>git</b>, a program that stores every version of every file. Two uses matter for this course: you download an exact version (so your line numbers match ours), and after you edit the code, git shows <b>exactly</b> what you changed and can undo it. Every code change in Part E will be shown this way.</p>
    <div class="tbl wide"><table>${head(['Command', 'What it does'])}
${row(['<code>git clone -b stable_22Jul2025_update6 --depth 1 https://github.com/lammps/lammps.git</code>', 'download one exact version (page 0.4)'])}
${row(['<code>git status</code>', 'list the files you have changed'])}
${row(['<code>git diff</code>', 'show every changed line'])}
${row(['<code>git diff --stat</code>', 'just count the changes per file'])}
${row(['<code>git restore src/pair_lj_cut.cpp</code>', 'undo your changes to one file'])}
${row(['<code>git stash</code> / <code>git stash pop</code>', 'put all changes aside / bring them back'])}</table></div>
  </section>

  <section class="block prose">
    <h2>2 · A real change, read as a diff</h2>
    <p>We changed one character in <code>src/pair_lj_cut.cpp</code>: the 4.0 in line 511 became 2.0. This is what git reported, copied from our terminal:</p>
    <div class="calc">$ git status
Not currently on any branch.
Changes not staged for commit:
	modified:   src/pair_lj_cut.cpp

$ git diff --stat
 src/pair_lj_cut.cpp | 2 +-
 1 file changed, 1 insertion(+), 1 deletion(-)</div>
    <div class="diff"><span class="dh">diff --git a/src/pair_lj_cut.cpp b/src/pair_lj_cut.cpp</span>
<span class="dh">index 55279f6..7e4a443 100644</span>
<span class="dh">--- a/src/pair_lj_cut.cpp</span>
<span class="dh">+++ b/src/pair_lj_cut.cpp</span>
<span class="dat">@@ -508,7 +508,7 @@ double PairLJCut::init_one(int i, int j)</span>

   lj1[i][j] = 48.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);
   lj2[i][j] = 24.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);
<span class="del">-  lj3[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);</span>
<span class="add">+  lj3[i][j] = 2.0 * epsilon[i][j] * pow(sigma[i][j], 12.0);</span>
   lj4[i][j] = 4.0 * epsilon[i][j] * pow(sigma[i][j], 6.0);

   if (offset_flag && (cut[i][j] > 0.0)) {</div>
    <div class="tbl wide"><table>${head(['Line of the diff', 'Meaning'])}
${row(['<code>--- a/…</code>, <code>+++ b/…</code>', 'a is the file before, b the file after'])}
${row(['<code>@@ -508,7 +508,7 @@</code>', 'this piece starts at line 508 and shows 7 lines, before and after'])}
${row(['<code>double PairLJCut::init_one(…)</code>', 'the function the change is in (a helpful extra)'])}
${row(['lines starting with a space', 'unchanged, shown for context'])}
${row(['<span class="del">− red</span>', 'removed'])}
${row(['<span class="add">+ green</span>', 'added: a changed line appears as one removal and one addition'])}
${row(['"Not currently on any branch"', 'normal: you downloaded a release tag, not a line of development'])}</table></div>
    <p>Notice the "Not currently on any branch" line: it is not an error, only a reminder that you are on a fixed release.</p>
  </section>

  <section class="block prose">
    <h2>3 · Predict by hand, then observe the real run</h2>
    <p>lj3 is used only in the <b>energy</b> line (page 0.1, line 129). The force lines use lj1 and lj2, which did not change. So before running anything we can predict:</p>
    <div class="calc">energy per atom:  ½ Σ n (lj3 r⁻¹² − lj4 r⁻⁶) = ½(2 x² S6 − 4 x S3) = x² S6 − 2x S3
                = 1.54016742 − 2 × 4.92685145 = <b>−8.3135355</b>   (was −6.7733681)
TotEng          = −8.31353548 + 4.49887500    = <b>−3.8146605</b>
Press           = <b>−3.7033504</b>  unchanged: it comes from the forces
Temp            unchanged at <b>every</b> step: same forces, so the atoms move exactly the same</div>
    <p>We recompiled (16 s: two files, pair_lj_cut.cpp and lammps.cpp, which records the git state in the version information) and ran melt. Real output, original and changed:</p>
    <div class="tbl"><table>${head(['Step', 'Temp', 'E_pair original', 'E_pair changed', 'Press (both)'])}
${row(['0', '3', '−6.7733681', '<b>−8.3135355</b>', '−3.7033504'])}
${row(['50', '1.6842865', '−4.8082494', '−8.4398955', '5.5666131'])}
${row(['100', '1.6712577', '−4.7875609', '−8.4401865', '5.6613913'])}
${row(['150', '1.6444751', '−4.7471034', '−8.4454307', '5.8614211'])}
${row(['200', '1.6471542', '−4.7509053', '−8.4561267', '5.8805431'])}
${row(['250', '1.6645597', '−4.7774327', '−8.4536814', '5.7526089'])}</table></div>
    <p>Every prediction holds: step 0 matches to all digits, and Temp and Press are identical at every step. The change made the <b>reported</b> energy wrong without changing the motion at all. That is a warning for Part E: energy and force are written in separate lines, and a code change must keep them consistent (F = −dE/dr), or the program will run happily and print nonsense.</p>
    <div class="calc">$ git restore src/pair_lj_cut.cpp       <b>undo</b>
$ make -j 2                              <b>13 s</b>
$ ./lmp -in in.melt
         0   3             -6.7733681      0             -2.2744931     -3.7033504      <b>back to the original</b></div>
  </section>

  <section class="block">
    <h2>4 · Your own edit: see the diff and the prediction</h2>
    <p class="muted">Change the four numbers of lines 509–512. The page writes the diff git would show, and predicts the melt step-0 line by hand. To run your edit for real you would compile it (page 0.4); the run above shows that such predictions hold.</p>
    <div class="lmp" id="mydiff">
      <div class="row">${['lj1', 'lj2', 'lj3', 'lj4'].map((n, i) => `<label class="rlab">${n} <input class="r-in coef" value="${ORIG[i]}"></label>`).join('')}<button class="btn sm ghost" data-k="orig">Original</button></div>
      <div class="diff mine"></div>
      <div class="tbl pred"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Change only lj1 from 48.0 to 24.0. What happens to E_pair, Press and the motion?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>E_pair at step 0 is unchanged (−6.7733681): lj1 appears only in the force. The repulsive force halves, so Press = 2.53196685 + (0.8442/6)(24 x² S6 − 24 x S3) = <b>−8.9041878</b>, and the atoms move differently from step 1 on. Energy and force no longer belong together: energy is not conserved. Try it in the box above.</div></div>
    <div class="ex"><p><b>2.</b> A diff begins with <code>@@ -116,1 +116,1 @@</code>. What does it tell you?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The change is at line 116 and covers one line before and one line after: a single line was replaced, with no lines added or removed. Line 116 of pair_lj_cut.cpp is the forcelj line.</div></div>
    <div class="ex"><p><b>3.</b> You edited three files and want to start again from the release. Which one command undoes everything?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>git restore .</code> (the dot means "everything here"). If you might want the changes later, <code>git stash</code> puts them aside instead, and <code>git stash pop</code> brings them back. Then run make again.</div></div>
  </section>`,

  mount(root) {
    const box = root.querySelector('#mydiff'), ins = [...box.querySelectorAll('.coef')];
    function draw() {
      const n = ins.map(i => i.value.trim() || '0'), A = LINES(ORIG), B = LINES(n);
      let html = `<span class="dh">--- a/src/pair_lj_cut.cpp</span>\n<span class="dh">+++ b/src/pair_lj_cut.cpp</span>\n<span class="dat">@@ -508,7 +508,7 @@ double PairLJCut::init_one(int i, int j)</span>\n`;
      const changed = A.some((l, i) => l !== B[i]);
      A.forEach((l, i) => { html += l === B[i] ? ` ${esc(l)}\n` : `<span class="del">-${esc(l)}</span>\n<span class="add">+${esc(B[i])}</span>\n`; });
      box.querySelector('.mine').innerHTML = changed ? html : '<span class="muted">No change: git diff prints nothing.</span>';
      const v = n.map(Number);
      if (v.some(x => !isFinite(x))) { box.querySelector('.pred').innerHTML = '<p class="lmp-err">Use plain numbers, such as 2.0.</p>'; return; }
      const o = step0(48, 24, 4, 4), p = step0(...v);
      const cell = (a, b) => g8(a) === g8(b) ? g8(b) : `<b>${g8(b)}</b>`;
      box.querySelector('.pred').innerHTML = `<table><tr><th>melt step 0</th><th>original</th><th>your edit (by hand)</th></tr>
        <tr><td>Temp</td><td>3</td><td>3</td></tr>
        <tr><td>E_pair</td><td>${g8(o.E)}</td><td>${cell(o.E, p.E)}</td></tr>
        <tr><td>TotEng</td><td>${g8(o.TotEng)}</td><td>${cell(o.TotEng, p.TotEng)}</td></tr>
        <tr><td>Press</td><td>${g8(o.Press)}</td><td>${cell(o.Press, p.Press)}</td></tr></table>
        <p class="muted small">${v[0] / v[2] === 12 && v[1] / v[3] === 6 ? 'Energy and force still belong together (lj1 = 12·lj3, lj2 = 6·lj4): a consistent potential.' : '<span class="m-bad">Energy and force no longer belong together: lj1 must be 12 × lj3 and lj2 must be 6 × lj4 for F = −dE/dr.</span>'}</p>`;
    }
    ins.forEach(i => i.oninput = draw);
    box.querySelector('[data-k=orig]').onclick = () => { ins.forEach((i, k) => { i.value = ORIG[k]; }); draw(); };
    draw();
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
