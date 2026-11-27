// Research · F.11: tools and contributing. What tools/ offers, how a contribution reaches LAMMPS, and the checks CI
// runs on it — with one of those checks caught missing every tab in the tree because of a one-character regex typo.
import { NATIVE } from './tools-native.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const show = s => esc(s).replace(/\t/g, '<span style="color:var(--accent,#c60)">→</span>\t').replace(/ +$/gm, m => `<span style="color:var(--accent,#c60)">${'·'.repeat(m.length)}</span>`);

const TOOLS = [
  ['tools/python/', 'log2txt.py, logplot.py, dump2xyz/pdb/cfg.py, dumpsort.py, neb_combine/neb_final.py; the Pizza.py log and dump readers they use'],
  ['tools/regression-tests/', 'run_tests.py: rerun examples/ and compare logs (page F.3)'],
  ['tools/coding_standard/', 'whitespace.py, permissions.py, homepage.py, errordocs.py, versiontags.py, fmtlib.py: the style checks CI runs'],
  ['tools/msi2lmp, ch2lmp, moltemplate (external), …', 'converters from other programs\' force-field and structure files to data files'],
  ['tools/replica/', 'reorder_remd_traj.py: sort temper trajectories by temperature (page F.8)'],
  ['tools/lammps-gui/', 'a graphical editor/runner for inputs (also distributed as a binary)'],
  ['unittest/', 'GoogleTest suites and the force-style YAML tests (page E.11)'],
];
const FLOW = `# the path of a contribution
git clone https://github.com/<you>/lammps.git     # your fork
git checkout -b pair-lj-mine develop               # branch from develop, not stable
#   add src/EXTRA-PAIR/pair_lj_mine.{cpp,h}, doc/src/pair_lj_mine.rst,
#   unittest/force-styles/tests/mol-pair-lj_mine.yaml, an example; update the style lists in the docs
cd src && make check-whitespace check-permissions check-homepage check-errordocs check-fmtlib
clang-format -i EXTRA-PAIR/pair_lj_mine.cpp        # src/.clang-format holds the LAMMPS style
git commit -am "add pair style lj/mine" && git push origin pair-lj-mine
#   open a pull request against lammps/lammps:develop; fill in the template's checklist
#   CI: style check, builds (GNU make + CMake, several compilers/platforms), unit tests, docs; then review`;
const CHECKLIST = ['the feature is complete', 'licensing information is complete', 'corresponding author information is complete', 'the source code follows the LAMMPS formatting guidelines', 'suitable new or updated documentation is included', 'the documentation is integrated and builds with the documentation system', 'verified with the conventional (make) build', 'verified with the CMake build', 'suitable tests are added to unittest/', 'a package README is included or updated', 'one or more example inputs are included'];
const ISSUE = `Title: coding_standard/whitespace.py misses TAB characters on lines with an "n" after the tab

tools/coding_standard/whitespace.py, check_tabs() (line 81):

    pattern = re.compile(r'[^\\n]*\\t+[^n]*\\n$')

The second character class is [^n] ("not the letter n") instead of [^\\n]. A line with a tab
is flagged only if no letter "n" follows the last tab: "\\tx = 1;" is reported, "\\tint n = 0;"
is not. With [^\\n]* the check reports 48 lines with TABs in 11 files of stable_22Jul2025_update6
(also present in develop as of 3002014), none of which "make check-whitespace" finds today.
Because fix_file() only fixes reported lines, "make fix-whitespace" leaves them in place too.

Fix: [^n]* -> [^\\n]* (one character). The TABs found could be removed in the same PR, or
excluded where they are intentional (unittest/commands/test_variables.cpp tests TAB parsing?).`;
const CLANG = `cd src
clang-format --dry-run pair_lj_cut.cpp                    # how far is a file from the LAMMPS style?
clang-format pair_lj_cut.cpp | diff pair_lj_cut.cpp - | grep -c '^[<>]'
#   pair_lj_cut.cpp      3 changed lines   (clang-format 18 vs the version the style was tuned for)
#   fix_nve.cpp          0
#   pair_lj_mine.cpp   157   (page E.3's Mie style, written by hand)`;

export default {
  id: 'tools-contributing',
  title: 'Tools and contributing',
  minutes: 12,
  summary: 'The last page of the LAMMPS part: what the tools/ directory gives you, and how code you write becomes part of LAMMPS. You use log2txt.py and check its output to the last bit. Then you run the style checks that every pull request must pass, on an eight-line test file. One of them, the TAB check, misses lines it should catch. You find the cause, a one-character typo in a regular expression, and count what it hid in the source tree: 48 lines in 11 files. The page ends with the issue you would file.',
  learn: [
    'What is in tools/, and which tool for which job',
    'Fork, branch, checks, pull request: the contribution path and its checklist',
    'Running and reading the CI style checks yourself',
    'Testing the testers: a regex bug in check-whitespace, and a draft report',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: tools/ and the way in</h2>
    <div class="tbl tc-tools"></div>
    <pre class="lmp-mini">${esc(FLOW)}</pre>
    <div class="prose"><p>The style check that CI runs on every pull request is ${gh('.github/workflows/style-check.yml', 31)}: from src/, <code>make check-whitespace check-permissions check-homepage check-errordocs check-fmtlib</code>, which call the scripts in tools/coding_standard/. The pull-request template then asks you to tick eleven boxes:</p></div>
    <div class="prose tc-check"></div>
  </section>

  <section class="block">
    <h2>2 · Use a tool: log2txt.py</h2>
    <div class="prose"><p>Page F.10's 5000-step log, reduced to three columns:</p></div>
    <pre class="lmp-mini">LAMMPS_PYTHON_TOOLS=tools/python/pizza python3 tools/python/log2txt.py log.lammps data.txt Step Temp PotEng -n</pre>
    <pre class="lmp-mini">${esc(NATIVE.data.split('\n').slice(0, 4).join('\n'))}\n…</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Check data.txt against the log</button></div>
    <div class="prose tc-asum"></div>
  </section>

  <section class="block">
    <h2>3 · Run the CI style check yourself</h2>
    <div class="prose"><p>An eight-line test file (shown with → for TAB and · for trailing blanks), checked with <code>python3 tools/coding_standard/whitespace.py .</code>:</p></div>
    <pre class="lmp-mini">${NATIVE.demo.split('\n').map((l, i) => `${String(i + 1).padStart(2)}  ${show(l)}`).join('\n')}</pre>
    <pre class="lmp-mini">${esc(NATIVE.orig)}</pre>
    <div class="prose"><p>Lines 3 and 5 are reported, lines 4 and 6 are not, though all four start with or contain a TAB. Exit code 1 here; but on the LAMMPS tree itself the check passes. Why? ${gh('tools/coding_standard/whitespace.py', 81)}:</p></div>
    <pre class="lmp-mini">pattern = re.compile(r'[^\\n]*\\t+[^n]*\\n$')      # check_tabs()</pre>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Test the pattern, and the fixed one</button></div>
    <div class="prose tc-bsum"></div>
    <div class="tbl tc-btab"></div>
  </section>

  <section class="block">
    <h2>4 · Formatting with clang-format</h2>
    <pre class="lmp-mini">${esc(CLANG)}</pre>
    <div class="prose"><p>Formatting is not enforced by CI the way whitespace is, but reviewers will ask for it, and it is easier to run clang-format once than to discuss braces. Do it in a separate commit, so the reviewer can see your real changes apart from the reformatting.</p></div>
  </section>

  <section class="block">
    <h2>5 · Report it</h2>
    <div class="prose"><p>A bug in a tool is reported like a bug in the code: what, where, a minimal example, the effect, the fix. The issue you would file (the second small finding, that tools/python/pizza/log.py writes regexes like "\\s" in ordinary strings and so prints SyntaxWarnings under Python 3.12+, would go in its own issue):</p></div>
    <pre class="lmp-mini">${esc(ISSUE)}</pre>
    <div class="prose"><p>With this page the LAMMPS part ends. You can read the source (part C), change it (E), extend it as a plugin or through Python (E, F), test, debug, profile and reproduce it (F), and send what you find back upstream.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> If the check was broken, why does the tree have only 48 TAB lines and not thousands?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Because the check half works. Every TAB line with no "n" after the last TAB was reported, and contributors fixed it (or <code>make fix-whitespace</code> did). What survived is exactly the TABs the regex cannot see: all 48 remaining lines have an "n" after the TAB. A test that catches some cases leaves behind a population selected to be the ones it misses. That is a good reason to test your tests on examples where you know the answer, as section 3 did.</div></div>
    <div class="ex"><p><b>2.</b> log2txt.py printed −6.77336805325309, the log −6.7733680532530904. Did the tool lose precision?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. Python prints the shortest decimal that reads back as the same double; %.17g always prints 17 digits, some of them redundant. Both strings are the same double, as section 2 checked for all 153 values. It would lose precision if the log had been written with the default %g (6 digits): no tool can restore digits that were never printed (page F.7).</div></div>
    <div class="ex"><p><b>3.</b> Your new pair style passes all checks and unit tests. What will a reviewer still ask?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The things no script checks. Is it physically sound and documented with a reference? Does it duplicate an existing style (lj/cut with Mie exponents exists as mie/cut)? Does it support everything the docs claim (single(), restart, mixing, r-RESPA)? Is the example small and fast? And is there someone, you, who will answer issues about it later?</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('.tc-tools').innerHTML = `<table>${head(['where', 'what'])}${TOOLS.map(row).join('')}</table>`;
    $('.tc-check').innerHTML = `<ol>${CHECKLIST.map(c => `<li>${c}</li>`).join('')}</ol>`;

    $('[data-k=a]').onclick = () => {
      const L = NATIVE.log, hdr = L[0].trim().split(/\s+/), cols = ['Step', 'Temp', 'PotEng'].map(c => hdr.indexOf(c));
      const logRows = L.slice(1).map(l => l.trim().split(/\s+/)), txt = NATIVE.data.trim().split('\n').filter(l => !l.startsWith('#')).map(l => l.trim().split(/\s+/));
      let same = 0, n = 0, shorter = 0;
      logRows.forEach((r, i) => cols.forEach((c, k) => { n++; if (Number(r[c]) === Number(txt[i][k])) same++; if (txt[i][k].replace(/[-.]/g, '').length < r[c].replace(/[-.]/g, '').length) shorter++; }));
      $('.tc-asum').innerHTML = `<p>${logRows.length} rows × 3 columns: ${same} of ${n} values the same double as in the log ${same === n ? OK : BAD}, although ${shorter} of them are printed with fewer digits. Python writes the shortest decimal that rounds back to the same double (−6.77336805325309), while the log's %.17g writes all 17 digits (−6.7733680532530904). JavaScript's Number() reads both to the same bits. Step comes out as 0.0, 100.0: log.py stores every column as a float.</p>`;
    };

    $('[data-k=b]').onclick = () => {
      const bad = /^[^\n]*\t+[^n]*\n$/, good = /^[^\n]*\t+[^\n]*\n$/;
      const D = NATIVE.demo.split('\n').slice(0, -1).map((l, i) => [i + 1, l, bad.test(l + '\n'), good.test(l + '\n'), /\t/.test(l)]);
      const T = NATIVE.tree; let hb = 0, hg = 0, withN = 0;
      for (const [, , l] of T) { if (bad.test(l + '\n')) hb++; if (good.test(l + '\n')) hg++; if (/n/.test(l.slice(l.lastIndexOf('\t') + 1))) withN++; }
      const files = new Set(T.map(t => t[0])).size;
      $('.tc-btab').innerHTML = `<table>${head(['line', 'text', 'has a TAB', "LAMMPS's pattern", 'fixed pattern'])}${D.map(([i, l, b, g, t]) => row([i, `<code>${show(l)}</code>`, t ? 'yes' : '', b ? 'reported' : '', g ? 'reported' : ''])).join('')}</table>`;
      $('.tc-bsum').innerHTML = `<p>The second character class is <code>[^n]</code>, "anything but the letter n", where <code>[^\\n]</code>, "anything but a newline", was meant. The regex matches a TAB line only if the rest of the line after the last TAB contains no letter n. Line 4 (<code>int n</code>) and line 6 (<code>nlocal</code>) have one, so they pass. The pattern applied here in JavaScript reproduces the native report exactly ${D.filter(d => d[2]).map(d => d[0]).join(',') === '3,5' ? OK : BAD}, and the fixed pattern reports all four TAB lines ${D.filter(d => d[3]).map(d => d[0]).join(',') === '3,4,5,6' ? OK : BAD}.</p>
      <p>On the whole stable_22Jul2025_update6 tree, with the script's own include/exclude lists: the original reports 0 TAB lines, the one-character fix reports ${T.length} in ${files} files (25 in lib/colvars/colvarmodule_refs.h, the rest in KOKKOS, LATBOLTZ, READMEs, a CMake module and a unit test). Applied here to those ${T.length} lines: LAMMPS's pattern matches ${hb} ${hb === 0 ? OK : BAD}, the fixed one ${hg} ${hg === T.length ? OK : BAD}, and all ${withN} have an "n" after their last TAB. The pattern is unchanged in the develop branch as of 24 November 2026.</p>`;
    };

    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
