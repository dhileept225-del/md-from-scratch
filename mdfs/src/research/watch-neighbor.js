// Research · E.8: watch the neighbour list decide. Prints added to Neighbor::decide and check_distance; every printed
// largest displacement recomputed from the dumped positions (bit for bit); the blind spot of the "dangerous builds"
// counter; and the triclinic corner rule (neighbor.cpp:2414) caught in the act: delta2 is 0 at 999 of 1000 checks.
const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';
const s17 = x => (x === 0 ? '0' : x.toPrecision(17).replace(/\.?0+(e|$)/, '$1'));

const LIQ = nm => `units           lj
atom_style      atomic
lattice         fcc 0.8442
region          box block 0 3 0 3 0 3
create_box      1 box
create_atoms    1 box
mass            1 1.0
velocity        all create 1.44 87287 loop geom
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neigh_modify    ${nm}
fix             1 all nve
dump            d all custom 1 x.dump id x y z ix iy iz
dump_modify     d format float %.17g sort id
run             60`;
const BOX = `units           lj
atom_style      atomic
region          box prism 0 10 0 10 0 10 0 0 0 units box
create_box      1 box
create_atoms    1 single 9 9 9 units box
mass            1 1.0
pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5
neighbor        1.0 bin
neigh_modify    every 1 delay 0 check yes
fix             1 all nve
fix             2 all deform 1 x final 0 10.5 y final 0 11 z final 0 12 remap x
run             1000`;

const b64 = s => Uint8Array.from(atob(s), c => c.charCodeAt(0)).buffer;
const parseWatch = lines => lines.filter(l => /max_dr2/.test(l)).map(l => { const m = l.match(/step (\d+) ago (\d+) max_dr2 (\S+) trigger (\S+) (\w+)( DANGEROUS)?/); return { s: +m[1], ago: +m[2], m: +m[3], tr: +m[4], build: m[5] === 'build', danger: !!m[6] }; });

// largest |x − xhold|² by hand. At a build step the dump shows positions after the periodic wrap, but decide() ran
// before it; atoms whose image flag changed are put back with ± L (exact when the wrap was exact).
function handMax(N, L, X, I, frameOf, W) {
  const at = (f, i, k) => X[(f * N + i) * 3 + k], im = (f, i, k) => I[(f * N + i) * 3 + k];
  let hold = frameOf(0); const out = [];
  for (const w of W) {
    const f = frameOf(w.s); let mx = 0;
    for (let i = 0; i < N; i++) {
      // same association as LAMMPS: delx*delx + dely*dely + delz*delz
      const dx = [0, 1, 2].map(k => { let x = at(f, i, k); const di = im(f, i, k) - im(hold, i, k); if (di === 1) x = x + L; else if (di === -1) x = x - L; return x - at(hold, i, k); });
      const r = dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2];
      if (r > mx) mx = r;
    }
    out.push(mx); if (w.build) hold = f;
  }
  return out;
}

export default {
  id: 'watch-neighbor',
  title: 'Watch the neighbour list decide',
  minutes: 18,
  summary: 'Three print lines inside Neighbor::decide and check_distance make every rebuild decision visible: the step, steps since the last build, the largest displacement and the trigger. You recompute each printed displacement from the dumped positions, bit for bit, and predict every build; with checks only every 10 steps you watch atoms move 60 % past half the skin without a single "dangerous build" being counted; and in a growing triclinic box you catch the corner rule from page D.3 in the act: its second-largest corner move is 0 at 999 of 1000 checks.',
  learn: [
    'Instrumenting a decision: print what the code compares, not what it concludes',
    'Reproducing max |x − xhold|² across periodic wraps',
    'What "Dangerous builds" counts, and what it cannot see',
    'Evidence for a bug report, straight from the code',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: one comparison per step</h2>
    <div class="prose"><p>Every step Verlet::run asks ${gh('src/neighbor.cpp', 2363)} whether to rebuild (page D.3). With <code>neigh_modify every N delay D check yes</code>: count steps since the last build (ago); if ago ≥ D and ago is a multiple of N, call check_distance, which compares each atom's displacement since the last build with half the skin, stopping at the first atom over the limit. "Dangerous" is a build triggered at the first allowed check (ago = max(N, D)): the limit may have been crossed before anyone looked.</p>
<p>The loop stops early, so to see the numbers the change adds a separate loop for the largest displacement, prints it with the trigger and the decision, and (for changing triclinic boxes) prints the two corner moves the code picked next to the two largest.</p></div>
  </section>

  <section class="block">
    <h2>2 · Change</h2>
    <pre class="lmp-mini gd wn-diff">loading…</pre>
    <div class="prose"><p>Rebuild 11 s (neighbor.cpp and the version file). Afterwards <code>git checkout src/neighbor.cpp</code> and a rebuild restore the original; the build counts below are the same with and without the prints.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: a liquid checked every step</h2>
    <div class="prose"><p>108 LJ atoms at T = 1.44, default skin 0.3σ, so the trigger is (0.3/2)² = 0.0225σ². Positions and image flags dumped every step with 17 digits (native run).</p></div>
    <pre class="lmp-mini">${esc(LIQ('every 1 delay 0 check yes'))}</pre>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Recompute every check</button><span class="lmp-status muted wn-sta"></span></div>
    <div class="tbl wn-a"></div>
    <div class="prose wn-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: checking every 10 steps</h2>
    <div class="prose"><p>Same run with <code>neigh_modify every 10 delay 0 check yes</code>.</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ Show the checks</button></div>
    <div class="tbl wn-b"></div>
    <div class="prose wn-bsum"></div>
  </section>

  <section class="block">
    <h2>5 · Observe: the corner rule, caught</h2>
    <div class="prose"><p>The triclinic branch of check_distance (${gh('src/neighbor.cpp', 2414)}) should subtract the two largest corner moves from the skin; page D.3 showed that it keeps a maximum but does not pass the old maximum down to second place. The prints show it directly. One atom, a prism box growing from 10³ to 10.5 × 11 × 12 with the lower corner fixed, skin 1.0σ:</p></div>
    <pre class="lmp-mini">${esc(BOX)}</pre>
    <div class="row"><button class="btn sm primary" data-k="c">▶ Show the corner moves</button></div>
    <div class="tbl wn-c"></div>
    <div class="prose wn-csum"></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why does the instrumented code compute the maximum in a separate loop instead of printing inside the existing one?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The existing loop breaks at the first atom over the limit, so its last rsq is not the maximum, and on "keep" steps it never sees all atoms at all. A print must not change the decision it reports: a separate read-only loop leaves the original logic untouched.</div></div>
    <div class="ex"><p><b>2.</b> With every 10 and a build at step 30 with max_dr2 = 0.047, why is it not counted as dangerous, and could pairs have been missed?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The counter only fires when ago = max(every, delay), i.e. at the first check after a build; at step 30, ago = 20. At step 20 everything was within 0.126σ; by step 30 one atom had moved 0.217σ. Two atoms approaching each other can close a gap by up to twice that between checks, more than the 0.3σ skin, so yes, a pair could have come within the cutoff unlisted. Use every 1 (the default) unless you have measured that it is safe.</div></div>
    <div class="ex"><p><b>3.</b> What is the one-line fix for the corner rule?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden><code>if (delta &gt; delta1) { delta2 = delta1; delta1 = delta; } else if (delta &gt; delta2) delta2 = delta;</code> — exactly what the WATCH line computes as "two largest".</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    import('./watch-native.js').then(({ NATIVE }) => this.wire(root, $, NATIVE));
    return () => {};
  },

  wire(root, $, NATIVE) {
    $('.wn-diff').textContent = NATIVE.diff.replace(/\n(?=diff --git)/g, '\n');
    const N = NATIVE.n, L = NATIVE.L;
    const A1 = parseWatch(NATIVE.a1w), A2 = parseWatch(NATIVE.a2w);
    $('[data-k=a]').onclick = () => {
      const X = new Float64Array(b64(NATIVE.a1x)), I = new Int8Array(b64(NATIVE.a1i)), H = handMax(N, L, X, I, s => s, A1);
      let ok = 0, dec = 0; const trig = 0.25 * 0.3 * 0.3, rows = [];
      A1.forEach((w, k) => { const same = H[k] === w.m, d = (H[k] > trig) === w.build && trig === w.tr; ok += same; dec += d; if (w.build || k < 3 || A1[k + 1]?.build) rows.push(row([w.s, w.ago, `<span class="mono">${s17(w.m)}</span>`, w.build ? '<b>build</b>' : 'keep', same && d ? OK : BAD])); });
      $('.wn-a').innerHTML = `<table>${head(['step', 'ago', 'max_dr2 (LAMMPS)', 'decision', 'hand'])}${rows.join('')}</table>`;
      const builds = A1.filter(w => w.build);
      $('.wn-asum').innerHTML = `<p>All ${A1.length} checks: max_dr2 identical to the last bit at ${ok}, decision and trigger (${s17(trig)}) right at ${dec}. The table shows the builds and the step before each. Builds at steps ${builds.map(w => w.s).join(', ')}: every ${builds.slice(1).map((w, k) => w.s - builds[k].s).join(', ')} steps, each time the fastest atom has just passed 0.15σ (√0.0225). The hand code had one subtlety: decide() runs before the periodic wrap of a build step, but the dump is written after it, so atoms whose image flag changed are moved back by ±L.</p>`;
    };
    $('[data-k=b]').onclick = () => {
      const X = new Float64Array(b64(NATIVE.a2x)), I = new Int8Array(b64(NATIVE.a2i)), H = handMax(N, L, X, I, s => s / 10, A2);
      const rows = A2.map((w, k) => row([w.s, w.ago, `<span class="mono">${s17(w.m)}</span>`, Math.sqrt(w.m).toFixed(3), w.build ? (w.danger ? '<b>build, DANGEROUS</b>' : '<b>build</b>') : 'keep', H[k] === w.m ? OK : BAD]));
      $('.wn-b').innerHTML = `<table>${head(['step', 'ago', 'max_dr2', 'max |dr| (σ)', 'decision', 'hand'])}${rows.join('')}</table>`;
      const nd = A2.filter(w => w.danger).length, hidden = A2.filter(w => w.build && !w.danger);
      $('.wn-bsum').innerHTML = `<p>${A2.length} checks (the other 54 steps print "not checked"). LAMMPS reports ${nd} dangerous builds, at the steps where the first check after a build already found an atom past 0.15σ. The builds at steps ${hidden.map(w => w.s).join(' and ')} are not counted, although the fastest atom had moved ${hidden.map(w => Math.sqrt(w.m).toFixed(3)).join(' and ')}σ, up to ${Math.round(100 * (Math.max(...hidden.map(w => Math.sqrt(w.m))) / 0.15 - 1))} % past half the skin: the check 10 steps earlier said "keep", so the counter assumes the crossing happened late. "Dangerous builds = 0" with every &gt; 1 is a lower bound, not a guarantee.</p>`;
    };
    $('[data-k=c]').onclick = () => {
      const rows = NATIVE.dw.filter((r, k) => k < 4 || r[4] === 'build').slice(0, 9).map(r => row([r[0], r[1], Number(r[5]).toFixed(5), Number(r[6]).toFixed(5), Number(r[7]).toFixed(5), Number(r[8]).toFixed(5), Number(r[3]).toFixed(5), (Math.max(0, 0.5 * (1.0 - (Number(r[7]) + Number(r[8])))) ** 2).toFixed(5), r[4]]));
      $('.wn-c').innerHTML = `<table>${head(['step', 'ago', 'delta1 (code)', 'delta2 (code)', 'largest', '2nd largest', 'trigger (code)', 'trigger (2 largest)', 'decision'])}${rows.join('')}</table>`;
      const fb = NATIVE.dw.find(r => r[4] === 'build');
      $('.wn-csum').innerHTML = `<p>At ${NATIVE.dzero} of ${NATIVE.dn} checks the code's delta2 is 0 while the true second-largest corner move is almost as large as the largest: the corners are visited (0,0,0), (1,0,0), (0,1,0), (1,1,0), (0,0,1), …, (1,1,1), their moves grow in that order here, and each new maximum overwrites delta1 without moving the old one down. The trigger comes out too large, so LAMMPS rebuilt 6 times (steps 157, 313, …, 937) where the documented rule asks for 8 (page D.3's hand model: every 116 steps). At the first build the atom had moved ${Math.sqrt(Number(fb[2])).toFixed(3)}σ, against ${Math.sqrt(Number(fb[3])).toFixed(3)}σ allowed by the code and ${Math.max(0, 0.5 * (1.0 - (Number(fb[7]) + Number(fb[8])))).toFixed(3)}σ by the two-largest rule. A missed rebuild only matters if a pair comes within the cutoff unlisted, which needs fast atoms and a fast box change at the same time; but it is a bug, with a one-line fix (exercise 3), and these prints are the evidence a report needs.</p>`;
    };
  },
};
