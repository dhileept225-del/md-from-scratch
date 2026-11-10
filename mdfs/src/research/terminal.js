// A small practice terminal for the research pages. Files live in memory; `lmp` runs the real LAMMPS (WebAssembly).
import { runLammps } from './lammps-runner.js';

export const IN_MELT = `# 3d Lennard-Jones melt

units           lj
atom_style      atomic

lattice         fcc 0.8442
region          box block 0 10 0 10 0 10
create_box      1 box
create_atoms    1 box
mass            1 1.0

velocity        all create 3.0 87287 loop geom

pair_style      lj/cut 2.5
pair_coeff      1 1 1.0 1.0 2.5

neighbor        0.3 bin
neigh_modify    every 20 delay 0 check no

fix             1 all nve

#dump           id all atom 50 dump.melt

#dump           2 all image 25 image.*.jpg type type &
#               axes yes 0.8 0.02 view 60 -30
#dump_modify    2 pad 3

#dump           3 all movie 25 movie.mpg type type &
#               axes yes 0.8 0.02 view 60 -30
#dump_modify    3 pad 3

thermo          50
run             250
`;

const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');

function tokenize(line) {
  const out = []; let cur = '', q = null, has = false;
  for (const ch of line) {
    if (q) { if (ch === q) q = null; else cur += ch; continue; }
    if (ch === '"' || ch === "'") { q = ch; has = true; continue; }
    if (/\s/.test(ch)) { if (cur || has) { out.push(cur); cur = ''; has = false; } continue; }
    cur += ch;
  }
  if (cur || has) out.push(cur);
  return out;
}

// mountTerminal(el, { files, start, onRun }) ; files: { 'path/to/file': 'content' } relative to /home/you
export function mountTerminal(el, opts = {}) {
  const HOME = ['home', 'you'];
  const root = { type: 'dir', kids: {} };
  const mk = (parts) => { let d = root; for (const p of parts) { d.kids[p] = d.kids[p] || { type: 'dir', kids: {} }; d = d.kids[p]; } return d; };
  mk(HOME);
  for (const [p, c] of Object.entries(opts.files || {})) { const parts = [...HOME, ...p.split('/')]; const name = parts.pop(); mk(parts).kids[name] = { type: 'file', text: c }; }
  let cwd = [...HOME, ...((opts.start || '').split('/').filter(Boolean))];
  const hist = []; let hi = 0, busy = false;

  el.classList.add('term');
  el.innerHTML = `<div class="term-out" aria-live="polite"></div>
    <div class="term-line"><span class="term-prompt"></span><input class="term-in" autocomplete="off" autocapitalize="off" spellcheck="false" aria-label="terminal command"></div>
    <div class="term-editor" hidden><div class="term-ed-head"></div><textarea spellcheck="false" wrap="off"></textarea><div class="row"><button class="btn sm primary" data-ed="save">Save and exit (Ctrl-O, Ctrl-X)</button><button class="btn sm ghost" data-ed="quit">Exit without saving</button></div></div>`;
  const out = el.querySelector('.term-out'), input = el.querySelector('.term-in'), promptEl = el.querySelector('.term-prompt');
  const ed = el.querySelector('.term-editor'), edTa = ed.querySelector('textarea');

  const pathStr = parts => '/' + parts.join('/');
  const prompt = () => { const p = pathStr(cwd).replace(pathStr(HOME), '~'); promptEl.textContent = `you@mdfs:${p}$ `; };
  const print = (s, cls = '') => { const d = document.createElement('div'); if (cls) d.className = cls; d.innerHTML = s; out.appendChild(d); out.scrollTop = out.scrollHeight; };
  const resolve = p => {
    if (!p || p === '~') return [...HOME];
    let parts = p.startsWith('/') ? [] : p.startsWith('~/') ? [...HOME] : [...cwd];
    for (const s of p.replace(/^~\//, '').split('/')) { if (!s || s === '.') continue; if (s === '..') parts.pop(); else parts.push(s); }
    return parts;
  };
  const node = parts => { let d = root; for (const p of parts) { if (!d || d.type !== 'dir') return null; d = d.kids[p]; } return d || null; };
  const parentOf = parts => node(parts.slice(0, -1));
  const readFile = p => { const n = node(resolve(p)); if (!n) return { err: `${p}: No such file or directory` }; if (n.type === 'dir') return { err: `${p}: Is a directory` }; return { text: n.text }; };

  const cmds = {
    help: () => print(esc('Commands here: pwd  ls  cd  cat  less  grep  mkdir  cp  nano  lmp (or ./lmp)  mpirun  clear')),
    pwd: () => print(esc(pathStr(cwd))),
    ls: a => { const t = a.filter(x => !x.startsWith('-'))[0]; const n = node(resolve(t)); if (!n) return print(esc(`ls: cannot access '${t}': No such file or directory`), 'err');
      if (n.type === 'file') return print(esc(t)); print(Object.keys(n.kids).sort().map(k => n.kids[k].type === 'dir' ? `<span class="t-dir">${esc(k)}/</span>` : esc(k)).join('    ') || ''); },
    cd: a => { const p = resolve(a[0]); const n = node(p); if (!n) return print(esc(`cd: ${a[0]}: No such file or directory`), 'err'); if (n.type !== 'dir') return print(esc(`cd: ${a[0]}: Not a directory`), 'err'); cwd = p; },
    cat: a => { if (!a.length) return print('cat: give a file name', 'err'); for (const f of a) { const r = readFile(f); if (r.err) print(esc('cat: ' + r.err), 'err'); else print(esc(r.text.replace(/\n$/, ''))); } },
    less: a => { if (!a.length) return print('less: give a file name', 'err'); const r = readFile(a[0]); if (r.err) return print(esc('less: ' + r.err), 'err');
      print(esc(r.text.replace(/\n$/, ''))); print('(less shows one screen at a time: Space for the next page, q to quit. Here the whole file is shown.)', 'note'); },
    grep: a => {
      const flags = a.filter(x => /^-[a-zA-Z]+$/.test(x)).join(''), rest = a.filter(x => !/^-[a-zA-Z]+$/.test(x));
      if (rest.length < 2) return print('usage: grep [-n] [-c] [-i] PATTERN FILE', 'err');
      let re; try { re = new RegExp(rest[0], flags.includes('i') ? 'i' : ''); } catch { return print(esc(`grep: invalid pattern ${rest[0]}`), 'err'); }
      const r = readFile(rest[1]); if (r.err) return print(esc('grep: ' + r.err), 'err');
      const lines = r.text.replace(/\n$/, '').split('\n'), hits = lines.map((l, i) => [i + 1, l]).filter(([, l]) => re.test(l));
      if (flags.includes('c')) return print(String(hits.length));
      print(hits.map(([n, l]) => (flags.includes('n') ? `<span class="t-num">${n}:</span>` : '') + esc(l)).join('\n'));
    },
    mkdir: a => { if (!a.length) return print('mkdir: give a folder name', 'err'); const p = resolve(a[0]); if (node(p)) return print(esc(`mkdir: cannot create directory '${a[0]}': File exists`), 'err');
      const par = parentOf(p); if (!par || par.type !== 'dir') return print(esc(`mkdir: cannot create directory '${a[0]}': No such file or directory`), 'err'); par.kids[p.at(-1)] = { type: 'dir', kids: {} }; },
    cp: a => { if (a.length < 2) return print('usage: cp SOURCE DESTINATION', 'err'); const r = readFile(a[0]); if (r.err) return print(esc('cp: ' + r.err), 'err');
      let p = resolve(a[1]); const n = node(p); if (n && n.type === 'dir') p = [...p, resolve(a[0]).at(-1)];
      const par = parentOf(p); if (!par || par.type !== 'dir') return print(esc(`cp: cannot create '${a[1]}': No such file or directory`), 'err'); par.kids[p.at(-1)] = { type: 'file', text: r.text }; },
    nano: a => { if (!a.length) return print('nano: give a file name', 'err'); const p = resolve(a[0]); const n = node(p);
      if (n && n.type === 'dir') return print(esc(`nano: ${a[0]} is a directory`), 'err');
      const par = parentOf(p); if (!par || par.type !== 'dir') return print(esc(`nano: cannot open ${a[0]}`), 'err');
      ed.hidden = false; ed.dataset.path = JSON.stringify(p); ed.querySelector('.term-ed-head').textContent = `  GNU nano   ${a[0]}${n ? '' : '   [ New File ]'}`; edTa.value = n ? n.text : ''; edTa.focus(); return 'wait'; },
    clear: () => { out.innerHTML = ''; },
  };
  async function lmp(a, np = 1) {
    const i = a.indexOf('-in'); const f = i >= 0 ? a[i + 1] : null;
    if (!f) return print('LAMMPS needs an input script: lmp -in FILE', 'err');
    const r = readFile(f); if (r.err) return print(esc(`ERROR: Cannot open input script ${f}: No such file or directory`), 'err');
    if (np > 1) print(esc(`(This practice terminal has one processor, so LAMMPS runs on 1 instead of ${np}. On a real computer, mpirun -np ${np} starts ${np} copies that share the work.)`), 'note');
    print('<span class="note">running the real LAMMPS in your browser…</span>');
    const res = await runLammps(r.text);
    if (res.version) print(esc(res.version));
    print(esc(res.lines.join('\n')));
    if (res.error) print(esc(res.error), 'err');
    node(cwd).kids['log.lammps'] = { type: 'file', text: [res.version, ...res.lines].join('\n') + '\n' };
    if (opts.onRun) opts.onRun({ script: r.text, res });
  }
  async function exec(line) {
    print(`<span class="term-prompt">${esc(promptEl.textContent)}</span>${esc(line)}`);
    const t = tokenize(line); if (!t.length) return;
    let [c, ...a] = t;
    if (c === 'mpirun' || c === 'mpiexec') { const k = a.indexOf('-np'); const np = k >= 0 ? parseInt(a[k + 1]) || 1 : 1; const rest = k >= 0 ? a.slice(k + 2) : a; if (!/^(\.\/)?lmp/.test(rest[0] || '')) return print('usage: mpirun -np 4 lmp -in FILE', 'err'); return lmp(rest.slice(1), np); }
    if (/^(\.\/)?lmp(_serial|_mpi)?$/.test(c)) return lmp(a);
    if (cmds[c]) return cmds[c](a);
    print(esc(`bash: ${c}: command not found`), 'err');
  }
  ed.querySelector('[data-ed=save]').onclick = () => { const p = JSON.parse(ed.dataset.path); parentOf(p).kids[p.at(-1)] = { type: 'file', text: edTa.value.endsWith('\n') ? edTa.value : edTa.value + '\n' }; ed.hidden = true; print('<span class="note">[ saved ]</span>'); input.focus(); };
  ed.querySelector('[data-ed=quit]').onclick = () => { ed.hidden = true; input.focus(); };
  input.addEventListener('keydown', async e => {
    if (e.key === 'Enter') {
      if (busy) return; const line = input.value; input.value = ''; if (line.trim()) { hist.push(line); hi = hist.length; }
      busy = true; input.disabled = true; try { await exec(line); } finally { busy = false; input.disabled = false; prompt(); if (ed.hidden) input.focus(); }
    } else if (e.key === 'ArrowUp') { if (hi > 0) { hi--; input.value = hist[hi]; e.preventDefault(); } }
    else if (e.key === 'ArrowDown') { if (hi < hist.length) { hi++; input.value = hist[hi] || ''; e.preventDefault(); } }
  });
  el.addEventListener('click', e => { if (!e.target.closest('.term-editor') && !e.target.closest('button')) input.focus({ preventScroll: true }); });
  prompt();
  print('<span class="note">Practice terminal. Type a command and press Enter, or tap a command in the list. Type help for the list.</span>');
  return {
    type: async line => { if (busy) return; input.value = ''; busy = true; input.disabled = true; try { await exec(line); } finally { busy = false; input.disabled = false; prompt(); } },
    file: p => { const n = node(resolve(p)); return n && n.type === 'file' ? n.text : null; },
  };
}
