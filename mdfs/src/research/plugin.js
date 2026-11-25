// Research · E.10: same style, no rebuild: a plugin. fix temp/andersen (page E.6) removed from LAMMPS and loaded at
// run time from a shared library instead; byte-identical output, checked again by hand; what the PLUGIN package costs
// once, what a plugin costs each time, and what LAMMPS does and does not check when it loads one.
import { runLammps } from './lammps-runner.js';
import { NATIVE } from './plugin-native.js';
import { NATIVE as E6 } from './andersen-native.js';
import { RanMars } from './ranmars.js';
import { crlog } from './crpow.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

const BUILD = `# once: LAMMPS with the PLUGIN package (exports its symbols to loaded libraries)
cd lammps/build
cmake -D PKG_PLUGIN=yes .
make -j2                                   # every file recompiles: ${NATIVE.rebuild.toFixed(0)} s here

# each time the style changes: compile two files into a shared library
F="-std=c++17 -O3 -DNDEBUG -DLAMMPS_SMALLBIG -fPIC -I lammps/src"
mpicxx $F -c fix_temp_andersen.cpp
mpicxx $F -c andersenplugin.cpp
mpicxx -shared -rdynamic -o andersenplugin.so andersenplugin.o fix_temp_andersen.o   # ${NATIVE.compile.toFixed(1)} s`;
const frames = t => { const L = t.split('\n'), F = []; for (let i = 0; i < L.length; i++) if (L[i].startsWith('ITEM: TIMESTEP')) { const n = +L[i + 3]; F.push(L.slice(i + 9, i + 9 + n).map(l => l.trim().split(/\s+/).map(Number))); i += 8 + n; } return F; };

export default {
  id: 'plugin',
  title: 'Same style, no rebuild: a plugin',
  minutes: 15,
  summary: 'Every change so far meant rebuilding LAMMPS. A plugin is a style compiled on its own into a shared library and loaded by the input script: plugin load andersenplugin.so. You take the Andersen thermostat of page E.6 out of LAMMPS, wrap it in a 25-line registration function, compile it in under a second, load it, and get output byte-identical to the built-in version (checked again by hand). Then the fine print: the PLUGIN package costs one full rebuild, a version mismatch is only a message, compiler settings must match, and a plugin can silently override a built-in style.',
  learn: [
    'lammpsplugin_init: registering a style name and a creator at run time',
    'The PLUGIN package: exported symbols, and the one-time rebuild',
    'Building a plugin: same headers, defines and MPI as the host',
    'plugin load / list / unload, and what LAMMPS checks (not much)',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: registering at run time</h2>
    <div class="prose"><p>A built-in style is found through the map that style_pair.h, style_fix.h … fill at compile time (page E.2). The PLUGIN package adds the <code>plugin</code> command (${gh('src/PLUGIN/plugin.cpp')}): it opens a shared library with dlopen, looks up one function, <code>lammpsplugin_init</code>, and calls it with a callback. The library uses the callback to hand over a style type ("fix"), a name ("temp/andersen") and a creator function, and LAMMPS adds them to the same map. From then on <code>fix … temp/andersen</code> works exactly as if it had been compiled in.</p>
<p>For this to link, the plugin's code calls LAMMPS's own functions (Fix::Fix, RanMars, error->all …) which live in the <code>lmp</code> executable; the PLUGIN package builds it with exported symbols (CMake ENABLE_EXPORTS, i.e. -rdynamic) so the dynamic linker can find them.</p></div>
  </section>

  <section class="block">
    <h2>2 · Change: out of LAMMPS, into a library</h2>
    <div class="prose"><p>fix_temp_andersen.{h,cpp} from page E.6 are moved out of src/ (so LAMMPS no longer has the style: <code>lmp -h</code> does not list it) and compiled with one extra file:</p></div>
    <pre class="lmp-mini">${esc(NATIVE.wrapper)}</pre>
    <pre class="lmp-mini">${esc(BUILD)}</pre>
    <div class="prose"><p>The flags matter: the same C++ standard and, above all, the same <code>-DLAMMPS_SMALLBIG</code> (it sets the integer sizes of atom IDs and step numbers, so class layouts must agree), and mpicxx because this LAMMPS uses MPI (a serial LAMMPS needs <code>-I src/STUBS</code> instead). The library is ${(NATIVE.sosize / 1024).toFixed(0)} kB.</p></div>
  </section>

  <section class="block">
    <h2>3 · Match: the E.6 gas, from a plugin</h2>
    <pre class="lmp-mini">${esc(NATIVE.gasin)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.out.join('\n'))}</pre>
    <div class="prose"><p>The dump of this run is byte for byte the file the built-in fix wrote on page E.6 (compared with <code>cmp</code>). So the hand check of E.6 applies unchanged; rerun it on that dump:</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ Rerun the random stream by hand</button></div>
    <div class="prose pl-asum"></div>
  </section>

  <section class="block">
    <h2>4 · Observe: the fine print</h2>
    <div class="prose"><p><b>A program without the PLUGIN package</b> does not know the command at all. Your browser's LAMMPS:</p></div>
    <div class="row"><button class="btn sm primary" data-k="b">▶ plugin load here</button><span class="lmp-status muted pl-stb"></span></div>
    <pre class="lmp-mini pl-b" hidden></pre>
    <div class="prose"><p><b>Versions.</b> The same library rebuilt with <code>plugin.version = "2 Aug 2023"</code> loads with one informational line and runs:</p>
<pre class="lmp-mini">Loading plugin: Andersen thermostat (teaching example) v1.0 by MD from scratch
${esc(NATIVE.oldver)}</pre>
<p>Nothing stops it (${gh('src/PLUGIN/plugin.cpp', 162)}). Here it worked because the classes it uses did not change, but a real version mismatch (a member added to Fix, a changed function signature) gives a crash or silently wrong numbers. Rebuild plugins for every LAMMPS you load them into.</p>
<p><b>Overriding.</b> A plugin may register a name LAMMPS already has (<code>nve</code>, say). For pair, bond, angle … styles LAMMPS prints a warning ("Overriding built-in … style from plugin"), and from then on the name means the plugin. Loading the same plugin twice is ignored with a message; <code>plugin unload fix temp/andersen</code> deletes any fix of that style and forgets the name.</p></div>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Why did enabling the PLUGIN package recompile all ~450 files?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>CMake adds the definition -DLMP_PLUGIN to every compilation of the library, and make rebuilds any object whose compile command changed. 193 s once; afterwards each plugin change costs one second instead of a ~12 s rebuild and relink.</div></div>
    <div class="ex"><p><b>2.</b> Why -fPIC?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A shared library is mapped at an address chosen at load time; position-independent code reaches its data and functions through tables instead of fixed addresses. Without it the linker refuses to build the .so on x86-64.</div></div>
    <div class="ex"><p><b>3.</b> When is a plugin the right choice, and when a normal build?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>A plugin suits code you share with people who use their own LAMMPS build (one file to give them, no need to rebuild their LAMMPS), or a style you iterate on quickly. For production on a cluster, compiling it in is simpler: no library path, no ABI worries, and it shows in <code>lmp -h</code>.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('[data-k=a]').onclick = () => {
      const F = frames(E6.dump), mass = { 1: 1.0, 2: 3.0 }, rng = new RanMars(58923 + 0), typ = F[0].map(a => a[1]);
      let v = F[0].map(a => a.slice(2, 5)), ok = 0, n = 0, nc = 0;
      for (let s = 1; s <= 30; s++) {
        for (let i = 0; i < v.length; i++) if (rng.uniform() < 20.0 * 0.005) { const sg = Math.sqrt(1.0 * 1.5 / (mass[typ[i]] * 1.0)); v[i] = [sg * rng.gaussian(crlog), sg * rng.gaussian(crlog), sg * rng.gaussian(crlog)]; nc++; }
        for (let i = 0; i < v.length; i++) for (let k = 0; k < 3; k++) { n++; if (F[s][i][2 + k] === v[i][k]) ok++; }
      }
      $('.pl-asum').innerHTML = `<p>${ok} of ${n} velocity components identical to the last bit ${ok === n ? OK : BAD}, ${nc} collisions, as f_2 printed. The same source was compiled twice, once into lmp (page E.6) and once into a library, with the same flags, and loading it at run time changed nothing.</p>`;
    };
    $('[data-k=b]').onclick = async e => {
      const b = e.currentTarget, st = $('.pl-stb'); b.disabled = true; st.textContent = 'Running…';
      const res = await runLammps('plugin load andersenplugin.so'); b.disabled = false; st.textContent = '';
      $('.pl-b').hidden = false; $('.pl-b').textContent = res.lines.filter(l => /ERROR|Last input/.test(l)).join('\n') || res.error || '(no error)';
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
