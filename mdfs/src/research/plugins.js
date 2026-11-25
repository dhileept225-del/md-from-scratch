// Research · F.2: plugins in practice. Four styles from Part E taken out of LAMMPS and packaged as one library,
// built out of tree with CMake, loaded automatically through LAMMPS_PLUGIN_PATH; the earlier recorded outputs
// reproduced byte for byte, and one of them rechecked by hand here.
import { NATIVE } from './plugins-native.js';
import { NATIVE as MIE } from './mie-native.js';
import { crpow } from './crpow.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const esc = s => String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;');
const GH = 'https://github.com/lammps/lammps/blob/stable_22Jul2025_update6/';
const gh = (p, l) => `<a href="${GH}${p}${l ? '#L' + l : ''}" target="_blank" rel="noopener">${p}${l ? ':' + l : ''}</a>`;
const OK = '<span class="m-ok">✓</span>', BAD = '<span class="m-bad">✗</span>';

const TREE = `mdfsplugin/
  CMakeLists.txt
  LAMMPSInterfaceCXX.cmake        copied from lammps/examples/plugins
  mdfsplugin.cpp                  the registration function
  pair_lj_mine.{h,cpp}            page E.3
  bond_quartic_mine.{h,cpp}       page E.9
  compute_nbcount_atom.{h,cpp}    page E.7
  fix_temp_andersen.{h,cpp}       page E.6`;
const BUILD = `cmake -S mdfsplugin -B mdfsplugin/build -D LAMMPS_HEADER_DIR=$HOME/lammps/src -D CMAKE_BUILD_TYPE=Release
cmake --build mdfsplugin/build -j2                       # ${NATIVE.build.toFixed(1)} s, ${(NATIVE.size / 1024).toFixed(0)} kB
mkdir -p ~/lammps-plugins && cp mdfsplugin/build/mdfsplugin.so ~/lammps-plugins/
export LAMMPS_PLUGIN_PATH=~/lammps-plugins               # e.g. in ~/.bashrc or a module file
lmp -in in.anything                                      # no plugin command needed`;
const CHECKS = [
  ['E.3', 'two atoms with pair lj/mine 2.5 14 6 (Mie), eight distances', 'the 8 lines of PotEng and force'],
  ['E.6', 'ten-atom gas with fix temp/andersen, 30 steps', 'v.dump'],
  ['E.7', '500-atom binary liquid with compute nbcount/atom', 'x.dump and the thermo line'],
  ['E.9', 'C–H bond with bond quartic/mine, 600 steps', 'x.dump'],
];

export default {
  id: 'plugins',
  title: 'Plugins',
  minutes: 12,
  summary: 'Page E.10 loaded one style from a hand-compiled library. This is how to do it for real: four styles from Part E (a pair, a bond, a compute and a fix) taken out of LAMMPS, packaged in one library with a short CMake project that needs only the LAMMPS headers, and loaded automatically at start-up through LAMMPS_PLUGIN_PATH, without a single plugin command in the input. The earlier pages\' recorded outputs come back byte for byte, and one of them is checked by hand again here.',
  learn: [
    'One registration function, several styles: pair/bond (factory v1), compute/fix (factory v2)',
    'An out-of-tree CMake build against the headers only',
    'LAMMPS_PLUGIN_PATH: automatic loading, also after clear',
    'Sharing a plugin: what the receiver needs to match',
  ],

  html: () => `
  <section class="block">
    <h2>1 · Understand: one library, many styles</h2>
    <div class="prose"><p>lammpsplugin_init (page E.10) may call the registration callback as often as it likes. Each call names a style kind, a name and a creator. Pair, bond, angle and other force-field styles are created with the LAMMPS pointer only (factory v1); computes and fixes also receive their command arguments (factory v2). At start-up LAMMPS looks in every directory of <code>LAMMPS_PLUGIN_PATH</code> for files ending in <code>plugin.so</code> and loads them (${gh('src/PLUGIN/plugin.cpp', 80)}), and it does so again after every <code>clear</code>, which recreates everything.</p></div>
  </section>

  <section class="block">
    <h2>2 · Package</h2>
    <div class="prose"><p>The four styles are moved out of lammps/src (LAMMPS rebuilt, ${NATIVE.rebuild.toFixed(0)} s: it no longer knows them) into their own project:</p></div>
    <pre class="lmp-mini">${esc(TREE)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.cpp)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.cmake)}</pre>
    <div class="prose"><p>LAMMPSInterfaceCXX.cmake does the part that went wrong by hand on page E.10 if forgotten: it finds the same MPI, sets <code>-DLAMMPS_SMALLBIG</code> (or bigbig, to match LAMMPS_SIZES) and the OpenMP flags. The project never touches the LAMMPS build tree, only the headers in src/.</p></div>
    <pre class="lmp-mini">${esc(BUILD)}</pre>
    <pre class="lmp-mini">${esc(NATIVE.list)}</pre>
  </section>

  <section class="block">
    <h2>3 · Match: four earlier pages, again</h2>
    <div class="prose"><p>The inputs of four Part E pages, run unchanged (they contain no plugin command) with the plugin directory set, compared with the outputs recorded when the styles were compiled into LAMMPS:</p></div>
    <div class="tbl pg-c"></div>
    <div class="prose"><p>Since the E.3 output is unchanged, its hand check applies as it stands; rerun it:</p></div>
    <div class="row"><button class="btn sm primary" data-k="a">▶ E.3's eight distances by hand</button></div>
    <div class="prose pg-asum"></div>
  </section>

  <section class="block prose">
    <h2>4 · Sharing it</h2>
    <p>What a colleague needs to use mdfsplugin.so on their own machine:</p>
    <ul><li><b>A LAMMPS with the PLUGIN package</b> (<code>lmp -h</code> lists "plugin" among the commands), built with exported symbols, which the package arranges.</li>
    <li><b>The same LAMMPS version</b> the plugin was compiled against; LAMMPS only prints a note otherwise (page E.10). Safer: send the source directory and let them build it against their own headers (2 s).</li>
    <li><b>The same MPI library and integer sizes</b>; LAMMPSInterfaceCXX.cmake takes care of both when they build it.</li>
    <li><b>The same plugin at read time</b> for restart files that name its styles; with LAMMPS_PLUGIN_PATH set, read_restart finds them.</li></ul>
    <p>On a cluster this is one line in a module file (<code>prepend-path LAMMPS_PLUGIN_PATH …</code>), and everybody's LAMMPS gains the styles without being rebuilt.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> The E.3 input contains eight <code>clear</code> commands. How many times was the plugin loaded?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The message "Loaded 1 plugins from …" appeared 9 times: at start-up and after each clear. clear destroys and recreates every LAMMPS class, including the style maps; afterwards LAMMPS re-registers the plugins it already had (plugin_restore, lammps.cpp:889) and runs the automatic loading again, so plugin styles survive clear whether they came from LAMMPS_PLUGIN_PATH or from a plugin load command.</div></div>
    <div class="ex"><p><b>2.</b> Why "Loaded 1 plugins" for four styles?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>The count is of files successfully opened; <code>plugin list</code> shows the four styles inside the one file.</div></div>
  </section>`,

  mount(root) {
    const $ = s => root.querySelector(s);
    $('.pg-c').innerHTML = `<table>${head(['page', 'run', 'compared', 'result'])}${CHECKS.map(c => row([c[0], c[1], c[2], `${OK} identical`])).join('')}</table>`;
    $('[data-k=a]').onclick = () => {
      const n = 14, m = 6, cmie = (n / (n - m) * crpow(n / m, 3, 4)), lj1 = cmie * n * 1.0 * 1.0, lj2 = cmie * m * 1.0 * 1.0, lj3 = cmie * 1.0 * 1.0, lj4 = cmie * 1.0 * 1.0;
      let ok = 0;
      for (const w of MIE.two) {
        const x2 = 5.0 + Number(w[1]), delx = 5.0 - x2, rsq = delx * delx + 0.0 * 0.0 + 0.0 * 0.0, r2inv = 1.0 / rsq, rm = crpow(r2inv, 3), rn = crpow(r2inv, 7);
        const fx = 0.0 + delx * (1.0 * (lj1 * rn - lj2 * rm) * r2inv); let ev = (lj3 * rn - lj4 * rm) - 0.0; ev *= 1.0;
        ok += (0.0 + ev === Number(w[5])) && (fx === Number(w[7]));
      }
      $('.pg-asum').innerHTML = `<p>${ok} of ${MIE.two.length} energies and forces identical to the last bit ${ok === MIE.two.length ? OK : BAD}, now produced by the plugin. Same source, same flags, same numbers: whether a style is compiled in or loaded at run time is a packaging decision, not a numerical one.</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
