// Research · Part 0.4: building LAMMPS from source, with real output from our build.
import { g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;

export default {
  id: 'build-lammps',
  title: 'Building LAMMPS from source',
  minutes: 25,
  summary: 'Compiling LAMMPS yourself takes four commands. We did it on a small computer and show the real output of every step: how long it took, what each line means, how to add a package, and the three errors beginners meet most. The proof that it worked is, again, the melt step-0 line matched by hand.',
  learn: [
    'The four commands that download, configure, compile and test LAMMPS',
    'How to read the configuration summary and the list of installed packages',
    'Why adding a package means compiling almost everything again, measured',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: what you need</h2>
    <p>Building means running the compiler on the LAMMPS source (page 0.1). Four tools do the work. Check that each one answers in your terminal (page 0.2):</p>
    <div class="tbl wide"><table>${head(['Tool', 'Check with', 'What it does', 'Ours'])}
${row(['git', '<code>git --version</code>', 'downloads the source', '2.43.0'])}
${row(['CMake, at least 3.16', '<code>cmake --version</code>', 'reads your choices and writes the build instructions', '3.28.3'])}
${row(['a C++ compiler', '<code>g++ --version</code>', 'turns the source into a program (C++17 preferred, C++11 the minimum)', 'GCC 13.3.0'])}
${row(['MPI (optional)', '<code>mpirun --version</code>', 'lets one run use many processors', 'MPICH'])}</table></div>
    <p>On Ubuntu one command installs all four: <code>sudo apt install git cmake g++ mpich</code>. On macOS, <code>xcode-select --install</code> and then <code>brew install cmake open-mpi</code>. On Windows, do this inside Ubuntu for Windows (WSL).</p>
  </section>

  <section class="block prose">
    <h2>2 · The four commands, with our real output</h2>
    <p><b>Step 1: download</b> exactly the version this course quotes (the newest is <code>stable_30Sep2026</code>):</p>
    <div class="calc">git clone -b stable_22Jul2025_update6 --depth 1 https://github.com/lammps/lammps.git
cd lammps</div>
    <p><b>Step 2: configure.</b> Make a separate folder for the build, then tell CMake where the instructions are (the <code>cmake</code> folder) and which packages you want:</p>
    <div class="calc">mkdir build
cd build
cmake ../cmake -D PKG_MOLECULE=yes -D PKG_KSPACE=yes -D PKG_RIGID=yes</div>
    <p>CMake checks your computer and prints a summary. These are the lines worth reading in ours:</p>
    <div class="calc">-- Found MPI: TRUE (found version "5.0") found components: CXX      <b>MPI found: parallel runs possible</b>
   LAMMPS Version:   2025.7.22.6 stable_22Jul2025_update6        <b>the version you asked for</b>
   Build type:       Release                                     <b>full optimisation</b>
-- Enabled packages: KSPACE;MOLECULE;RIGID                       <b>your package choice</b>
      C++ Standard:  17                                          <b>the compiler is modern enough</b>
      C++ Flags:     -O3 -DNDEBUG                                <b>-O3: optimise for speed</b></div>
    <p>A ready-made choice also exists: <code>cmake -C ../cmake/presets/basic.cmake ../cmake</code> switches on KSPACE, MANYBODY, MOLECULE and RIGID.</p>
    <p><b>Step 3: compile.</b> <code>-j 2</code> uses 2 processors at once; use the number your computer has.</p>
    <div class="calc">make -j 2

[  8%] Building CXX object CMakeFiles/lammps.dir/…/src/angle_deprecated.cpp.o
[  9%] Building CXX object CMakeFiles/lammps.dir/…/src/angle.cpp.o
…                                                         <b>451 files, one line each</b>
[100%] Linking CXX static library liblammps.a
[100%] Linking CXX executable lmp
[100%] Built target lmp                                   <b>done: 266 s on 2 processors</b></div>
    <p><b>Step 4: test.</b> The new program is the file <code>lmp</code> in the build folder, so run it as <code>./lmp</code>:</p>
    <div class="calc">./lmp -h | head -3
./lmp -in ../examples/melt/in.melt

   Step          Temp          E_pair         E_mol          TotEng         Press
         0   3             -6.7733681      0             -2.2744931     -3.7033504</div>
  </section>

  <section class="block prose">
    <h2>3 · Match: is the build right?</h2>
    <p>A program that compiles can still be wrong (a bad compiler flag, a broken library). The step-0 line is a test with a known answer, worked out by hand on the Start page. Paste your build's line:</p>
    <div class="lmp paste">
      <input class="paste-in mono" value="         0   3             -6.7733681      0             -2.2744931     -3.7033504" aria-label="your step 0 line">
      <div class="row"><button class="btn sm primary" data-k="check">Check my build</button><span class="muted small">(the box holds our line; replace it with yours)</span></div>
      <div class="paste-out"></div>
    </div>
  </section>

  <section class="block prose">
    <h2>4 · Modify and observe: add a package</h2>
    <p>Our build stopped at <code>pair_style tersoff</code>: that force belongs to the MANYBODY package. Adding a package means configuring again (the dot means "this build folder") and compiling again:</p>
    <div class="calc">cmake -D PKG_MANYBODY=yes .
make -j 2</div>
    <div class="tbl wide"><table>${head(['', 'before', 'after adding MANYBODY'])}
${row(['installed packages', 'KSPACE MOLECULE RIGID', 'KSPACE MANYBODY MOLECULE RIGID'])}
${row(['pair styles available', '52', '91'])}
${row(['size of lmp', '11.6 MB', '12.8 MB'])}
${row(['files compiled', '451 (from nothing)', '489 (almost everything again)'])}
${row(['time on 2 processors', '266 s', '274 s'])}
${row(['melt step 0', '−6.7733681', '−6.7733681 (unchanged)'])}</table></div>
    <p>Why compile almost everything again, when only 39 new files were added? Many core files read small generated lists of "which styles exist" and "which packages are installed". Adding a package changes those lists, so every file that reads them must be compiled again. Changing one existing file (page 0.1: 3.8 s) is cheap; changing the package list is not. Choose your packages before starting long work.</p>
  </section>

  <section class="block prose">
    <h2>5 · The three errors beginners meet most</h2>
    <div class="tbl wide"><table>${head(['You see', 'Why', 'Fix'])}
${row(['<code>CMake Error: The source directory "…/lammps" does not appear to contain CMakeLists.txt.</code>', 'cmake was pointed at the lammps folder; its instructions are in lammps/cmake', '<code>cmake ../cmake</code> from inside lammps/build'])}
${row(['<code>make: *** No targets specified and no makefile found.  Stop.</code>', 'make was run before cmake, or in the wrong folder', 'run cmake first, then make, both inside build'])}
${row(['<code>ERROR: Unrecognized pair style \'tersoff\' is part of the MANYBODY package which is not enabled in this LAMMPS binary.</code>', 'the input uses a package you did not build', '<code>cmake -D PKG_MANYBODY=yes .</code> then <code>make</code>'])}</table></div>
    <p>All three messages above are copied from our own terminal.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Your computer has 8 processors. Roughly how long should the build of section 2 take, if 2 processors took 266 s?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Compiling 451 independent files splits well, so about 266 × 2/8 ≈ 67 s, a little more in practice because linking at the end uses one processor. Use <code>make -j 8</code>.</div></div>
    <div class="ex"><p><b>2.</b> You edit one line in <code>src/pair_lj_cut.cpp</code> and type <code>make -j 2</code>. Will it take 266 s again?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. make compiles only the changed file and links the program again: 3.8 s on our computer (page 0.1). Only changing the package list or the compiler settings forces a near-complete rebuild.</div></div>
    <div class="ex"><p><b>3.</b> After adding MANYBODY, why must melt still print E_pair = −6.7733681?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>melt uses only <code>pair_style lj/cut</code>, which did not change. Adding a package adds new styles but must not change existing ones. If step 0 changed, something in the build would be broken.</div></div>
  </section>`,

  mount(root) {
    const HAND = { Temp: 3, E_pair: -6.7733681, TotEng: -2.2744931, Press: -3.7033504 };
    const inp = root.querySelector('.paste-in'), out = root.querySelector('.paste-out');
    root.querySelector('[data-k=check]').onclick = () => {
      const w = inp.value.trim().split(/\s+/);
      if (w.length < 6 || w[0] !== '0') { out.innerHTML = '<p class="lmp-err">Paste the line that starts with 0 (step 0): six numbers, from Step to Press.</p>'; return; }
      const idx = { Temp: 1, E_pair: 2, TotEng: 4, Press: 5 }; let all = true;
      const rows = Object.keys(idx).map(k => { const m = matchMarkup(g8(HAND[k]), w[idx[k]]); all = all && m.full; return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${w[idx[k]]}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('');
      out.innerHTML = `<div class="tbl"><table><tr><th>step 0</th><th>by hand</th><th>your build</th><th></th></tr>${rows}</table></div><p>${all ? '<span class="m-ok">Every digit matches: your build computes forces and energies correctly.</span>' : '<span class="m-bad">Something differs. If the input is the unmodified in.melt, the build is suspect.</span>'}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
