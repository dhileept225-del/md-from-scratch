// Research · Part 0.3: getting LAMMPS without compiling, and checking the installation by hand.
import { g8, matchMarkup } from './lammps-runner.js';

const row = c => `<tr>${c.map(x => `<td>${x}</td>`).join('')}</tr>`;
const head = c => `<tr>${c.map(x => `<th>${x}</th>`).join('')}</tr>`;
const HAND = { Step: 0, Temp: 3, E_pair: -6.7733681, E_mol: 0, TotEng: -2.2744931, Press: -3.7033504 };

const OS = {
  windows: `<p>Download an installer from the official LAMMPS site and run it. It installs the command-line program <code>lmp.exe</code>, LAMMPS-GUI, example inputs and potential files.</p>
<div class="calc">1. open  https://packages.lammps.org/windows.html
   (or the LAMMPS-GUI installer on the LAMMPS releases page on GitHub)
2. run the .exe installer
   (Windows may warn that it is unsigned: allow it only for files from lammps.org)
3. open PowerShell and type:   lmp -h</div>
<p>To run on several processors you also need Microsoft MPI (MS-MPI), as the LAMMPS Windows instructions explain.</p>`,
  mac: `<p>Two easy ways. LAMMPS-GUI comes as a disk image with the command-line program inside; Homebrew installs the command-line program only.</p>
<div class="calc"># way 1: LAMMPS-GUI
open the LAMMPS-macOS-multiarch-GUI-&lt;version&gt;.dmg from the releases page on GitHub,
drag the app into Applications

# way 2: Homebrew (install Homebrew first from brew.sh)
brew install lammps
ls $(brew --prefix)/bin | grep lmp     # shows the program's name (lmp, lmp_serial or lmp_mpi)</div>`,
  linux: `<p>LAMMPS-GUI is available as a flatpak bundle that also contains the command-line program.</p>
<div class="calc">cd ~/Downloads
flatpak install --user LAMMPS-GUI-Linux-x86_64-GUI-&lt;version&gt;.flatpak
flatpak run --command=lmp org.lammps.lammps-gui -in in.melt</div>
<p>Your distribution may also offer a package (for example <code>sudo apt install lammps</code> on Ubuntu), but it is often an older version: check with <code>lmp -h</code>.</p>`,
  conda: `<p>Works on Linux and macOS once Miniconda is installed. These are the commands from the LAMMPS manual:</p>
<div class="calc">conda config --add channels conda-forge
conda create -n my-lammps-env
conda activate my-lammps-env
conda install lammps
lmp -h</div>`,
};

export default {
  id: 'get-lammps',
  title: 'Getting LAMMPS without compiling',
  minutes: 15,
  summary: 'Ready-made LAMMPS programs exist for Windows, macOS and Linux. You install one, ask it which version and packages it has, run the melt example, and prove the installation works by matching its first line to the hand calculation.',
  learn: [
    'The ready-made ways to install LAMMPS on each system, and what they limit',
    'How to read lmp -h: the version and the installed packages',
    'Check an installation by matching the melt step-0 line to a hand calculation',
  ],

  html: () => `
  <section class="block prose">
    <h2>1 · Understand: two ways to get LAMMPS</h2>
    <p>Someone has to compile LAMMPS (page 0.1). Either the LAMMPS developers do it and you download the program (this page), or you compile it yourself (page 0.4). Downloading is quicker; compiling lets you choose the version and packages, and is the only way to change the code.</p>
    <div class="tbl wide"><table>${head(['', 'Download a program', 'Compile it yourself'])}
${row(['time to first run', 'minutes', '15–30 minutes the first time, including installing the tools'])}
${row(['version', 'whatever was packaged', 'any release, including the newest'])}
${row(['packages', 'chosen by the packager', 'exactly the ones you need'])}
${row(['change the C++ code', 'not possible', 'yes (Part E)'])}</table></div>
    <p>No installation at all is also possible: the <a href="/lammps/lammps-start">Start page</a> runs LAMMPS inside your browser.</p>
  </section>

  <section class="block">
    <h2>2 · Install for your system</h2>
    <div class="seg os-tabs" role="tablist"><button class="on" data-os="windows">Windows</button><button data-os="mac">macOS</button><button data-os="linux">Linux</button><button data-os="conda">conda</button></div>
    <div class="prose os-body"></div>
    <p class="muted small">Steps follow the official LAMMPS installation pages (sources at the end). The program may be called <code>lmp</code>, <code>lmp_serial</code>, <code>lmp_mpi</code> or <code>lmp.exe</code> depending on the installer; use your name wherever this course writes lmp.</p>
  </section>

  <section class="block prose">
    <h2>3 · Observe: ask LAMMPS about itself</h2>
    <p><code>lmp -h</code> prints help. Its first lines say which version you have, and further down it lists the installed packages. Here is the start of ours, on the 2-processor computer used for this course:</p>
    <div class="calc">$ lmp -h

Large-scale Atomic/Molecular Massively Parallel Simulator - 22 Jul 2025 - Update 6
Git info (HEAD / stable_22Jul2025_update6)
…
Installed packages:

KSPACE MOLECULE RIGID</div>
    <div class="tbl wide"><table>${head(['Version line says', 'What it means'])}
${row(['30 Sep 2026', 'the newest stable release (the one running in your browser on these pages)'])}
${row(['22 Jul 2025', 'the previous stable release (the source code quoted in this course)'])}
${row(['a date with no "stable" tag', 'a feature release: newer features, less tested'])}</table></div>
    <p>A command from a package that is not installed stops LAMMPS with a clear message. Ours, asked for <code>pair_style tersoff</code>:</p>
    <div class="calc">ERROR: Unrecognized pair style 'tersoff' is part of the MANYBODY package which is not enabled in this LAMMPS binary.</div>
  </section>

  <section class="block prose">
    <h2>4 · Match: prove your installation works</h2>
    <p>Run the melt example (it comes with every installer, in <code>examples/melt/</code>, or copy it from the <a href="/lammps/lammps-start">Start page</a>):</p>
    <div class="calc">lmp -in in.melt</div>
    <p>Its first thermo line is fixed by the crystal and the starting temperature, so every correct installation of every version must print exactly the numbers we computed by hand on the Start page. Paste your line below.</p>
    <div class="lmp paste">
      <input class="paste-in mono" value="         0   3             -6.7733681      0             -2.2744931     -3.7033504" aria-label="your step 0 line">
      <div class="row"><button class="btn sm primary" data-k="check">Check my line</button><span class="muted small">(the box holds our own line; replace it with yours)</span></div>
      <div class="paste-out"></div>
    </div>
    <p>Later lines (step 50, 100, …) may differ slightly between versions, compilers and numbers of processors: tiny rounding differences grow in a chaotic liquid. Step 0 has no such freedom.</p>
  </section>

  <section class="block prose">
    <h2>Your turn</h2>
    <div class="ex"><p><b>1.</b> Your <code>lmp -h</code> lists <code>KSPACE MANYBODY MOLECULE RIGID</code>. Can you use <code>pair_style tersoff</code>? And <code>pair_style reaxff</code>?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>tersoff yes: it belongs to MANYBODY. reaxff no: it belongs to the REAXFF package, so LAMMPS stops with "Unrecognized pair style 'reaxff' is part of the REAXFF package which is not enabled". You need another installer or your own build (page 0.4).</div></div>
    <div class="ex"><p><b>2.</b> A student's hand calculation gives TotEng = −2.2733681 instead of −2.2744931. What did they forget?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>Their kinetic energy per atom is 4.5 = 1.5 × 3, without the factor (N−1)/N. LAMMPS removes the 3 degrees of freedom of the whole crystal's motion, so KE/N = 1.5 × 3 × 3999/4000 = 4.4988750, and −6.7733681 + 4.4988750 = −2.2744931.</div></div>
    <div class="ex"><p><b>3.</b> Two correct installations print the same step 0 but different numbers at step 250. Is one of them broken?</p><button class="btn sm ghost" data-ex>Show answer</button><div class="ans" hidden>No. Different compilers, versions or processor counts add numbers in a slightly different order, so they round differently in the last digit. In a liquid such tiny differences grow step by step until the trajectories differ, while averages (temperature, energy) stay the same. Step 0 involves no time stepping, so it must agree exactly.</div></div>
  </section>

  <section class="block prose">
    <h2>Sources</h2>
    <ul class="small">
      <li><a href="https://www.lammps.org/download/" target="_blank" rel="noopener">LAMMPS download page</a> (stable releases 30 Sep 2026 and 22 Jul 2025)</li>
      <li><a href="https://docs.lammps.org/Install_windows.html" target="_blank" rel="noopener">LAMMPS manual: Windows</a>, <a href="https://docs.lammps.org/Install_mac.html" target="_blank" rel="noopener">macOS</a>, <a href="https://docs.lammps.org/Install_conda.html" target="_blank" rel="noopener">conda</a></li>
      <li><a href="https://lammps-gui.lammps.org/installation.html" target="_blank" rel="noopener">LAMMPS-GUI installation</a></li>
    </ul>
  </section>`,

  mount(root) {
    const body = root.querySelector('.os-body');
    const show = k => { body.innerHTML = OS[k]; root.querySelectorAll('[data-os]').forEach(b => b.classList.toggle('on', b.dataset.os === k)); };
    root.querySelectorAll('[data-os]').forEach(b => b.onclick = () => show(b.dataset.os));
    const ua = navigator.userAgent; show(/Windows/.test(ua) ? 'windows' : /Mac OS X|Macintosh/.test(ua) && !/iPhone|iPad/.test(ua) ? 'mac' : /Linux/.test(ua) && !/Android/.test(ua) ? 'linux' : 'windows');
    const inp = root.querySelector('.paste-in'), out = root.querySelector('.paste-out');
    root.querySelector('[data-k=check]').onclick = () => {
      const w = inp.value.trim().split(/\s+/);
      if (w.length < 6 || w[0] !== '0') { out.innerHTML = '<p class="lmp-err">Paste the line that starts with 0 (step 0): six numbers, from Step to Press.</p>'; return; }
      const keys = ['Temp', 'E_pair', 'TotEng', 'Press'], idx = { Temp: 1, E_pair: 2, TotEng: 4, Press: 5 };
      let all = true;
      const rows = keys.map(k => { const m = matchMarkup(g8(HAND[k]), w[idx[k]]); all = all && m.full; return `<tr><td>${k}</td><td class="mono">${m.html}</td><td class="mono">${w[idx[k]]}</td><td>${m.full ? '<span class="m-ok">✓</span>' : '<span class="m-bad">✗</span>'}</td></tr>`; }).join('');
      out.innerHTML = `<div class="tbl"><table><tr><th>step 0</th><th>by hand</th><th>yours</th><th></th></tr>${rows}</table></div><p>${all ? '<span class="m-ok">Every digit matches: your LAMMPS works.</span>' : '<span class="m-bad">Something differs. Check that you ran the unmodified in.melt; if you did, the installation is suspect.</span>'}</p>`;
    };
    root.querySelectorAll('[data-ex]').forEach(b => b.onclick = () => { const a = b.nextElementSibling; a.hidden = !a.hidden; b.textContent = a.hidden ? 'Show answer' : 'Hide answer'; });
    return () => {};
  },
};
