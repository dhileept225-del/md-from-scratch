# MD from Scratch

A free, hands-on course in molecular dynamics: from school maths to running and reading LAMMPS.

**Live site: https://md-from-scratch.vercel.app**

- **Lessons and labs**: every idea of an MD code (forces, integrators, thermostats, Ewald, analysis) built step by step, with live simulations and every number checked by hand.
- **LAMMPS from scratch**: the real LAMMPS runs in your browser (WebAssembly). Each page explains part of the program, lets you change the input or the code, and matches the result with a hand calculation.

Made by Thummuru Dhileep (Okayama University). Corrections and suggestions are welcome as GitHub issues.

## Run it locally

```
cd mdfs
npm install
npm run dev      # development server
npm run build    # static site in mdfs/dist, one HTML page per lesson
```
