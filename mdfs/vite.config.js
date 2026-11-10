import { defineConfig } from 'vite';

// base '/' because pages live at nested paths (/lessons/<id>); every Vercel deploy serves from its domain root.
// worker format 'es': the LAMMPS WebAssembly worker uses ES module syntax.
export default defineConfig({ base: '/', worker: { format: 'es' } });
