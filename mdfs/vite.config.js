import { defineConfig } from 'vite';

// base '/' because pages live at nested paths (/lessons/<id>); every Vercel deploy serves from its domain root.
export default defineConfig({ base: '/' });
