import { defineConfig } from 'vite';

// base './' keeps every link relative, so the built site works on Vercel and in previews.
export default defineConfig({ base: './' });
