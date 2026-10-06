import { defineConfig } from 'vitest/config';
import react from '@vitejs/plugin-react';

export default defineConfig({
  plugins: [react()],
  test: {
    environment: 'jsdom',
    include: ['src/**/*.test.{ts,tsx}'],
    globals: true,
    // jsdom has no layout engine; the editor measures the document. See the
    // setup file for exactly which gaps this fills.
    setupFiles: ['./vitest.setup.ts'],
  },
});
