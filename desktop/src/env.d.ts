/// <reference types="vite/client" />

import type { DesktopApi } from './preload';

declare global {
  interface Window {
    /** Injected by the preload script via contextBridge. */
    readonly desktop: DesktopApi;
  }
}

export {};
