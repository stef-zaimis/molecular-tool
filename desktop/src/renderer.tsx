import { StrictMode } from 'react';
import { createRoot } from 'react-dom/client';

// JetBrains Mono, bundled as woff2 via @fontsource so the app renders
// identically on Windows, macOS and Linux without a locally installed copy.
import '@fontsource/jetbrains-mono/400.css';
import '@fontsource/jetbrains-mono/500.css';
import '@fontsource/jetbrains-mono/700.css';

import './styles/tokens.css';
import './styles/global.css';

import { App } from './app/App';
import { ProjectProvider } from './app/state/ProjectContext';

const container = document.getElementById('root');
if (!container) throw new Error('Root container #root was not found.');

createRoot(container).render(
  <StrictMode>
    <ProjectProvider>
      <App />
    </ProjectProvider>
  </StrictMode>,
);
