import { HomeIcon } from '../icons/Icons';
import './TabStrip.css';

export interface TabDescriptor {
  readonly id: string;
  readonly label: string;
  readonly enabled: boolean;
}

interface TabStripProps {
  readonly tabs?: readonly TabDescriptor[];
  readonly activeId?: string | null;
  /** True on the project-creation screen, where home itself is the current view. */
  readonly homeActive?: boolean;
  readonly onHome: () => void;
  readonly onSelect?: (id: string) => void;
  readonly homeLabel?: string;
}

/**
 * The strip below the title bar.
 *
 * In the mockups the home button and the analysis tabs share one visual
 * language: the current item takes the body colour and appears to merge into
 * the page beneath it, while the others sit dark and recede. Disabled tabs use
 * the same treatment as inactive ones — in the reference, unselected analyses
 * are rendered with body-coloured text on the dark fill, which is why they read
 * as unavailable rather than merely unselected.
 */
export function TabStrip({
  tabs = [],
  activeId = null,
  homeActive = false,
  onHome,
  onSelect,
  homeLabel = 'Home',
}: TabStripProps): JSX.Element {
  return (
    <nav className="tabstrip" aria-label="Analyses">
      <button
        type="button"
        className={`tabstrip__home${homeActive ? ' is-active' : ''}`}
        onClick={onHome}
        aria-label={homeLabel}
        aria-current={homeActive ? 'page' : undefined}
      >
        <HomeIcon size={36} />
      </button>

      {tabs.length > 0 && (
        <div className="tabstrip__tabs" role="tablist">
          {tabs.map((tab) => {
            const active = tab.id === activeId;
            return (
              <button
                key={tab.id}
                type="button"
                role="tab"
                className={`tabstrip__tab${active ? ' is-active' : ''}`}
                aria-selected={active}
                disabled={!tab.enabled}
                title={tab.enabled ? undefined : `${tab.label} was not selected for this project`}
                onClick={() => onSelect?.(tab.id)}
              >
                {tab.label}
              </button>
            );
          })}
        </div>
      )}
    </nav>
  );
}
