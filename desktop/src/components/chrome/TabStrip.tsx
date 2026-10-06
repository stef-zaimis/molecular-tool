import { HomeIcon } from '../icons/Icons';
import './TabStrip.css';

export interface TabDescriptor {
  readonly id: string;
  readonly label: string;
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
 * Home is the leftmost item and represents the project page itself, not a
 * "back" action. The analysis tabs are whatever the caller passes in — the
 * caller derives them from the analyses selected for the project, so an
 * unselected analysis has no tab at all rather than a disabled one.
 *
 * Visual language: the current item takes the body colour and merges into the
 * page beneath it; the others sit dark and recede.
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
        // Already on the project page: activating Home would be a no-op, so it
        // is inert rather than pretending to navigate.
        onClick={homeActive ? undefined : onHome}
        aria-label={homeLabel}
        aria-current={homeActive ? 'page' : undefined}
        aria-disabled={homeActive || undefined}
      >
        <HomeIcon size={34} />
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
