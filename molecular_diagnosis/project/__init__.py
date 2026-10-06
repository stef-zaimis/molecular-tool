"""
Project persistence: one SQLite database per project, owned by this service.

Nothing above this package may execute SQL. The Electron renderer and the
main process talk to `molecular_diagnosis.service`, which talks to
`ProjectService`, which is the only thing that opens the database.
"""

from molecular_diagnosis.project.service import ProjectService
from molecular_diagnosis.project.db import (
    Capabilities,
    connect,
    detect_capabilities,
    migrate,
    open_project_db,
)

__all__ = [
    "Capabilities",
    "ProjectService",
    "connect",
    "detect_capabilities",
    "migrate",
    "open_project_db",
]
