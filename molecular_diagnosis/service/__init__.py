"""
UI-independent service boundary.

`molecular_diagnosis.service` is the only surface the desktop frontend talks
to. It exposes plain-data operations over a newline-delimited JSON protocol on
stdin/stdout (see `__main__.py`) and knows nothing about Electron, Tkinter, or
any UI toolkit.

Run it with:  python -m molecular_diagnosis.service
"""

from molecular_diagnosis.service.errors import ErrorCode, ServiceError
from molecular_diagnosis.service.handlers import METHODS, dispatch

__all__ = ["METHODS", "ErrorCode", "ServiceError", "dispatch"]
