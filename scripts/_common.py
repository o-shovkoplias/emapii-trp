"""Shared bootstrap for the numbered analysis scripts: make ``src/`` importable when run
from the repository root without installation."""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
if str(REPO / "src") not in sys.path:
    sys.path.insert(0, str(REPO / "src"))
