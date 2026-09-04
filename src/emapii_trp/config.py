"""Load ``config.yaml`` from the repository root."""
from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from . import REPO_ROOT


def load_config(path: Path | None = None) -> dict[str, Any]:
    """Return the parsed configuration dictionary (default: ``<repo>/config.yaml``)."""
    cfg_path = path or REPO_ROOT / "config.yaml"
    with cfg_path.open() as fh:
        return yaml.safe_load(fh)
