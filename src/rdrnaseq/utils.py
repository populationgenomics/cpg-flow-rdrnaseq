"""
suggested location for any utility methods or constants used across multiple stages
"""

import logging
from datetime import datetime
from os.path import exists
from pathlib import Path

from cpg_utils import config

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005


def can_reuse(
    path: list[Path | str] | Path | str | None,
    overwrite: bool = False,
) -> bool:
    """
    Checks if the object at `path` is good to reuse:
    * overwrite has the default value of False,
    * check_intermediates has the default value of True,
    * object exists.

    If `path` is a collection, it requires all paths to exist.
    """
    if overwrite:
        return False

    if not config.config_retrieve(['workflow', 'check_intermediates'], True):
        return False

    if not path:
        return False

    paths = path if isinstance(path, list) else [path]
    if not all(exists(fp) for fp in paths):
        return False

    logging.debug(f'Reusing existing {path}')
    return True
