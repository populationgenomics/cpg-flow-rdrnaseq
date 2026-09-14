from dataclasses import dataclass
from typing import Any

from cpg_utils import config


@dataclass(frozen=True)
class QcFlag:
    flag: str
    value: float
    comparison: str
    threshold: float
    section: str
    date: str
    ar_guid: str
    severity: str = 'fail'
    resolved: bool = False
    resolution_date: str | None = None


DIRECTIONS: dict[str, tuple[str, Any]] = {
    'under': ('<', lambda val, thresh: val < thresh),
    'over': ('>', lambda val, thresh: val > thresh),
}
SEVERITIES: tuple[str, ...] = ('fail', 'warn')


def load_thresholds(seq_type: str) -> dict[str, dict[str, dict[str, float]]]:
    """Read the nested qc_thresholds config into {direction: {metric: {severity: threshold}}}."""
    thresholds: dict[str, dict[str, dict[str, float]]] = {d: {} for d in DIRECTIONS}
    for severity in SEVERITIES:
        for direction in DIRECTIONS:
            configured = config.config_retrieve(['qc_thresholds', seq_type, severity, direction], {})
            for metric, threshold in configured.items():
                thresholds[direction].setdefault(metric, {})[severity] = threshold
    return thresholds


def worst_breach(val: float, tiers: dict[str, float], direction: str) -> tuple[str, float] | None:
    """Return (severity, threshold) of the most severe tier breached, else None."""
    _, breaches = DIRECTIONS[direction]
    for severity in SEVERITIES:
        if severity in tiers and breaches(val, tiers[severity]):
            return severity, tiers[severity]
    return None
