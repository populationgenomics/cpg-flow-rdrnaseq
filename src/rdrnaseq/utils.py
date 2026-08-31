from dataclasses import dataclass


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
