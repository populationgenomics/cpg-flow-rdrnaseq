"""
suggested location for any utility methods or constants used across multiple stages
"""

from datetime import datetime

DATE_STRING: str = datetime.now().strftime('%y-%m')  # noqa: DTZ005
