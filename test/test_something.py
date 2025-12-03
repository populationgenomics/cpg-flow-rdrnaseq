"""
add relevant tests to this directory
"""

import pytest


def test_something():
    """
    add a test here
    """
    with pytest.raises(ZeroDivisionError):
        _ = 1 / 0
