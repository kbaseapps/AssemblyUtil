'''
Configure pytest fixtures and helper functions for this directory.
'''

import time
import traceback

def assert_exception_correct(got: Exception, expected: Exception):
    err = "".join(traceback.TracebackException.from_exception(got).format())
    assert got.args == expected.args, err
    assert type(got) == type(expected)


def assert_close_to_time(time_, seconds_slop=1):
    """
    Checks that a timestamp in seconds since the epoch is within a given amount of the current
    time.
    """
    now_ms = time.time()
    assert now_ms + seconds_slop > time_

    assert now_ms - seconds_slop < time_