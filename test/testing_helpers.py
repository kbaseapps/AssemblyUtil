import traceback

def assert_exception_correct(got: Exception, expected: Exception):
    err = "".join(traceback.TracebackException.from_exception(got).format())
    assert got.args == expected.args, err
    assert type(got) == type(expected)