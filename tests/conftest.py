import pytest


def pytest_addoption(parser):
    parser.addoption(
        "--ebtel_idl_path", action="store", default=None, help="Path to EBTEL IDL code"
    )


@pytest.fixture
def ebtel_idl_path(request):
    return request.config.getoption("--ebtel_idl_path")
