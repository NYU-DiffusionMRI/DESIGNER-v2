def pytest_addoption(parser):
    parser.addoption(
        "--no-cleanup",
        action="store_true",
        default=False,
        help="Skip cleanup of DESIGNER's scratch (processing) directory and TMI's params directory.",
    )
