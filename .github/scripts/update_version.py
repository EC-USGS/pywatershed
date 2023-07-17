import argparse
import textwrap
from datetime import datetime
from pathlib import Path

from filelock import FileLock
from packaging.version import Version

_project_name = "pywatershed"
_project_root_path = Path(__file__).parent.parent.parent
_version_txt_path = _project_root_path / "version.txt"
_version_py_path = _project_root_path / _project_name / "version.py"
_initial_version = Version("0.0.1")
_current_version = Version(_version_txt_path.read_text().strip())


def update_version_txt(version: Version):
    with open(_version_txt_path, "w") as f:
        f.write(str(version))
    print(f"Updated {_version_txt_path} to version {version}")


def update_version_py(timestamp: datetime, version: Version):
    with open(_version_py_path, "w") as f:
        f.write(
            f"# {_project_name} version file automatically created using "
            f"{Path(__file__).name} on {timestamp:%B %d, %Y %H:%M:%S}\n\n"
        )
        f.write(f'__version__ = "{version}"\n')
        f.close()
    print(f"Updated {_version_py_path} to version {version}")


def update_version(
    timestamp: datetime = datetime.now(),
    version: Version = None,
):
    lock_path = Path(_version_txt_path.name + ".lock")
    try:
        lock = FileLock(lock_path)
        previous = Version(_version_txt_path.read_text().strip())
        version = (
            version
            if version
            else Version(previous.major, previous.minor, previous.micro)
        )

        with lock:
            update_version_txt(version)
            update_version_py(timestamp, version)
    finally:
        try:
            lock_path.unlink()
        except:
            pass


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        prog=f"Update {_project_name} version",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=textwrap.dedent(
            """\
            Update version information stored in version.txt in the project root,
            as well as several other files in the repository. If --version is not
            provided, the version number will not be changed. A file lock is held
            to synchronize file access. The version tag must comply with standard
            '<major>.<minor>.<patch>' format conventions for semantic versioning.
            """
        ),
    )
    parser.add_argument(
        "-v",
        "--version",
        required=False,
        help="Specify the release version",
    )
    parser.add_argument(
        "-g",
        "--get",
        required=False,
        action="store_true",
        help="Just get the current version number, don't update anything (defaults to false)",
    )
    args = parser.parse_args()

    if args.get:
        print(Version((_project_root_path / "version.txt").read_text().strip()))
    else:
        update_version(
            timestamp=datetime.now(),
            version=Version(args.version) if args.version else _current_version,
        )
