# Release guide
This document describes release procedures, conventions, and utilities for
`pywatershed`.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Table of Contents**

  - [Conventions](#conventions)
  - [Releasing `pywatershed`](#releasing-pywatershed)
  - [Utility scripts](#utility-scripts)
    - [Updating version numbers](#updating-version-numbers)
- [A '+' is appended to the version number in `version.txt` (if it was not already
there) to indicate that the repository is in a preliminary/development state. To
omit the '+' for release-ready versions, use the `--approve` (short `-a`)
flag. To set the version number, use the `--version` (short `-v`) option.](#a--is-appended-to-the-version-number-in-versiontxt-if-it-was-not-already%0Athere-to-indicate-that-the-repository-is-in-a-preliminarydevelopment-state-to%0Aomit-the--for-release-ready-versions-use-the---approve-short--a%0Aflag-to-set-the-version-number-use-the---version-short--v-option)
- [<<<<<<< HEAD
To get the current version number without writing any changes to the
repository's files, use the `--get` (short `-g`) flag:](#-head%0Ato-get-the-current-version-number-without-writing-any-changes-to-the%0Arepositorys-files-use-the---get-short--g-flag)
    - [Preparing for PRs](#preparing-for-prs)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->



## Conventions

- Releases follow the [git
  flow](https://nvie.com/posts/a-successful-git-branching-model/).
- Release numbers follow [semantic version](https://semver.org/) conventions.
- Minor and major releases branch from `develop`. Patches branch from `main`.

## Releasing `pywatershed`

The release procedure is mostly automated. The workflow is defined in
`.github/workflows/release.yaml` and triggers when a release or patch branch is
pushed to this repo.

Prior to release:
1. Run asv tests without asv, perform an asv regression test:
1. On develop update the `what's new.rst` to include the date of the release


To release a new version:

0. Perform ASV performance benchmarks against previous release on Windows,
   Linux, and Mac. Gather results.

1. On your local machine, create a release branch from `develop` or a patch
   branch from `main`.  The branch's name must follow format
   `v{major}.{minor}.{patch}` ([semantic version](https://semver.org/) number
   with a leading 'v'). For instance, for a minor release, if this repo is an
   `upstream` remote and one's local `develop` is up to date with upstream
   `develop`, then from `develop` run `git switch -c vx.y.z`.

2. If this is a patch release, make changes/fixes locally. If this is a major or
   minor release, no changes are needed.

3. Push the branch to this repo. For instance, if this repo is an `upstream`
   remote: `git push -u upstream vx.y.z`. This starts a job to:

    - Check out the release branch
    - Update version number in `version.txt` and `pywatershed/version.py` to
      match the version in the branch name
    - Build and check the Python package
    - Generate a changelog since the last release
    - Prepend the changelog to the cumulative `HISTORY.md`
    - Upload the package and changelog as artifacts
    - Draft a PR against `main` with the updated version files and cumulative
      changelog. The cumulative `HISTORY.md` is version-controlled, release
      changelogs are not.

3. Inspect the package and changelog. If they look good, merge the PR to `main`.

    **Note**: it is critical to *merge* the PR to `main`, not squash as is
    conventional for development PRs. Squashing causes `develop` and `main` to
    diverge. Merging to `main` preserves commit history and ensures `develop`
    and `main` don't diverge.

    Merging the PR to `main` will trigger another job to draft a [GitHub
    release](https://github.com/EC-USGS/pywatershed/releases). The release is
    not yet publicly visible at this point. The release notes are autofilled as
    the changelog since the last release.

4. Manually add the following to the Github release.
    - Frozen environments for all platforms.
    - Attach ASV results for all platforms.

   Inspect the GitHub release. If needed, make any manual edits to the release
   notes. If the release looks good, publish it via GitHub UI or CLI. This tags
   the head of `main` with the release version number (**Note**: release tags,
   unlike branches, don't include an initial `v`, as is common in some projects)
   and triggers jobs to:

    - Publish the package to PyPI
    - Check out `main`
    - Run `.github/scripts/update_version.py -v x.y+1.0.dev0` to update `version.txt` and `pywatershed/version.py` with the minor version number incremented. The `.dev0` suffix indicates preliminary development status.
    - Draft a PR against `develop` with the updated version files and the updates previously merged to `main`.

5. Merge the PR to `develop`. As above, it is important to *merge* the PR, not
   squash, to preserve history and keep `develop` and `main` from diverging.

6. Manually update `main` image on WholeTale to have the current release.

## Utility scripts

The automated release procedure uses a few scripts, located in
`.github/scripts`:

- `update_version.py`
- `pull_request_prepare.py`

The former should never need to be run manually. The latter is convenient for
formatting source files before opening PRs.

### Updating version numbers

<<<<<<< HEAD
The `update_version.py` script can be used to update version numbers embedded in
the repository. The script acquires a file lock to make sure only one process
edits version files at a given time.

If the script is run with no arguments, updated timestamp comments are written
but the version number is not changed.

A '+' is appended to the version number in `version.txt` (if it was not already
there) to indicate that the repository is in a preliminary/development state. To
omit the '+' for release-ready versions, use the `--approve` (short `-a`)
flag. To set the version number, use the `--version` (short `-v`) option.
=======
The `update_version.py` script can be used to update version numbers embedded in the repository. The script acquires a file lock to make sure only one process edits version files at a given time. If the script is run with no arguments, updated timestamp comments are written but the version number is not changed. To set the version number, use the `--version` (short `-v`) option.
>>>>>>> upstream/develop

For instance, to set the version number before a release:

```shell
python .github/scripts/update_version.py -a -v 0.1.3
```

<<<<<<< HEAD
To get the current version number without writing any changes to the
repository's files, use the `--get` (short `-g`) flag:
=======
Or to set the version number on `develop` following a release:

```shell
python .github/scripts/update_version.py -a -v 0.2.0.dev0
```

To get the current version number without writing any changes to the repository's files, use the `--get` (short `-g`) flag:
>>>>>>> upstream/develop

```shell
python .github/scripts/update_version.py -g
```

### Preparing for PRs

The `pull_request_prepare.py` script lints Python source code files by running
`black` and `isort` on the `pywatershed` subdirectory. This script should be run
before opening a pull request, as CI will fail if the code is not properly
formatted. For instance, from the project root:

```shell
python .github/scripts/pull_request_prepare.py
```
