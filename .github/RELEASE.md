# Release guide This document describes release procedures, conventions, and
utilities for `pywatershed`.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->
**Contents**

- [Conventions](#conventions)
- [Releasing `pywatershed`](#releasing-pywatershed)
- [Utility scripts](#utility-scripts)
  - [Updating version numbers](#updating-version-numbers)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Conventions

- Releases follow the [git
  flow](https://nvie.com/posts/a-successful-git-branching-model/).  - Release
  numbers follow [semantic version](https://semver.org/) conventions.  - Minor
  and major releases branch from `develop`. Patches branch from `main`.

## Releasing `pywatershed`

The release procedure is mostly automated. The workflow is defined in
`.github/workflows/release.yaml` and triggers when a release or patch branch is
pushed to this repo.

To release a new version:

1. Test asv benchmarking with the `-q` flag to ensure it is working (multiple
   platformas a bonus). You can merge the ref you are testing, but you cant
   test the release ref on main until after the release is made, so that step
   will come below and those details will merge to develop and not the
   release.

1. If there are extended release notes to put in the gh-pages branch, once
   you have approval, change the dates of the appropriate pages to the
   desired release date. Change any links to the extended documentation in
   develop branch to match. Push gh-pages branch to confirm in develop before
   releasing.

1. On your local machine, create a release branch from `develop` or a patch
   branch from `main`.  The branch's name must follow format
   `v{major}.{minor}.{patch}` ([semantic version](https://semver.org/) number
   with a leading 'v'). For instance, for a minor release, if this repo is an
   `upstream` remote and one's local `develop` is up to date with upstream
   `develop`, then from `develop` run `git switch -c vx.y.z`. 

1. If this is a patch release, make changes/fixes locally. If this is a major or
   minor release, no changes may be needed. In either case, add the release version
   and date to the top of `doc/whats-new.rst`. If a patch, put it below the
   pending minor release. Also update the version in `doc/index.rst`. Update the
   CITATION.cff file. If a major release, get the provisional new DOI from USGS
   and add it to CITATION.cff and the main README.md.

1. Push the branch to this repo. For instance, if this repo is an `upstream`
   remote: `git push -u upstream vx.y.z`. This starts a job to:
    - Check out the release branch and update version number in `version.txt` and
      `pywatershed/version.py` to match the version in the branch name.
    - Build and check the Python package.
    - Upload the package as an artifact.
    - Draft a PR against `main` with the updated version files.

1. On all platforms, pull the release from upstream and perform ASV performance
   benchmarks against previous release , e.g., ```asv run --verbose
   -show-stderr HASHFILE:pws_refs_for_asv.txt``` after editing the file to
   contain the previous and current release. Collect performance results from
   various machines into a single report and use `asv publish` to generate
   the static webpages to be included with the release as artifacts in that
   step below.

1. Inspect the package. If it looks good, merge the PR to `main`.

   **Note**: it is critical to *merge* the PR to `main`, not squash as is
   conventional for development PRs. Squashing causes `develop` and `main` to
   diverge. Merging to `main` preserves commit history and ensures `develop`
   and `main` don't diverge.

   Merging the PR to `main` will trigger another job to draft a [GitHub
   release](https://github.com/EC-USGS/pywatershed/releases). The release is
   not yet publicly visible at this point.

1. Inspect the GitHub release. If needed, make any manual edits to the release
   notes. If the release looks good, publish it via GitHub UI or CLI. Manually
   add the asv static web pages and frozen conda dependencies for each platform.

   Publishing the release on GitHub automatically tags the head of `main` with
   the release version number (**Note**: release tags, unlike branches, don't
   include an initial `v`, as is common in some projects) and triggers jobs to
   publish the package to PyPI.

1. Delete the release branch on upstream.

1. On all platforms, pull the release from upstream and perform ASV performance
   benchmarks against previous release. First delete `.asv/` entirely (unless
   there's some good reason to keep `.asv/envs`), then update
   `pws_refs_for_asv.txt` to have the previous release and the above release,
   next `asv run --verbose --show-stderr HASHFILE:pws_refs_for_asv.txt`. Collect
   performance results from various machines into a single `.asv/results/`, then
   use `asv publish` to generate the static webpages in `.asv/html`. When these
   are to your liking, add them to the gh-pages branch in the `asv/vM-m-p`
   directory.

1. In the case of a minor or major release,  (a patch would be applied to both
   main and develop?) a couple of manual steps:
   - Check out `main` as a new branch to get it back into develop,
     eg feat_main_to_dev.
   - Run `.github/scripts/update_version.py -v x.y+1.0.dev0` to update
     `version.txt` and `pywatershed/version.py` with the minor version number
     incremented. The `.dev0` suffix indicates preliminary development status.
   - Add a new minor release to the top of `doc/whats-new.rst`
   - Incorporate the ASV details from the previous section above.
   - Open a PR against `develop` with the updated version files and the
     updates previously merged to `main`.

1. Merge the PR to `develop`. As above, it is important to *merge* the PR, not
   squash, to preserve history and keep `develop` and `main` from diverging.

1. Update `main` image on WholeTale to have the current release.


## Utility scripts

The automated release procedure uses a few scripts, located in
`.github/scripts`.

### Updating version numbers

The `update_version.py` script can be used to update version numbers embedded in
the repository. The script acquires a file lock to make sure only one process
edits version files at a given time. If the script is run with no arguments,
updated timestamp comments are written but the version number is not changed. To
set the version number, use the `--version` (short `-v`) option.

For instance, to set the version number before a release:

```shell python .github/scripts/update_version.py -a -v 0.1.3 ```

Or to set the version number on `develop` following a release:

```shell python .github/scripts/update_version.py -a -v 0.2.0.dev0 ```

To get the current version number without writing any changes to the
repository's files, use the `--get` (short `-g`) flag:

```shell python .github/scripts/update_version.py -g ```

**Note**: this script should not need to be run manually, as it is run automatically in the release automation.
