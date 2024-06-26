# Contributing

We welcome community development! We ask  contributors to follow a few guidelines, outlined below.

<!-- START doctoc generated TOC please keep comment here to allow auto update -->
<!-- DON'T EDIT THIS SECTION, INSTEAD RE-RUN doctoc TO UPDATE -->


- [Found a bug?](#found-a-bug)
- [Missing a feature?](#missing-a-feature)
- [Submission guidelines](#submission-guidelines)
  - [Issues](#issues)
  - [Pull requests](#pull-requests)

<!-- END doctoc generated TOC please keep comment here to allow auto update -->

## Found a bug?

If you find a bug in the source code, you can help us by [submitting an issue](https://github.com/modflowpy/flopy/blob/develop/CONTRIBUTING.md#submit-issue). Even better, you can [submit a Pull Request with a fix](https://github.com/modflowpy/flopy/blob/develop/CONTRIBUTING.md#submit-pr).

## Missing a feature?

You can request a new feature by [submitting an issue](https://github.com/modflowpy/flopy/blob/develop/CONTRIBUTING.md#submit-issue). If you would like to implement a new feature, please submit an issue with a proposal for your work first, to be sure that we can use it. Please consider what kind of change it is:

- For a major feature, first open an issue and outline your proposal so that it can be discussed. This will allow us to better coordinate our efforts, prevent duplication of work, and help you craft the change so that it is successfully accepted into the project.
- Small features can be crafted and directly submitted as a Pull Request.

## Submission guidelines

### Issues

Before you submit an issue, please search the issue tracker, maybe an issue for your problem already exists and the discussion might inform you of workarounds readily available.

We want to fix all the issues as soon as possible, but before fixing a bug we need to reproduce and confirm it. In order to reproduce bugs, we will systematically ask you to provide a minimal, complete, and verifiable example. Having a minimal, complete, and verifiable example gives us a wealth of important information without going back & forth to you with additional questions like:

- version of `pywatershed` used
- and most importantly &mdash; a use-case that fails

We will be insisting on a minimal, complete, and verifiable example in order to save maintainers time and ultimately be able to fix more bugs. We understand that sometimes it might be hard to extract essentials bits of code from a larger code-base but we really need to isolate the problem before we can fix it.

Unfortunately, we are not able to investigate / fix bugs without a minimal, complete, and verifiable example, so if we don't hear back from you we are going to close an issue that doesn't have enough info to be reproduced.

### Pull requests

To submit a pull request (PR) please follow these steps:

1. Search GitHub for an open or closed PR that relates to your submission. You don't want to duplicate effort.

2. Fork the repository.

3. Make your changes in a new branch. Be sure to include test cases.

4. Run `isort` and `black` on the `pywatershed` module.

    Python source code files should be formatted by running `black` and `isort` on the
    `pywatershed` subdirectory before opening a pull request, as CI will fail if the code
    is not properly formatted. For instance, from the project root:

    ```shell
    black --check --diff pywatershed
    isort --check --diff pywatershed
    ```

    **Note**:  PRs must pass format checks in GitHub Actions before they can be accepted.

5. Run the full test suite and make sure all tests pass.

    **Note**: PRs must pass all tests in GitHub Actions before they can be accepted.

6. Commit your changes with a descriptive message that follows [conventional commit](https://github.com/modflowpy/flopy/blob/develop/CONTRIBUTING.md#commit) guidelines.

7. Push your branch to your fork.

8. Open a pull request against the upstream repository's `develop` branch.

9. If we suggest changes:

    a. Make the required updates.
    b. Rerun the tests, ensuring they still pass.
    c. If the upstream `develop` has moved on since you opened the PR, rebase your branch and force push to your PR branch. For instance, assuming you have updated your fork's `develop` from upstream on GitHub, and your fork's remote is called `origin` as usual:
    
        ```
        git fetch origin
        git rebase develop -i
        git push -f origin <your PR branch>
        ```

That's it! Thanks for your contribution.
