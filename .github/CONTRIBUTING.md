# Contributing to LazyModeler

Thanks for your interest in contributing to LazyModeler! We appreciate bug
reports, feature ideas, documentation improvements, tests, examples, and code
contributions.

You do not need to make a huge contribution for it to matter. Small fixes,
clear bug reports, better examples, and improved documentation are all very
welcome.

## Table of contents  
  
- [Code of Conduct](#code-of-conduct)  
- [Ways to contribute](#ways-to-contribute)  
- [Reporting bugs](#reporting-bugs)  
- [Suggesting features](#suggesting-features)  
- [Development workflow](#development-workflow)
	- [1. Fork and clone the repository](#1-fork-and-clone-the-repository)  
	- [2. Create a branch](#2-create-a-branch)  
	- [3. Install development dependencies](#3-install-development-dependencies)  
	- [4. Make your changes](#4-make-your-changes)  
	- [5. Add or update tests](#5-add-or-update-tests)  
	- [6. Check the package](#6-check-the-package)  
	- [7. Check style](#7-check-style)  
	- [8. Commit your changes](#8-commit-your-changes)  
	- [9. Keep your branch up to date](#9-keep-your-branch-up-to-date)  
	- [10. Open a pull request](#10-open-a-pull-request)  
		- [Pull request checklist](#pull-request-checklist)  
- [Documentation changes](#documentation-changes)  
- [Tests and examples](#tests-and-examples)  
- [Versioning and releases](#versioning-and-releases)  
- [Questions](#questions)

## Code of Conduct

This project follows a [`CODE_OF_CONDUCT`](CODE_OF_CONDUCT.md). By participating in this project, you
agree to follow it.

If you need to report unacceptable behaviour, contact:

- `Lara Koesters`: `lkoesters@bgc-jena.mpg.de`

## Ways to contribute

You can help by:

- reporting bugs or unexpected behaviour
- suggesting new features or improvements
- improving documentation, examples, or tutorials
- adding or improving tests
- fixing bugs
- improving code quality, performance, or maintainability

If you are unsure whether something would be useful, feel free to open an issue
first and describe your idea.

## Reporting bugs

When encountering a bug, please open a new issue [`here`](https://github.com/LMKoesters/LazyModeler/issues). Before opening an issue, please check whether you are using the latest version of LazyModeler and whether the problem has already been reported.

Bug reports should include as much of the following as possible:

- a short, clear description of the problem
- a minimal reproducible example
- the expected behaviour
- the actual behaviour
- any error messages or warnings
- your R version
- your LazyModeler version
- your operating system

A good bug report makes it much easier to understand and fix the problem. If the issue involves model output, plots, or data, please try to use a small example dataset whenever possible. Avoid sharing sensitive or private data.
## Suggesting features

Feature requests are welcome. Please make sure the feature has not been proposed before and include:

- what problem the feature would solve
- how you would expect the feature to work
- a short example of the desired behaviour
- whether you would be interested in helping implement it

Not every feature request can be accepted, but clear use cases are very helpful
for deciding what belongs in the package.
## Development workflow

> ### Legal Notice
> When contributing to this project, you must agree that you have authored 100% of the content, that you have the necessary rights to the content and that the content you contribute may be provided under the project licence.

This section describes the recommended workflow for contributing code.
### 1. Fork and clone the repository

Fork the repository on GitHub, then clone your fork:

```bash
git clone https://github.com/LMKoesters/LazyModeler.git
cd LazyModeler
```

Add the original repository as upstream:

```bash
git remote add upstream https://github.com/LMKoesters/LazyModeler.git
```
### 2. Create a branch

Create a new branch for your change:

```bash
git checkout -b my-short-description
```

Use a short, descriptive branch name.
### 3. Install development dependencies

Open R in the package directory and install the package dependencies:

```r
install.packages("devtools")
devtools::install_dev_deps()
```

Then load the package locally:

```r
devtools::load_all()
```
### 4. Make your changes

Please keep pull requests focused. A pull request that fixes one bug or adds one
well-scoped feature is much easier to review than a large mixed change.

If you change user-facing behaviour, please update the relevant documentation
and examples.

If you add or change exported functions, update the roxygen comments and run:
```r
devtools::document()
```
### 5. Add or update tests

If your change fixes a bug, please add a test that would have failed before the
fix.

If your change adds new functionality, please add tests for the main expected
behaviour and relevant edge cases.

Run the tests with:

```r
devtools::test()
```

For a single test file, you can use:

```r
testthat::test_file("tests/testthat/test-file-name.R")
```
### 6. Check the package

Before opening a pull request, please run:

```r
devtools::check()
```

If your change affects test coverage, you can also check coverage with:

```r
covr::package_coverage()
```
### 7. Check style

Please aim for readable, consistent R code.

You can run:
```r
lintr::lint_package()
```

Some lint messages may be stylistic rather than blocking, but please address
issues that indicate possible bugs, unclear code, unused objects, or namespace
problems.

### 8. Commit your changes

Commit your changes with a clear message:

```bash
git add -A
git commit -m "Your short and descriptive commit message"
```
### 9. Keep your branch up to date

If your branch gets out of date, update it from the main repository:

```bash
git fetch upstream
git checkout main
git merge upstream/main
git checkout my-short-description
git merge main
```

If you are comfortable with rebasing, you may use git rebase instead. If not,
merging is totally fine.

### 10. Open a pull request

Push your branch:

```bash
git push origin my-short-description
```

Then open a pull request on GitHub.

In your pull request, please include:

- what changed
- why the change is needed
- how you tested it
- any known limitations or follow-up work

If your pull request fixes an issue, mention it in the description.
#### Pull request checklist

Before submitting a pull request, please check:

 - [ ] The change is focused and reasonably small.
 - [ ] The code runs locally.
 - [ ] Relevant tests were added or updated.
 - [ ] devtools::test() passes.
 - [ ] devtools::check() passes, or remaining notes are explained.
 - [ ] Documentation was updated if needed.
 - [ ] The pull request description explains the motivation and the test results.
## Documentation changes

Documentation improvements are very welcome.

This includes:

- fixing typos
- clarifying examples
- improving function documentation
- adding explanations to the README
- improving vignettes or tutorials

For roxygen documentation changes, run:

```r
devtools::document()
```
## Tests and examples

Please keep tests reasonably small and deterministic.

When adding tests:

- use small synthetic datasets where possible
- set seeds for random data generation
- avoid depending on local files, local paths, or interactive state
- avoid long-running tests unless they are clearly necessary
## Versioning and releases

Maintainers handle version bumps and releases. Contributors generally do not need to update the package version unless asked.
## Questions

If you are stuck or unsure how to approach a contribution, feel free to open an
issue and ask. It is completely fine to discuss an idea before writing code.

Thanks again for helping improve LazyModeler!