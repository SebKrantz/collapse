# Contributing to *collapse*

- If you found a problem or have a feature suggestion, please [open an issue](https://github.com/SebKrantz/collapse/issues/new/choose) using one of the templates.
- For broader proposals start a [discussion](https://github.com/SebKrantz/collapse/discussions).
- To contribute directly, read [How to contribute code](#how-to-contribute-code) section below.
- Contributors will be mentioned in the `DESCRIPTION` file as `"ctb"` **only if** the contribution is a substantial improvement or new functionality. 

## How to contribute code

For code contributions, follow the five simple steps below.

1. [Open an issue](https://github.com/SebKrantz/collapse/issues/new/choose) to let others know about the problem/feature you aim to address.
You may skip this step if your contribution is minor.
2. [Fork the repository](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/working-with-forks/fork-a-repo).
Make sure to uncheck "Copy the main branch only" box to include all branches.
3. [Clone the repository](https://docs.github.com/en/repositories/creating-and-managing-repositories/cloning-a-repository) locally.
4. Switch to `development` branch or create a new one off `development` and make your changes in it.
5. Once your code is ready, push it to your fork and 
[open a pull request](https://docs.github.com/en/pull-requests/collaborating-with-pull-requests/proposing-changes-to-your-work-with-pull-requests/creating-a-pull-request) 
to `development` branch in the source repository.

> [!IMPORTANT]  
> When opening a pull request, make sure to select `development` branch as the target (base) branch.

## How to develop locally

The project uses `renv` to allow you to easily install the dependencies in an isolated manner.

> [!NOTE]  
> You may need to [install or update `renv`](https://rstudio.github.io/renv/index.html#installation) package on your system before proceeding.

1. Open `collapse.Rproj`. This will automatically activate the environment.
2. Run `renv::install()` in the console to install the packages required for development and testing.
3. For interactive testing, you may want to use `devtools::load_all()` to load the complete `collapse` package.
4. Once your code is ready, run tests with `devtools::test()` and run checks with `devtools::check()`.
