# Scripts expect the repository root as the current working directory.
if (!dir.exists("code") || !dir.exists("data")) {
  stop("Open the project root containing code/ and data/. See CODE_SETUP.md.")
}
dir.create("results", showWarnings = FALSE)
