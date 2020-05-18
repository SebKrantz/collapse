
.onLoad <- function(libname, pkgname) {

  res <- .Call(C_collapse_init, "init.success")
  if(!is.character(res) || res != "init.success") stop("collapse not succesfully loaded!")

  options("collapse_unused_arg_action" = "warning") # error, warning, message or none

  invisible(res)
}

.onAttach <- function(libname, pkgname) {
  packageStartupMessage(paste0("collapse ",packageVersion("collapse"),", see ?`collapse-package` or ?`collapse-documentation`"))
}

.onUnload <- function (libpath) {
  library.dynam.unload("collapse", libpath)
}

release_questions <- function() {
  c(
    "Have you updated the version number in DESCRIPTION, NEWS.md, NEWS.Rd, cran.comments and .onAttach?",
    "Updated Readme?",
    "Spell check ?",
    "built vignettes properly with Sys.setenv(NCRAN = TRUE)?",
    "Have you updated all help files with code changes, even if it's only documenting arguments or links?",
    "All function in global_macros.R?"
  )
}
