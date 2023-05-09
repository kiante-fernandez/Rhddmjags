#' Get an example
#'
#' @param name The name of the example
#' @param type The name of sampler template
#' @param filename What filename you want to save (defaults to the name of the exercise in the working directory)
#'
#' @return Saves a file to the working directory (or path from filename)
#' @export
#'
#' @examples
#' \dontrun{
#' example("regression", "JAGS") # get example for the regression jags
#' example("nolapse", "JAGS") # get example for the nolapse stan
#' example("nolapse", "STAN", "exercises/simulation_study.Rmd") # save into exercises directory
#' }
example <- function(name = c("regression", "nolapse", "simple", "CCP", "blocked"), type = c("JAGS", "STAN"), filename = NULL) {
  fname <- sprintf("stubs/%s-stub.Rmd", match.arg(name))
  f <- system.file(fname, package = "Rhddmjags")

  if (f == "") stop("Example ", name, " doesn't exist")

  if (is.null(filename)) {
    filename <- gsub("^stubs/", "", fname)
  }

  file.copy(f, filename)
  utils::browseURL(filename)
}
