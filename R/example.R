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
#' example("regression", "jags") # get example for the regression jags
#' example("nolapse", "jags") # get example for the nolapse stan
#' example("nolapse", "stan", "exercises/simulation_study.Rmd") # save into exercises directory
#' }
example <- function(name = c("blocked",
                             "CPP_sim",
                             "CPP",
                             "nolapse",
                             "recovery",
                             "regression",
                             "simple"), type = c(NULL, "jags", "stan"), filename = NULL) {
  if (is.null(type)) {
    fname <- sprintf("templates/%s-example.R", match.arg(name))
  } else {
    fname <- sprintf("templates/%s-example-%s.R", match.arg(name), match.arg(type))
  }

  f <- system.file(fname, package = "Rhddmjags")

  if (f == "") stop("Example ", name, " doesn't exist")

  if (is.null(filename)) {
    filename <- gsub("^templates/", "", fname)
  }

  file.copy(f, filename)
  utils::browseURL(filename)
}
