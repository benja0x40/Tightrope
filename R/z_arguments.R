# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Provide default values to unspecified arguments
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
DefaultArgs <- function(default, ignore = NULL, from = NULL, to = NULL) {

  lst <- names(default)

  if(is.null(from)) from <- parent.frame()
  if(is.null(to)) to <- parent.frame()

  if(is.function(from)) {
    lst <- methods::formalArgs(from)
    from <- parent.frame()
  }

  lst <- setdiff(lst, ignore)

  for(a in lst) {
    if(! (is.null(from[[a]]) | identical(from, to))) {
      to[[a]] <- from[[a]]
    }
    if(is.null(to[[a]]) & ! is.null(default[[a]])) {
      to[[a]] <- default[[a]]
    }
  }
}
