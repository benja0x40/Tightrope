# COMMON #######################################################################

# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Provide default values to unspecified arguments.
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
    if(is.null(to[[a]]) & ! (is.null(from[[a]]) | identical(from, to))) {
      to[[a]] <- from[[a]]
    }
    if(is.null(to[[a]]) & ! is.null(default[[a]])) {
      to[[a]] <- default[[a]]
    }
  }
}

# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Standardize the length of vector arguments.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
VectorArgs <- function(lst, from = NULL, size = NULL) {

  if(is.null(from)) from <- parent.frame()
  if(is.null(size)) {
    size <- 0
    for(x in lst) size <- max(size, length(from[[x]]))
  }

  for(x in lst) from[[x]] <- rep(from[[x]], length.out = size)

  if(! is.environment(from)) from
}

# =============================================================================.
#' ** RESERVED FOR INTERNAL USE **
# -----------------------------------------------------------------------------.
#' @description
#' Standardize the value of clonal arguments.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
ClonalArg <- function(u, a, d) { # user value, arg names, default value

  n <- length(a)
  r <- rep(list(d), n)
  names(r) <- a

  if(is.null(names(u))) {
    d[] <- rep(u, length.out = length(d))
    r[] <- rep(list(d), n)
  } else {
    u <- lapply(u, rep, length.out = length(d))
    for(k in names(u)) r[[k]][] <- u[[k]]
  }

  r
}

# Tightrope ####################################################################
