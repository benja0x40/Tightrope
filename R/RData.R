# =============================================================================.
#' Save an R object as RData with semi-automated path and file name
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadObj},
#' \link{AvailableObj},
#' \link{saveRDS},
#' \link{save.image},
#' \link{reassign}
# -----------------------------------------------------------------------------.
#' @param obj
#' an R object to be saved as RData.
#'
#' @param path
#' folder location used to save the R object.
#' When none is specified it is the current working directory.
#'
#' @param name
#' a file name. When omitted (default) this name is automatically set to the
#' name of the R object being passed as argument.
#'
#' @param ...
#' optional parameters passed to the \link{saveRDS} function.
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @examples
#' # R objects are automatically saved in the current working directory
#' a <- rep("a", 10)
#' SaveObj(a)
#' rm(a)
#' LoadObj(a)
#' print(a)
#' file.remove("a.rdata") # delete the file generated for this example
#'
#' # R objects can be saved and loaded from any path
#' b <- rep("b", 10)
#' SaveObj(b, path = "MyRData") # path is automatically created if necessary
#' b <- NULL
#' LoadObj(b, path = "MyRData")
#'
#' # Not yet defined R objects may be created only when no saved version exists
#' if(! LoadObj(a, "MyRData")) {
#'   a <- 0
#'   SaveObj(a, "MyRData")
#' }
#' print(a)
#'
#' # Saved object can be updated by overwriting
#' a <- 1
#' SaveObj(a, path = "MyRData")
#' rm(a)
#' LoadObj(a)
#' print(a)
# -----------------------------------------------------------------------------.
#' @export
SaveObj <- function(obj, path = NULL, name = NULL, ...) {
  obj.name <- name
  if(is.null(obj.name)) obj.name <- deparse(substitute(obj))
  fpath <- MakePath(path, obj.name, ext = ".rdata")
  path  <- dirname(fpath)
  if(path!="" & ! file.exists(path)) {
    message("[creating] ", path)
    system(paste("mkdir -p", path))
  }
  if(! file.exists(fpath)) {
    message("[saving] ", fpath)
  } else {
    message("[updating] ", fpath)
  }
  saveRDS(obj, fpath, ...)
}

# =============================================================================.
#' Load an R object from RData with semi-automated path and file name
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{SaveObj},
#' \link{AvailableObj},
#' \link{readRDS},
#' \link{reassign}
# -----------------------------------------------------------------------------.
#' @param obj
#' an R object previously saved using \link{SaveObj}.
#'
#' @param path
#' folder location used to load the R object.
#' When none is specified it is the current working directory.
#'
#' @param name
#' a file name. When omitted (default) this name is automatically set to the
#' name of the R object being passed as argument.
#'
#' @param need
#' logical, when \code{TRUE} LoadObj raises an error if object loading fails
#' (default = F).
#'
#' @param pos
#' the target \link{environment} where the object should be loaded
#' (default = .GlobalEnv).
#' See \link{assign} for documentation on the different ways to specify
#' environments with the \code{pos} paramter.
#'
#' @param overload
#' logical. When the R object to be loaded is already defined in the target
#' environment LoadObj can avoid (overload = F) or force (overload = T) the
#' reloading of this object.
#'
#' @param ...
#' optional parameters passed to the \link{readRDS} function.
#'
#' @return
#' logical, \code{TRUE} when object was loaded successfully,
#' \code{FALSE} otherwise.
# -----------------------------------------------------------------------------.
#' @examples
#' # R objects are automatically saved in the current working directory
#' a <- rep("a", 10)
#' SaveObj(a)
#' rm(a)
#' LoadObj(a)
#' print(a)
#' file.remove("a.rdata") # delete the file generated for this example
#'
#' # R objects can be saved and loaded from any path
#' b <- rep("b", 10)
#' SaveObj(b, path = "MyRData") # path is automatically created if necessary
#' b <- NULL
#' LoadObj(b, path = "MyRData")
#'
#' # Not yet defined R objects may be created only when no saved version exists
#' if(! LoadObj(a, "MyRData")) {
#'   a <- 0
#'   SaveObj(a, "MyRData")
#' }
#' print(a)
#'
#' # Saved object can be updated by overwriting
#' a <- 1
#' SaveObj(a, path = "MyRData")
#' rm(a)
#' LoadObj(a)
#' print(a)
# -----------------------------------------------------------------------------.
#' @export
LoadObj <- function(
  obj, path = NULL, name = NULL, need = F, pos = .GlobalEnv, overload = F, ...
) {
  need <- need # Force evaluation of need to avoid overloading issues
  obj.name <- name
  if(is.null(obj.name) & ! missing(obj)) {
    obj.name <- deparse(substitute(obj))
  }
  if(is.null(obj.name)) {
    stop("either obj or name must be specified")
  }
  fpath <- MakePath(path, obj.name, ext = ".rdata")
  chk <- exists(x = obj.name, where = pos)
  msg <- ifelse(chk, "[overloading] ", "[loading] ")
  res <- NULL
  if(file.exists(fpath)) {
    if(chk & ! overload) {
      message("[skipped] ", fpath)
      res <- T
    } else {
      message(msg, fpath)
      res <- assign(obj.name,readRDS(fpath, ...), pos = pos)
    }
  }
  if(need & is.null(res)) {
    stop("failed to load ", fpath)
  }
  ! is.null(res)
}

# HIDDEN #######################################################################

# =============================================================================.
#' Concatenate several strings to form a filesystem path
# -----------------------------------------------------------------------------.
#' @param ... character strings forming a file path
#' @param ext file name extension (default = none)
#'
#' @return
#' \code{MakePath} returns a \code{character} value
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
MakePath <- function(..., ext = "") {
  clean_slashes <-function(x) {
    x <- gsub("[/]+", "/", x)
    x <- gsub("^/", "", x)
    x <- gsub("/$", "", x)
    if(x != "") x <- paste0(x, "/")
    x
  }
  path <- list(...)
  path <- path[! sapply(path, is.null)]
  if(length(path) > 0) {
    path <- lapply(path, as.character)
    root <- ""
    if(substr(path[1], 1, 1) == "/") root <- "/"
    n <- length(path)
    path[-n] <- sapply(path[-n], clean_slashes)
    path <- paste0(root, paste(path, collapse = ""), ext)
    path <- gsub("[/]+", "/", path)
    # path <- gsub("[\\.]+([^\\.]+)$", ".\\1", path)
  } else {
    path <- NULL
  }
  path
}

# =============================================================================.
#' Check if an R object is available as RData
# -----------------------------------------------------------------------------.
#' @seealso
#'   \link{LoadObj},
#'   \link{SaveObj}
# -----------------------------------------------------------------------------.
#' @inheritParams LoadObj
#'
#' @return
#' logical, \code{TRUE} when RData is available and \code{FALSE} otherwise.
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
AvailableObj <- function(obj, path = NULL, name = NULL) {
  obj.name <- name
  if(is.null(obj.name) & ! missing(obj)) {
    obj.name <- deparse(substitute(obj))
  }
  if(is.null(obj.name)) {
    stop("either obj or name must be specified")
  }
  fpath <- MakePath(path, obj.name, ext = ".rdata")

  file.exists(fpath)
}

# =============================================================================.
#' Reassign object to a different environment
# -----------------------------------------------------------------------------.
#' @seealso
#' \link{LoadObj},
#' \link{SaveObj},
#' \link{AvailableObj},
#' \link{assign},
#' \link{environment}
# -----------------------------------------------------------------------------.
#' @param obj
#' an R object
#'
#' @param pos
#' the target \link{environment} where the object should be reassigned.
#' See \link{assign} for documentation on the different ways to specify
#' environments with the \code{pos} paramter.
#'
#' @param src
#' the source \link{environment} where the object is currently located
#' (default = .GlobalEnv).
#' See \link{assign} for documentation on the different ways to specify
#' environments with the \code{pos} paramter.
#'
#' @param keep
#' logical indicating if the original object should be removed or preserved
#' at its current location (defautl = F, remove).
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
reassign <- function(obj, pos, src = .GlobalEnv, keep = F) {
  obj.name <- deparse(substitute(obj))
  assign(obj.name, obj, pos = pos)
  if(! keep) {
    rm(list = obj.name, pos = src)
  }
}

# NOT EXPORTED #################################################################

# =============================================================================.
#' Convert RData file to tab delimited file
# -----------------------------------------------------------------------------.
#' @param path
#' character
#'
#' @param name
#' character
#'
#' @param col.names
#' logical
#'
#' @param row.names
#' logical
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
rdata2tsv <- function(path, name, col.names = T, row.names = F) {
  rdp <- paste(path, name, ".rdata", sep = "")
  tdp <- paste(path, name, ".txt", sep = "")
  x <- readRDS(rdp)
  write.table(
    x, file = tdp,
    quote = F, sep = "\t", col.names = col.names, row.names = row.names
  )
}

# =============================================================================.
#' Convert tab delimited file to RData
# -----------------------------------------------------------------------------.
#' @param path
#' character
#'
#' @param name
#' character
#'
#' @return NULL
# -----------------------------------------------------------------------------.
#' @keywords internal
#' @export
tsv2rdata <- function(path, name) {
  rdp <- paste(path, name, ".rdata", sep = "")
  tdp <- paste(path, name, ".txt", sep = "")
  x <- read.delim(tdp, stringsAsFactors = F)
  saveRDS(x, rdp)
}
