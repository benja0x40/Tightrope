# =============================================================================.
#
# -----------------------------------------------------------------------------.
.onAttach <- function(...) {

  # Initialize global options of the Tightrope package
  Tightrope::ResetOptions()

}

# =============================================================================.
#
# -----------------------------------------------------------------------------.
.onDetach <- function(...) {

  # Remove global options of the Tightrope package from the R environment
  Tightrope::RemoveOptions()

}
