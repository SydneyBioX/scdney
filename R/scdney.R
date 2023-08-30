# Adapted from tidyverse codes

core <- c("scMerge", "scClassify", "CiteFuse", "scDC", "scClust", "Cepo", "scHOT", "scFeatures")

core_unloaded <- function() {
  search <- paste0("package:", core)
  core[!search %in% search()]
}

same_library <- function(pkg) {
  loc <- if (pkg %in% loadedNamespaces()) dirname(getNamespaceInfo(pkg, "path"))
  do.call(
    "library",
    list(pkg, lib.loc = loc, character.only = TRUE, warn.conflicts = FALSE)
  )
}

scdney_attach <- function() {
  to_load <- core_unloaded()
  if (length(to_load) == 0)
    return(invisible())
  
  msg(
    cli::rule(
      left = "Attaching packages",
      right = paste0("scdney ", package_version("scdney"))
    ),
    startup = TRUE
  )
  
  versions <- vapply(to_load, package_version, character(1))
  packages <- paste0(
    crayon::green(cli::symbol$tick), " ", 
    crayon::magenta(crayon::italic(format(to_load))), " ",
    crayon::col_align(versions, max(crayon::col_nchar(versions)))
  )
  
  if (length(packages) %% 2 == 1) {
    packages <- append(packages, "")
  }
  col1 <- seq_len(length(packages) / 2)
  info <- paste0(packages[col1], "     ", packages[-col1])
  
  msg(paste(info, collapse = "\n"), startup = TRUE)
  
  suppressPackageStartupMessages(
    lapply(to_load, same_library)
  )
  
  invisible()
}

package_version <- function(x) {
  version <- as.character(unclass(utils::packageVersion(x))[[1]])
  
  if (length(version) > 3) {
    version[4:length(version)] <- crayon::red(as.character(version[4:length(version)]))
  }
  paste0(version, collapse = ".")
}

is_attached <- function(x) {
  paste0("package:", x) %in% search()
}

.onAttach <- function(...) {
  needed <- core[!is_attached(core)]
  if (length(needed) == 0)
    return()
  
  crayon::num_colors(TRUE)
  scdney_attach()

  
}


msg <- function(..., startup = FALSE) {
  if (startup) {
    packageStartupMessage(text_col(...))
  } else {
    message(text_col(...))
  }
}

text_col <- function(x) {
  # If RStudio not available, messages already printed in black
  if (!rstudioapi::isAvailable()) {
    return(x)
  }
  
  if (!rstudioapi::hasFun("getThemeInfo")) {
    return(x)
  }
  
  theme <- rstudioapi::getThemeInfo()
  
  if (isTRUE(theme$dark)) crayon::white(x) else crayon::black(x)
  
}


