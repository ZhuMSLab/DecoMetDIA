DeIsotope <- function(spec, num.isotope = 2, ppm = 30, prop = c(0.9, 0.5)) {
  if (length(prop) != num.isotope) {
    stop('Parameters num.isotope and prop do not match!!')
  }
  spec <- spec[order(spec[, 1]), , drop = FALSE]
  idx.isotope <- which(findIsotope(spec, num.isotope, ppm, prop) == 1)
  if (length(idx.isotope) > 0) {
    spec[-idx.isotope, , drop = FALSE]
  } else {
    spec
  }
}