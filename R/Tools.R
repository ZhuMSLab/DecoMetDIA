setGeneric('CatLineSeparator', function(info) {
  info.sep <- paste(rep('-', 38), collapse  = '')
  cat('', info.sep, info, info.sep, '', sep = '\n')
})

setGeneric('CheckSkip', function(fn, rerun) {
  # check if use the existed results. if rerun = FALSE, ignore this check
  a <- file.exists(fn) & {!rerun}
  if (a) {
    cat('Using previous results: ', fn, '\n')
  }
  a
})

setGeneric('CreateResultDir', function(d.out, is.overwrite = TRUE) {
  if (file.exists(d.out) & !is.overwrite) {
    d.out.it <- 1
    while (TRUE) {
      d.out.tmp <- sprintf('%s.%s', d.out, d.out.it)
      if (!file.exists(d.out.tmp)) {
        break
      } else {
        d.out.it <- d.out.it + 1
      }
    }
    d.out <- d.out.tmp
  }

  if (!file.exists(d.out)) {
    dir.create(d.out, recursive = TRUE)
  }

  return(d.out)
})

setGeneric('OutputSpec', function(spec, idx, dt.peaktable, fn.save = NULL, grp = NULL) {
  if (!is.null(fn.save)) {
    sink(fn.save)
  }
  if (length(spec) == 0) {
    return(NULL)
  }
  ft.nm <- dt.peaktable[idx, 'name']
  mz <- dt.peaktable[idx, 'mzmed']
  rt <- dt.peaktable[idx, 'rtmed']

  cat('NAME: ', ft.nm, '\n', sep = '')
  if (!is.null(grp)) {
    cat('SAMPLEGROUP: ', grp, '\n', sep = '')
  }
  cat('PRECURSORMZ: ', mz, '\n', sep = '')
  cat('RETENTIONTIME: ', rt, '\n', sep = '')
  cat('Num Peaks: ', nrow(spec), '\n', sep = '')
  for (nr in seq(nrow(spec))) {
    cat(paste(spec[nr, ], collapse = ' '), '\n', sep = '')
  }
  cat('\n')
  if (!is.null(fn.save)) {
    sink()
  }
})

setGeneric('Gauss', function(x, h, mu, sigma){
  h*exp(-(x-mu)^2/(2*sigma^2))
})

setGeneric('GetCombinations', function(n) {
  res <- c()
  for(i in seq(n-1)) {
    res <- rbind(res, cbind(i, seq(i+1, n)))
  }
  colnames(res) <- NULL
  return(res)
})
