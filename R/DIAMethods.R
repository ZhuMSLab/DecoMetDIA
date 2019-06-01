setGeneric('MatchScanidx', function(fp, info.swath){
  # get swath setup
  seq.window <- seq(info.swath$nSWATH)
  names(seq.window) <- seq.window

  # read raw data with ms2 info
  xr <- xcmsRaw(fp, includeMSn = TRUE)
  xr.ms2 <- msn2xcmsRaw(xr)
  scantime.ms1 <- xr@scantime
  scantime.ms2 <- xr.ms2@scantime

  # preparing for ms1 and ms2 scan assignment
  scantime.range <- cbind(scantime.ms1, c(tail(scantime.ms1, -1),
                                          max(scantime.ms2) + 0.1))
  colnames(scantime.range) <- c('start', 'end')
  seq.scantime <- seq(nrow(scantime.range))
  seq.scantime.ms2 <- seq_along(scantime.ms2)
  # matching scantime for ms1 and ms2
  scantime.match <- lapply(seq.scantime, function(nr) {
    scantime.range.r <- scantime.range[nr, ]
    which(scantime.ms2 > scantime.range.r[1] & scantime.ms2 < scantime.range.r[2])
  })

  scantime.ms2.interval <- median(diff(scantime.ms2))
  scantime.wait <- median(scantime.ms2[sapply(scantime.match, `[`, 1)] - scantime.ms1) - scantime.ms2.interval
  scantime.ms2.theory <- seq(info.swath$nSWATH) * scantime.ms2.interval + scantime.wait/2

  # assigning ms2 scans to related ms1
  scidx.match <- lapply(seq.scantime, function(idx) {
    idx.st.match <- scantime.match[[idx]]
    if (length(idx.st.match) == info.swath$nSWATH) {
      idx.match <- seq.scantime.ms2[idx.st.match]
      names(idx.match) <- seq.window
    } else {
      scantime.theory <- scantime.range[idx, 1] + scantime.ms2.theory
      scantime.ms2.match <- scantime.ms2[idx.st.match]
      idx.match.theory <- sapply(scantime.ms2.match, function(st) {
        tmp <- st - scantime.theory
        which.min(tmp[tmp >= 0])
      })
      idx.match <- rep(NA, info.swath$nSWATH)
      names(idx.match) <- seq.window
      idx.match[idx.match.theory] <- seq.scantime.ms2[idx.st.match[idx.match.theory]]
    }
    idx.match
  })
  names(scidx.match) <- seq.scantime
  scidx.match
})


setGeneric('GetSWATHinfo', function(fp = '../SWATHsetup.csv') {
  infoSWATH <- list()
  if (file.exists(fp)) {
    swath.window <- read.csv(fp)
    nSWATH <- nrow(swath.window)
    overlapped <- 1 + which(swath.window[-1, 1] - swath.window[-nSWATH, 2] < 0)
    swath.window[overlapped, 1] <- swath.window[overlapped - 1, 2]
  } else {
    stop('No SWATH setup found!!!')
  }
  infoSWATH$nSWATH <- nSWATH
  infoSWATH$window <- swath.window
  return(infoSWATH)
})

# smoothing peaks with loess method
setGeneric('SmoothLoess', function(pre.data, span = 0.3, degree = 1L,
                                   is.span.constant = FALSE) {
  n.spec <- nrow(pre.data)
  spec.smooth <- data.frame(pre.data)
  if (is.span.constant) {
    pre.inf <- data.frame(idx = seq(nrow(pre.data)), intensity = pre.data[, 'intensity'])
    pre.fit <- loess(intensity~idx, data = pre.inf, span = span, degree = degree)
    pre.predict <- predict(pre.fit, data.frame(idx = seq(nrow(pre.data))))
    spec.smooth$intensity.s <- pre.predict
  } else {
    if (n.spec > 5) {
      tmp <- 4/n.spec
      span <- ifelse(tmp > span, tmp, span)
      pre.inf <- data.frame(idx = seq(nrow(pre.data)), intensity = pre.data[, 'intensity'])
      pre.fit <- loess(intensity~idx, data = pre.inf, span = span, degree = degree)
      pre.predict <- predict(pre.fit, data.frame(idx = seq(nrow(pre.data))))
      spec.smooth$intensity.s <- pre.predict
    } else{
      spec.smooth$intensity.s <- pre.data[, 'intensity']
    }
  }
  return(as.matrix(spec.smooth))
})

# get the index of the given data located in a range defined by row of matrix
setGeneric('GetIntervalIndex', function(dataPoint, dataMatrix) {
  isDone <- FALSE
  idx <- 1
  while (!isDone) {
    if ((dataPoint >= dataMatrix[idx, 1] & dataPoint < dataMatrix[idx, 2]) | idx == length(dataMatrix[,2])) {
      isDone <- TRUE
    } else {
      idx <- idx + 1
    }
  }
  return(idx)
})

# output a data range with given reference data and ppm
setGeneric('PpmRange', function(ref, ppm, spliter = 400) {
  dev <- ppm * 1e-6
  ref + c(-1, 1) * max(ref, spliter) * dev
})
