setGeneric(
  'DetectPeaks',
  function(eic, peakwidth, num.scantime, idx.apex.eic,
           snthr = 3, is.smooth = TRUE, n.skip.max = 0) {
    is.peak <- FALSE
    eic <- data.frame(SmoothLoess(eic, span = 0.1))
    colnames(eic) <- c('s', 'rt', 'mz', 'i', 'i.s')
    eic.int <- eic[, ifelse(is.smooth, 'i.s', 'i')]
    is.roi.skip <- is.roi <- getContinuousPtsAboveThrIdx(eic[, 'i'], thr=0,
                                                         num=peakwidth[1],
                                                         iStart = 0, nSkipMax = 1)
    idx.fr.roi <- GetRoi(is.roi, idx.apex.eic)
    idx.fr.roi.skip <- GetRoi(is.roi.skip, idx.apex.eic)
    if (is.null(idx.fr.roi)) {
      return(NULL)
    }
    noise <- EstimateChromNoise(eic[, 'i'], trim = 0.05, min.pts = 12)

    noise.local <- GetLocalNoiseEstimate(eic[, 'i'],
                                         td = eic[, 's'],
                                         ftd = eic[idx.fr.roi, 's'],
                                         noiserange = c(12, 42),
                                         Nscantime = num.scantime,
                                         threshold = noise, num = 4,
                                         n.skip.max = 0)
    ## Final baseline & Noise estimate
    baseline <- max(1, min(noise.local[1], noise))
    sdnoise <- max(1, noise.local[2])
    sdnoise <- ifelse(is.na(sdnoise), 0, sdnoise)
    intthr <-  sdnoise * snthr

    not.noise <- getContinuousPtsAboveThrIdx(eic[, 'i'],
                                             thr=baseline, num=4, iStart = 0,
                                             nSkipMax = n.skip.max)
    idx.not.noise <- GetRoi(not.noise, idx.apex.eic)

    if (is.null(idx.not.noise)) {
      return(NULL)
    }
    idx.fr.roi.denoise <- idx.fr.roi.skip[idx.fr.roi.skip %in% idx.not.noise]
    if(length(idx.fr.roi.denoise) == 0) {
      return(NULL)
    }
    ## is there any data above S/N * threshold ?
    if (!(any(eic[, 'i'] - baseline >= intthr))) {
      lb <- 1
      rb <- nrow(eic)
      idx.apex.frroi <- idx.apex.eic
    } else {
      idx.apex.roi <- FindLocalMax(eic$i.s[idx.fr.roi.denoise], v = baseline)
      idx.bnd.roi <- FindLocalMin(eic$i.s[idx.fr.roi.denoise])

      if (is.null(idx.apex.roi) | length(idx.apex.roi) == 0) {
        lb <- 1
        rb <- nrow(eic)
        idx.apex.frroi <- idx.apex.eic
      } else {
        if(is.null(idx.bnd.roi)) {
          idx.bnd.frroi <- range(idx.fr.roi.denoise)
        }
        idx.apex.frroi <- idx.fr.roi.denoise[idx.apex.roi]
        idx.bnd.frroi <- idx.fr.roi.denoise[idx.bnd.roi]

        idx.bnd.can <- sort(unique(c(range(idx.fr.roi.denoise), idx.bnd.frroi)))


        lb <- c()
        rb <- c()
        sharpness <- c()
        i.rm <- c()
        for (idx in seq_along(idx.apex.frroi)) {
          idx.diff <- idx.bnd.can - idx.apex.frroi[idx]
          s.idx.diff <- sign(idx.diff)
          idx.sign.change <- which(diff(s.idx.diff) > 0)
          lb[idx] <- idx.bnd.can[idx.sign.change]
          if (idx.sign.change == length(idx.bnd.can)) {
            rb[idx] <- idx.bnd.can[idx.sign.change]
          } else {
            rb[idx] <- idx.bnd.can[idx.sign.change + 1]
          }

          if (idx > 1) if (lb[idx] < idx.apex.frroi[idx - 1]) {
            lb[idx] <- idx.apex.frroi[idx - 1] + which.min(eic.int[idx.apex.frroi[idx - 1] : idx.apex.frroi[idx]]) - 1
          }
          if (idx < length(idx.apex.frroi)) if(rb[idx] > idx.apex.frroi[idx + 1]) {
            rb[idx] <- idx.apex.frroi[idx] + which.min(eic.int[idx.apex.frroi[idx] : idx.apex.frroi[idx+1]]) - 1
          }
          if (lb[idx] == idx.apex.frroi[idx] | rb[idx] == idx.apex.frroi[idx]) {
            i.rm <- append(i.rm, idx)
            next
          }
        }

        if (length(i.rm) > 0) {
          lb <- lb[-i.rm]
          rb <- rb[-i.rm]
          idx.apex.frroi <- idx.apex.frroi[-i.rm]
        }
        if (length(idx.apex.frroi) > 0) {
          is.peak <- TRUE
        }
      }
    }

    info.pk <-  data.frame('lb' = lb,
                           'apex' = idx.apex.frroi,
                           'rb' = rb,
                           'bl' = baseline,
                           'is.peak' = is.peak)
    return(list('info.pk' = info.pk,
                'info.eic' = eic))
  })


setGeneric('FindLocalMax', function(x, m = 3, v = 200){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape < 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] <= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  unname(pks[x[pks] >= v])
})


setGeneric('FindLocalMin', function(x, m = 3){
  shape <- diff(sign(diff(x, na.pad = FALSE)))
  pks <- sapply(which(shape > 0), FUN = function(i){
    z <- i - m + 1
    z <- ifelse(z > 0, z, 1)
    w <- i + m + 1
    w <- ifelse(w < length(x), w, length(x))
    if(all(x[c(z : i, (i + 2) : w)] >= x[i + 1])) return(i + 1) else return(numeric(0))
  })
  pks <- unlist(pks)
  unname(pks)
})

GetRoi <- function(is.roi, idx.apex.eic) {
  if (all(is.roi == 0)) {
    return(NULL)
  }
  i.s <- which(diff(is.roi) == 1)
  i.e <- which(diff(is.roi) == -1)
  if (is.roi[1] == 1) i.s <- c(0, i.s)
  if (is.roi[length(is.roi)] == 1) i.e <- c(i.e, length(is.roi))
  i.s <- ifelse1(length(i.s) == 0, 1, i.s+1)
  i.e <- ifelse1(length(i.e) == 0, length(is.roi), i.e)
  range.roi <- cbind(i.s, i.e)

  idx.fr.roi <- NULL
  for (idx in seq(nrow(range.roi))) {
    if (range.roi[idx, 1] <= idx.apex.eic &
        range.roi[idx, 2] >= idx.apex.eic) {
      idx.fr.roi <- range.roi[idx, 1] : range.roi[idx, 2]
      break
    }
  }
  return(idx.fr.roi)
}

EstimateChromNoise <- function(x, trim=0.05, min.pts=20) {
  # modified from xcms
  gz <- which(x > 0)
  if (length(gz) < min.pts)
    return(mean(x))

  mean(x[gz], trim=trim)
}

GetLocalNoiseEstimate <- function(d, td, ftd, noiserange, Nscantime, threshold, num, n.skip.max = 0) {
  # modified from xcms
  drange <- which(td %in% ftd)

  CalculateBL <- function(d, drange, threshold, num, n.skip.max, noiserange) {
    n1 <- d[-drange] ## region outside the detected ROI (wide)
    n1.cp <-try(getContinuousPtsAboveThrIdx(n1, thr=threshold, num=num, iStart = 0, nSkipMax = n.skip.max)) ## continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n1 <- n1[!n1.cp]
    if (length(n1) > 1)  {
      baseline1 <- mean(n1)
      sdnoise1 <- sd(n1)
    } else
      baseline1 <- sdnoise1 <- 1

    ## noiserange[1]
    d1 <- drange[1]
    d2 <- drange[length(drange)]
    nrange2 <- c(max(1,d1 - noiserange[1]) : d1, d2 : min(length(d),d2 + noiserange[1]))
    # region outside the detected ROI (narrow)
    n2 <- d[nrange2]
    # continousPtsAboveThreshold (probably peak) are subtracted from data for local noise estimation
    n2.cp <- getContinuousPtsAboveThrIdx(n2, thr=threshold, num=num, iStart = 0,
                                         nSkipMax = n.skip.max)
    n2 <- n2[!n2.cp]
    if (length(n2) > 1)  {
      baseline2 <- mean(n2)
      sdnoise2 <- sd(n2)
    } else
      baseline2 <- sdnoise2 <- 1

    c(min(baseline1,baseline2),min(sdnoise1,sdnoise2))
  }

  if (length(drange) < Nscantime) {
    res <- CalculateBL(d, drange, threshold, num, n.skip.max, noiserange)
  } else {
    mid <- ceiling(length(d) / 2)
    drange <- seq(max(1, mid - noiserange[1]), min(mid + noiserange[1], length(d)))
    res <- res1 <- CalculateBL(d, drange, threshold, num, n.skip.max, noiserange)
    drange <- seq(max(1, mid - noiserange[2]), min(mid + noiserange[2], length(d)))

    if (length(drange) < Nscantime) {
      res2 <- CalculateBL(d, drange, threshold, num, n.skip.max, noiserange)
      res <- colMins(rbind(res1, res2))
    }
  }

  return(res)
}

trimm <- function(x, trim=c(0.05,0.95)) {
  a <- sort(x[x>0])
  Na <- length(a)
  quant <- round((Na*trim[1])+1):round(Na*trim[2])
  a[quant]
}
