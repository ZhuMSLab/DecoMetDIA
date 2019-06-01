setGeneric('ConsensusSpec', function(spectra.all, smpclass, num.spec,
                                     minfrac.vote = 0.5,
                                     int.thr.abs = 0,
                                     is.norm.consensus = FALSE) {
  spec.final <- lapply(seq(num.spec), function(idx) {
    spec.cons.grp <- lapply(unique(smpclass), function(nm.smpclass) {
      idx.smp <- which(smpclass %in% nm.smpclass)
      num.smp <- length(idx.smp)
      spec.pk.grp <- lapply(idx.smp, function(idx.smp.1) {
        spec <- spectra.all[[idx.smp.1]][[idx]]$spec[, 1:2, drop = FALSE]
        DeNoise(spec, int.thr.abs = int.thr.abs)
      })
      if (all(sapply(spec.pk.grp, is.null))) {
        return(NULL)
      }
      mz.merged <- MergeMZ(spec.pk.grp)
      spec.extract <- lapply(spec.pk.grp, ExtractSpectra, mz.merged)
      is.contained <- matrix(sapply(spec.extract, function(spec) {
        attributes(spec)$is.found
      }), ncol = num.smp)
      num.contained <- rowSums(is.contained)
      if (max(num.contained) < min(2, num.smp)) {
        return(NULL)
      }
      int.fragment <- matrix(sapply(spec.extract, function(spec) {
        attributes(spec)$intensity
      }), ncol = num.smp)
      idx.keep <- which(num.contained == max(num.contained))
      if (length(idx.keep) == 0) {
        return(NULL)
      }
      idx.norm <- idx.keep[which.max(rowSums(int.fragment[idx.keep, , drop = FALSE]))]
      col.keep <- which(is.contained[idx.norm, ])
      num.contained <- rowSums(is.contained[, col.keep, drop = FALSE])
      if (is.norm.consensus) {
        factor.norm <- int.fragment[idx.norm, ]
      } else {
        factor.norm <- rep(1, ncol(int.fragment))
      }
      int.norm <- t(t(int.fragment[, col.keep, drop = FALSE])/factor.norm[col.keep])
      idx.keep <- which(num.contained/length(col.keep) >= minfrac.vote)
      if (minfrac.vote == 0) {
        idx.keep <- which(num.contained/length(col.keep) > minfrac.vote)
      }
      if (length(idx.keep) == 0) {
        return(NULL)
      }
      mz.consensus <- mz.merged[idx.keep]
      int.consensus <- sapply(idx.keep, function(idx) {
        median(int.norm[idx, is.contained[idx, col.keep]])
      })
      spec.consensus <- as.matrix(data.frame('mz' = mz.consensus,
                                             'intensity' = int.consensus))
      res <- RemoveRingEffect(spec.consensus)
      return(res)
    })
    names(spec.cons.grp) <- unique(smpclass)
    idx.keep <- which(!sapply(spec.cons.grp, is.null))
    if (length(idx.keep) == 0) {
      return(NULL)
    } else {
      return(spec.cons.grp)
    }
  })
  return(spec.final)
})

setGeneric('DeNoise', function(spec, is.norm = FALSE, int.thr.rel = 0.0, int.thr.abs = 200,
                               ppm = 30, num.isotope = 2, prop = c(0.9, 0.5)) {
  if (is.null(spec)) return(NULL)
  spec <- RemoveRingEffect(spec)
  spec <- DeIsotope(spec, num.isotope = num.isotope, ppm = ppm, prop = prop)
  spec <- spec[spec[, 2] >= int.thr.abs, , drop = FALSE]
  if (nrow(spec) == 0) return(NULL)
  idx.int.max <- which.max(spec[, 2])
  if (is.norm) {
    spec[, 2] <- spec[, 2] / spec[idx.int.max, 2]
  }

  spec <- spec[spec[, 2] > spec[idx.int.max, 2] * int.thr.rel, , drop = FALSE]
  if (nrow(spec) == 0) {
    return(NULL)
  } else {
    return(spec)
  }
})

setGeneric('MergeMZ', function(ms2.all, ppm = 30, mz.bin.min = 0.004) {
  GetWeightedMZ <- function(mz.info) {
    sum(mz.info[, 'mz'] * mz.info[, 'intensity']) / sum(mz.info[, 'intensity'])
  }
  spec.all <- do.call(rbind, ms2.all)
  spec.all <- spec.all[order(spec.all[, 'mz']), , drop = FALSE]

  idx.left <- seq(nrow(spec.all))

  mz.merged <- {}
  while (length(idx.left) > 0 ) {
    idx <- tail(idx.left, 1)
    mz <- spec.all[idx, 'mz']
    mz.range <- c(-1, 0) * max(prod(mz, ppm, 1e-6), mz.bin.min) + mz
    idx.range <- idx.left[spec.all[idx.left, 'mz'] >= mz.range[1] &
                            spec.all[idx.left, 'mz'] <= mz.range[2]]
    mz.tmp <- GetWeightedMZ(spec.all[idx.left[idx.range], , drop = FALSE])
    mz.merged <- rbind(mz.merged, mz.tmp)
    idx.left <- idx.left[-idx.range]
  }
  mz.merged <- sort(mz.merged)

  return(mz.merged)
})

setGeneric('ExtractSpectra', function(spec, mz, ppm = 30, mz.bin.min = 0.004) {
  spec.tmp <- lapply(mz, function(mz.1) {
    mz.range <- mz.1 + max(prod(mz.1, ppm, 1e-6), mz.bin.min) * c(-1, 1)
    idx <- which(spec[, 'mz'] >= mz.range[1] & spec[, 'mz'] <= mz.range[2])
    if (length(idx) == 0) {
      return(NULL)
    } else if (length(idx) > 1) {
      idx <- idx[which.min(abs(spec[idx, 'mz'] - mz.1))]
    }
    spec[idx, ]
  })
  res <- do.call(rbind, spec.tmp)

  attributes(res)$is.found <- !sapply(spec.tmp, is.null)
  attributes(res)$intensity <- sapply(spec.tmp, function(spec) {
    ifelse(is.null(spec), 0, spec[2])
  })
  return(res)
})

setGeneric('RemoveRingEffect', function(spec, mz.diff.thr = 0.3, int.rel.thr = 0.2) {
  nr.ring <- nrow(spec) + 1
  mz <- spec[, 'mz']

  mz.diff <- diff(mz)
  idx.mzdiff <- which(mz.diff <= mz.diff.thr)
  if (length(idx.mzdiff) == 0) {
    return(spec)
  }

  nr.ring.possible <- unique(c(idx.mzdiff, idx.mzdiff + 1))
  while (TRUE) {

    idx.int.max <- which.max(spec[nr.ring.possible, 2])
    nr.int.max <- nr.ring.possible[idx.int.max]
    int.thr <- spec[nr.int.max, 2] * int.rel.thr

    mz.diff <- abs(mz[nr.ring.possible[-idx.int.max]] - mz[nr.int.max])
    int <- spec[nr.ring.possible[-idx.int.max], 2]
    nr.ring <- append(nr.ring, nr.ring.possible[-idx.int.max][which(mz.diff <= mz.diff.thr & int <= int.thr)])
    nr.ring.possible <- nr.ring.possible[!nr.ring.possible %in% c(nr.ring, nr.int.max)]
    if (length(nr.ring.possible) == 0) {
      break
    }
  }

  return(spec[-nr.ring, , drop = FALSE])
})
