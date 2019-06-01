setGeneric('DecoSpectra', function(idx.pg, nm.smp, spectra.eics, pk.ms1, num.scantime,
                                    info.pk.ms1, peakwidth,
                                    idx.apex.eic, snthr = 5,
                                    is.PlotEIC.MS2 = TRUE,
                                    is.smooth = TRUE, isFWHM = FALSE, is.simple.mpk = TRUE,
                                    is.dec.all = FALSE, is.dec.smoothed = FALSE) {

  peaks.ms2.info <- lapply(seq_along(spectra.eics), function(idx.eic) {
    DetectPeaks(spectra.eics[[idx.eic]], peakwidth, num.scantime,
                idx.apex.eic, snthr = snthr, n.skip.max = 1)
  })
  if (all(sapply(peaks.ms2.info, is.null))) {
    return(NULL)
  }
  names(peaks.ms2.info) <- seq_along(spectra.eics)
  peaks.ms2.info <- peaks.ms2.info[which(!sapply(peaks.ms2.info, is.null))]

  res <- lapply(names(peaks.ms2.info), function(nm.info) {
    info <- peaks.ms2.info[[nm.info]]
    idx.a <- info$info.pk$apex
    idx.s <- info$info.pk$lb
    idx.e <- info$info.pk$rb
    is.peak <- info$info.pk$is.peak
    bl <- info$info.pk$bl

    info.adj <- data.frame(t(sapply(seq_along(idx.a), function(idx) {
      idx.eic.pk <- idx.s[idx]:idx.e[idx]
      eic.pk <- info$info.eic[idx.eic.pk, ]
      eic.int <- info$info.eic$i.s
      idx.eic.pk.keep <- idx.eic.pk[which(eic.pk$i > 0)]

      idx.l <- idx.s[idx] : (idx.a[idx] - 1)
      idx.r <- (idx.a[idx] + 1) : idx.e[idx]
      sharpness.l <- median((eic.int[idx.a[idx]] - eic.int[idx.l]) %/% (idx.a[idx] - idx.l))
      sharpness.r <- median((eic.int[idx.a[idx]] - eic.int[idx.r]) %/% (idx.r - idx.a[idx]))
      sharpness <- mean(c(sharpness.l, sharpness.r))

      sd.smooth <- sd(info$info.eic$i[idx.eic.pk.keep] - info$info.eic$i.s[idx.eic.pk.keep])/diff(range(info$info.eic$i.s[idx.eic.pk.keep]))

      c('shrp' = sharpness,
        'sd.smooth' = sd.smooth)
    })))

    data.frame(info.adj,
               'lb' = idx.s,
               'apex' = idx.a,
               'apex.t' = info$info.eic$rt[idx.a],
               'rb' = idx.e,
               'bl' = bl,
               'is.peak' = is.peak,
               stringsAsFactors = FALSE)
  })
  names(res) <- names(peaks.ms2.info)
  info.pk.ms2.raw <- info.pk.ms2 <- do.call(rbind, res)

  nr.rm <- which(info.pk.ms2$sd.smooth > 0.35 | !info.pk.ms2$is.peak)
  if (length(nr.rm)) {
    info.pk.ms2 <- info.pk.ms2[-nr.rm, , drop = FALSE]
  }
  if (nrow(info.pk.ms2) == 0) {
    return(NULL)
  }
  nm.split <- sapply(strsplit(rownames(info.pk.ms2), split = '\\.'), `[[`, 1)
  nm.dup <- nm.split[duplicated(nm.split)]
  info.pk.ms2$simp <- !nm.split %in% nm.dup
  info.pk.ms1.md <- data.frame(matrix(NA, ncol = ncol(info.pk.ms2)), row.names = '0')
  names(info.pk.ms1.md) <- names(info.pk.ms2)
  info.pk.ms1.md[, names(info.pk.ms1)] <- info.pk.ms1
  info.pk.ms1.md$is.peak <- TRUE
  info.pk.all <- rbind(info.pk.ms2, info.pk.ms1.md)

  idx.pk <- floor(as.numeric(rownames(info.pk.all)))
  rnm.info.pk <- rownames(info.pk.all)

  rt <- info.pk.all$apex.t
  names(rt) <- rownames(info.pk.all)
  cl.rt <- hclust(dist(rt), method = 'centroid')
  lb.cl.rt <- 15
  h.cutree <- max(lb.cl.rt, floor(max(cl.rt$height)/4))
  nc <- nrow(cl.rt$merge) + 2L -
    apply(outer(c(cl.rt$height, Inf), h.cutree, ">"), 2, which.max)
  r.cl.rt <- cutree(cl.rt, k = nc)
  cl.rt.nm <- lapply(unique(r.cl.rt), function(cl) {
    names(which(r.cl.rt == cl))
  })

  cl.rt.nm <- cl.rt.nm[sapply(cl.rt.nm, length) > 1]
  if (length(cl.rt.nm) == 0) {
    return(NULL)
  }
  idx.cl.rt.main <- which(sapply(cl.rt.nm, function(nm.pk) {
    '0' %in% nm.pk
  }))

  if (length(idx.cl.rt.main) == 0) {
    return(NULL)
  }

  eic.peak <- lapply(rnm.info.pk, function(nm.pk) {
    idx.eic <- info.pk.all[nm.pk, 'lb']:info.pk.all[nm.pk, 'rb']
    if (nm.pk == '0') {
      cbind('rt' = idx.eic,
            'i.s' = pk.ms1[, 'intensity.s'])
    } else {
      cbind('rt' = idx.eic,
            'i.s' = peaks.ms2.info[[gsub('\\.\\d+', '', nm.pk)]]$info.eic$i.s[idx.eic])
    }
  })
  names(eic.peak) <- rnm.info.pk
  cl.components <- list()
  dist.main <- NULL
  for(nm.pk in cl.rt.nm) {
    if (length(nm.pk) < 2) {
      next
    }
    idx.paires <- GetCombinations(length(nm.pk))
    distant.mtx <- matrix(0, ncol = length(nm.pk), nrow = length(nm.pk))
    colnames(distant.mtx) <- rownames(distant.mtx) <- nm.pk
    for (nr in seq(nrow(idx.paires))) {
      idx.pair <- idx.paires[nr, ]
      distant <- GetDistantP(eic.peak[[nm.pk[idx.pair[1]]]],
                             eic.peak[[nm.pk[idx.pair[2]]]])
      distant.mtx[idx.pair[1], idx.pair[2]] <- distant.mtx[idx.pair[2], idx.pair[1]] <- distant
    }
    distant.mtx[distant.mtx > 1] <- 1
    cl.pk <- hclust(as.dist(distant.mtx), method = 'centroid')
    h.cutree <- 0.25
    nc <- nrow(cl.pk$merge) + 2L -
      apply(outer(c(cl.pk$height, Inf), h.cutree, ">"), 2, which.max)
    r.cl.pk <- cutree(cl.pk, k = nc)
    cl.pk.nm <- lapply(unique(r.cl.pk), function(cl) {
      names(which(r.cl.pk == cl))
    })
    cl.pk.nm <- cl.pk.nm[sapply(cl.pk.nm, length) > 1]
    if (length(cl.rt.nm) == 0) {
      next
    }
    cl.components <- append(cl.components, cl.pk.nm)
  }

  if (length(cl.components) == 0) {
    return(NULL)
  }

  nm.mpk.main <- NULL
  nm.pk.main <- NULL
  nm.mpk <- c()
  for (nm.pk in cl.components){
    if ('0' %in% nm.pk) {
      if (!all(nm.pk ==  '0')) {
        nm.pk.main <- nm.pk <- nm.pk[-which(nm.pk == '0')]
        dist.main <- sapply(eic.peak[nm.pk], GetDistantP, eic.peak[['0']])
        nm.mpk.main <- nm.pk[which.min(dist.main[nm.pk])]
      } else {
        nm.mpk.main <- nm.pk.main <- nm.pk
      }
      nm.mpk <- append(nm.mpk, nm.mpk.main)
    } else {
      info.pk <- info.pk.ms2[nm.pk, ]
      if (any(info.pk$simp)) {
        nm.pk <- nm.pk[info.pk$simp]
        info.pk <- info.pk.ms2[nm.pk, ]
      }
      nm.mpk <- append(nm.mpk, nm.pk[which.max(info.pk$shrp)])
    }
  }
  names(cl.components) <- nm.mpk
  if (is.null(nm.pk.main)) {
    nm.pk.main <- nm.mpk.main <- '0'
    cl.components <- c(cl.components, '0' = list('0'))
    nm.mpk <- c(nm.mpk, '0')
  }
  nm.pk.simp <- nm.pk.main[info.pk.all[nm.pk.main, 'simp']]
  if (length(nm.pk.simp) > 0) {
    spec.smp <- lapply(nm.pk.simp, function(nm.pk) {
      nm.eic <- gsub('\\.\\d+', '', nm.pk)
      c('mz' = peaks.ms2.info[[nm.eic]]$info.eic[idx.apex.eic, 'mz'],
        'intensity' = peaks.ms2.info[[nm.eic]]$info.eic[idx.apex.eic, 'i'] - peaks.ms2.info[[nm.eic]]$info.pk[1, 'bl'],
        'dist' = unname(dist.main[nm.pk]))
    })
    spec.smp <- do.call(rbind, spec.smp[!is.na(spec.smp)])
  } else {
    spec.smp <- NULL
  }
  eic.model.pk <- lapply(nm.mpk, function(nm.pk) {
    if (nm.pk == '0') {
      if (!is.na(info.pk.all[nm.pk.main, 'bl'])) {
        eic <- rep(0, nrow(peaks.ms2.info[[1]]$info.eic))
        eic[info.pk.all[nm.pk.main, 'lb']:info.pk.all[nm.pk.main, 'rb']] <- pk.ms1[, 'intensity.s'] - info.pk.all[nm.pk.main, 'bl']
      } else {
        return(NULL)
      }
    } else {
      idx.eic <- gsub('\\.\\d+', '', nm.pk)
      eic <- peaks.ms2.info[[idx.eic]]$info.eic[, 'i.s']
      info.mpk <- info.pk.all[nm.pk, ]
      eic[info.mpk$lb:info.mpk$rb] <- eic[info.mpk$lb:info.mpk$rb] - info.mpk$bl
      eic[-(info.mpk$lb:info.mpk$rb)] <- 0
    }
    eic[eic < 0] <- 0
    if (max(eic) == 0) {
      return(NULL)
    }
    eic/max(eic)
  })
  names(eic.model.pk) <- nm.mpk
  if (any(sapply(eic.model.pk, is.null))) {
    return(NULL)
  }
  mpk.mtx <- do.call(cbind, eic.model.pk)
  if (is.simple.mpk) {
    idx.sep <- which(rowSums(mpk.mtx) == 0)
    if (length(idx.sep) > 0) {
      idx.sep <- unique(c(0, idx.sep, nrow(mpk.mtx)+1))
      i.start <- max(idx.sep[idx.sep < idx.apex.eic]) + 1
      i.end <- min(idx.sep[idx.sep > idx.apex.eic]) - 1
    } else {
      i.start <- 1
      i.end <- nrow(mpk.mtx)
    }
    mpk.mtx <- mpk.mtx[i.start:i.end, , drop = FALSE]
    idx.keep <- which(apply(mpk.mtx, 2, function(dc) any(dc > 0)))
    mpk.mtx <- mpk.mtx[, idx.keep, drop = FALSE]
  } else {
    i.start <- 1
    i.end <- nrow(mpk.mtx)
  }
  if (ncol(mpk.mtx) <= 1) {
    return(list('spec' = spec.smp,
                'spec.smp' = spec.smp,
                'spec.cpl' = NULL))
  }
  idx.mpk.main <- which(colnames(mpk.mtx) == nm.mpk.main)
  info.pk.ms1$apex <- info.pk.ms1$apex - i.start + 1
  if (is.dec.all) {
    nm.pk.dec <- rownames(info.pk.ms2.raw)[!(rownames(info.pk.ms2.raw) %in% nm.pk.simp)]
  } else {
    nm.pk.dec <- cl.rt.nm[[idx.cl.rt.main]]
    nm.pk.dec <- nm.pk.dec[-which(nm.pk.dec == '0')]
    nm.pk.dec <- nm.pk.dec[!info.pk.ms2[nm.pk.dec, 'simp']]
    nm.pk.dec <- rownames(info.pk.ms2)[!info.pk.ms2[, 'simp']]
  }
  idx.pk.dec <- floor(as.numeric(rownames(info.pk.ms2.raw)))

  nm.int.dec <- ifelse(is.dec.smoothed, 'i.s', 'i')

  if (length(nm.pk.dec) > 0) {
    nm.eic.dec <- unique(gsub('\\.\\d+', '', nm.pk.dec))
    eic.pre.dec <- lapply(nm.eic.dec, function(idx.eic) {
      info.pk.i <- info.pk.ms2.raw[which(idx.pk.dec == idx.eic), , drop = FALSE]
      i.scan <- min(info.pk.i$lb):max(info.pk.i$rb)
      eic <- peaks.ms2.info[[as.character(idx.eic)]]$info.eic[, nm.int.dec]
      eic[i.scan] <- eic[i.scan] - info.pk.i[1, ]$bl
      eic[-i.scan] <- 0
      eic[eic < 0] <- 0
      eic[i.start:i.end]
    })
    names(eic.pre.dec) <- nm.eic.dec
    eic.mtx <- do.call(cbind, eic.pre.dec)
    info.decomposite <- lapply(names(eic.pre.dec), function(nm.eic) {
      eic <- eic.pre.dec[[nm.eic]]
      a <- optim(par = rep(0, ncol(mpk.mtx)), fn = funOptim, gr = NULL,
                 eic, mpk.mtx,
                 lower = 0,
                 method="L-BFGS-B")$par
      tmp <- sapply(seq_along(a), function(i) {
        a[i] * mpk.mtx[, i]
      })
      eic.fit.mpk <- apply(tmp, 1, sum)
      resdual <- eic - eic.fit.mpk
      list('a' = a,
           'eic.dcps.mpk' = tmp[, idx.mpk.main],
           'mz' = unname(peaks.ms2.info[[nm.eic]]$info.eic[info.pk.ms1$apex, 'mz']),
           'into' = unname(peaks.ms2.info[[nm.eic]]$info.eic[info.pk.ms1$apex, 'i']),
           'intd' = tmp[info.pk.ms1$apex, idx.mpk.main],
           'bl' = unname(peaks.ms2.info[[nm.eic]]$info.pk$bl[1]))
    })
    spec.cpl <- t(sapply(info.decomposite, function(info) {
      c('mz' = info$mz, 'intensity' = info$'intd', 'dist' = 0)
    }))
    is.keep.cpl <- spec.cpl[, 'intensity'] > 0
    if (any(is.keep.cpl)) {
      spec.cpl <- spec.cpl[is.keep.cpl, , drop = FALSE]
    } else {
      spec.cpl <- NULL
    }
  } else {
    spec.cpl <- NULL
  }

  t.mpk <- peaks.ms2.info[[1]][[2]]$rt[i.start:i.end]
  if (is.PlotEIC.MS2) {
    d.plot <- file.path('components_hclust', nm.smp)
    if (!dir.exists(d.plot)) {
      dir.create(d.plot, recursive = TRUE)
    }
    colnames(pk.ms1) <- c('s', 'rt', 'mz', 'i', 'i.s')

    pdf(paste0(d.plot, '/', idx.pg, '.pdf'))
    rg.pk <- which(spectra.eics[[1]][, 'rt'] >= pk.ms1[1, 'rt'] &
                     spectra.eics[[1]][, 'rt'] <= tail(pk.ms1[, 'rt'], 1))
    xlm <- range(spectra.eics[[1]][, 'rt'])
    rg.int <- c(0, max(sapply(spectra.eics, function(eic) {
      eic[rg.pk, 'intensity']
    })))
    ylm <- range(c(rg.int, pk.ms1[, 'i']))
    plot(0,0,type = 'n', col = 'red', lwd = 2,
         xlim = xlm, ylim = ylm, main = idx.pg,
         xlab = 'Retention time (s)',
         ylab = 'Intensity')
    idx.eic.spec <- floor(as.numeric(c(nm.pk.main, nm.pk.dec)))
    sapply(seq_along(spectra.eics), function(idx.eic) {
      eic <- spectra.eics[[idx.eic]]
      if (idx.eic %in% idx.eic.spec) {
        lwd.plot <- 2
        clr <- 'red'
      } else {
        lwd.plot <- 1
        clr <- 'blue'
      }
      lines(eic[, c('rt', 'intensity')], lwd = lwd.plot, col = clr)
    })
    lines(pk.ms1[, c('rt', 'i')],type = 'l', col = 'black', lwd = 3, lty = 2)
    cls <- rainbow(length(cl.rt.nm))
    plot(pk.ms1[, 'rt'], pk.ms1[, 'i.s']/max(pk.ms1[, 'i.s']),
         col = 'red', lwd = 3, type = 'n',
         xlim = range(spectra.eics[[1]][,'rt']),
         ylim = c(0, 1),
         main = idx.pg,
         xlab = 'Retention time (s)',
         ylab = 'Relative intensity')

    for (idx in seq_along(cl.rt.nm)) {
      nms.pk <- cl.rt.nm[[idx]]
      for (nm.pk in nms.pk) {
        nm.eic <- strsplit(nm.pk, split = '\\.')[[1]][1]
        if (nm.pk == '0') {
          abline(v = info.pk.all['0', 'apex.t'], col = cls[idx], lwd = 2)
        } else {
          dr <- info.pk.all[nm.pk, ]
          idx.plot <- dr$lb:dr$rb

          lines(peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'rt'],
                peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'i.s']/max(peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'i.s']),
                col = cls[idx],
                lwd = 1)
        }
      }
    }
    lines(pk.ms1[, 'rt'], pk.ms1[, 'i.s']/max(pk.ms1[, 'i.s']),
          col = 'black', lwd = 2, lty = 2)

    cls <- rainbow(length(cl.components))
    plot(pk.ms1[, 'rt'], pk.ms1[, 'i.s']/max(pk.ms1[, 'i.s']),
         col = 'red', lwd = 3, type = 'n',
         xlim = range(spectra.eics[[1]][,'rt']),
         ylim = c(0, 1),
         main = idx.pg,
         xlab = 'Retention time (s)',
         ylab = 'Relative intensity')

    for (idx in seq_along(cl.components)) {
      nms.pk <- cl.components[[idx]]
      for (nm.pk in nms.pk) {
        nm.eic <- strsplit(nm.pk, split = '\\.')[[1]][1]
        if (nm.pk == '0') {
          abline(v = info.pk.all['0', 'apex.t'], col = cls[idx], lwd = 2)
        } else {
          dr <- info.pk.all[nm.pk, ]
          idx.plot <- dr$lb:dr$rb

          lines(peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'rt'],
                peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'i.s']/max(peaks.ms2.info[[nm.eic]]$info.eic[idx.plot, 'i.s']),
                col = cls[idx],
                lwd = 1)
        }
      }
    }
    lines(pk.ms1[, 'rt'], pk.ms1[, 'i.s']/max(pk.ms1[, 'i.s']),
          col = 'black', lwd = 2, lty = 2)
    cls <- rainbow(length(eic.model.pk))
    plot(0,0, type = 'n', xlim = range(t.mpk), ylim = c(0, 1),
         xlab = 'Retention time (s)',
         ylab = 'Relative intensity')
    sapply(seq_along(eic.model.pk), function(idx) {
      eic <- eic.model.pk[[idx]][i.start:i.end]
      rt
      lines(t.mpk, eic/max(eic), col = cls[idx], lwd = 1)
    })
    lines(pk.ms1[, 'rt'], pk.ms1[, 'i.s']/max(pk.ms1[, 'i.s']),
          col = 'black', lwd = 2, lty = 2)
    dev.off()
  }

  if(is.null(ncol(spec.smp)) & is.null(ncol(spec.cpl))) {
    return(NULL)
  }
  spec <- na.omit(rbind(spec.smp, spec.cpl))
  spec <- spec[order(spec[, 'mz']), , drop = FALSE]
  rownames(spec) <- NULL
  if (is.PlotEIC.MS2) {
    if (!dir.exists('decomposite')) dir.create('decomposite')
    d.plot <- file.path('decomposite', nm.smp)
    if (!dir.exists(d.plot)) {
      dir.create(d.plot, recursive = TRUE)
    }
    pdf(paste0(d.plot, '/', idx.pg, '.pdf'))
    plot(0, 0, type = 'n',
         xlim = range(spectra.eics[[1]][,'rt']),
         ylim = c(0, max(spec[, 2])),
         main = idx.pg,
         xlab = 'Retention time (s)',
         ylab = 'Intensity')
    i.s <- which(spectra.eics[[1]][,'rt'] >= pk.ms1[1, 'rt'])[1]
    i.e <- max(which(spectra.eics[[1]][,'rt'] <= tail(pk.ms1[, 'rt'], 1)))
    if (length(nm.pk.dec) > 0) {
      eic.decomposite <- lapply(info.decomposite, `[[`, 'eic.dcps.mpk')
      sapply(eic.decomposite, function(eic.1) {
        lines(spectra.eics[[1]][i.start:i.end,'rt'], eic.1)
      })
    }
    if (length(nm.pk.simp) > 0) {
      sapply(nm.pk.simp, function(nm.pk) {
        nm.eic <- gsub('\\.\\d+', '', nm.pk)
        lines(peaks.ms2.info[[nm.eic]]$info.eic[i.s:i.e, 'rt'],
              peaks.ms2.info[[nm.eic]]$info.eic[i.s:i.e, nm.int.dec] - peaks.ms2.info[[nm.eic]]$info.pk[1, 'bl'],
              col = 'blue')
      })
    }
    dev.off()
  }
  return(list('spec' = spec,
              'spec.smp' = spec.smp,
              'spec.cpl' = spec.cpl))
})


GetDistantP <- function(x, y) {
  dt <- mergeEIC(x, y)
  c.g <- Gauss(x = which.max(dt[, 2]) -
                 which.max(dt[, 3]),
               h = 1,
               mu = 0,
               sigma = 12)
  if (c.g < 0.7) return (1)
  c.p <- cor(dt[, 2], dt[, 3]) * c.g
  if (is.na(c.p) | is.nan(c.p) | is.null(c.p)) {
    c.p <- 0
  }
  1 - c.p
}
