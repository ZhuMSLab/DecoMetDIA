setMethod(
  'GetDIASpectra',
  'diaSpectra',
  function(obj, fp.smp, dia.feature, dt.peaktable,
           d.spec = 'SpecDIA',
           peakwidth = c(5, 30), snthr = 5,
           int.filter = 100, ppm.ms2.mtp = 15,
           is.PlotEIC.MS2 = TRUE,
           span = 0.1,
           rerun = FALSE,
           is.span.constant = FALSE,
           fp.swath = '../SWATHsetup.csv',
           rt.extend = 2.5,
           is.dec.all = FALSE,is.dec.smoothed = FALSE) {
    pg.raw <- dia.feature@peakgroup$raw
    seq.pg <- seq_along(pg.raw)
    idx.smp <- which(dia.feature@filepaths %in% fp.smp)
    if (length(idx.smp) == 0) {
      if (file.exists(fp.smp)) {
        idx.smp <- which(basename(dia.feature@filepaths) %in% basename(fp.smp))
      } else {
        stop(paste0("file not exist: \n  ", fp.smp))
      }
    }

    xr <- xcmsRaw(fp.smp, includeMSn = TRUE)
    scantime.ms1 <- xr@scantime
    mz.range <- xr@mzrange

    scan.ms1 <- lapply(seq_along(xr@scantime) , function(idx) {
      getScan(xr, idx)
    })
    xr.ms2 <- msn2xcmsRaw(xr)
    scantime.ms2 <- xr.ms2@scantime
    scan.ms2 <- lapply(seq_along(xr.ms2@scantime) , function(idx) {
      getScan(xr.ms2, idx)
    })
    rm(list = c('xr', 'xr.ms2'))

    info.swath <- GetSWATHinfo(fp = fp.swath)
    scanidx.match <- MatchScanidx(fp.smp, info.swath)

    rt.info <- do.call(rbind, lapply(seq.pg, function(idx.pg) {
      pk <- pg.raw[[idx.pg]]
      data.frame('idx' = idx.pg, pk[pk$sample == idx.smp, c('rt', 'rtmin', 'rtmax')])
    }))

    if (!is.null(rt.extend)) {
      rt.info$rtmin.ext <- rt.info$rt - rt.extend * (rt.info$rtmax - rt.info$rtmin)
      rt.info$rtmax.ext <- rt.info$rt + rt.extend * (rt.info$rtmax - rt.info$rtmin)
      rt.info$rtmin.ext[rt.info$rtmin.ext < 0] <- 0
      rt.info$rtmax.ext[rt.info$rtmax.ext > max(scantime.ms1)] <- max(scantime.ms1)
    } else {
      rt.info$rtmin.ext <- rt.info$rtmin
      rt.info$rtmax.ext <- rt.info$rtmax
    }
    info.all <- lapply(seq.pg, function(idx.pg) {
      idx.ms1.raw <- idx.ms1 <- which(scantime.ms1 >= rt.info[idx.pg, 'rtmin'] &
                                        scantime.ms1 <= rt.info[idx.pg, 'rtmax'])
      idx.ms1.ext.raw <- idx.ms1.ext <- which(scantime.ms1 >= rt.info[idx.pg, 'rtmin.ext'] &
                                                scantime.ms1 <= rt.info[idx.pg, 'rtmax.ext'])
      mz <- pg.raw[[idx.pg]][idx.smp, 'mz']
      mz.max <- PpmRange(ref = mz, ppm = 25)[2]

      mz.range <- as.numeric(pg.raw[[idx.pg]][idx.smp, c('mzmin', 'mzmax')])
      idx.swath <- GetIntervalIndex(mz, info.swath$window)
      idx.ms2 <- unname(sapply(scanidx.match[as.character(idx.ms1)], `[`, idx.swath))
      idx.ms2.ext <- unname(sapply(scanidx.match[as.character(idx.ms1.ext)], `[`, idx.swath))

      eic.ms1.pre <- extractEIC(scan.ms1[idx.ms1.ext],
                                matrix(mz.range, ncol = 2),
                                mz)[[1]]
      peaks.ms1.ext <- cbind('scanidx' = idx.ms1.ext,
                             'rt' = scantime.ms1[idx.ms1.ext],
                             'mz' = eic.ms1.pre$mz,
                             'intensity' = eic.ms1.pre$intensity)
      is.keep <- idx.ms1.ext %in% idx.ms1
      peak.ms1 <- peaks.ms1.ext[is.keep, , drop = FALSE]

      peak.ms1.ext.smooth.raw <- peak.ms1.ext.smooth <- SmoothLoess(peaks.ms1.ext, span = 0.1)
      peak.ms1.smooth.raw <- peak.ms1.smooth <- peak.ms1.ext.smooth[is.keep, , drop = FALSE]

      idx.na <- which(is.na(idx.ms2.ext))
      if (length(idx.na) > 0) {
        idx.ms2.ext <- idx.ms2.ext[-idx.na]
        if (length(idx.ms2.ext) == 0) {
          return(NULL)
        }
        idx.ms1.na <- idx.ms1[idx.ms1 == idx.ms1.ext[idx.na]]
        if (length(idx.ms1.na) > 0) {
          peak.ms1.ext.smooth <- peak.ms1.ext.smooth[-idx.ms1.na, , drop = FALSE]
          peaks.ms1.ext <- peaks.ms1.ext[-idx.ms1.na, , drop = FALSE]
          idx.ms1 <- idx.ms1.raw[-idx.ms1.na]
          idx.ms1.ext <-idx.ms1.ext[-idx.ms1.na]
        }
      }

      if (all(peak.ms1.smooth[, 'intensity.s'] == 0)) {
        return(NULL)
      }

      idx.apex.ms1 <- which(idx.ms1.ext == idx.ms1[which.max(peak.ms1.smooth[,'intensity.s'])])

      rt.apex.ms1 <- peak.ms1.ext.smooth[idx.apex.ms1, 'rt']
      if (nrow(peak.ms1.smooth) < 4) {
        return(NULL)
      }

      spec.exp.ext <- scan.ms2[idx.ms2.ext]
      spec.apex <- spec.exp.ext[[idx.apex.ms1]]
      idx.ms2.intfilter <- which(spec.apex[,2] >= int.filter &
                                   spec.apex[, 1] <= mz.max)
      mz.ms2 <- spec.exp.ext[[idx.apex.ms1]][idx.ms2.intfilter, 1]
      if (length(mz.ms2) == 0) {
        return(NULL)
      }
      mz.ms2.range <- t(sapply(mz.ms2, function(mz) {
        PpmRange(mz, ppm.ms2.mtp)
      }))

      eic.ms2.pre <- extractEIC(spec.exp.ext, mz.ms2.range, mz.ms2)
      ms2.eic.ext <- lapply(eic.ms2.pre, function(spec) {
        cbind('scanidx' = idx.ms2.ext,
              'rt' = scantime.ms2[idx.ms2.ext],
              do.call(cbind, spec))
      })

      range.pk.ms1 <- range(which(is.keep))
      info.pk.ms1 <- data.frame('lb' = range.pk.ms1[1],
                                'apex' = idx.apex.ms1,
                                'apex.t' = rt.apex.ms1,
                                'rb' = range.pk.ms1[2],
                                'bl' = pg.raw[[idx.pg]][idx.smp, ]$maxo/pg.raw[[idx.pg]][idx.smp, ]$sn,
                                'simp' = FALSE)
      rownames(info.pk.ms1) <- 0
      nm.smp <- gsub('(?i)\\.mzxml', '', basename(fp.smp))
      spec.decon <- DecoSpectra(idx.pg,nm.smp, ms2.eic.ext, peak.ms1.smooth,
                                 length(idx.ms2.ext),
                                 idx.apex.eic = idx.apex.ms1,
                                 info.pk.ms1 = info.pk.ms1,
                                 peakwidth = peakwidth,
                                 is.PlotEIC.MS2 = is.PlotEIC.MS2,
                                 snthr = snthr, isFWHM = isFWHM,
                                 is.dec.all = is.dec.all,
                                 is.dec.smoothed = is.dec.smoothed)


      spec.decon
    })
    obj@spectra <- info.all
  })
