#' @title DecoMetDIA: Deconvolution of Multiplexed MS/MS Spectra for Metabolite
#'     Identification in SWATH-MS based Untargeted Metabolomics
#' @description DecoMetDIA picks peaks from SWATH acquired data with XCMS and
#'     extracts MS2 EICs under related SWATH window for each feature of different
#'     sample. Deconvolution algorithms are applied to MS2 EICs for each feature
#'     and thus MS2 spectra are generated. The features are identified by matching
#'     the obtained spectra with librarial ones.
#' @param d.in Working directory, where the well orgnized data are stored
#' @param d.out Directory for saving results
#' @param nSlaves No. of threads for data processing
#' @param peakwidth Chromatographic peak width, given as range (min,max) in seconds
#'     (\code{xcms})
#' @param sn Signal to noise ratio cutoff for peak detection using \code{xcms}
#' @param ppm.pd Maxmial tolerated m/z deviation in consecutive scans, in ppm
#'     for MS1 peak detection
#' @param polarity Polarity setup for ionization
#' @param lc LC column used for chromatographic seperation
#' @param minfrac minimum fraction of samples necessary in at least one of the
#'     sample groups for it to be a valid group when group peaks with 'denisty'
#'     method in \code{xcms}
#' @param bw bandwidth (standard deviation or half width at half maximum)
#'     of gaussian smoothing kernel to apply to the peak density chromatogram when
#'     group peaks with 'denisty' method in \code{xcms}
#' @param mzwid width of overlapping m/z slices to use for creating peak density
#'     chromatograms and grouping peaks across samples when group peaks with
#'     'denisty' method in \code{xcms}
#' @param is.ms1.filter logical value for indicating if filter isotopic or low
#'     abundant peaks
#' @param thr.ms1.filter MS1 peak filter threshold for filtering low abundant peaks
#' @param filed.ms1.filter MS1 peak filter method (maxo:intensity, into:raw peak area or
#'     intb:base-line corrected peak area)
#' @param is.plot.eic.feature logical value for indicating if plot EICs of features
#' @param is.check.eic logical value for indicating if plot peak boundaries
#'     when plotting feature EICs
#' @param is.PlotEIC.MS2 logical value for indicating if plot EICs of fragment ions
#' @param snthr.ms2 signal to noise threshold for MS2 peak detection
#' @param fp.swath.setup file path for SWATH setup file (.csv format recording swath windows)
#' @param ppm.ms2.mtp ppm threshold for EIC extraction of fragment ions in multiplexed spectrum
#' @param int.filter intensity threshold for removing low abundant fragment ions
#'     in multiplexed spectrum
#' @param minfrac.vote = 0.5,
#' @param is.output.spec.smp logical value for indicating if output deconvoluted spectra
#'     for each sample
#' @param is.span.constant logical value for indicating if using constant span
#'     when smoothing EICs
#' @param rerun logical value for indicating if rerun the already processed peak
#'     detection results
#' @param rerun.ms2 logical value for indicating if rerun the already prcoessed
#'     MS2 data analysis process
#' @param is.ms1.only logical value for indicating if only detect peaks
#' @param is.dec.all logical value for indicating if deconvolute all EICs, if not,
#'     the 'simple' peaks are ignored
#' @param is.dec.smoothed logical value for indicating if using smoothed peaks as
#'     model peaks.
#' @param is.norm.consensus logical value for indicating if normalize intensities
#'     to 1 when consensusing spectra cross samples.

DecoMetDIA <- function(d.in = '.',
                       d.out='DecoMetDIA_Result',
                       nSlaves = 6,
                       # peak detection setup
                       peakwidth = c(5, 30),
                       sn = 6,
                       ppm.pd = 25,
                       polarity = c('positive', 'negative'),
                       lc = c('HILIC', 'RP'),
                       minfrac = 0.5,
                       bw = 5,
                       mzwid = 0.015,
                       is.ms1.filter = FALSE,
                       thr.ms1.filter = NULL,
                       filed.ms1.filter = c('maxo', 'into', 'intb'),
                       is.plot.eic.feature = F,
                       is.check.eic = FALSE,
                       # SWATH MS2 setup
                       is.PlotEIC.MS2 = FALSE,
                       snthr.ms2 = 5,
                       fp.swath.setup ='../SWATHsetup.csv',
                       ppm.ms2.mtp = 15,
                       int.filter = 50,
                       minfrac.vote = 0.5,
                       is.plot.eic.spec = F,
                       is.output.spec.smp = TRUE,
                       is.span.constant = FALSE,
                       rerun = FALSE,
                       rerun.ms2 = FALSE,
                       is.ms1.only = FALSE,
                       is.dec.all = FALSE,
                       is.dec.smoothed = FALSE,
                       is.norm.consensus = TRUE) {
  span <- 0.3
  mz.diff <- 0.01
  wd0 <- getwd()
  setwd(d.in)
  polarity <- match.arg(polarity)
  lc <- match.arg(lc)
  filed.ms1.filter <- match.arg(filed.ms1.filter)

  d.out <- CreateResultDir(d.out, is.overwrite = FALSE)
  files_mzxml <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzxml$')
  files_mzml <- list.files(d.in, recursive = TRUE, full.names = TRUE, pattern = '(?i)mzml$')
  files <- c(files_mzxml, files_mzml)

  nSlaves <- min(parallel::detectCores() - 1, nSlaves, length(files))
  bpparam <- SnowParam(workers = nSlaves, type = 'SOCK')
  CatLineSeparator('Detecting and aligning features ...')
  fn.skip <- 'dia.feature.RData'
  if ((!rerun) & file.exists(fn.skip)) {
    cat('using existing results:', fn.skip, '...\n')
    load(fn.skip)
  } else {
    dia.feature <- GetFeatures(new('diaFeature'), files,
                               polarity = polarity, lc = lc,
                               ppm.pd = ppm.pd, sn = sn, peakwidth = peakwidth,
                               mzdiff = mz.diff, minfrac = minfrac,
                               method.rector = 'obiwarp',
                               is.ms1.filter = is.ms1.filter,
                               thr.ms1.filter = thr.ms1.filter,
                               filed.ms1.filter = filed.ms1.filter,
                               is.plot.eic.feature = is.plot.eic.feature,
                               is.check.eic = is.check.eic,
                               nSlaves = nSlaves, bw = bw, mzwid = mzwid)

    save(dia.feature, file = fn.skip)
  }

  if(is.ms1.only) {
    return(NULL)
  }
  dt.peaktable <- read.csv('PeakTable-annotated_filtered.csv', stringsAsFactors = FALSE)

  cat('ALL peak detection work are finished!!\n')
  t.start <- Sys.time()

  CatLineSeparator('Start processing SWATH data...')
  spectra.all <- bplapply(files, function(fp.smp) {
    suppressMessages(require(DecoMetDIA))
    spec.dia <- GetDIASpectra(new('diaSpectra'), fp.smp, dia.feature,
                              dt.peaktable=dt.peaktable,
                              d.spec = file.path(d.out, 'SpecDIA'),
                              peakwidth = peakwidth, snthr = snthr.ms2,
                              int.filter = int.filter, ppm.ms2.mtp = ppm.ms2.mtp,
                              is.PlotEIC.MS2 = is.PlotEIC.MS2,
                              span = span,
                              rerun = rerun.ms2,
                              is.span.constant = is.span.constant,
                              fp.swath = fp.swath.setup,
                              is.dec.all = is.dec.all,
                              is.dec.smoothed = is.dec.smoothed)
    return(spec.dia)
  },BPPARAM = bpparam)
  save(spectra.all, file = 'dataSpectraDIA.Rda')


  spec.consensus <- ConsensusSpec(spectra.all, dia.feature@smpclass,
                                  nrow(dt.peaktable), minfrac.vote,
                                  is.norm.consensus = is.norm.consensus)
  names(spec.consensus) <- dt.peaktable$name
  spec.consensus <- spec.consensus[!sapply(spec.consensus, is.null)]
  sink(file.path(d.out, 'spec_consensus.msp'))
  for (nm.ft in names(spec.consensus)) {
    for (nm.grp in unique(dia.feature@smpclass)) {
      spec <- spec.consensus[[nm.ft]][[nm.grp]]
      if (!is.null(spec)) {
        OutputSpec(spec, which(dt.peaktable$name == nm.ft), dt.peaktable, grp = nm.grp)
      }
    }
  }
  sink()

  if (is.output.spec.smp) {
    names(spectra.all) <- files
    for (fp.smp in files) {
      idx.smp <- which(dia.feature@filepaths %in% fp.smp)
      if (length(idx.smp) == 0) {
        if (file.exists(fp.smp)) {
          idx.smp <- which(basename(dia.feature@filepaths) %in% basename(fp.smp))
        } else {
          stop(paste0("file not exist: \n  ", fp.smp))
        }
      }

      d.spec <- file.path(file.path(d.out, 'SpecDIA'), dia.feature@smpclass[idx.smp],  gsub('(?i).mzxml', '', basename(fp.smp)))
      if (!dir.exists(d.spec)) {
        dir.create(d.spec, recursive = TRUE)
      }
      info.all <- spectra.all[[fp.smp]]
      for(nm.spec in do.call(rbind, lapply(info.all, names))[1,]) {
        fn.spec <- file.path(d.spec, paste0(nm.spec ,'.msp'))
        sink(fn.spec)
        for (idx.pg in seq_along(info.all)) {
          if (!is.null(info.all[[idx.pg]])) {
            spec.output <- info.all[[idx.pg]][[nm.spec]]
            if (!is.null(spec.output)) {
            }
            OutputSpec(spec.output, idx.pg, dt.peaktable)
          }
        }
        sink()
      }
    }
  }

  t.end <- Sys.time()
  cat('ALL WORK DONE!!')
  cat(t.start, '\n', t.end, '\n')
  setwd(wd0)
  return(file.path(d.out, 'SpecDIA'))
}

