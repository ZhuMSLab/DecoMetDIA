setClass('diaFeature', slots = list(peaks = 'data.frame', peakgroup = 'list',
                                    rt = 'list', samplenames = 'character',
                                    filepaths = 'character', smpclass = 'character'))
setClass('diaSpectra', slots = list(featurePeaks = 'list', featurePeaksSmooth = 'list',
                                    spectraPeaks = 'list', spectraPeakArea = 'list',
                                    ppcScore = 'list', spectra = 'list'))

setGeneric('GetFeatures', function(obj, ...) standardGeneric('GetFeatures'))
setGeneric('GetDIASpectra', function(obj, ...) standardGeneric('GetDIASpectra'))
