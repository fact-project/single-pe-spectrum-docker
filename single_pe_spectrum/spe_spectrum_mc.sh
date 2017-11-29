#!/bin/bash

# extract singles for mc
drsfile=${drsfile:-"stdMcDrsFile.drs.fits.gz"}
suffix=${suffix:-".001_P_MonteCarlo000_Events.fits.gz"}
infolder=${infolder:-"/input_data/"}

first_run=${first_run:-0}
last_run=${last_run:-first_run}

outpath_spectra=${outpath_spectra:-"/output_data/single_pe_spectra.root"}
outpath_fit=${outpath_fit:-"/output_data/single_pe_spectra_fit.root"}

ROOT_HIST=0 root -q -b -l extract_singles_mc.C+\(\""$infolder\"",\""$suffix\"",${first_run},${last_run},\""$drsfile\"",\""$outpath_spectra\""\)

ROOT_HIST=0 root -q -b -l fit_spectra_mc.C+\(\""$outpath_spectra\"",\""$outpath_fit\""\)
