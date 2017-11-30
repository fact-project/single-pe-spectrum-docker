#!/bin/bash

# extract singles for mc
drsfile=${drsfile:-"stdMcDrsFile.drs.fits.gz"}
infile=${infile:-"/input_data/20120723_"}
outpath_spectra=${outpath_spectra:-"/output_data/single_pe_spectra.root"}
outpath_fit=${outpath_fit:-"/output_data/single_pe_spectra_fit.root"}

first_run=${first_run:-9}
last_run=${last_run:-$first_run}

ROOT_HIST=0 root -q -b -l extract_singles.C+\(\""$infile\"",$first_run,$last_run,\""$drsfile\"",\""${outpath_spectra}\""\)

ROOT_HIST=0 root -q -b -l fit_spectra.C+\(\""$outpath_spectra\"",\""$outpath_fit\""\)
