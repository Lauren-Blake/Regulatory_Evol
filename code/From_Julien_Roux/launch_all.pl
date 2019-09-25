#!/usr/bin/perl
use strict;
use warnings;

## Julien Roux
## Nov 28, 2013

## Launch script to find species DMRs
#system("echo \"R CMD BATCH --no-restore --no-save find_speciesDMRs.R find_speciesDMRs.Rout\" | qsub -l h_vmem=128g -v PATH -cwd -N speciesDMRs -o dmr.out -e dmr.err");
system("echo \"R CMD BATCH --no-restore --no-save find_speciesSpecificDMRs.R find_speciesSpecificDMRs.Rout\" | qsub -l h_vmem=110g -v PATH -cwd -N speciesDMRs -o dmr.out -e dmr.err");

## Launch script to find tissue DMRs
#system("echo \"R CMD BATCH --no-restore --no-save find_tissueDMRs.R find_tissueDMRs.Rout\" | qsub -l h_vmem=128g -v PATH -cwd -N tissueDMRs -o dmr.out -e dmr.err");

## Launch script to find individual-specific DMRs
#system("echo \"R CMD BATCH --no-restore --no-save find_individualDMRs.R find_individualDMRs.Rout\" | qsub -l h_vmem=128g -v PATH -cwd -N individualDMRs -o dmr.out -e dmr.err");

exit;
