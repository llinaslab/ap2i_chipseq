from os.path import join, basename, dirname
from os import getcwd
from subprocess import check_output

###############################################################################
# Globals
###############################################################################

ROOT = "/Users/llinaslab/projects/ap2i_chipseq"

LOGS = join(ROOT, "logs")

###############################################################################
# Functions
###############################################################################

# source R script
def source_r(path, file_name):
  return 'Rscript --vanilla --default-packages=methods,stats,utils -e \'setwd("{0}")\' -e \'source("{1}")\''.format(path, file_name)

# source Rmarkdown document
def source_rmd(path, file_name, outdir):
  return 'Rscript --vanilla --default-packages=methods,stats,utils,knitr -e \'setwd("{0}")\' -e \'rmarkdown::render("{1}",output_dir="{2}")\''.format(path, file_name, outdir)

# remove suffix from file path
# don't include period
def rmsuffix(text, suffix):
  if not text.endswith(suffix):
    out = text
  else:
    out = text[:len(text)-len(suffix)]
  return out

# change suffix of file path
# don't include period
def chgsuffix(file, old, new):
  out = rmsuffix(file, old)
  return out + new

# create an empty file
# don't modify it if it already exists
def touch(path):
  with open(path, 'a'):
    os.utime(path, None)
