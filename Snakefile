# vim: set ft=python:
shell.prefix("set -eo pipefail; ")

from os.path import join, basename, dirname, exists
from os import getcwd, makedirs
from subprocess import check_output
from glob import glob

from snakemake.utils import min_version

min_version("3.4.1")

# Globals ---------------------------------------------------------------------

include: "config.py"

# Job Handlers ----------------------------------------------------------------

onsuccess:
  print("Workflow finished, without any errors!")

#onerror:
#  shell("send_email.py -t philippross369@gmail.com -s 'Snakefile error' -b {log}")

# Rules -----------------------------------------------------------------------

rule all:
  input:
    join(ROOT, "results/final/report.html")


# Subworkflows
subworkflow alignment:
  workdir: join(ROOT, "workflows/alignment")

subworkflow correlations:
  workdir: join(ROOT, "workflows/correlations")

subworkflow peaks:
  workdir: join(ROOT, "workflows/peaks")

subworkflow annotation:
  workdir: join(ROOT, "workflows/annotation")

# Final report
rule report:
  input:
    alignment("../../results/alignment/report.html"),
    correlations("../../results/correlations/report.html"),
    peaks("../../results/peaks/report.html"),
    annotation("../../results/annotation/report.html"),
    script = "report.Rmd"
  output:
    join(ROOT, "results/final/report.html")
  run:
    shell(source_rmd(getcwd(), input.script, join(ROOT, "results/final/")))

# Clean up final report
rule clean:
  run:
    shell("rm -rf "+join(ROOT, "results/final/report.html"))
