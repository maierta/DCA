#!/bin/bash

# Run analysis on a local workstation (no scheduler) using MPI.
# The threaded solver is configured for 1 walker + 1 accumulator per MPI rank.

RUNNER="mpiexec -n 4"
RUN_DCA="${RUNNER} ../../applications/analysis/main_analysis"

date

$RUN_DCA ./T=0.1/input_tp.json
$RUN_DCA ./T=0.09/input_tp.json
$RUN_DCA ./T=0.08/input_tp.json
$RUN_DCA ./T=0.07/input_tp.json

date
