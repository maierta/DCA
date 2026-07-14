#!/bin/bash

# Run all tc calculations on a local workstation (no scheduler) using MPI.
# The threaded solver is configured for 1 walker + 1 accumulator per MPI rank.

RUNNER="mpiexec -n 4"
RUN_DCA="${RUNNER} ../../applications/dca/main_dca"
date

$RUN_DCA ./T=1/input_sp.json
$RUN_DCA ./T=0.75/input_sp.json
$RUN_DCA ./T=0.5/input_sp.json
$RUN_DCA ./T=0.25/input_sp.json
$RUN_DCA ./T=0.125/input_sp.json
$RUN_DCA ./T=0.1/input_sp.json
$RUN_DCA ./T=0.1/input_tp.json
$RUN_DCA ./T=0.09/input_sp.json
$RUN_DCA ./T=0.09/input_tp.json
$RUN_DCA ./T=0.08/input_sp.json
$RUN_DCA ./T=0.08/input_tp.json
$RUN_DCA ./T=0.07/input_sp.json
$RUN_DCA ./T=0.07/input_tp.json

date
