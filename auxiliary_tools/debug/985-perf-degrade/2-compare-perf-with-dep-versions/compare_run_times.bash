#!/bin/bash

# Usage: bash /gpfs/fs1/home/ac.tvo/E3SM-Project/e3sm_diags/auxiliary_tools/debug/985-perf-degrade/2-compare-perf-with-dep-versions/compare_run_times.bash

set -e

SCRIPT="auxiliary_tools/debug/985-perf-degrade/2-compare-perf-with-dep-versions/mvce.py"

# Function to run the script in a given conda env and extract run time
run_and_get_time() {
    local ENV_NAME=$1
    echo "Activating $ENV_NAME..."
    source $(conda info --base)/etc/profile.d/conda.sh
    conda activate "$ENV_NAME"
    echo "Running $SCRIPT in $ENV_NAME..."
    # Capture output
    OUTPUT=$(python "$SCRIPT")
    # Extract run time from output (assumes line: "Total run time: X.XX seconds")
    RUNTIME=$(echo "$OUTPUT" | grep "Total run time:" | awk '{print $(NF-1)}')
    echo "$ENV_NAME run time: $RUNTIME seconds"
    conda deactivate
    echo "$RUNTIME"
}

# Run in e3sm_diags_dev_985
TIME1=$(run_and_get_time "e3sm_diags_dev_985")

# Run in e3sm_diags_dev_985_pinned
TIME2=$(run_and_get_time "e3sm_diags_dev_985_pinned")

echo
echo "Summary of run times:"
echo "e3sm_diags_dev_985:         $TIME1 seconds"
echo "e3sm_diags_dev_985_pinned:  $TIME2 seconds"