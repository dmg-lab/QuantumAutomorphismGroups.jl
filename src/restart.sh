#!/bin/bash

# This script is used to restart a computation that was interrupted


if [ $# -lt 1 ]; then
    echo "Usage: $0 <computation name>"
    exit 1
fi


filename="$1"

julia -e "using Pkg; Pkg.activate(\"../\"); include(\"./quantumMatroid.jl\"); restartComputation(\"$filename\")"


