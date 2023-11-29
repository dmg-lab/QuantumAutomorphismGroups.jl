#!/bin/bash

# Hardcoded strings
path_to_csv="../data/data_table.csv"
path_to_tex="../../matroid_quantum_automorphism.tex"

# Path to the Julia script
julia_script="create_table.jl"


julia "$julia_script" "$string1" "$string2"









