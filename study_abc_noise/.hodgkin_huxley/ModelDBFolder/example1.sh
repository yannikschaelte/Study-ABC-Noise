#!/bin/bash

# EXAMPLE 1:
# Run all models for a constant input (strength of DC input is 7 micro
# amp / cm^2)
# Write out first 30 interspike intervals (in ms)

END=6
for ((i=0;i<=END;i++)); do
    echo running method $i
    ./HH_run $i 100 1E7 0.01 30 7. 0. 0. 0. 0 2 123 > ${i}_output.txt
done
