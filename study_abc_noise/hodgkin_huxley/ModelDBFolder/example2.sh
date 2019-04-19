#!/bin/bash

# EXAMPLE 2:
# Run all models in voltage clamp (voltage clamp to 20 mV) for 100 ms
# (1E4 time steps with 0.01ms step size)
# Write out time, voltage, proportion of open Na channels and
# proportion of open K channels

END=6
for ((i=0;i<=END;i++)); do
    echo running method $i
    ./HH_run $i 100 1E4 0.01 100 20. 0. 0. 0. 1 1 123 > ${i}_output.txt
done
