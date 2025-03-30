#!/bin/bash

echo "Enter filename"
read input_file

bond1=(23 19) # first atom is fixed
bond2=(23 11) # first atom is fixed

fragment1="15,16,17,18,20,21,22,32" # atoms to move in first fragment (excluding bonded atom)
fragment2="9,10,12,13,14,33" # atoms to move in second fragment (excluding bonded atom)

bond1range=(0 2) # range to diplace first fragement in Angstrom, can include negative (shorter than reference)
bond2range=(0 2) # range to diplace second fragement in Angstrom, can include negative (shorter than reference)

numfiles1=10 # how many files will be generated for bond 1, increments will be evenly distributed based on displacement
numfiles2=10 # how many files will be generated for bond 2, increments will be evenly distributed based on displacement

modred_maker.py "${input_file}" "${bond1[0]}" "${bond1[1]}" "${bond1range[0]}" "${bond1range[1]}" "${numfiles1}" B1 "[${fragment1}]"
for i in B1*; do modred_maker.py "$i" "${bond2[0]}" "${bond2[1]}" "${bond2range[0]}" "${bond2range[1]}" "${numfiles2}" "${i%.gjf}_B2" "[${fragment2}]"; done

