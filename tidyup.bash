#!/bin/bash

#remove the unnecessary middle of the file name
#REMOVE echo BEFORE EXECUTION

for f in *.fastq; do echo mv "$f" "${f/220329_2022_03_29_Miseq_Run8_Nikola_Manon_Naina.220415.MiSeq.FCB.lane1.gcap_20_01./}"; done

#remove GC121272_ from the beginning of the file name
#REMOVE echo BEFORE EXECUTION

for f in *.fastq; do echo mv "$f" "${f/GC121272_/}"; done #The GC number may need to be changed depending on the number of the run
