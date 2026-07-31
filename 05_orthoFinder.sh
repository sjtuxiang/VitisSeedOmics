#!/bin/bash

orthofinder -f data -S diamond -M msa -A mafft -T iqtree3 -t 64 -a 32 -o results
