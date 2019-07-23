# sample_bias

This folder contains scripts that we used to process data for the paper "The effect of sample bias and experimental artifacts on statistical phylogenetic analysis of picornaviruses."

## fetch_date.py

Gets modification and collection dates from GenBank files defined in script. Plots graph and histogram of dates frequency distribution for each GenBank file.

## random_sample.py

Generates random samples with number of sequences *n_seq_max* the following way:

* All sequences in the reference alignment are partitioned into groups by
the first *threshold* characters of the GenBank ID.
* Then the random group is chosen, and all sequences from this group are
added to alignment. This step is repeated until the number of sequences reaches
*n_seq_max*.
