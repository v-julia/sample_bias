# sample_bias

This folder contains scripts that we used to process data for the paper "The effect of sample bias and experimental artifacts on statistical phylogenetic analysis of picornaviruses."

## date_distribution.py

Gets modification and collection dates from GenBank files defined in script. Plots graph and histogram of dates frequency distribution for each GenBank file.

## fasta2nex.py

Converts fasta file to nexus format.

### Usage

```
usage: fasta2nex.py [-h] -input INPUT_FILE

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file in fasta format
```


## random_sample.py

Generates random sample of sequences from file in fasta-format.

Three algorithms are provided:
1) *single_picking* - random picking *n_seq_max* sequences

2) *group_picking* - generates random samples with number of sequences *n_seq_max* the following way:

    * All sequences in the reference alignment are partitioned into groups by
    the first *threshold* characters of the GenBank ID.
    * Then the random group is chosen, and all sequences from this group are
    added to alignment. This step is repeated until the number of sequences reaches
    *n_seq_max*.
3) *smart_picking* -  divides sequences into groups by the first *threshold*
 characters in GenBank Accession. Randomly removes k percent sequences in groups which size exceed m.
    
### Usage

```
usage: random_sample.py [-h] -input INPUT_FILE [-out_dir OUTPUT_DIR] -n_samp
                        N_SAMPLES -n_seq_max N_SEQ_MAX [-threshold THRESHOLD]
                        -alg ALGORITHM

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file
  -out_dir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory to save generated alignments
  -n_samp N_SAMPLES, --n_samples N_SAMPLES
                        Number of random alignments to generate
  -n_seq_max N_SEQ_MAX, --n_seq_max N_SEQ_MAX
                        Number of sequences in generated alignment
  -threshold THRESHOLD, --threshold THRESHOLD
                        Threshold for the first characters to generate random
                        groups
  -alg ALGORITHM, --algorithm ALGORITHM
                        Algorithm of generating random alignment
                        'single_picking' - picking n_seq_max random sequences;
                        'group_picking' - divides sequences into groups by the
                        first ~threshold~ characters in GenBank Accession.
                        Then picks random groups and adds sequences from them
                        to the new alignment till the total number of
                        sequences becomes n_seq_max; 'smart_picking' - divides
                        sequences into groups by the first ~threshold~
                        characters in GenBank Accession. Randomly removes k
                        percent sequences in groups which size exceed m.
```


## add_mut.py

### Usage

```
usage: add_mut.py [-h] -input INPUT_FILE -sn SEQ_NAME -m MUTATIONS -t TYPE_MUT
                  [-pout PATH_OUT]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input xml-file generated in BEAUTi
  -sn SEQ_NAME, --seq_name SEQ_NAME
                        The id of sequence which sequence will be mutated
  -m MUTATIONS, --mutations MUTATIONS
                        List with numbers of mutations which should be added
                        to the nucleotide sequence of seq_name, should be
                        separated by spaces. Example: '1,5,10'
  -t TYPE_MUT, --type_mut TYPE_MUT
                        The type of mutations to introduce: 1 if synonymous, 0
                        if non-synonymous, -1 if any type
  -pout PATH_OUT, --path_out PATH_OUT
                        Output directory. If not defined the output files will
                        be saved in 'mutations' folder in the directory of
                        input file
```

## change_year.py

### Usage

```
usage: change_year.py [-h] -input INPUT_FILE -sn SEQ_NAME -y YEARS
                      [-pout PATH_OUT]

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input xml-file generated in BEAUTi
  -sn SEQ_NAME, --seq_name SEQ_NAME
                        The id of sequence which collection years will be
                        changed
  -y YEARS, --years YEARS
                        The numbers which should be added to the collection
                        date of seq_name, should be separated by spaces
  -pout PATH_OUT, --path_out PATH_OUT
                        Output directory. If not defined the output files will
                        be saved in 'years' folder in the directory of input
                        file
```

## get_rate.py

Derives tree's branch substitution rates from tree-file in nexus format (output of TreeAnnotator program). 
Saves derived rates and logarithms of rates into "rates.txt" file located in output_dir.
 Plots histograms of rates and log rates showing mean rate and mean + (1,2,3)sd as vertical lines.

### Usage

```
usage: get_rate.py [-h] -input INPUT_FILE [-out_dir OUTPUT_DIR] -t TITLE -f
                   FORMAT

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file with tree-file generated by TreeAnnotator
  -out_dir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory to save generated alignments
  -t TITLE, --title TITLE
                        Title for plot
  -f FORMAT, --format FORMAT
                        format of figure
```

## plot_hist_rates.py

Plots histogram of rates, shows min rate, standard deviations as vertical lines on hist. Saves the figure with histogram in the *filename*

### Usage

```
usage: plot_hist_rates.py [-h] -input INPUT_FILE [-out_dir OUTPUT_DIR] -t
                          TITLE -f FORMAT

optional arguments:
  -h, --help            show this help message and exit
  -input INPUT_FILE, --input_file INPUT_FILE
                        Input file with tree-file generated by TreeAnnotator
  -out_dir OUTPUT_DIR, --output_dir OUTPUT_DIR
                        Output directory to save generated figures
  -t TITLE, --title TITLE
                        Title for plot
  -f FORMAT, --format FORMAT
                        format of figure
```

## rates_stat.py


### Usage

```
usage: rates_stat.py [-h] -orig ORIG_TREE_PATH -changed FOLDER_CHANGED_PATH -m
                     MODE

optional arguments:
  -h, --help            show this help message and exit
  -orig ORIG_TREE_PATH, --orig_tree_path ORIG_TREE_PATH
                        Path to original MCC tree generated using
                        TreeAnnotator
  -changed FOLDER_CHANGED_PATH, --folder_changed_path FOLDER_CHANGED_PATH
                        Path to the folder with trees were sequences'
                        collection years or nucleotide sequences are changed
  -m MODE, --mode MODE  Type of sequences' change:'mutations' if mutations
                        were introduced into nucleotide sequence of seq_name,
                        'years' if some values were added to the collection
                        date of seq_name
```

## analyse_log.R

Derives median and HPDs from log-files...