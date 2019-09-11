import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import seaborn as sns
#from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO



def gen_random_sample3(input_file, output_dir, N_align, N_seq_max, threshold):

    '''
    Divides all sequences into groups by the first threshold characters in GenBank accession number.
    Then randomly chooses the groups and adds all sequences from the group to the new alignment 
    till the number of sequences exceeds the maximum number defined by user.

    Input:
        input_file - file with sequences in fasta-format used to g
        output_dir - output directory for generated files
        N_align - number of alignments to generate
        N_seq_max - number of sequences in generated alignment
        threshold - the number of letters that should coincide in 
        GB Accessions of two sequences to join them into one group

    Ouput:
        Saves new file with sequences in fasta-format in output directory
    '''

    input_file = '/'.join(input_file.split('\\')) # change path to linux style
    aln = list(SeqIO.parse(open(input_file), 'fasta'))

    #Dictionary with groups which gb ids start with the same string
    ids_dict = {}

    for rec in aln:
        #print(rec)
        group_id = rec.id[0:threshold] #id of new group
        if group_id in ids_dict.keys():
            ids_dict[group_id].append(rec.id)
        else:
            ids_dict[group_id] = [rec.id]
    
    '''
    #this code is for printing length of group and Genbank accessions of sequences that belong to each group
    l_gsize = [] #list with group sizes
    for key in sorted(ids_dict.keys()):
        #print(key + ','+ str(len(ids_dict[key])))
        l_gsize.append(len(ids_dict[key]))
        #for entr in ids_dict[key]:
        #    print('\t'+entr)
    #print(len(ids_dict.keys()))

    #print(sorted(l_gsize))
    sns.set()
    plt.hist(l_gsize)
    plt.title("Distribution of group sizes")
    plt.xlabel("Size of group")
    plt.ylabel("Number of sequences")
    plt.show()
    '''

    for j in range(N_align):
        
        print('Generating {} file'.format(str(j+1)))
        
        #current size of random sample
        Sample_size = 0 
        #list with ids in random sample
        sample_ids = [] 
        #list with random groups that has been chosen in loop
        random_groups = []
        #alignment objects of sequences that will be added to random sample
        list_aln_obj = []
        
        while Sample_size<N_seq_max:
            print("Current size of alignment {}".format(Sample_size))
            #choose random group
            random_group = random.choice(list(ids_dict.keys()))

            #If the group has already been chosen we continue our search
            if random_group in random_groups:
                continue
            else:
                print(random_group, len(ids_dict[random_group]))
                random_groups.append(random_group)
                #add all sequences from the group to random sample
                
                for id in ids_dict[random_group]:
                    if len(sample_ids) < N_seq_max:
                        sample_ids.append(id)
                    else:
                        break
            Sample_size = len(sample_ids)
        print("Final sample size {}".format(Sample_size))


        for rec in aln:

            if rec.id  in sample_ids:
                list_aln_obj.append(rec)
                

        output_file_n = ".".join(input_file.split("/")[-1].split(".")[:-1])+"_3random_"+str(j)+".fasta"
        if output_dir == None:
            output_dir = "/".join(input_file.split("/")[:-1]) + "/"
        output = output_dir + output_file_n
        AlignIO.write(AlignIO.MultipleSeqAlignment(list_aln_obj), output, "fasta")




def gen_random_single_pick(input_file, output_dir, N_align, N_seq_max):

    '''
    This function generates N_align random alignments consisting of N_seq_max sequences from input_file fasta-file,
    then it saves the new alignments in ouput_dir

    Input arguments:
        input_file, type str - path to the file with sequences in fasta-format used to generate random samples
        output_dir, type str - output directory for generated files
        N_align, type int - number of alignments to generate
        N_seq_max, type int - number of sequences in generated alignment
    Ouput:
        Saves new file with sequences in fasta-format in output directory
    '''

    input_file = '/'.join(input_file.split('\\')) # change path to linux style
    aln = list(SeqIO.parse(open(input_file), 'fasta')) # list with records


    for j in range(N_align):
        
        print('Generating {} file'.format(str(j+1)))
        
        #current size of random sample
        Sample_size = 0 
        #list with ids in random sample
        sample_ids = [] 
        #list with random groups that has been chosen in loop
        
        #list with random records
        list_aln_obj = random.sample(aln, k = N_seq_max)
        #name of output file
        output_file_n = ".".join(input_file.split("/")[-1].split(".")[:-1])+"_singlerandom_"+str(j)+".fasta"
        if output_dir == None:
            output_dir = "/".join(input_file.split("/")[:-1]) + "/"
        output = output_dir + output_file_n
        AlignIO.write(AlignIO.MultipleSeqAlignment(list_aln_obj), output, "fasta")


def gen_rnd_smart(input_file, min_num_seq, N_align, output_dir):
    '''
    Divides all sequences into groups by the first 5 characters in GenBank accession number
    Randomly removes k% sequences in groups which size exceed m. k and m are defined by user.
    Saves the resulting sequences in fasta format
    fasta - name of file with sequences in fasta_format
    Input:
        input_file -str - path to file with sequences in fasta-format
        min_num_seq - int - minimal number of sequences to add from the groups which size exceeds m
        N_align - int - number of alignments to generate
        output_dir - str - path to output directory
    Ouput:
        Saves new file with sequences in fasta-format in output directory

    '''
    
    input_file = '/'.join(input_file.split('\\')) # change path to linux style

    #creates list of sequence names (the order of sequences is significant)
    
    records = list(SeqIO.parse(open(input_file), "fasta"))
    seq_names_ordered = []
    for rec in records:
        seq_names_ordered.append(rec.id)


    #dictionary with sequences
    record_dict = SeqIO.to_dict(records)

    #creating dictionary with names of groups (first five chars of GB ids) as keys and list of sequences' names in this
    #group as values
    GB_groups = dict()
    #this iteration includes ref_name
    for key in record_dict.keys():
            if key[:5] not in GB_groups:
                    GB_groups[key[:5]] = []
            GB_groups[key[:5]].append(key)

    #creating dictionary and list of number of seqs in each groups for hist
    GB_groups_num = dict()
    nums_for_hist = []
    for key in GB_groups.keys():
            GB_groups_num[key] = len(GB_groups[key])
            nums_for_hist.append(len(GB_groups[key]))

    #creating histogram
    #plt.figure(figsize=(15,10))
    plt.hist(nums_for_hist, bins = range(min(nums_for_hist),max(nums_for_hist),1))
    plt.title('Group size distribution', size = "25")
    plt.xlabel('Size of the group ', size = "20")
    plt.ylabel('Number of groups', size = "20")
    plt.show()


    print("Type in the maximal group size; then press Enter")
    max_size = int(input())
    print("Type in the percent of sequences to remove from the groups with size larger than adjusted; then press Enter")
    per_of_seq = float(input())

    if max_size+1 < min_num_seq:
        print("The minimal number of sequences to get from large groups is higher than \
        the maximal group size. The minimal number of sequences will be set equal to \
        the (maximal group size + 1)")
        min_num_seq = max_size+1

    # generating N_align alignments
    for i in N_align:
        # creating list with names of sequences in final sample
        final_sample_ids = []
        # number of sequences in the new alignment
        num = 0
        for value in GB_groups.values():
                if len(value) <= max_size:
                    num = num + len(value)
                    for i in value:
                            final_sample_ids.append(i)
                else:
                    
                    random_sample = random.sample(value, max(min_num_seq,int(len(value)*per_of_seq/100)))
                    num = num + (len(random_sample))
                    for seq in random_sample:
                        final_sample_ids.append(seq)
        print('Number of sequences in resulting alignment {}'.format(num))



        #creating list with final sequences
        final_seq_list = []
        for name in seq_names_ordered:
            if name in final_sample_ids:
                #print(record)
                final_seq_list.append(record_dict[name])

        #name of output file
        output_file_n = os.path.splitext(os.path.split(input_file)[1])[0] + "_random_" + str(max_size)+ '_' + str(per_of_seq) + '_' + str(i) + ".fasta"
        
        if output_dir == None:
            output_dir = os.path.split(input_file)[0]
        else:
            output_dir = '/'.join(input_file.split('\\'))
        output = output_dir + output_file_n
        #writing sequences into file
        SeqIO.write(final_seq_list, output, "fasta")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input file", required=True)
    parser.add_argument("-out_dir", "--output_dir", type=str,
                        help="Output directory to save generated alignments")
    parser.add_argument("-n_samp", "--n_samples", type=int,
                        help="Number of random alignments to generate", required=True)
    parser.add_argument("-n_seq_max", "--n_seq_max", type=int,
                        help="Number of sequences in generated alignment", required=True)
    parser.add_argument("-threshold", "--threshold", type=int,
                        help="Threshold for the first characters to generate random groups")
    parser.add_argument("-alg", "--algorithm", type=str,
                        help="Algorithm of generating random alignment\n \
                        'single_picking' - picking n_seq_max random sequences; \
                        'group_picking' - divides sequences into groups by the first ~threshold~ \
                        characters in GenBank Accession. Then picks random groups and \
                        adds sequences from them to the new alignment till the total \
                        number of sequences becomes n_seq_max \n \
                        'smart_picking' - divides sequences into groups by the first ~threshold~ \
                        characters in GenBank Accession. Randomly removes k% sequences in groups \
                        which size exceed m. ", required=True)
    args = parser.parse_args()
    if args.algorithm == "group_picking":
        if args.threshold == None:
            print("Please, define the threshold")
        else:
            gen_random_sample3(args.input_file, args.output_dir, args.n_samples, args.n_seq_max, args.threshold)
    else:
        if args.algorithm == "single_picking":
            gen_random_single_pick(args.input_file, args.output_dir, args.n_samples, args.n_seq_max)
        else:
            if args.algorithm == "smart_picking":
                print("Type in the minimal number of sequences to be added from large groups; then press Enter")
                min_num_seq = int(input())
                gen_rnd_smart(args.input_file, min_num_seq, args.n_samples, args.output_dir)