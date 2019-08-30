import argparse
import matplotlib.pyplot as plt
import numpy as np
import os
import random
import seaborn as sns
#from Bio.Phylo.TreeConstruction import DistanceCalculator
from Bio import AlignIO
from Bio import SeqIO

'''
Эта программа создает выборку следующим образом:
Все последовательности разбиваются на группы по первым 4 знакам genbank ID. 
Затем случайно выбираем группу. Из этой группы все последовательности добавляем в новую выборку, пока последовательностей не станет, сколько нужно.
Если последовательностей меньше, чем 
'''
# input_file - file with sequences in fasta-format used to generate random samples
# output_dir - output directory for generated files
# N_align - number of alignments to generate
# N_seq_max - number of sequences in generated alignment
# threshold - the number of letters that should coincide in 
# GB Accessions of two sequences to join them into one group

def gen_random_sample3(input_file, output_dir, N_align, N_seq_max, threshold):

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
                        number of sequences becomes n_seq_max", required=True)
    args = parser.parse_args()
    if args.algorithm == "group_picking":
        if args.threshold == None:
            print("Please, define the threshold")
        else:
            gen_random_sample3(args.input_file, args.output_dir, args.n_samples, args.n_seq_max, args.threshold)
    else:
        gen_random_single_pick(args.input_file, args.output_dir, args.n_samples, args.n_seq_max)


