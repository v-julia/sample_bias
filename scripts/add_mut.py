import argparse
import copy
from lxml import etree
import os
import random
import sys
import statistics

global aa_dict

aa_dict = {

"TTT" : "F", "TTC" : "F", "TTA" : "L", "TTG" : "L",
"CTT" : "L", "CTC" : "L", "CTA" : "L", "CTG" : "L",
"ATT" : "I", "ATC" : "I", "ATA" : "I", "ATG" : "M",
"GTT" : "V", "GTC" : "V", "GTA" : "V", "GTG" : "V",
"TCT" : "S", "TCC" : "S", "TCA" : "S", "TCG" : "S",
"CCT" : "P", "CCC" : "P", "CCA" : "P", "CCG" : "P",
"ACT" : "T", "ACC" : "T", "ACA" : "T", "ACG" : "T",
"GCT" : "A", "GCC" : "A", "GCA" : "A", "GCG" : "A",
"TAT" : "Y", "TAC" : "Y", "TAA" : "X", "TAG" : "X",
"CAT" : "H", "CAC" : "H", "CAA" : "Q", "CAG" : "Q",
"AAT" : "N", "AAC" : "N", "AAA" : "K", "AAG" : "K",
"GAT" : "D", "GAC" : "D", "GAA" : "E", "GAG" : "E",
"TGT" : "C", "TGC" : "C", "TGA" : "X", "TGG" : "W",
"CGT" : "R", "CGC" : "R", "CGA" : "R", "CGG" : "R",
"AGT" : "S", "AGC" : "S", "AGA" : "R", "AGG" : "R",
"GGT" : "G", "GGC" : "G", "GGA" : "G", "GGG" : "G"
}


'''
Introduces N synonymous or non-synonymous mutations into sequence seq

Input:
    seq - str - nucleotide sequence
    N - int - number of mutations to add to the sequence seq
    syn - int (0 or 1) - the type of mutations (synonymous if syn==1)
Output:
    new_seq - str- sequence seq with N mutations
'''

def add_mut(seq, N, syn='0'):

    # nucleotides
    nuc = ['A', 'T', 'G', 'C']

    # list of codons of seq
    codons = []
    i = 0
    while i<len(seq):
        codons.append(seq[i:i+3])
        i+=3
    print(codons)
    
    if syn =='1':
        #list with positions of codons that have been used
        used_list = []
        i=0
        while i<N:
            n = random.randint(0,len(codons)-1)
            if '-' in codons[n]:
                print(codons[n])
                used_list.append(n)
            if n not in used_list:
                #print(n)
                #print(codons[n])
                nuc1 = nuc[:]
                nuc1.remove(codons[n][2])
                #random nucleotide
                rn = random.choice(nuc1)
                #print(rn)
                new_cod = codons[n][0:2]+rn
                #print(new_cod)
                if aa_dict[new_cod] == aa_dict[codons[n]]:
                    print(n)
                    print(codons[n])
                    codons[n] = new_cod
                    print(codons[n])
                    used_list.append(n)
                    i+=1
                else:
                    continue
        new_seq = ''.join(codons)
    else:
        print('nonsynonymous mutations')
        #numbers of positions of seq that have already been changed
        used_pos = []
        new_seq = seq[:]
        i=0

        while i<N:
            
            n = random.randint(0,len(new_seq)-1)
            print(n)
            if new_seq[n] == '-':
                used_pos.append(n)
            if n not in used_pos:
                nuc1 = nuc[:]
                nuc1.remove(new_seq[n])
                old_codon = seq[(n//3)*3:(n//3)*3+3]
                old_codon1 = codons[(n//3)]
                print(old_codon, old_codon1)
                rn = random.randint(0,len(nuc1)-1)
                #print(n%3)
                new_codon = old_codon[:n%3]+nuc1[rn] + old_codon[(n%3)+1:]
                #print(new_seq[n])
                #print(nuc1[rn])                
                if aa_dict[new_codon] != 'X' and aa_dict[new_codon]!=aa_dict[old_codon]:
                    new_seq = new_seq[:n] + nuc1[rn]+new_seq[n+1:]
                    used_pos.append(n)
                    #print(new_seq)
                    print(new_codon)
                    i+=1
                else:
                    continue

    print(seq)
    print(new_seq)
    return(new_seq)

'''
Introduces any N  mutations into sequence seq

Input:
    seq - str - nucleotide sequence
    N - int - number of mutations to add to the sequence seq
Output:
    new_seq - str- sequence seq with N mutations
'''
#N - number of mutations
#seq - sequence as str object
def add_any_mut(seq, N):
    
    nuc = ['A', 'T', 'G', 'C']
    #serial numbers of positions of seq that have already been changed
    used_pos = []
    seq = seq.upper()
    new_seq = seq[:]
    
    codons = []
    i = 0
    while i<len(seq):
        codons.append(seq[i:i+3])
        i+=3
    print(codons)
    
    i=0
    while i<N:
        # position to introduce the mutation
        n = random.randint(0,len(new_seq)-1)
        print(n)
        # skip if gap
        if new_seq[n] == '-':
            used_pos.append(n)
        if n not in used_pos:
            nuc1 = nuc[:]
            nuc1.remove(new_seq[n])
            old_codon = seq[(n//3)*3:(n//3)*3+3]
            old_codon1 = codons[(n//3)]
            print(old_codon, old_codon1)
            # random nucleotide
            rn = random.randint(0,len(nuc1)-1)
            #print(n%3)
            new_codon = old_codon[:n%3]+nuc1[rn] + old_codon[(n%3)+1:]
            #print(new_seq[n])
            #print(nuc1[rn])                
            if aa_dict[new_codon] != 'x':
                new_seq = new_seq[:n] + nuc1[rn]+new_seq[n+1:]
                used_pos.append(n)
                print(new_seq)
                print(new_codon)
                i+=1
            else:
                continue
    return(new_seq)

'''
Input:
    input_xml - str - path to xml-file generated using BEAUYi
    seq_name - str - id of sequence]
    list_mutations - list with str object - list with number of mutations to introduce
    syn - type of mutations (syn = 1 - synonymous, syn = 0 - non-synonymous, -1 - any)
    path_out - path to output directory
Output:
    Creates new xml-files where sequence seq_name has N (from list_mutations)
    synonymous(syn=1), non-syn (syn=0) or any (syn=-1) mutations.
'''

def write_xml(input_xml, seq_name, list_mutations, syn, path_out):

    # xml tree of original file generated in beauty
    tree = etree.parse(input_xml)
    # root element of xml tree
    root = tree.getroot()
    
    # alignment section
    alignment = root.findall("alignment")[0]
    
    #seial sequence number counter
    k = 0
    #iteration over sequences in alignment section
    for seqel in alignment.iter('sequence'):
        #print(seqel.tag)
        #print((seqel.itertext()))
        #sequence id
        cur_seq_name = (seqel.findall('taxon')[0]).attrib['idref']
        print(seq_name)
        if seq_name == cur_seq_name:
            t = etree.tostring(seqel, encoding='unicode')
            #print(t)
            #nucleotide sequence of seq_name
            seq = (t.split('\n')[2]).strip('\t\n')
            for N in list_mutations:
                # adds mutations to the sequence according to mutation type
                if syn == '-1':
                    new_seq = add_any_mut(seq, int(N))
                else:
                    new_seq = add_mut(seq, int(N), syn)

                # sequence element in alignment section with mutated sequence of seq_name 
                new_t = '\t\t<sequence>\n\t\t\t<taxon idref="'+seq_name+'"/>\n\t\t\t'+new_seq+'\n\t\t</sequence>'
                print(new_t)
                # convert string to etree object
                new_seqel = etree.fromstring(new_t)
                #print(new_seqel)
                #copy of xml-tree
                root1 = copy.deepcopy(root)

                # adding new seq element to alignment section
                alignment1 = root1.findall("alignment")[0]
                alignment1[k] = new_seqel

                # new xml in string format
                new_tree = etree.tostring(root1, encoding='unicode')

                # creating the directory to save output xml file
                if sys.platform == 'win32' or sys.platform == 'cygwin':
                    if path_out != '':
                        path_out = path_out + '\\'
                    if not os.path.exists(path_out + "mutations"):
                        try:
                            os.system('mkdir ' + path_out + "mutations")
                        except:
                            print('Can\'t create output folder')
                            continue
                    if not os.path.exists(path_out + "mutations\\" + N):
                        try:
                            print(os.path.exists(path_out + "mutations\\" + N))
                            os.system('mkdir ' + path_out + "mutations\\" + N)
                        except:
                            print('Can\'t create output folder')
                            continue
                    if not os.path.exists(path_out + "mutations\\" + N + "\\" + seq_name):
                        try:
                            os.system('mkdir ' + path_out + "mutations\\" + N + "\\" + seq_name)
                        except:
                            print('Can\'t create output folder')
                            continue
                    path_out_l = path_out + 'mutations\\' + N + '\\' + seq_name + '\\'
                else:
                    if path_out != '':
                        path_out = path_out + '/'
                    if not os.path.exists(path_out + "mutations"):
                        try:
                            os.system('mkdir ' + path_out + "mutations")
                        except:
                            print('Can\'t create output folder')
                            continue
                    if not os.path.exists(path_out + "mutations/" + N):
                        try:
                            print(os.path.exists(path_out + "mutations/" + N))
                            os.system('mkdir ' + path_out + "mutations/" + N)
                        except:
                            print('Can\'t create output folder')
                            continue
                    if not os.path.exists(path_out + "mutations/" + N + "/" + seq_name):
                        try:
                            os.system('mkdir ' + path_out + "mutations/" + N + "/" + seq_name)
                        except:
                            print('Can\'t create output folder')
                            continue
                    path_out_l = path_out + 'mutations/' + N + '/' + seq_name + '/'

                # creating new xml-file
                temp =  os.path.splitext(os.path.split(input_xml)[1])[0] # template for output file name
                new_file = open(path_out_l + temp + '-' +N + '.xml', 'w')
                new_file.write(new_tree)
                new_file.close()
        k+=1

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input xml-file generated in BEAUTi", required=True)
    parser.add_argument("-sn", "--seq_name", type=str,
                        help="The id of sequence which sequence will be mutated", required=True)
    parser.add_argument("-m", "--mutations", type=str,
                        help="List with numbers of mutations which should be added to \
                        the nucleotide sequence of seq_name, should be separated by spaces. \
                        Example: '1,5,10'", required=True)
    parser.add_argument("-t", "--type_mut", type=str,
                        help="The type of mutations to introduce: 1 if synonymous, \
                        0 if non-synonymous, -1 if any type", required=True)
    parser.add_argument("-pout", "--path_out", type=str,
                        help="Output directory. If not defined the output files will be saved \
                        in 'mutations' folder in the directory of input file")
    args = parser.parse_args()

    if not args.path_out:
        args.path_out = os.path.split(args.input_file)[0]

    args.mutations = args.mutations.split(" ")

    write_xml(args.input_file, args.seq_name, args.mutations, args.type_mut, args.path_out)
