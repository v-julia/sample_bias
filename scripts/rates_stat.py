import argparse
import numpy as np
import matplotlib.pyplot as plt
import os
import re
import statistics
import sys

'''
Plots frequency distribution of rates in output directory
Input:
    rates - list or np.array of float
    output_dir - str - path to output directory
    name_fig - str - name of file with figure
    mode - type of sequence change (collection date or mutation)

'''
def plot_hist_rates(rates, output_dir, seq_name, name_fig, mode='', num_mut=''):

    # mean rate
    mu = round(statistics.mean(rates),4)
    # standart deviation of rate
    sigma = round(statistics.stdev(rates),4)

    # bins for histogram
    bins = []
    # step in histogram
    step = (max(rates) - min(rates))/50
    i = min(rates) - step * 3
    while i < max(rates) + step * 3:
        bins.append(i)
        i += step
    # plots histogram
    plt.hist(rates, bins)
    plt.xlim(min(rates) - step * 3, max(rates) + step * 3)
    # plots mean, mean+-(1,2,3)sd
    xposition = [mu-3*sigma, mu-2*sigma, mu-sigma,  mu+sigma, mu+2*sigma, mu+3*sigma]
    for xc in xposition:
        plt.axvline(x=xc, color='red', linestyle='--')
    plt.axvline(x=mu, color='red', linestyle='-')
    
    plt.title(' '.join([seq_name, mode, str(num_mut)]))
    plt.text(0.5,0.5, '$\mu='+str(mu)+', \sigma$='+str(sigma), \
            fontsize=12, horizontalalignment='center', verticalalignment='center')
    plt.xlabel("Substitution rate")
    plt.ylabel("Frequency")
    plt.savefig(Path(output_dir,name_fig), transparent = True)
    #plt.show()
    plt.close()



'''
Calculates mean, sd of substitution rates and their logs. Plots frequency distribution 
of rates/ their logs, shows mean rate (solid red line) and mean+(1,2,3)sd (dashed red lines)

Input:
    tree_path - str - path to file with MCC tree generated by TreeAnnotator
    seq_name -str - id of sequence
    mode - str - type of sequence change
    num_ch - str - the number of mutations introduced into sequences or years added
    to sequence colletion date
Output:
    rate_num - float - subsitution rate of branch leading to seq_name
    mu - float - mean subsitution rate in tree
    sigma - float - standard deviation of subsitution rates in tree
'''
def get_stat(tree_path, seq_name, mode='', num_ch = ''):
    print(seq_name)
    # output directory
    output_dir = os.path.split(tree_path)[0]

    #open file with tree
    print(tree_path)
    file = open(tree_path, 'r')
    # flag , k=1 if Translation block has started
    k=0
    for line in file:
        if line == '\tTranslate\n':
            k = 1
            continue
        if k == 1:
            # line contains sequences' ids and their serial number in parentheses tree
            line = ((line.strip('\n')).strip('\t')).strip(',')
            if line.split(' ')[1] == seq_name:
                #number of seq in parentheses tree
                num = line.split(' ')[0]
                k = 0
        # parentheses tree
        if line[0:4] == 'tree':
            tree_par = line
            # all rates in MCC tree
            m = re.findall("rate=[0-9]+\.[0-9E\-]+", line)
    file.close()
    #print(m)

    # get rate of sequence seq_name
    tree_par = tree_par.replace('%', '')
    r = re.search( r'[[(,]%s\[&[a-zA-z,_0-9\-\.=\{\}]*\]' %str(num), tree_par)
    
    # substring with parameters for seq number num
    numsubstr = r.group()
    ratesubstr = re.search(r'rate=[0-9]+\.[0-9E\-]+',numsubstr)
    rate_num = round(float(ratesubstr.group()[5:]),4)
    rate_num_log = np.log10(rate_num)

    # lists with rates and log rates
    rates = []
    rates_log = []
    

    # output text file with rates retrieved from tree
    file_out = open(Path(output_dir,'rates.txt'),'w')

    # converts all rates found in tree into int
    for rate in m:
        rate = rate[5:]
        #print(rate)
        rates.append(round(float(rate),4))
        rate_log = np.log10(float(rate))
        rates_log.append(rate_log)
        file_out.write(rate+'\t'+str(rate_log)+'\n')

    # mean subsitution rate
    mu = round(statistics.mean(rates),4)
    # sd of subsitution rate
    sigma = round(statistics.stdev(rates),4)
    
    # mean log subsitution rate
    mu_log = round(statistics.mean(rates_log),4)
    # sd of subsitution rate
    sigma_log = round(statistics.stdev(rates_log),4)

    file_out.write("mean="+str(mu)+",sd="+str(sigma)+"\nmean_log="+str(mu_log)+",sigma_log="+str(sigma_log))
    file_out.close()


            
    plot_hist_rates(rates, output_dir, seq_name, "hist_rates.svg", mode, num_ch)
    plot_hist_rates(rates_log, output_dir, seq_name, "hist_log_rates.svg", mode, num_ch)

    return rate_num, mu, sigma, rate_num_log, mu_log, sigma_log
'''
Input:
    path_fold_changed - str- path to folder with results of changing xml files
    path_orig_tree -str - path to original tree
    mode - str - type of changes introduces to sequence
Ouput:
    dict_seq - dict - dictionary with rates for sequences after changing isolation year
    dict_seq[seq_name][year][rate] = rate
    dict_seq[seq_name][year][mean] = mean rate
    dict_seq[seq_name][year][sigma] = sd rate
'''
def get_rates(path_orig_tree, path_fold_changed, mode):

    # list of values of added mutations or years

    list_param = [int(x) for x in os.listdir(path_fold_changed)]
    print(list_param)

    # dictionaries with rates/logrates of sequences after changing year/seq 
    # dict_seq[seq_name][year][rate] = rate - rate of changed sequence
    # dict_seq[seq_name][year][mean] = mean rate - mean rate in tree with changed seq
    # dict_seq[seq_name][year][sigma] = sd rate - sd of rate in tree with changed seq
    dict_seq = {}
    dict_seq_log = {}

    # list of sequences with changes collection years or mutated sequence
    list_seq = os.listdir(path_fold_changed + str(list_param[0]))

    print(list_seq)

    # template for name of filwe with MCC tree
    tree_name_temp = os.path.split(path_orig_tree)[1]

    #print(sample_path)
    for i in range(len(list_param)):
        for seq_name in list_seq:
            tree_path = Path(path_fold_changed+ str(list_param[i]),seq_name,tree_name_temp)
            if os.path.exists(tree_path):
                if seq_name not in dict_seq.keys():
                    dict_seq[seq_name] = {}
                    #dict_seq[seq_name][list_param[i]] = {}
                    dict_seq_log[seq_name] = {}
                    #dict_seq_log[seq_name][list_param[i]] = {}
                
                if list_param[i] not in dict_seq[seq_name].keys():
                    dict_seq[seq_name][list_param[i]] = {}
                    dict_seq_log[seq_name][list_param[i]] = {}

                rate, mean, sigma, rate_log, mean_log, sigma_log = get_stat(tree_path, seq_name, mode, list_param[i])
                # writes values to to dictionaries
                dict_seq[seq_name][list_param[i]]['rate'] = round(rate,4)
                dict_seq[seq_name][list_param[i]]['mean'] = mean
                dict_seq[seq_name][list_param[i]]['sigma'] = sigma
                dict_seq_log[seq_name][list_param[i]]['rate'] = round(rate_log,4)
                dict_seq_log[seq_name][list_param[i]]['mean'] = mean_log
                dict_seq_log[seq_name][list_param[i]]['sigma'] = sigma_log
            else:
                continue

    #sample_tree_path = 'sample/sample_al.tree'

    # add sequences' rates from the original file
    for seq_name in list_seq:
        tree_path = Path(path_fold_changed+ str(list_param[0]),seq_name,tree_name_temp)
        print("this is weird")
        print(tree_path)
        if os.path.exists(tree_path):
            rate, mean, sigma, rate_log, mean_log, sigma_log = get_stat(path_orig_tree, seq_name)
            dict_seq[seq_name][0] = {}
            dict_seq[seq_name][0]['rate'] = round(rate,4)
            dict_seq[seq_name][0]['mean'] = mean
            dict_seq[seq_name][0]['sigma'] = sigma
            dict_seq_log[seq_name][0] = {}
            dict_seq_log[seq_name][0]['rate'] = round(rate_log,4)
            dict_seq_log[seq_name][0]['mean'] = mean_log
            dict_seq_log[seq_name][0]['sigma'] = sigma_log
            
    return dict_seq, dict_seq_log

'''
Input:
    dict_seq - dict - dictionaries with rates/logrates of sequences after changing year/seq 
    dict_seq[seq_name][year][rate] = rate - rate of changed sequence with id seq_name
    dict_seq[seq_name][year][mean] = mean rate - mean rate in tree with changed seq
    dict_seq[seq_name][year][sigma] = sd rate - sd of rate in tree with changed seq

    path_fold_changed - str - path to folder with results of changing xml files
    mode - str - type of changes introduces to sequence
Output:
    Writes dictionary to file with name path_fold_changed + 'table_'+mod+'.txt'
'''
def print_table(dict_seq, path_fold_changed, mod):

    # list with line for future filw with table
    table_str = []

    # first line
    str1 = ' '
    # list of changed sequences' ids
    list_seq_avail = list(dict_seq.keys())
    #print(list_seq_avail)
    # list of values of added years or mutations
    list_param = list(dict_seq[list_seq_avail[0]].keys())
    list_param.sort()
    #print(list_param)
    for par in list_param:
        str1 = str1 +'\t'+str(par)
    str1+='\n'
    table_str.append(str1)

    # strings with rates' values
    for seq_name in list_seq_avail:
        stri = seq_name
        #print(seq_name)
        for param in list_param:
            if param in dict_seq[seq_name].keys():
                stri = stri+'\t'+str(dict_seq[seq_name][param]['rate'])
            else:
                stri = stri+'\tNone'
        stri+='\n'
        table_str.append(stri)

    file_out = open(Path(path_fold_changed,'table_'+mod+'.txt'), 'w')
    file_out.writelines(table_str)
    file_out.close()

'''
Input:
    dict_seq - dict - dictionaries with rates/logrates of sequences after changing year/seq 
    dict_seq[seq_name][year][rate] = rate - rate of changed sequence with id seq_name
    dict_seq[seq_name][year][mean] = mean rate - mean rate in tree with changed seq
    dict_seq[seq_name][year][sigma] = sd rate - sd of rate in tree with changed seq

    path_fold_changed - str - path to folder with results of changing xml files
    x_axis_title - str - type of changes introduces to sequence
    format - str - format of figure
Output:
    Draws linear plot: x-axis - values that were added to collection date of sequence or 
    number of mutations introduced to sequence with id seq_name,
    y-axis - rates of sequence seq_name in mcc trees where the collection date of 
    nucleotide sequences of seq_name were changed, red solid line - mean rate in trees,
    grey dashed lines - mean +- (1,2,3)sd
'''

def plot_rates(dict_seq, path_fold_changed, x_axis_title, format):
    # ids of changed sequences
    seq_names = list(dict_seq.keys())
    # list of values of added years or mutations
    params = [int(x) for x in list(dict_seq[seq_names[0]].keys())]
    print(params)
    params.sort()
    params = np.array(params)

    # draws plot for each sequence
    for seq in seq_names:
        rates = np.empty(1)
        means = np.empty(1)
        sds = np.empty(1)
        #print(rates)
        for param in params:
            rates = np.append(rates, dict_seq[seq][param]['rate'])
            means = np.append(means, dict_seq[seq][param]['mean'])
            sds = np.append(sds, dict_seq[seq][param]['sigma'])
        print(rates)
        print(sds)
        print(means)
        print(params)
        plt.plot(params, rates[1:], linestyle = '-', c='black')
        plt.plot(params, means[1:], linestyle= '-', c = 'red')
        plt.plot(params, means[1:] + sds[1:] , linestyle = '--', c='grey')
        plt.plot(params, means[1:]- sds[1:], linestyle = '--', c='grey')
        plt.plot(params, means[1:] + 2*sds[1:] , linestyle = '--', c='grey')
        plt.plot(params, means[1:]- 2*sds[1:], linestyle = '--', c='grey')
        plt.plot(params, means[1:] + 3*sds[1:] , linestyle = '--', c='grey')
        plt.plot(params, means[1:]- 3*sds[1:], linestyle = '--', c='grey')
        plt.title(seq)
        plt.xlabel(x_axis_title)
        plt.ylabel("Substitution rate")
        plt.savefig(Path(path_fold_changed,seq + format), transparent = True)
        #plt.show()
        plt.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-orig", "--orig_tree_path", type=str,
                        help="Path to original MCC tree generated using TreeAnnotator", required = True)
    parser.add_argument("-changed", "--folder_changed_path", type=str,
                        help="Path to the folder with trees were sequences' collection years or \
                        nucleotide sequences are changed", required = True)
    parser.add_argument("-m", "--mode", type=str,
                        help="Type of sequences' change:'mutations' if mutations were introduced \
                        into nucleotide sequence of seq_name, 'years' if some values were added to \
                        the collection date of seq_name", required = True)
    args = parser.parse_args()

    # creates dictionary with individual rate, mean rate of tree, sd of rates in tree for each changed sequence
    dict_rates, dict_rates_log  = get_rates(args.orig_tree_path, args.folder_changed_path, args.mode)
    # draws the dependence of seq_name rate on the value of parameter change
    plot_rates(dict_rates, args.folder_changed_path, '# of added {}'.format(args.mode), '.svg')
    plot_rates(dict_rates_log, args.folder_changed_path, '# of added {}'.format(args.mode), '_log.svg')
    # prins table
    print_table(dict_rates, args.folder_changed_path, 'years')