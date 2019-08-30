#this script gets rate of isolates with changed isoltaion date or sequence
#then creates table with rates
import os
import re
import numpy as np
import statistics
import matplotlib.pyplot as plt

#gets rate of seq_name in tree from tree_path file
def get_rate(tree_path, seq_name):
	file_tree = open(tree_path, 'r')
	k = 0
	for line in file_tree:
		if line == '\tTranslate\n':
			k =1
			continue
		if k == 1:
			
			line = ((line.strip('\n')).strip('\t')).strip(',')
			#print(line)
			#number of seq in parenthesis tree
			if line.split(' ')[1] == seq_name:
				num = line.split(' ')[0]
				k = 0
		if line[0:4] == 'tree':
			tree_par = line
			break
	#print(num)
	tree_par = tree_par.replace('%', '')
	r = re.search( r'[[(,]%s\[&[a-zA-z,_0-9\-\.=\{\}]*\]' %str(num), tree_par)
	#substring with parameters for seq number num
	numsubstr = r.group()
	print(numsubstr)
	ratesubstr = re.search(r'rate=[0-9\.]+',numsubstr)
	print(ratesubstr.group())
	rate = round(float(ratesubstr.group()[5:]),5)
	print(rate)
	file_tree.close()
	return rate

#return mean rate and stdev from tree in tree_path
def get_stat(tree_path):
	#open file with tree
	print(tree_path)
	file = open(tree_path, 'r')
	for line in file:
		if line[0:4] == 'tree':
			m = re.findall("rate=[0-9]+\.[0-9E\-]+", line)
	file.close()
	print(m)
	
	#lists with rates and log rates
	rates = []
	rates_log = []
	
	output_dir = '/'.join(tree_path.split('/')[:-1])+'/'
	file_out = open(output_dir+'rates.txt','w')
	
	for rate in m:
		rate = rate[5:]
		#print(rate)
		rates.append(float(rate))
		rate_log = np.log10(float(rate))
		rates_log.append(rate_log)
		file_out.write(rate+'\t'+str(rate_log)+'\n')
	file_out.close()
	
	mu = round(statistics.mean(rates),4)
	sigma = round(statistics.stdev(rates),4)
	
	mu_log = round(statistics.mean(rates_log),4)
	sigma_log = round(statistics.stdev(rates_log),4)
	
	print('mean ',mu)
	print('sigma ',sigma)
	

	list_file = tree_path.split('/')
	if len(list_file) == 3:
		return mu, sigma, mu_log, sigma_log
	else:

		seq_name = tree_path.split('/')[-2]
		mode = tree_path.split('/')[-4]
		print('mode ', mode)
		num_mut = tree_path.split('/')[-3]


		bins = []
		i = 0
		while i<0.09:
			bins.append(i)
			i+=0.001

		plt.hist(rates, bins)
		plt.xlim(0, 0.01)
		plt.title(seq_name +' '+ mode+' '+num_mut)
		plt.text(0.006,25, '$\mu='+str(mu)+', \sigma$='+str(sigma), fontsize=12)
		plt.xlabel("Substitution rate")
		#plt.yscale('log')
		#plt.xscale('log')
		plt.ylabel("Frequency")
		plt.savefig(output_dir+ 'rates_hist.png', transparent = True)
		#plt.show()
		plt.close()


		bins_log = []
		i = -3.9
		while i<-2:
			bins_log.append(i)
			i+=0.01

		plt.hist(rates_log, bins_log)
		#plt.xlim(0, 0.01)
		plt.title(seq_name +' '+ mode+' '+num_mut)
		plt.text(-3.5,2.6, '$\mu='+str(mu_log)+', \sigma$='+str(sigma_log), fontsize=12)
		plt.xlabel("Substitution rate, log")
		#plt.yscale('log')
		plt.ylabel("Frequency")
		plt.savefig(output_dir +'rates_hist_log.png', transparent = True)
		#plt.show()
		plt.close()
		
		return mu, sigma, mu_log, sigma_log


#dictionary with rates for sequences after changing isolation year
#dict_seq[seq_name][syn][num_mut] = rate
dict_seq = {}

#list_seq = os.listdir('years/'+ list_years_str[0])
list_seq = os.listdir('mutations/syn/5/')
print(list_seq)

list_mut = os.listdir('mutations/syn/')
for i in range(len(list_mut)):
	mut = int(list_mut[i])
	list_mut[i] = mut
print(list_mut)
list_seq_avail = []
for i in range(len(list_mut)):
	for seq_name in list_seq:
		tree_path = 'mutations/syn/'+ str(list_mut[i])+'/'+seq_name + '/sample_al_VP1.tree'
		if os.path.exists(tree_path):
			print('tree_path ',tree_path)
			list_seq_avail.append(seq_name)
			if seq_name not in dict_seq.keys():
				dict_seq[seq_name] = {}
				dict_seq[seq_name]['syn'] = {}
				dict_seq[seq_name]['syn'][0] = {}
			if 'syn' not in dict_seq[seq_name].keys():
				dict_seq[seq_name]['syn'] = {}
				dict_seq[seq_name]['syn'][0] = {}
			if list_mut[i] not in dict_seq[seq_name]['syn'].keys():
					dict_seq[seq_name]['syn'][list_mut[i]] = {}
			
			
			dict_seq[seq_name]['syn'][list_mut[i]]['rate'] = get_rate(tree_path, seq_name)
			dict_seq[seq_name]['syn'][list_mut[i]]['rate_log'] = round(np.log10(dict_seq[seq_name]['syn'][list_mut[i]]['rate']),4)
			mu, sigma, mu_log, sigma_log = get_stat(tree_path)
			dict_seq[seq_name]['syn'][list_mut[i]]['mean'] = mu
			dict_seq[seq_name]['syn'][list_mut[i]]['sigma'] = sigma
			dict_seq[seq_name]['syn'][list_mut[i]]['mean_log'] = mu_log
			dict_seq[seq_name]['syn'][list_mut[i]]['sigma_log'] = sigma_log
		else:
			continue
		tree_path = 'mutations/nonsyn/'+ str(list_mut[i])+'/'+seq_name + '/sample_al_VP1.tree'
		if os.path.exists(tree_path):
			print('tree_path ',tree_path)
			list_seq_avail.append(seq_name)
			if seq_name not in dict_seq.keys():
				dict_seq[seq_name] = {}
				dict_seq[seq_name]['nonsyn'] = {}
				dict_seq[seq_name]['nonsyn'][0] = {}
				print(dict_seq[seq_name]['nonsyn'].keys())
			if 'nonsyn' not in dict_seq[seq_name].keys():
				dict_seq[seq_name]['nonsyn'] = {}
				dict_seq[seq_name]['nonsyn'][0] = {}
			if list_mut[i] not in dict_seq[seq_name]['nonsyn'].keys():
				dict_seq[seq_name]['nonsyn'][list_mut[i]] = {}
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['rate'] = get_rate(tree_path, seq_name)
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['rate_log'] = round(np.log10(dict_seq[seq_name]['nonsyn'][list_mut[i]]['rate']),4)
			mu, sigma, mu_log, sigma_log = get_stat(tree_path)
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['mean'] = mu
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['sigma'] = sigma
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['mean_log'] = mu_log
			dict_seq[seq_name]['nonsyn'][list_mut[i]]['sigma_log'] = sigma_log
		else:
			continue

#sample_tree_path = 'sample/sample_al.tree'
sample_tree_path = 'sample/VP1_BrCr_less/sample_al_VP1.tree'
print(dict_seq)

print('zero mutations')

list_seq_avail = set(list_seq_avail)

for seq_name in list_seq_avail:
	print(seq_name)
	if 'syn' in dict_seq[seq_name].keys():
		print('syn')
		dict_seq[seq_name]['syn'][0]['rate'] = get_rate(sample_tree_path, seq_name)
		mu, sigma, mu_log, sigma_log = get_stat(sample_tree_path)
		dict_seq[seq_name]['syn'][0]['mean'] = mu
		dict_seq[seq_name]['syn'][0]['sigma'] = sigma
		dict_seq[seq_name]['syn'][0]['rate_log'] = round(np.log10(dict_seq[seq_name]['syn'][0]['rate']),4)
		dict_seq[seq_name]['syn'][0]['mean_log'] = mu_log
		dict_seq[seq_name]['syn'][0]['sigma_log'] = sigma_log

	if 'nonsyn' in dict_seq[seq_name].keys():
		print('nonsyn')
		dict_seq[seq_name]['nonsyn'][0]['rate'] = get_rate(sample_tree_path, seq_name)
		mu, sigma, mu_log, sigma_log = get_stat(sample_tree_path)
		dict_seq[seq_name]['nonsyn'][0]['mean'] = round(mu,4)
		dict_seq[seq_name]['nonsyn'][0]['sigma'] = sigma
		dict_seq[seq_name]['nonsyn'][0]['rate_log'] = round(np.log10(dict_seq[seq_name]['nonsyn'][0]['rate']),4)
		dict_seq[seq_name]['nonsyn'][0]['mean_log'] = mu_log
		dict_seq[seq_name]['nonsyn'][0]['sigma_log'] = sigma_log



print(dict_seq)

table_str = []
table_str.append('synonymous\n')
str1 = ' '
list_mut.append(0)
list_mut.sort()
print(list_mut)
for num in list_mut:
	str1 = str1 +'\t'+str(num) + '\tmean\tstdev\t-3sigma\t+3sigma\trate_log\tmean_log\tstdev_log\t-3sigma_log\t+3sigma_log'
str1+='\n'
table_str.append(str1)



for seq_name in list_seq_avail:
	stri = seq_name
	if 'syn' in dict_seq[seq_name].keys():
		print(seq_name)
		for num in list_mut:
			if num in dict_seq[seq_name]['syn'].keys():
				mean =dict_seq[seq_name]['syn'][num]['mean']
				stdev = dict_seq[seq_name]['syn'][num]['sigma']
				mean_log =dict_seq[seq_name]['syn'][num]['mean_log']
				stdev_log = dict_seq[seq_name]['syn'][num]['sigma_log']
				stri = stri+'\t'+str(dict_seq[seq_name]['syn'][num]['rate']) + '\t'+str(mean)+'\t'+str(stdev)+'\t'+str(round(mean- 3*stdev,4))+'\t'+str(round(mean+ 3*stdev,4))+'\t'+str(dict_seq[seq_name]['syn'][num]['rate_log'])+ '\t'+str(mean_log)+'\t'+str(stdev_log)+'\t'+str(round(mean_log- 3*stdev_log,4))+'\t'+str(round(mean_log+ 3*stdev_log,4))
			else:
				stri = stri+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'
		stri+='\n'
		table_str.append(stri)

table_str.append('nonsynonymous\n')

table_str.append(str1)
for seq_name in list_seq_avail:
	stri = seq_name
	if 'nonsyn' in dict_seq[seq_name].keys():
		print(seq_name)
		for num in list_mut:
			if num in dict_seq[seq_name]['nonsyn'].keys():
				mean =dict_seq[seq_name]['nonsyn'][num]['mean']
				stdev = dict_seq[seq_name]['nonsyn'][num]['sigma']
				mean_log =dict_seq[seq_name]['nonsyn'][num]['mean_log']
				stdev_log = dict_seq[seq_name]['nonsyn'][num]['sigma_log']
				stri = stri+'\t'+str(dict_seq[seq_name]['nonsyn'][num]['rate']) + '\t'+str(mean)+'\t'+str(stdev)+'\t'+str(round(mean- 3*stdev,4))+'\t'+str(round(mean+ 3*stdev,4))+'\t'+str(dict_seq[seq_name]['nonsyn'][num]['rate_log'])+ '\t'+str(mean_log)+'\t'+str(stdev_log)+'\t'+str(round(mean_log- 3*stdev_log,4))+'\t'+str(round(mean_log+ 3*stdev_log,4))
			else:
				stri = stri+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'+'\tNone'
		stri+='\n'
		table_str.append(stri)



file_out = open('table_mut.txt', 'w')
file_out.writelines(table_str)
file_out.close()
	
