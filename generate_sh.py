import os

file = open('beast_mut.sh', 'w')

#names of folders
#years = ['20-','10-','5-','3-','1+','3+','5+','10+','20+']
times = ['5','10','20']

#list_seqs = os.listdir('years/'+years[0])

list_seqs = os.listdir('mutations/syn/'+times[0])

print(list_seqs)


cur_path = os.getcwd()
print(cur_path)

#list with strings of sh script
file_str = []

k = 0
for seq in list_seqs:
	#for year in years:
	for time in times:
		k+=1
		#string = 'cd '+ cur_path + '/years/'+year+'/'+seq+'\n'
		string = 'cd '+ cur_path + '/mutations/syn/'+time+'/'+seq+'\n'
		string1 = 'cd '+ cur_path + '/mutations/nonsyn/'+time+'/'+seq+'\n'
		print(string)
		file_str.append(string)
		file_str.append('beast\n')
		file_str.append(string1)
		file_str.append('beast\n')


file.writelines(file_str)
file.close()
