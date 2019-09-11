import argparse
from lxml import etree
import copy
import os
import sys

'''
Input:
    file_name -str - path to xml-file with alignment and BEAST specifications
    seq_name - str - id of sequence (according to xml-file) which collection year will be changed
    years - list with str objects - list with numbers to add or subtract from seq_name collection date
    path_out - output directory
Output:
    Creates new xml files with changed collection year of chosen sequence. Saves files to the directory:
    path_out + 'years' + number_of_added_years
'''
#changes isolation date in sequence seq_name
def change_year(file_name, seq_name, years, path_out):
    #xml tree of original file generated in beauty
    tree = etree.parse(file_name)
    #root element of xml tree
    root = tree.getroot()

    #find sequences' names
    taxa = root.findall("taxa")[0]

    for year in years:
        print(year)
        if sys.platform == 'win32' or sys.platform == 'cygwin':
            #output path
            if path_out != '':
                path_out = path_out + "\\"
            
            print(os.path.exists(path_out + "years"))

            if not os.path.exists(path_out + "years"):
                try:
                    os.system('mkdir ' + path_out + "years")
                except:
                    print('Can\'t create output folder')
                    continue
            if not os.path.exists(path_out + "years\\" + year):
                try:
                    print(os.path.exists(path_out + "years\\" + year))
                    os.system('mkdir ' + path_out + "years\\" + year)
                except:
                    print('Can\'t create output folder')
                    continue
            if not os.path.exists(path_out + "years\\" + year + "\\" + seq_name):
                try:
                    os.system('mkdir ' + path_out + "years\\" + year + "\\" + seq_name)
                except:
                    print('Can\'t create output folder')
                    continue
            path_out_l = path_out + 'years\\' + year + '\\' + seq_name
        else:
            #output path
            if path_out != '':
                path_out = path_out + "/"
            if not os.path.exists(path_out + "years"):
                try:
                    os.system('mkdir ' + path_out + "years")
                except:
                    print('Can\'t create output folder')
                    continue
            if not os.path.exists(path_out + "years/" + year):
                try:
                    os.system('mkdir ' + path_out + "years/" + year)
                except:
                    print('Can\'t create output folder')
                    continue
            if not os.path.exists(path_out + "years/" + year + "/" + seq_name):
                try:
                    os.system('mkdir ' + path_out + "years/" + year + "/" + seq_name)
                except:
                    print('Can\'t create output folder')
                    continue
            path_out_l = path_out + 'years/' + year + '/' + seq_name
        
        print(path_out_l)
        #deep copy of xml tree root element
        root1 = copy.deepcopy(root)
        print(root1)

        #for child in root:
        #    print(child)
        '''
        mcmc = root1.findall("mcmc")[0]
        print(mcmc.attrib)

        logTree = mcmc.findall("logTree")[0]
        logTree.attrib['fileName'] = (logTree.attrib['fileName'])[:-6] + '.trees'
        print(logTree.attrib)

        log = mcmc.findall("log")[1]
        log.attrib['fileName'] = (log.attrib['fileName']).strip('.log')  + '.log'
        print(log.attrib)
        '''
        taxa = root1.findall("taxa")[0]
        print(taxa)
        # finds seq_names and changes its collection date
        for taxon in taxa:
            if taxon.attrib['id'] == seq_name:
                print(taxon.attrib)
                date = taxon.findall('date')[0]
                print(date.attrib)
                date.attrib['value'] = str(float(date.attrib['value']) + int(year))
                print(date.attrib)

        s =etree.tounicode(root1)
        #print(s)
        # creating new xml-file

        temp = os.path.splitext(os.path.split(file_name)[1])[0] # template for output file name
        
        new_file = open(path_out_l + '/' + temp + '_' + year + '.xml', 'w')
        new_file.write(s)
        new_file.close()


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-input", "--input_file", type=str,
                        help="Input xml-file generated in BEAUTi", required=True)
    parser.add_argument("-sn", "--seq_name", type=str,
                        help="The id of sequence which collection years will be changed", required=True)
    parser.add_argument("-y", "--years", type=str,
                        help="The numbers which should be added to \
                        the collection date of seq_name, should be separated by spaces", required=True)
    parser.add_argument("-pout", "--path_out", type=str,
                        help="Output directory. If not defined the output files will be saved \
                        in 'years' folder in the directory of input file")
    args = parser.parse_args()

    '''
    if sys.platform == 'win32' or sys.platform == 'cygwin':
        args.input_file = '/'.join(args.input_file.split('\\'))
        if  args.path_out:
            args.path_out = '/'.join(args.path_out.split('\\'))
    '''
    if not args.path_out:
        args.path_out = os.path.split(args.input_file)[0]

    years = args.years.split(" ")
    change_year(args.input_file, args.seq_name, years, args.path_out)