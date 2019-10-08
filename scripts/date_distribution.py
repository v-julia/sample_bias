import re
import matplotlib.pyplot as plt
import numpy as np



def get_modif_years(file_gb_name):
    '''
    Gets modification date of each entry in file_gb_name file in Genbank format
    file_gb_name - str, name of GenBank files
    '''

    modif_years = [] # list with modification years

    with open(file_gb_name) as input_file:
        for line in input_file:
            if line.startswith("LOCUS"):
                modif_date = re.search(r"[0-9]{1,2}-[a-zA-Z]{3}-[0-9]{4}", line).group()
                modif_year = int(re.search(r"[0-9]{4}",modif_date).group())
                modif_years.append(modif_year)
    input_file.close()
    
    return modif_years

def get_collection_years(file_gb_name):
    '''
     Gets collection date of each entry in in file_gb_name file in Genbank format
     collection date is retrived from source - collection date qualifier
     file_gb_name - str, name of GenBank files
     '''
    collection_years = [] # list with collection dates
    
    with open(file_gb_name) as input_file:
        for line in input_file:
            m = re.search(r"\s+/collection_date=\"",line)
            
            if m:
                
                year = re.search(r"[0-9]{4}", line)
                if year:
                    collection_years.append(int(year.group()))
    input_file.close()
    
    return(collection_years)



def plot_years(input_dir, viral_groups, type):
    '''
    Input:
        input_dir - str - directory with genbank-files
        viral_groups - list - list with groups (names of GenBank files without ext)
        type - str - "collect" if collection dates should be retrieved from gb-file,
                     "modif" if modification dates -""-
    Ouput:
        Plots collection/modification dates for different viral groups in one graph
    '''

    years = {} # dictionary with dats for each viral group
    # years[group][type_of_dates] = list with dates

    for vir_group in viral_groups:
        years[vir_group] = {}

        # get modification dates from gb file with sequences
        if type == "modif":
            years[vir_group]["years_all"] = get_modif_years(input_dir+vir_group+'.gb')
            
        # get collection dates from gb file with sequences
        if type == "collect":
            years[vir_group]["years_all"] = get_collection_years(input_dir+vir_group+'.gb')
        #print(years[vir_group]["years_all"])
        # list with years without replicates
        years[vir_group]["years"] = list(set(years[vir_group]["years_all"]))
        years[vir_group]["years"].sort()
        # counts of years
        years[vir_group]["counts"] = [years[vir_group]["years_all"].count(year) for year in years[vir_group]["years"]]
        
        plt.plot(years[vir_group]["years"], (years[vir_group]["counts"]), label = vir_group)
        
    if type == "modif":
        plt.title("Modification years")
    if type == "collect":
        plt.title("Collection dates")
    plt.xlabel("Year")
    plt.ylabel("Frequency")
    plt.legend()
    plt.savefig(input_dir+type+".svg")
    plt.savefig(input_dir+type+".png")
    plt.show()


def plot_hist_years(input_dir, viral_groups, type, fig_type):
    '''
    Input:
        input_dir - str - directory with genbank-files
        viral_groups - list - list with groups (names of GenBank files without ext)
        type - str - "collect" if collection dates should be retrieved from gb-file,
                     "modif" if modification dates -""-
    Ouput:
        Plots collection/modification dates for different viral groups in one graph
    
    plots histogram of collection/modification dates for different viral groups,
    four graphs in one figure
    viral_groups - dictionary
    viral_groups[name_of_group] = [list of dates]
    type - type of dates to get from gb file,
    "collect" for collection dates, "modif" for modification years
    '''

    print(fig_type)
    if fig_type == "all":
        fig, subplots = plt.subplots(nrows=2, ncols=2, sharex=True, sharey=True)
        i = 0
        if type == "modif":
            fig.suptitle("Distribution of modification years")
        if type == "collect":
            fig.suptitle("Distribution of collection dates")
    years = {} # dictionary with dates for each viral group

    
    for vir_group in viral_groups:
        print(vir_group)
        years[vir_group] = {}
        # get modification dates from gb file with sequences
        if type == "modif":
            years[vir_group]["years_all"] = get_modif_years(input_dir+vir_group+'.gb')
        # get collection dates from gb file with sequences
        if type == "collect":
            years[vir_group]["years_all"] = get_collection_years(input_dir+vir_group+'.gb')
        
        # list with years without replicates
        years[vir_group]["years"] = list(set(years[vir_group]["years_all"]))
        years[vir_group]["years"].sort()
        if fig_type == "all":
            fig.axes[i].hist(years[vir_group]["years_all"], bins = list(range(min(years[vir_group]["years"]), max(years[vir_group]["years"]))), log=True)
            fig.axes[i].set_xlim(1924,2018)
            fig.axes[i].set_title(vir_group)
            
            if i ==0 or i == 2:
                fig.axes[i].set_ylabel('Frequency')
            if i ==2 or i ==3:
                fig.axes[i].set_xlabel('Year')
            i = i + 1
        else:
            plt.hist(years[vir_group]["years_all"], bins = list(range(min(years[vir_group]["years"]), max(years[vir_group]["years"]))), log=True)
            plt.title(vir_group)
            plt.xlabel("Year")
            plt.xlim(min(years[vir_group]["years"]),max(years[vir_group]["years"]))
            plt.ylabel("Frequency")
            plt.savefig(input_dir + vir_group + ".png")
            plt.savefig(input_dir + vir_group + ".svg")
            #plt.show()
            
            
        
    
    if fig_type == "all":
        fig.savefig(input_dir + type + "_hist.svg")
        fig.savefig(input_dir + type + "_hist.png")


if __name__ == "__main__":

    input_dir = ''

    viral_groups = ["EV-A",  "HepatoA","FMDV", "ParechoA"]

    #plot_years(viral_groups, "modif")
    #plot_years(viral_groups,"collect")

    #plot_hist_years(viral_groups, "modif")
    #plot_hist_years(viral_groups,"collect", "single")
    plot_hist_years(input_dir, viral_groups, "collect", "all")