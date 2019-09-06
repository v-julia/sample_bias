library("coda")
library("reshape2")
library("ggplot2")

# input:
#   filename - path to file with mcmc log
#   burnin_pt - burnin pertentage
# ouput:
#   preprocessed dataframe with mcmc sampling results with sufficient columns

load_data = function(filename, burnin_pt) {
  # read log file, skips the first 4 lines with comments
  data = read.csv(filename, skip=4, header=T, sep="\t")
  
  # number of lines in mcmc output
  n = nrow(data)

  # changes col name
  names(data)[names(data) == 'treeModel.rootHeight'] = 'treeHeight'
  
  
  # Select relevant and discard burnin
  data = data[(burnin_pt*n/100):n,c("meanRate", "treeHeight", "ucld.mean", "ucld.stdev", 
                                      "age.root.", "treeLength")]
  
  # Post-processing  
  colnames(data) = c("Clock rate", "Tree height", "ucld.mean", "ucld.stdev", "age(root)", "Tree length")
  #data = melt(data, measure.vars = 1:ncol(data))
  return (data)
}

# input:
#   data - dataframe created from log file using load_data function
# output:
#   df - dataframe with 95%HPD_lower,95%HPD_higher,median for each parameter (column) in data

basic_stat = function(data){
  
  # create mcmc object from dataframe
  mcmc_obj = mcmc(mcmc_table)
  
  # hpd intervals
  hpd = HPDinterval(mcmc_obj)
  
  # medians
  m = apply(mcmc_obj, 2, median)
  
  df = as.data.frame(t(cbind(hpd, m)))
  rownames(df) = c('95%HPD_lower','95%HPD_higher','median')

  return(df)
}

# input: 
#   input_dir - the directory with log files for the model
#   parameter - parameter estimated in mcmc (col name of log file)
# output:
#   df - dataframe with median, lower and higher HPD intervals of parameter

sample_table = function(input_dir){
  
  # parameters to plot
  parameters = c('Clock rate', 'Tree height', 'Tree length')
  
  # dataframe with hpd and median of rates for simulated samples
  rate_df = data.frame(row.names = c('95HPD%_lower','95HPD%_higher','median'))
  
  # dataframe with hpd and median of tree heights
  height_df = data.frame(row.names = c('95HPD%_lower','95HPD%_higher','median'))
  
  # dataframe with hpd and median of tree lengths
  length_df = data.frame(row.names = c('95HPD%_lower','95HPD%_higher','median'))
  
  # log files of simulated samples under specific model
  file_names = list.files(input_dir, pattern = "\\.log$")

  for (file_name in file_names){
    # table with lower and higher95% hpds and median values for each parameter 
    table_stat = basic_stat(load_data(file.path(input_dir, file_name), 15))
    print(file.path(input_dir, file_name))
    prn
    
    rate_df = cbind(rate_df, table_stat$`Clock rate`)
    height_df = cbind(height_df, table_stat$`Tree height`)
    length_df = cbind(length_df, table_stat$`Tree length`)
  }
  colnames(rate_df) = seq(1:length(file_names))
  colnames(height_df) = seq(1:length(height_df))
  colnames(length_df) = seq(1:length(length_df))
  
  return(list(rate_df,height_df,length_df))
}

l = sample_table(input_dir)


setwd("D:/DATA/samplebias/biases/EV-A71_alignments/")

# folders with beast results for different models
model_folders = c("random_single", "random_groups")

list_rates = list()
list_tee_heights = list()


input_dir = "D:/DATA/samplebias/biases/EV-A71_alignments/random_single"

mcmc_table = load_data('EV71_VP1_all_cut_dates_no_chin_singlerandom_3.log',15)
mcmc_obj = mcmc(mcmc_table)
q = HPDinterval(mcmc_obj)
m = as.data.frame(apply(mcmc_obj, 2, median))
colnames(m) = 'median'
q = cbind(q,m)

file_names = list.files(input_dir, pattern = "\\.log$")


plot(mcmc_obj)

mcmc_table$meanRate
mcmc_table$treeModel.rootHeight
mcmc_table$age.root.