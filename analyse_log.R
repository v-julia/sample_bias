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
  data = read.csv(filename, skip=2, header=T, sep="\t")
  
  # number of lines in mcmc output
  n = nrow(data)

  # changes col name
  names(data)[names(data) == 'treeModel.rootHeight'] = 'treeHeight'
  
  
  # Select relevant and discard burnin
  #data = data[(burnin_pt*n/100):n,c("meanRate", "treeHeight", "ucld.mean", "ucld.stdev", 
  #                                    "age.root.", "treeLength")]
  data = data[(burnin_pt*n/100):n,c("meanRate", "treeHeight", "ucld.mean", "ucld.stdev")]  

  # Post-processing  
  #colnames(data) = c("Clock rate", "Tree height", "ucld.mean", "ucld.stdev", "age(root)", "Tree length")
  colnames(data) = c("Clock rate", "Tree height", "ucld.mean", "ucld.stdev")
  
  
  #data = melt(data, measure.vars = 1:ncol(data))
  return (data)
}

# input:
#   data - dataframe created from log file using load_data function
# output:
#   df - dataframe with 95%HPD_lower,95%HPD_higher,median for each parameter (column) in data

basic_stat = function(data){
  
  # create mcmc object from dataframe
  mcmc_obj = mcmc(data)
  
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
  #length_df = data.frame(row.names = c('95HPD%_lower','95HPD%_higher','median'))
  
  # log files of simulated samples under specific model
  file_names = list.files(input_dir, pattern = "\\.log$")

  for (file_name in file_names){
    # table with lower and higher95% hpds and median values for each parameter 
    table_stat = basic_stat(load_data(file.path(input_dir, file_name), 15))
    print(file.path(input_dir, file_name))

    
    rate_df = cbind(rate_df, table_stat$`Clock rate`)
    height_df = cbind(height_df, table_stat$`Tree height`)
    #length_df = cbind(length_df, table_stat$`Tree length`)
  }
  colnames(rate_df) = seq(1:length(file_names))
  colnames(height_df) = seq(1:length(height_df))
  #colnames(length_df) = seq(1:length(length_df))
  
  return(list(rate_df,height_df)) #,length_df))
}

# prepares dataframe to be suitable for plotting
# input:
#   df - dataframe with lower and higher hpd
#       and median of parameter in log_files from simulations under model
#   x_st - start point for x axis value
#   model - model of generating samples
# output:
#   df - df with column x - x-axis coordinates, and 'model' - type of model
prepare_df = function(df, x_st, model){
  
  # positions of values on x-axis
  if (ncol(df)==1){
    df = rbind(df, x_st+0.5)
  }
  else{
    df = rbind(df, seq(x_st+0.4,x_st+0.6,length.out=ncol(df)))
  }
  # bind row with model name
  df1 = rbind(df, rep(model,ncol(df)))
  rownames(df1)[nrow(df1)-1]='x'
  rownames(df1)[nrow(df1)]='model'
  #print(df1)
  # eto ne rabotaet
  #df1[1,] = as.numeric(df1[1,])
  #df1[2,] = as.numeric(df1[2,])
  #df1[3,] = as.numeric(df1[3,])
  #print(df1[1,1])
  return(df1)
  
}


input_dir = "D:/MY_FILES/DATA/Lukashev/Enteroviruses/sample_bias/sample_biasev71/log_stat"

# dataframe with coorespondence of folder names to model names
# the name of folder coresponds to the model of sampling generation
models = data.frame(c('single','groups', '2.5threshold'), c('random single', 'random groups', '2.5% identity threshold'), stringsAsFactors = F)
colnames(models) = c('folder', 'model')

# creating the first row for tables with hpd and median
l = sample_table(file.path(input_dir, models[1,]$folder))

# dataframe with median and hpds of mean rate, x coordinate on plot, model name for 
# each simulation under model
df_rates = prepare_df(l[[1]],1,models[1,]$model)
df_rates

# dataframe with median and hpds of root height, x coordinate on plot, model name for 
# each simulation under model
df_heights  = prepare_df(l[[2]],1,models[1,]$model)
df_heights

# adds values from remaining models
for (i in 2:nrow(models)){
  l = sample_table(file.path(input_dir, models[i,]$folder))
  df_rates = cbind(df_rates,prepare_df(l[[1]],1,models[i,]$model))
  df_heights = cbind(df_heights,prepare_df(l[[2]],1,models[i,]$model))


}
colnames(df_rates) = seq(1, by=1, ncol(df_rates))
df_rates = as.data.frame(t(df_rates), stringsAsFactors = FALSE)

colnames(df_heights) = seq(1, by=1, ncol(df_heights))
df_heights = as.data.frame(t(df_heights), stringsAsFactors = FALSE)

df_rates
df_heights

ggplot(df_rates) + geom_point(aes(x=as.numeric(df_rates$x),y=as.numeric(df_rates$median))) + 
  geom_segment(aes(x=as.numeric(df_rates$x), y=as.numeric(df_rates$`95HPD%_lower`),
                   xend=as.numeric(df_rates$x), yend=as.numeric(df_rates$`95HPD%_higher`))) +
  facet_grid(~ model) + labs(y="mean Rate", x="") +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) +
  scale_y_continuous(labels = scales::scientific)
#  scale_y_discrete(breaks=seq(min(as.integer(df_rates$`95HPD%_lower`))-1,max(as.integer(df_rates$`95HPD%_higher`))+1,1))

ggplot(df_heights) + geom_point(aes(x=as.numeric(df_heights$x),y=as.numeric(df_heights$median))) + 
  geom_segment(aes(x=as.numeric(df_heights$x), y=as.numeric(df_heights$`95HPD%_lower`),
                   xend=as.numeric(df_heights$x), yend=as.numeric(df_heights$`95HPD%_higher`)))+ 
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + 
  facet_grid(~ model) + labs(y="treeHeight", x="") 
  
  #scale_y_continuous=c(61,63,67)


# folders with beast results for different models
#model_folders = c("random_single", "random_groups")
#input_dir = "D:/DATA/samplebias/biases/EV-A71_alignments/random_single"

