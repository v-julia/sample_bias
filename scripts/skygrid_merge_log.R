library(ggplot2)
library(grid)
library(gridExtra)



dir=setwd("D:\\R\\data\\skygrid")

# directory with folders that contain skygrid coordinates for different dataset preparation types
# outputs of each preparation type are located in distinct folders
input_dir = setwd("D:/DATA/samplebias/biases/EV-A71_alignments/skygrids_coord/")

# dataframe with coorespondence of folder names to model names
# the name of folder coresponds to the model of sampling generation
models = data.frame(c('random_single','random_groups', '3.0', '5.0', '2_1'),
                    c('random single', 'random groups', '3.0% identity threshold', '5% identity threshold', 
                      'max=2 seq, 1%'),
                    stringsAsFactors = F)
colnames(models) = c('folder', 'model')


colors = c('orange', 'red', 'blue', 'cyan', 'grey', 'brown', 'magenta', 'violet', 'khaki','green')


load_folder_data = function(folder_name, dir){
  dir = file.path(dir,folder_name)
  file_names <- list.files(dir)
  for (i in 1:length(file_names)){
    file_names[i] = file.path(dir, file_names[i])
    
  }
    
  print(file_names)
  collective_data_frame=lapply(file_names,read.table, sep = '\t',skip = 1,header = T)
  return(collective_data_frame)
}

all_data = lapply(models$folder, load_folder_data, input_dir)

plot_skygrid = function(col_df, title){
  a = ggplot()
  for (i in 1:length(col_df)){
    a = a +  geom_line(aes_string(x=col_df[[i]]$Time,y=log10(col_df[[i]]$Median)), color=colors[i] , alpha=0.8, size = 1)+
      geom_ribbon(aes_string(ymin = log10(col_df[[i]]$Lower), 
                             ymax = log10(col_df[[i]]$Upper), x=col_df[[i]]$Time), fill = colors[i], alpha=0.1) +
      ggtitle(title)
    
  }
  
  a=a+scale_x_continuous(name="Year")+
    scale_y_continuous(name="Population size")+   theme_bw() 
  return(a)
  
}

all_plots =  list()
for (i in 1:length(all_data)){
  all_plots[[i]] = plot_skygrid(all_data[[i]],models$model[i])
}

do.call(grid.arrange,all_plots)



