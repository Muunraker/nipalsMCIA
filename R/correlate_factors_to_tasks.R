library("ggplot2")
library("reshape2")

# load the factors 
load("results/factors/factorizations.RData")

# load the task data 
fn = "../cmi-pb-pertussis/output/tasks/tasks.tsv"
task_list = read.table(fn, sep='\t', header = T)

########## Correlate the tasks to factors ##########
##### Merge all the task datasets #####
for (i in seq(1, nrow(task_list))){
  
  # loading task data
  task_info = task_list[i,]
  base = "../cmi-pb-pertussis/output/database_dump/day_splits_curated/"
  fn = sprintf("%s.pma.day%s.proc.tsv", task_info$assay, task_info$day)
  fn = paste0(base, fn)

  # skip abiters for now 
  #if (task_info$assay == 'abtiters'){
  #  next
  #}
  
  # skip fim2/3
  if (task_info$fetchname == 'fim2_3'){
    next 
  }
  
  # extracting the specific column
  task_data = read.table(fn, header = T, row.names = 1)
  if (task_info$assay == 'abtiters'){
    col_name = paste0('IgG.', toupper(task_info$fetchname))
  }
  else{
    col_name = task_info$fetchname
  }
  
  print(paste(i, fn))
  print(paste0('col_name:', col_name))
  
  task_data = data.frame(task_data[, col_name])
  colnames(task_data) = paste0(task_info$shortname, '_day', task_info$day)
  
  # merging task datae
  if (i == 1){
    merged_tasks = task_data
  }
  else{
    merged_tasks = merge(merged_tasks, task_data, by=0)
    row.names(merged_tasks) = merged_tasks$Row.names
    merged_tasks = subset(merged_tasks, select=-c(Row.names))
  }
}

# correlate the tasks and factors
for (i in seq(1,5)){
  print(i)
  # harmonize the data
  method = out$method[[i]]
  factors = out$factorizations[[i]][[1]]
  row.names(factors) = gsub('X', '', row.names(factors))
  shared_samples = intersect(row.names(factors), row.names(merged_tasks))
  final_factors = factors[shared_samples, ]
  final_tasks = merged_tasks[shared_samples, ]

  # correlate the factors and task
  factor_corrs = cor(final_factors, final_tasks)
  factor_corrs = factor_corrs[ , colSums(is.na(factor_corrs)) == 0]
  factor_corrs = melt(factor_corrs)
  colnames(factor_corrs) = c('comp', 'task', 'comp_val')

  # make a heatmap of the correlations
  p = ggplot(factor_corrs, aes(comp, task, fill= comp_val)) + 
    geom_tile()
  p = p + ggtitle(method) 
  #p = p + scale_x_continuous(breaks=seq(1, 10, 1))
  p = p + scale_fill_gradient(low="blue", high="red")
  print(p)
}

