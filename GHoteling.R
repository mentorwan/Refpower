library(Hotelling)

get_pvalue <- function(case,control,case_vector,control_vector,d.l2.full,site_case,site_control) {

skin_case <- case[case_vector==1]
saliva_case <- case[case_vector==0]
saliva_control <- control[control_vector==0]
skin_control <- control[control_vector==1]

d.l2.case <- list()
if( length(saliva_case) > 0 )
{
  for(i in 1:length(saliva_case))
  {
      temp <- d.l2.full[((d.l2.full$id == saliva_case[i]) & (d.l2.full$site == site_control)),]
      d.l2.case <- rbind(d.l2.case,temp)
  }
}
if(length(skin_case) > 0 )
{
  for(i in 1:length(skin_case))
  {
    temp <- d.l2.full[((d.l2.full$id == skin_case[i]) & (d.l2.full$site == site_case)),]
    d.l2.case <- rbind(d.l2.case,temp)
  }
}

d.l2.control <- list()
for(i in 1:length(control))
{
  temp <- d.l2.full[((d.l2.full$id == saliva_control[i]) & (d.l2.full$site == site_control)),]
  d.l2.control <- rbind(d.l2.control,temp)
}


common_col_names <- c("id","case.dist.stool","case.dist.nasal","control.dist.stool","control.dist.nasal")

first <- cbind(d.l2.case[,c(1,5:6)],NA,NA)
names(first) <- common_col_names
fourth <- cbind(d.l2.control[,c(1)],NA,NA,d.l2.control[,c(5:6)])
names(fourth) <- common_col_names

#d.l2.final <- rbind(first,second,third,fourth)
d.l2.final <- rbind(first,fourth)
#row.names(d.l2.final) <- d.l2.final$id
d.l2.final <- d.l2.final[,-1]

#return(GHotelling(d.l2.final,nBoot=1000,print.details = F))

data <- d.l2.final
colnames(data) <- c('x1', 'y1', 'x2', 'y2')
n <- dim(data)[1]

t2 <- GHotelling.helper(data, print.details = F)
#result = ghotelling(data,nboot=1000)
#result$boot.summary[3,]       

mean.mx   <- matrix(rep(colMeans(data, na.rm = T), n), nrow = n, ncol = 4, byrow = T)
data.null <- data - mean.mx

boot.ix   <- get.boot.ix(data.null, 1000)
t2.boot   <- matrix(NA, nrow = 1000, ncol = 1)

for(i in 1:1000){
  t2.boot[i] <- GHotelling.helper(data = data.null[boot.ix[,i], ], print.details = F)
}

pval <- mean(t2.boot >= t2) # wrt original data

return(pval)

}



get_pvalue_standard_hotelling <- function(case,control,case_vector,control_vector,d.l2.full,site_case,site_control) {
  
  #rm(list=setdiff(ls(),c("case_id_1a","control_id_1a","case_id_1a_vector","control_id_1a_vector")))
  
  
  #case_vector <- case_vector
  #control_vector <- control_vector
  ##Normally skin is case, saliva is control
  
  skin_case <- case[case_vector==1]
  saliva_case <- case[case_vector==0]
  saliva_control <- control[control_vector==0]
  skin_control <- control[control_vector==1]
  
 
  d.l2.case <- list()
  if( length(saliva_case) > 0 )
  {
    for(i in 1:length(saliva_case))
    {
      temp <- d.l2.full[((d.l2.full$id == saliva_case[i]) & (d.l2.full$site == site_control)),]
      d.l2.case <- rbind(d.l2.case,temp)
    }
  }
  if(length(skin_case) > 0 )
  {
    for(i in 1:length(skin_case))
    {
      temp <- d.l2.full[((d.l2.full$id == skin_case[i]) & (d.l2.full$site == site_case)),]
      d.l2.case <- rbind(d.l2.case,temp)
    }
  }
  
  d.l2.control <- list()
  for(i in 1:length(control))
  {
    temp <- d.l2.full[((d.l2.full$id == saliva_control[i]) & (d.l2.full$site == site_control)),]
    d.l2.control <- rbind(d.l2.control,temp)
  }
  
  common_col_names <- c("id","dist.stool","dist.nasal")
  
  first <- cbind(d.l2.case[,c(1,5:6)])
  names(first) <- common_col_names
  first <- first[,-1]
  fourth <- cbind(d.l2.control[,c(1)],d.l2.control[,c(5:6)])
  names(fourth) <- common_col_names
  fourth <- fourth[,-1]
  
  t2_standard = hotelling.test(first,fourth)
  pval_standard = with(t2_standard$stats, 1-pf(m*statistic,df[1],df[2]))
  #return(hotelling.test(first,fourth)$pval)
  
  return(pval_standard)
}

get_pvalue_test <- function(case,control,case_vector,control_vector,d.l2.full,site_case,site_control) {
  
  #rm(list=setdiff(ls(),c("case_id_1a","control_id_1a","case_id_1a_vector","control_id_1a_vector")))
  
  
  #case_vector <- case_vector
  #control_vector <- control_vector
  ##Normally skin is case, saliva is control
  
  skin_case <- case[case_vector==1]
  saliva_case <- case[case_vector==0]
  saliva_control <- control[control_vector==0]
  skin_control <- control[control_vector==1]
  
  ###
  
  #load('mm-HMP_data_082816-5-sites-level2-nasal-ref-half-sample-mean.rdata')
  #require('repmis')
  #load('test-matrix.rdata')
  #test.mx[c(1:3,20:22,85:69),]
  
  
  #GHotelling(test.mx,nBoot=10000,print.details = F, seed = 1)
  #GHotelling(test.mx[-(7:72),], nBoot = 10000, print.details = F, seed = 1)
  
  #d.l2.hmp <- HMPdistance(tax.level = 'l2.phylum', d.new.filename = "HMP_L2_030317_phylum-example-input-file.csv", d.new.ix.col.not.rel.abu = 1:4, measure = 'bc')
  #d.l2.full <- cbind(d[,c(1:3,5,5:6)],d.l2.hmp)
  
  d.l2.case <- d.l2.full[((d.l2.full$id %in% skin_case) & (d.l2.full$site == site_case)) | ((d.l2.full$id %in% saliva_case) & (d.l2.full$site == site_control)) ,] 
  d.l2.control <- d.l2.full[(d.l2.full$id %in% control) & (d.l2.full$site == site_control),]
  d.l2.case <- d.l2.case[order(d.l2.case$id),]
  d.l2.control <- d.l2.control[order(d.l2.control$id),]
  
  #Common elements which have both case and control data
  #d.l2.case_control <- d.l2.full[((d.l2.full$id %in% skin_case) & (d.l2.full$site == site_control)) | ((d.l2.full$id %in% saliva_case) & (d.l2.full$site == site_case)),]
  #d.l2.control_case <- d.l2.full[(d.l2.full$id %in% control) & (d.l2.full$site == site_case),]
  #d.l2.case_control <- d.l2.case_control[order(d.l2.case_control$id),]
  #d.l2.control_case <- d.l2.control_case[order(d.l2.control_case$id),]
  
  #d.l2.case.common <- intersect(d.l2.case_control$id,d.l2.case$id)
  #d.l2.case.unique <- setdiff(d.l2.case$id, d.l2.case_control$id)
  #d.l2.control.common <- intersect(d.l2.control_case$id,d.l2.control$id)
  #d.l2.control.unique <- setdiff(d.l2.control$id, d.l2.control_case$id)
  
  common_col_names <- c("id","case.dist.stool","case.dist.nasal","control.dist.stool","control.dist.nasal")
  
  #first <- cbind(d.l2.case[(d.l2.case$id %in% d.l2.case.unique),c(1,5:6)],NA,NA)
  #names(first) <- common_col_names
  #second <- cbind(d.l2.case[(d.l2.case$id %in% d.l2.case.common),c(1,5:6)],d.l2.case_control[(d.l2.case_control$id %in% d.l2.case.common),c(5:6)])
  #names(second) <- common_col_names
  #third <- cbind(d.l2.control[(d.l2.control$id %in% d.l2.control.common),c(1,5:6)],d.l2.control_case[(d.l2.control_case$id %in% d.l2.control.common),c(5:6)])
  #names(third) <- common_col_names
  #fourth <- cbind(d.l2.control[(d.l2.control$id %in% d.l2.control.unique),c(1)],NA,NA,d.l2.control[(d.l2.control$id %in% d.l2.control.unique),c(5:6)])
  #names(fourth) <- common_col_names
  
  first <- cbind(d.l2.case[,c(1,5:6)],NA,NA)
  names(first) <- common_col_names
  fourth <- cbind(d.l2.control[,c(1)],NA,NA,d.l2.control[,c(5:6)])
  names(fourth) <- common_col_names
  
  #d.l2.final <- rbind(first,second,third,fourth)
  d.l2.final <- rbind(first,fourth)
  row.names(d.l2.final) <- d.l2.final$id
  d.l2.final <- d.l2.final[,-1]
  
  return(GHotelling(d.l2.final,nBoot=1000,print.details = F))
}
