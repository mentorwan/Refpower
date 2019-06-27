library(vegan)
library(reshape2)

library(MiRKAT)
library(GUniFrac)

#Bray_Curtis <- readDM("bray_curtis_HMP_YW_even1000.txt")


#Export report to csv file
write_report <- function(list_name, file_name, current)
{
  df <- data.frame(matrix(unlist(list_name),nrow=1000,byrow=T))
  colnames(df) <- names(current)
  write.table(df,file=file_name,col.names=TRUE, row.names=FALSE, quote=FALSE)
}


#Generate permanova pvalue and return R2,w2 and pvalue
get_permanova_value <- function(case_id,control_id,case_id_vector,control_id_vector,full,site_case,site_control,DM)
{
  case_sample_id = case_id
  for(j in 1:length(case_id))
  {
    if(case_id_vector[j] == 1)
    {
      case_sample_id[j] <- full[full$id == case_id[j] & full$site == site_case,]$sid 
    }
    else if(case_id_vector[j] == 0)
    {
      case_sample_id[j] <- full[full$id == case_id[j] & full$site == site_control,]$sid 
    }
  }
  
  
  control_sample_id = control_id
  for( j in 1:length(control_id))
  {
    if(control_id_vector[j] == 1)  ##skin
    {
      control_sample_id[j] <- full[full$id == control_id[j] & full$site == site_case,]$sid 
    }
    else if (control_id_vector[j] == 0 ) ##saliva
    {
      control_sample_id[j] <- full[full$id == control_id[j] & full$site == site_control,]$sid 
    }
  }
  
  all_sample_id <- c(case_sample_id,control_sample_id)
  all_sample_col_id <- c(paste0("g1s",1:length(case_id)),paste0("g2s",1:(length(control_id))))
 
  all_rarified_sample_id <- all_sample_id[all_sample_id %in% row.names(DM)]
  
  all_rarified_sample_col_id <- all_sample_col_id[all_sample_id %in% row.names(DM)]
  
  all_rarified_sample_id <- as.character(all_rarified_sample_id)
  
  
  
  DM_subset <- DM[all_rarified_sample_id,all_rarified_sample_id]
  
  #DM_subset <- as.matrix(DM_subset)
  
  ###change row and column labels to group and sample name.
  #colnames(DM_subset) <- all_sample_col_id[match(colnames(DM_subset),all_sample_id)]
  #rownames(DM_subset) <- all_sample_col_id[match(rownames(DM_subset),all_sample_id)]
  
  #colnames(DM_subset) <- all_sample_col_id[match(all_rarified_sample_id,all_sample_id)]
  #rownames(DM_subset) <- all_sample_col_id[match(all_rarified_sample_id,all_sample_id)]
  
  #####
  
  colnames(DM_subset) <-  all_rarified_sample_col_id
  rownames(DM_subset) <-  all_rarified_sample_col_id

  
  ######
  
  
  #PERMANOVA(DM_subset)
  R2=calcR2(DM_subset)
  w2=calcOmega2(DM_subset)
  pvalue=calcPERMANOVAp(DM_subset)
  
  #cat("R2=",R2,"\t")    #R2 value
  #cat("w2=",w2,"\t")    #w2 value
  #cat("pvalue=",pvalue,"\n")   #p value
  
  return(c(w2,pvalue))
  
}

#Generate Mirkat P-value.
get_mirkat_value <- function(case_id,control_id,case_id_vector,control_id_vector,full,site_case,site_control,DM)
{
  case_sample_id = case_id
  for(j in 1:length(case_id))
  {
    if(case_id_vector[j] == 1)
    {
      case_sample_id[j] <- full[full$id == case_id[j] & full$site == site_case,]$sid 
    }
    else if(case_id_vector[j] == 0)
    {
      case_sample_id[j] <- full[full$id == case_id[j] & full$site == site_control,]$sid 
    }
  }
  
  
  control_sample_id = control_id
  for( j in 1:length(control_id))
  {
    if(control_id_vector[j] == 1)  ##skin
    {
      control_sample_id[j] <- full[full$id == control_id[j] & full$site == site_case,]$sid 
    }
    else if (control_id_vector[j] == 0 ) ##saliva
    {
      control_sample_id[j] <- full[full$id == control_id[j] & full$site == site_control,]$sid 
    }
  }
  
  all_sample_id <- c(case_sample_id,control_sample_id)
  all_sample_col_id <- c(paste0("g1s",1:length(case_id)),paste0("g2s",1:(length(control_id))))
  
  all_rarified_sample_id <- all_sample_id[all_sample_id %in% row.names(DM)]
  all_rarified_sample_id <- as.character(all_rarified_sample_id)
  
  all_rarified_sample_col_id <- all_sample_col_id[all_sample_id %in% row.names(DM)]
  
  DM_subset <- DM[all_rarified_sample_id,all_rarified_sample_id]
  
  K.BC = D2K(DM_subset)  #Kernal from MiRKAT package
  
  Case_Control = (grepl("g1",all_rarified_sample_col_id))**2  # Case is skin and control is saliva.
  
  result=MiRKAT(y = Case_Control, Ks = K.BC, out_type = "D",  method = "permutation", nperm = 1000)

  return(result)
  
}



