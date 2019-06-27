run_test <- function(site_case,site_control,nsample,scenario,iter)
{
  nsample=as.numeric(nsample)
  skin_saliva_a <- vector("list",length=iter)
  for(i in 1:iter)
  {
    
    if(scenario == 6)   #0% cases
    {
      case_id_a <- sample(control_id,nsample,replace=TRUE)
      case_id_a_vector <- rep(0,nsample)
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
    }
    
    if(scenario == 5) #20% cases
    {
      control_sample = as.integer(0.8*nsample)
      case_sample = as.integer(0.2*nsample)
      temp <- sample(control_id,control_sample,replace=TRUE)
      case_id_a <- c(sample(case_id,case_sample,replace=TRUE),temp)
      case_id_a_vector <- c(rep(1,case_sample),rep(0,control_sample))
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
    }
    
    if(scenario == 4)  #40% cases
    {
      control_sample = as.integer(0.6*nsample)
      case_sample = as.integer(0.4*nsample)
      temp <- sample(control_id,control_sample,replace=TRUE)
      case_id_a <- c(sample(case_id,case_sample,replace=TRUE),temp)
      case_id_a_vector <- c(rep(1,case_sample),rep(0,control_sample))
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
      
    }
    
    if(scenario == 3)  #60% cases
    {
      control_sample = as.integer(0.4*nsample)
      case_sample = as.integer(0.6*nsample)
      temp <- sample(control_id,control_sample,replace=TRUE)
      case_id_a <- c(sample(case_id,case_sample,replace=TRUE),temp)
      case_id_a_vector <- c(rep(1,case_sample),rep(0,control_sample))
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
      
    }
    
    if(scenario == 2)   #80% cases
    {
      control_sample = as.integer(0.2*nsample)
      case_sample = as.integer(0.8*nsample)
      temp <- sample(control_id,control_sample,replace=TRUE)
      case_id_a <- c(sample(case_id,case_sample,replace=TRUE),temp)
      case_id_a_vector <- c(rep(1,case_sample),rep(0,control_sample))
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
      
    } 
    
    if(scenario == 1)   #100% cases
    {
      case_id_a <- sample(case_id,nsample,replace=TRUE)
      case_id_a_vector <- c(rep(1,nsample))
      control_id_a <- sample(control_id,nsample,replace=TRUE)
      control_id_a_vector <- rep(0,nsample)
      
    }
    
    results_standard <- get_pvalue(case_id_a,control_id_a,case_id_a_vector,control_id_a_vector,d.l2.full,site_case,site_control)
    results_standard_hotelling <- get_pvalue_standard_hotelling(case_id_a,control_id_a,case_id_a_vector,control_id_a_vector,d.l2.full,site_case,site_control)
    results_permanov <- get_permanova_value(case_id_a,control_id_a,case_id_a_vector,control_id_a_vector,full,site_case,site_control,Bray_Curtis)
    results_mirkat <- get_mirkat_value(case_id_a,control_id_a,case_id_a_vector,control_id_a_vector,full,site_case,site_control,Bray_Curtis)
    cat(i,":pval=",results_standard,"\t")
    cat(results_standard_hotelling,"\t")
    cat(results_permanov[2],"\t")
    cat(results_mirkat,"\t")
    cat("w2=",results_permanov[1],"\n")
    
    #current <- as.list(c(pval_g=results_standard, pval_s=results_standard_hotelling))
    current <- as.list(c(pval_g=results_standard, pval_s=results_standard_hotelling, pval_p=results_permanov[2],pval_m=results_mirkat,w2=results_permanov[1]))
    #current <- as.list(c(pval_Ghotelling=results_mirkat))
    skin_saliva_a[[i]] <- current
    
  }
  return(skin_saliva_a)
}

mcnemar_test <- function(string1,string2)
{
  df_1_test1 <- eval(parse(text=paste0("df_1$",string1)))
  df_1_test2 <- eval(parse(text=paste0("df_1$",string2)))
  df_2_test1 <- eval(parse(text=paste0("df_2$",string1)))
  df_2_test2 <- eval(parse(text=paste0("df_2$",string2)))
  df_3_test1 <- eval(parse(text=paste0("df_3$",string1)))
  df_3_test2 <- eval(parse(text=paste0("df_3$",string2)))
  
  
  A= sum(df_1_test1 < 0.05 & df_1_test2 < 0.05) + sum(df_2_test1 < 0.05 & df_2_test2 < 0.05) + sum(df_3_test1 < 0.05 & df_3_test2 < 0.05)
  B= sum(df_1_test1 < 0.05 & df_1_test2 >= 0.05) + sum(df_2_test1 < 0.05 & df_2_test2 >= 0.05) + sum(df_3_test1 < 0.05 & df_3_test2 >= 0.05)
  C = sum(df_1_test1 >= 0.05 & df_1_test2 < 0.05) + sum(df_2_test1 >= 0.05 & df_2_test2 < 0.05) + sum(df_3_test1 >= 0.05 & df_3_test2 < 0.05)
  D= sum(df_1_test1 >= 0.05 & df_1_test2 >= 0.05) + sum(df_2_test1 >= 0.05 & df_2_test2 >= 0.05) + sum(df_3_test1 >= 0.05 & df_3_test2 >= 0.05)
  
  matrix(c(A,B,C,D),2,2,byrow=T)
  
  mcmear= (B-C)^2/(B+C)
  
  p= pchisq(mcmear,df=1,ncp=0,lower.tail=FALSE,log.p=FALSE)
  
  return(p)
  
}
