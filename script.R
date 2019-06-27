#Script to run power analysis

rm(list=ls())
setwd("/data/wany/Mitch/Power/Power_code_101212/All_level/skin_nasal_update_062719/")

load("skin_nasal.RData")
rm(d.l2.full)
rm(full)

load("../all_level.RData")
load("../Bray_Curtis_L2_new.RData")
load("../Bray_Curtis_L3_new.RData")
load("../Bray_Curtis_L4_new.RData")
load("../Bray_Curtis_L5_new.RData")
load("../Bray_Curtis_L6_new.RData")

source("GHoteling.R")
source("microbiome-fixed-reference-functions.r")
source("micropower.R")
source("Permanova.R")
source("run_test.R")


args = commandArgs(trailingOnly = TRUE)

iter = 1000
split = args[1]
scenario = args[2]

#cat(args[1], "\t",args[2],"\t",args[3],"\n")

split = args[1]
scenario = args[2]
#output_file = args[3]

d.l2.full <- d.l3.full[d.l3.full$sid %in% rownames(Bray_Curtis_L2_new),]
#Bray_Curtis_old <- readDM("bray_curtis_HMP_YW_even1000.txt")
Bray_Curtis <- Bray_Curtis_L3_new

skin_nasal = run_test("skin","nasal",split, scenario, iter)
df_1 <- data.frame(matrix(unlist(skin_nasal),nrow=1000,byrow=T))

skin_nasal = run_test("skin","nasal",split, scenario, iter)
df_2 <- data.frame(matrix(unlist(skin_nasal),nrow=1000,byrow=T))

skin_nasal = run_test("skin","nasal",split, scenario, iter)
df_3 <- data.frame(matrix(unlist(skin_nasal),nrow=1000,byrow=T))

#matrix(c(sum(df$X1<0.01), sum(df$X1< 0.05), sum(df$X2 < 0.01), sum(df$X2 < 0.05), sum(df$X3 < 0.01), sum(df$X3 < 0.05), sum(df$X4 < 0.01), sum(df$X4 < 0.05)),4,2, byrow= TRUE)

##Compare power for two test

matrix(c(sum(df_1$X1< 0.05), sum(df_2$X1<0.05), sum(df_3$X1< 0.05),
         sum(df_1$X2< 0.05), sum(df_2$X2<0.05), sum(df_3$X2< 0.05),
         sum(df_1$X3< 0.05), sum(df_2$X3<0.05), sum(df_3$X3< 0.05),
         sum(df_1$X4< 0.05), sum(df_2$X4<0.05), sum(df_3$X4< 0.05)), 
        4,3, byrow= TRUE)

options(scipen=100,digits = 4)
result <- c(sum(df_1$X1< 0.05), sum(df_2$X1<0.05), sum(df_3$X1< 0.05),
            sum(df_1$X2< 0.05), sum(df_2$X2<0.05), sum(df_3$X2< 0.05),
            sum(df_1$X3< 0.05), sum(df_2$X3<0.05), sum(df_3$X3< 0.05),
            sum(df_1$X4< 0.05), sum(df_2$X4<0.05), sum(df_3$X4< 0.05))

result <- c(result,mean(df_1$X5),mean(df_2$X5),mean(df_3$X5))

Ghoteling_reject <- sum(result[1:3])/3000
Standard_hotelling_reject <- sum(result[4:6])/3000
PERMANOVA_reject <- sum(result[7:9])/3000
MIRKat_reject <- sum(result[10:12])/3000
w2_average <- mean(result[13:15])

output <- c(Ghoteling_reject,Standard_hotelling_reject,PERMANOVA_reject,MIRKat_reject,w2_average)

#save(result,output,)


#format(c(power_test1,power_test2,p),scientific = FALSE)

mcnemar_result <- c(mcnemar_test("X1","X3"),mcnemar_test("X2","X3"),
                    mcnemar_test("X1","X4"),mcnemar_test("X2","X4"),
                    mcnemar_test("X1","X2"),mcnemar_test("X3","X4"))

test <- c(result,output,mcnemar_result)
cat(test,"\n",file=args[3])

save(df_1,df_2,df_3,file=args[4])



  
  

