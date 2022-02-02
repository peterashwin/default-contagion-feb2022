#library(copula)
source("SG_CopulasLecceCourse.R")
copulas_list = c('amh','clayton','frank','gumbel','joe','normal','t','galambos','tawn','fgm','plackett') #huslerReiss does not work
copulas_to_evaluate = c(1:11) #1:11
optMeth_list = c("Nelder-Mead","BFGS","CG","L-BFGS-B","SANN","Brent")
optMeth_to_evaluate = c(1) #1:6, Brent works well with amh, NM with gumbel

#choose which data to model
data_type = c("a","s") #a - our model, #s - spiliopoulos
data_type_to_evaluate = c(1) #c(1:2)
data_number = c("01","02","03","04","05","06","07","08","09","10","11")
# Study multiple data sets in the same time - for instance coming from models with different parameter beta.
data_number_to_evaluate = c(1) #c(1:11)
# convert bt_et.m and bt_et_s.m generated in Matlab into csv files, each for different parameter set. 
# files a_01.csv corresponds to the first sample from bt_et.m whereas s_01.csv to the first sample from bt_et_s.m

#choose whether to evaluate with basic copula or survival
survival = TRUE #TRUE or FALSE
standarise = FALSE #TRUE or FALSE

cat("\014") #clear console

t0 = Sys.time()

df = data.frame( data_type=character(), beta_seq=character(), cop=character(), meth=character(),
                 survival=double(),     theta=double(),       pv=double(),     pv_over_0=double(),
                 tau=double(),          rho=double(),         stringsAsFactors=FALSE)

for(i in copulas_to_evaluate){
  i_iter = which(i==copulas_to_evaluate)
  for(j in optMeth_to_evaluate){
    j_iter = which(j==optMeth_to_evaluate)
    for(k in data_type_to_evaluate){
      k_iter = which(k==data_type_to_evaluate)
      for(l in data_number_to_evaluate){
        l_iter = which(l==data_number_to_evaluate)
        cat("\n\n\n")
        print("----------------------------------------------------")
        print(paste("Copula: ",i_iter,"/",length(copulas_to_evaluate),", data type: ",k_iter,"/",length(data_type_to_evaluate),", set no: ",l_iter,"/",length(data_number_to_evaluate),sep=""))
        data_file = paste0("data/",data_type[k],"_",data_number[l],".csv")
        
        results = GS_CopulaCoursePhD_R(copulas_list[i],optMeth_list[j],paste(data_type[k],data_number[l],copulas_list[i],optMeth_list[j],sep='_'),data_file,survival,standarise)
        
        work_done = (i_iter-1+(j_iter-1+(k_iter-1+l_iter/length(data_number_to_evaluate))/length(data_type_to_evaluate))/length(optMeth_to_evaluate))/length(copulas_to_evaluate)
        work_remaining = 1-work_done
        t1 = Sys.time()
        
        print(paste0("Time elapsed:   ",sprintf("%.3f",difftime(t1,t0,units="hours")),"hrs."))
        print(paste0("Time remaining: ",sprintf("%.3f",difftime(t1,t0,units="hours")*work_remaining/work_done),"hrs."))
        print("----------------------------------------------------")
        
        df[nrow(df)+1,] = list(data_type[k], data_number[l], copulas_list[i], optMeth_list[j],
                               ifelse(survival,1,0), results$theta, results$pv, ifelse(results$pv>=0.05,1,0),
                               results$tau, results$rho)
      }
    }
  }
}

clear_table = TRUE
if(!clear_table){
  df_old = read.csv(file="data/copulas_results.csv")
  df_old = df_old[,(ncol(df_old)-ncol(df)+1):ncol(df_old)]
  df = rbind(df_old,df)
  print('CSV file updated')
}
write.csv(df, file="data/copulas_results.csv")
