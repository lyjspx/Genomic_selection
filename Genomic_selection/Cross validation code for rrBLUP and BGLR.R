setwd("C:/Users/evan.salsman.2/Documents/AYT Data/Forward prediction rrBLUP/")

#load phenotype data
BLUP.AYT12_17=read.delim("2012-2017 agronomic and quality BLUPs with env fixed.txt", header=T, na.string="NA")
dim(BLUP.AYT12_17)

GBS50RF.AYT12_17=read.delim("AYT12-17_1413GBSna50AB_MAF5_LDKNNi_F-101_meanImpute.txt", header=T, sep=",", na.string="NA")
dim(GBS50RF.AYT12_17)

#GS model_3k 90vs10 cross-validation 12-16AYT
phenotype=as.matrix(BLUP.AYT12_17[1:1413,-c(1:2)])
dim(phenotype)
phenotype[1:3,]
Markers_impute=GBS50RF.AYT12_17[,-1]
dim(Markers_impute)
Markers_impute[1:3,1:5]

######Run this command in first time 
# install.packages("foreach")
# install.packages("doParallel")
library("foreach")
library("doParallel")

###rrBLUP cross-validation 12-16AYT

traits=9
cycles=100

# 
# for(r in 1:cycles)
# {
#   train=as.matrix(sample(1:1413,1271))
#   test=setdiff(1:1184, train)
#   Pheno_train=phenotype[train,]
#   m_train=Markers_impute[train,]
#   Pheno_valid=phenotype[test,]
#   m_valid=Markers_impute[test,]
#   for(i in 1:traits)
#   {
#     y_train=Pheno_train[,i]
#     y_answer=mixed.solve(y_train,Z=m_train, K=NULL, SE=FALSE, return.Hinv=FALSE)
#     summary(y_answer)
#     y_e=as.matrix(y_answer$u)
#     head(y_e)
#     pred_y_valid=m_valid %*%y_e
#     y_valid=Pheno_valid[,i]
#     #accuracy_rrBLUP[r,i]=cor(pred_y_valid, y_valid, use="pairwise.complete.obs")
#     cor(pred_y_valid, y_valid, use="pairwise.complete.obs")
#   }
# }


rrBLUP_cross_validation = function(num_traits,marker_rrBLUP,pheno_rrBLUP){
    one_cycle = c()
    train=as.matrix(sample(1:1413,1271))
    test=setdiff(1:1184, train)
    Pheno_train=pheno_rrBLUP[train,]
    m_train=marker_rrBLUP[train,]
    Pheno_valid=pheno_rrBLUP[test,]
    m_valid=marker_rrBLUP[test,]
    for(i in 1:num_traits)
    {
      y_train=Pheno_train[,i]
      y_answer=mixed.solve(y_train,Z=m_train, K=NULL, SE=FALSE, return.Hinv=FALSE)
      y_e=as.matrix(y_answer$u)
      pred_y_valid=m_valid %*%y_e
      y_valid=Pheno_valid[,i]
      
      one_cycle = c(one_cycle, cor(pred_y_valid, y_valid, use="pairwise.complete.obs"))
    }
    return(one_cycle)
}

num_cores = detectCores()
registerDoParallel(num_cores)

rrBLUP_result_list = foreach(m = seq(1,cycles),.packages = "rrBLUP" ,.combine = c, .inorder = T) %dopar%
  rrBLUP_cross_validation(num_traits = traits, marker_rrBLUP = Markers_impute, pheno_rrBLUP = phenotype)
  
rrBLUP_result = matrix(rrBLUP_result_list, ncol = traits,byrow = T)
stopImplicitCluster()

write.table(rrBLUP_result, file="accuracy_AYT12-17_GBSna50_rrBLUP_90vs10_agronomicT.txt", sep=",")

###BGLR
traits=9
cycles=100

#Wrap one cycle into a function
bays_cross_validation = function(pheno_bays, marker_bays,num_trait){
  test=as.matrix(sample(1:1413,142))
  Pheno_train=pheno_bays
  Pheno_train[test,1:9]=NA
  one_cycle = c()
  for(i in 1:num_trait)
  {
    y_train=Pheno_train[,i]
    fm.Ht_W.BA=BGLR(y=y_train, ETA=list(list(X=marker_bays, model='BayesA', saveEffects=TRUE)), nIter=12000, thin=2, burnIn=500,verbose = F)
    one_cycle=c(one_cycle,cor(fm.Ht_W.BA$yHat[test], phenotype[test,i], use="pairwise.complete.obs"))
  }
  return(one_cycle)
}

#####Parallel mode
num_cores = detectCores()
registerDoParallel(num_cores)

bays_result_list = foreach(m = seq(1,cycles),.packages = "BGLR" ,.combine = c, .inorder = T) %dopar%
  bays_cross_validation(pheno_bays = phenotype,marker_bays = Markers_impute,num_trait = traits)
  
bays_result = matrix(bays_result_list, ncol = traits,byrow = T) ###You may inspect result here
stopImplicitCluster()
######Parallel end

write.table(bays_result, file="accuracy_AYT12-17_GBSna50_rrBLUP_90vs10_agronomicT.txt", sep=",")




# for(r in 1:cycles)
# {
#   test=as.matrix(sample(1:1413,142))
#   Pheno_train=phenotype
#   Pheno_train[test,1:9]=NA
#   for(i in 1:traits)
#   {
#     y_train=Pheno_train[,i]
#     fm.Ht_W.BA=BGLR(y=y_train, ETA=list(list(X=Markers_impute, model='BayesA', saveEffects=TRUE)), nIter=12000, thin=2, burnIn=500,verbose = F)
#     accuracy_Bayes[r,i]=cor(fm.Ht_W.BA$yHat[test], phenotype[test,i], use="pairwise.complete.obs")
#   }
# }
# write.table(accuracy_Bayes, file="accuracy_AYT12-17_GBSna50_rrBLUP_90vs10_agronomicT.txt", sep=",")


