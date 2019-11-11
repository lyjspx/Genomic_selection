getwd()
setwd("Y:/Yuan/GSwithRelative")
list.files()

#####rrBLUP with relatives
library("rrBLUP")
library("snpReady")

phenotype=read.table("2012-2017 BLUPs for all traits(1).txt",header = T,stringsAsFactors = F,row.names = 1)
phenotype[1:5,1:5]

geneticInfo=read.table("AYT12-17_1413GBSna50AB_MAF5_LDKNNi_F-101_meanImpute(1).txt",header = T,row.names = 2,na.strings = "NA",sep = ",")
geneticInfo[1:5,1:5]

geneticInfo=geneticInfo[,2:13162]
geneticInfo[1:5,1:5]

intersect(row.names(phenotype),row.names(geneticInfo)) 
individualList=intersect(row.names(phenotype),row.names(geneticInfo))#get individual intersection

phenotype=phenotype[individualList,]#selection pheno and geno of individualList
geneticInfo=geneticInfo[individualList,]

#No missing in phenotypes, skip this part
#validIndividualInTrait=row.names(phenotype)[(complete.cases(phenotype))]
#length(validIndividualInTrait)

geneticInfoImputed=A.mat(geneticInfo,return.imputed = T)$imputed#impute missing

problemeticMarker=which(colMeans(geneticInfoImputed)<0.000001)
gMatrix=G.matrix(M=geneticInfoImputed[,-problemeticMarker],method = "VanRaden",format = "wide",plot = F)
gMatrixrrBLUP=A.mat(geneticInfo,return.imputed = T)$A
gMatrixAdditive=gMatrix

geneticInfo=geneticInfo[,-problemeticMarker]#delete three markers because of G.matrix()

predictBreedingValue=function(TPgeno,TPpheno,VPgeno){
  predictedRidgeRegression=mixed.solve(TPpheno,Z=TPgeno)
  predictedValue=VPgeno%*%predictedRidgeRegression$u
  return(predictedValue)
}

#17
testing=individualList[grep(pattern = "AYT17",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT17",individualList)]

result=c()
for(trait in colnames(phenotype)[2:8]){
  oneTrait=c(trait)
  for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
    oneTrait=c(oneTrait,numTopIndividual)
    for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
      cumulativeCor=0
      for(rep in 1:100){#define number of replicationss
        geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
        geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
        gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
        gMatrixAdditive=gMatrix$Ga
        averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
        averageRelativity=averageRelativity[training]
        individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
        predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
        corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
        cumulativeCor=cumulativeCor+corelation
      }
      oneTrait=c(oneTrait,cumulativeCor/100)
    }
    predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
    corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
    oneTrait=c(oneTrait,13158,corelation)
  }
  result=rbind(result,oneTrait)
}
write.csv(result,file = "2017 all traits with one reduced model.csv")

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=7
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[3:3], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
    {oneTrait=c(trait)
      for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
        oneTrait=c(oneTrait,numTopIndividual)
        for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
          cumulativeCor=0
          for(rep in 1:100){#define number of replicationss
            geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
            geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
            gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
            gMatrixAdditive=gMatrix$Ga
            averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
            averageRelativity=averageRelativity[training]
            individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
            predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
            corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
            cumulativeCor=cumulativeCor+corelation
          }
          oneTrait=c(oneTrait,cumulativeCor/100)
        }
        predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
        corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
        oneTrait=c(oneTrait,13158,corelation)
      }
      oneTrait
    }
  
###  
AYT17Result=rbind(test_para,AYT17Result)
write.csv(AYT17Result,file="AYT17_hundred_rep_all_traits.csv",row.names = F,quote = F)




#16
testing=individualList[grep(pattern = "AYT16",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT16",individualList)]

result=c()
for(trait in colnames(phenotype)[2:8]){
  oneTrait=c(trait)
  for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
    oneTrait=c(oneTrait,numTopIndividual)
    for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
      cumulativeCor=0
      for(rep in 1:100){#define number of replicationss
        geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
        geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
        gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
        gMatrixAdditive=gMatrix$Ga
        averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
        averageRelativity=averageRelativity[training]
        individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
        predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
        corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
        cumulativeCor=cumulativeCor+corelation
      }
      oneTrait=c(oneTrait,cumulativeCor/100)
    }
    predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
    corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
    oneTrait=c(oneTrait,13158,corelation)
  }
  result=rbind(result,oneTrait)
}

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=7
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[2:8], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
{oneTrait=c(trait)
for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
  oneTrait=c(oneTrait,numTopIndividual)
  for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
    cumulativeCor=0
    for(rep in 1:100){#define number of replicationss
      geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
      geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
      gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
      gMatrixAdditive=gMatrix$Ga
      averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
      averageRelativity=averageRelativity[training]
      individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
      predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
      corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
      cumulativeCor=cumulativeCor+corelation
    }
    oneTrait=c(oneTrait,cumulativeCor/100)
  }
  predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
  corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
  oneTrait=c(oneTrait,13158,corelation)
}
oneTrait
}
write.csv(test_para,file="AYT16_hundred_rep_all_traits.csv",row.names = F,quote = F)
###  

#15
testing=individualList[grep(pattern = "AYT15",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT15",individualList)]

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=20
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[2:8], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
{oneTrait=c(trait)
for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
  oneTrait=c(oneTrait,numTopIndividual)
  for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
    cumulativeCor=0
    for(rep in 1:100){#define number of replicationss
      geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
      geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
      gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
      gMatrixAdditive=gMatrix$Ga
      averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
      averageRelativity=averageRelativity[training]
      individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
      predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
      corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
      cumulativeCor=cumulativeCor+corelation
    }
    oneTrait=c(oneTrait,cumulativeCor/100)
  }
  predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
  corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
  oneTrait=c(oneTrait,13158,corelation)
}
oneTrait
}
write.csv(test_para,file="AYT15_hundred_rep_all_traits.csv",row.names = F,quote = F)
###  

#14
testing=individualList[grep(pattern = "AYT14",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT14",individualList)]

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=20
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[2:8], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
{oneTrait=c(trait)
for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
  oneTrait=c(oneTrait,numTopIndividual)
  for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
    cumulativeCor=0
    for(rep in 1:100){#define number of replicationss
      geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
      geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
      gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
      gMatrixAdditive=gMatrix$Ga
      averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
      averageRelativity=averageRelativity[training]
      individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
      predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
      corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
      cumulativeCor=cumulativeCor+corelation
    }
    oneTrait=c(oneTrait,cumulativeCor/100)
  }
  predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
  corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
  oneTrait=c(oneTrait,13158,corelation)
}
oneTrait
}
write.csv(test_para,file="AYT14_hundred_rep_all_traits.csv",row.names = F,quote = F)
###  


#13
testing=individualList[grep(pattern = "AYT13",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT13",individualList)]

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=20
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[2:8], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
{oneTrait=c(trait)
for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
  oneTrait=c(oneTrait,numTopIndividual)
  for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
    cumulativeCor=0
    for(rep in 1:100){#define number of replicationss
      geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
      geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
      gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
      gMatrixAdditive=gMatrix$Ga
      averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
      averageRelativity=averageRelativity[training]
      individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
      predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
      corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
      cumulativeCor=cumulativeCor+corelation
    }
    oneTrait=c(oneTrait,cumulativeCor/100)
  }
  predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
  corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
  oneTrait=c(oneTrait,13158,corelation)
}
oneTrait
}
write.csv(test_para,file="AYT13_hundred_rep_all_traits.csv",row.names = F,quote = F)
###  



#13
testing=individualList[grep(pattern = "AYT12",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT12",individualList)]

####Try to parallel
library("foreach")
library("doParallel")
numOfCore=20
registerDoParallel(cores=numOfCore)
test_para = foreach(trait = colnames(phenotype)[2:8], .combine =rbind,.packages = c("snpReady","rrBLUP","rgl"), .errorhandling = "pass", .inorder = T) %dopar%
{oneTrait=c(trait)
for(numTopIndividual in seq(from=50,to=length(training),by=100)){#define number of lines in trainging
  oneTrait=c(oneTrait,numTopIndividual)
  for(markerNum in c(100,500,1000,2000,3000,4000,5000,10000)){#define number of markers
    cumulativeCor=0
    for(rep in 1:100){#define number of replicationss
      geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
      geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
      gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
      gMatrixAdditive=gMatrix$Ga
      averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
      averageRelativity=averageRelativity[training]
      individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
      predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
      corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
      cumulativeCor=cumulativeCor+corelation
    }
    oneTrait=c(oneTrait,cumulativeCor/100)
  }
  predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
  corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")#calculate accuracy with all markers
  oneTrait=c(oneTrait,13158,corelation)
}
oneTrait
}
write.csv(test_para,file="AYT12_hundred_rep_all_traits.csv",row.names = F,quote = F)
###  

# #16
# testing=individualList[grep(pattern = "AYT16",individualList)]s
# training=individualList[-grep(pattern = "AYT16",individualList)]
# result=c()
# for(trait in colnames(phenotype)[2:8]){
#   oneTrait=c(trait)
#   for(numTopIndividual in seq(from=50,to=length(training),by=50)){
#     oneTrait=c(oneTrait,numTopIndividual)
#     for(markerNum in c(500,1000,2000,3000,4000)){
#       cumulativeCor=0
#       for(rep in 1:10){
#         geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
#         geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
#         #problemeticMarker=which(colMeans(geneticInfoImputed)<0.0001)
#         gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
#         gMatrixAdditive=gMatrix$Ga
#         averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
#         averageRelativity=averageRelativity[training]
#         individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
#         predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
#         corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#         cumulativeCor=cumulativeCor+corelation
#       }
#       oneTrait=c(oneTrait,cumulativeCor/10)
#     }
#     predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
#     corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#     oneTrait=c(oneTrait,13158,corelation)
#   }
#   result=rbind(result,oneTrait)
# }
# write.csv(result,file = "2016 all traits with one reduced model.csv")
# 
# #15
# testing=individualList[grep(pattern = "AYT15",individualList)]
# training=individualList[-grep(pattern = "AYT15",individualList)]
# result=c()
# for(trait in colnames(phenotype)[2:8]){
#   oneTrait=c(trait)
#   for(numTopIndividual in seq(from=50,to=length(training),by=50)){
#     oneTrait=c(oneTrait,numTopIndividual)
#     for(markerNum in c(500,1000,2000,3000,4000)){
#       cumulativeCor=0
#       for(rep in 1:10){
#         geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
#         geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
#         #problemeticMarker=which(colMeans(geneticInfoImputed)<0.0001)
#         gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
#         gMatrixAdditive=gMatrix$Ga
#         averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
#         averageRelativity=averageRelativity[training]
#         individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
#         predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
#         corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#         cumulativeCor=cumulativeCor+corelation
#       }
#       oneTrait=c(oneTrait,cumulativeCor/10)
#     }
#     predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
#     corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#     oneTrait=c(oneTrait,13158,corelation)
#   }
#   result=rbind(result,oneTrait)
# }
# write.csv(result,file = "2015 all traits with one reduced model.csv")
# #14
# testing=individualList[grep(pattern = "AYT14",individualList)]
# training=individualList[-grep(pattern = "AYT14",individualList)]
# result=c()
# for(trait in colnames(phenotype)[2:8]){
#   oneTrait=c(trait)
#   for(numTopIndividual in seq(from=50,to=length(training),by=50)){
#     oneTrait=c(oneTrait,numTopIndividual)
#     for(markerNum in c(500,1000,2000,3000,4000)){
#       cumulativeCor=0
#       for(rep in 1:10){
#         geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
#         geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
#         #problemeticMarker=which(colMeans(geneticInfoImputed)<0.0001)
#         gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
#         gMatrixAdditive=gMatrix$Ga
#         averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
#         averageRelativity=averageRelativity[training]
#         individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
#         predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
#         corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#         cumulativeCor=cumulativeCor+corelation
#       }
#       oneTrait=c(oneTrait,cumulativeCor/10)
#     }
#     predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
#     corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#     oneTrait=c(oneTrait,13158,corelation)
#   }
#   result=rbind(result,oneTrait)
# }
# write.csv(result,file = "2014 all traits with one reduced model.csv")
# #13
# testing=individualList[grep(pattern = "AYT13",individualList)]
# training=individualList[-grep(pattern = "AYT13",individualList)]
# result=c()
# for(trait in colnames(phenotype)[2:8]){
#   oneTrait=c(trait)
#   for(numTopIndividual in seq(from=50,to=length(training),by=50)){
#     oneTrait=c(oneTrait,numTopIndividual)
#     for(markerNum in c(500,1000,2000,3000,4000)){
#       cumulativeCor=0
#       for(rep in 1:10){
#         geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
#         geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
#         #problemeticMarker=which(colMeans(geneticInfoImputed)<0.0001)
#         gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
#         gMatrixAdditive=gMatrix$Ga
#         averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
#         averageRelativity=averageRelativity[training]
#         individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
#         predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
#         corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#         cumulativeCor=cumulativeCor+corelation
#       }
#       oneTrait=c(oneTrait,cumulativeCor/10)
#     }
#     predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
#     corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#     oneTrait=c(oneTrait,13158,corelation)
#   }
#   result=rbind(result,oneTrait)
# }
# write.csv(result,file = "2013 all traits with one reduced model.csv")
# #12
# testing=individualList[grep(pattern = "AYT12",individualList)]
# training=individualList[-grep(pattern = "AYT12",individualList)]
# result=c()
# for(trait in colnames(phenotype)[2:8]){
#   oneTrait=c(trait)
#   for(numTopIndividual in seq(from=50,to=length(training),by=50)){
#     oneTrait=c(oneTrait,numTopIndividual)
#     for(markerNum in c(500,1000,2000,3000,4000)){
#       cumulativeCor=0
#       for(rep in 1:10){
#         geneticInfoSelected=geneticInfo[,sort(sample(2:13158,markerNum))]
#         geneticInfoImputed=A.mat(geneticInfoSelected,return.imputed = T)$imputed
#         #problemeticMarker=which(colMeans(geneticInfoImputed)<0.0001)
#         gMatrix=G.matrix(M=geneticInfoImputed[,],method = "VanRaden",format = "wide",plot = F)
#         gMatrixAdditive=gMatrix$Ga
#         averageRelativity=apply(gMatrixAdditive[,testing],1,mean)
#         averageRelativity=averageRelativity[training]
#         individualsForTraining=names(sort(averageRelativity,decreasing = T)[1:numTopIndividual])
#         predictedBreedingValue=predictBreedingValue(geneticInfoImputed[individualsForTraining,],phenotype[individualsForTraining,][,trait],geneticInfoImputed[testing,])
#         corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#         cumulativeCor=cumulativeCor+corelation
#       }
#       oneTrait=c(oneTrait,cumulativeCor/10)
#     }
#     predictedBreedingValue=predictBreedingValue(geneticInfoImputed[training,],phenotype[training,][,trait],geneticInfoImputed[testing,])
#     corelation=cor(predictedBreedingValue,phenotype[testing,][,trait],use = "pairwise.complete.obs")
#     oneTrait=c(oneTrait,13158,corelation)
#   }
#   result=rbind(result,oneTrait)
# }
# write.csv(result,file = "2012 all traits with one reduced model.csv")

#############G matrix to a vector
dim(gMatrixAdditive$Ga)
length(individualList)
result=data.frame(numCol=3,numRow=997578)
i=1
for(numRow in 1:1412){
  for(numCol in (numRow+1):1413){
    # result=rbind(result,c(individualList[numRow],individualList[numCol],gMatrixAdditive$Ga[numRow,numCol]))
    result[i,1]=individualList[numRow]
    result[i,2]=individualList[numCol]
    result[i,3]=gMatrixAdditive$Ga[numRow,numCol]
    i=i+1
  }
}
write.csv(result,"gMatrixToVector.csv",row.names = F,quote = F)

effectiveTwoWayCross=read.csv("Effective two way cross.csv",header = T,row.names = 1,stringsAsFactors = F)

fullHalfSib=c()
for(line in 1:997578){
  if(result$numCol[line] %in% rownames(effectiveTwoWayCross)){
    if(result$numRow[line] %in% rownames(effectiveTwoWayCross)){
      l1P1=effectiveTwoWayCross[result$numCol[line],3]
      l1P2=effectiveTwoWayCross[result$numCol[line],4]
      l2P1=effectiveTwoWayCross[result$numRow[line],3]
      l2P2=effectiveTwoWayCross[result$numRow[line],4]
      line1Line2=c("line1"=result$numCol[line],"line2"=result$numRow[line],"accuracy"=result$V3[line],
                   "line1P1"=l1P1,"line1P2"=l1P2,"line2P1"=l2P1,"line2P2"=l2P2,
                   "numOfTotalPar"=length(intersect(c(l1P1,l1P2),c(l2P1,l2P2))))
      # if(line1Line2[8]>0){
      #   fullHalfSib=rbind(fullHalfSib,line1Line2)
      # }
      fullHalfSib=rbind(fullHalfSib,line1Line2)
    }
  }
}
intersect(c("1","2"),c("2","3"))
length(intersect(c("1","2"),c("2","3")))

result$numCol[1] %in% rownames(effectiveTwoWayCross)
effectiveTwoWayCross[result$numCol[1],3]
line1Line2[8]>-1

write.csv(fullHalfSib,file = "all sib info.csv",row.names = F,quote = F)
