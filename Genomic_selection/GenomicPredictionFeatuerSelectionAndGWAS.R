#Currently work in personal office desktop
getwd()
pheno = read.table("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/GenomicPredictionWithNeuralNetwork/2012-2017 BLUPs for all traits.txt",
                   header = T,row.names = 1)
pheno[1:5,]

geno = read.table("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/GenomicPredictionWithNeuralNetwork/numberic_geno.txt",
                  header = T,na.strings = NA,row.names = 1)
geno[1:5,1:5]

individualList = intersect(row.names(pheno),row.names(geno))

geno = geno[individualList,]
pheno = pheno[individualList,]

GWASResult12_16 = read.table("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/GenomicPredictionWithNeuralNetwork/GWAS_MixedModel_stats.txt",
                        header = T,stringsAsFactors = F)
GWASResult12_16[1:5,1:5]
Traits = unique(GWASResult$Trait)
Traits

library("rrBLUP")

genoImputed = A.mat(geno,return.imputed = T)$imputed
genoImputed[1:5,1:5]

write.csv(genoImputed,"numeric_geno_rrBLUP_Imputed.csv")

# pca_result = prcomp(genoImputed,center = T)
# pca_result$x[1:5,1:5]

predictBreedingValue=function(TPgeno,TPpheno,VPgeno){
  predictedRidgeRegression = mixed.solve(TPpheno,Z=TPgeno)
  predictedValue = VPgeno%*%predictedRidgeRegression$u
  return(predictedValue)
}

###########
#By GWAS result
GWASResult12_16[1:5,1:5]

markerInUse = GWASResult12_16[(GWASResult12_16$Trait=="SEblup") & (GWASResult12_16$p < 0.1),]$Marker

genoInUse = genoImputed[,markerInUse]

testing=individualList[grep(pattern = "AYT17",individualList)]#define training and testing
training=individualList[-grep(pattern = "AYT17",individualList)]
predictedBreedingValue=predictBreedingValue(genoInUse[training,],pheno[training,][,"SEMEXT"],genoInUse[testing,])
corelation=cor(predictedBreedingValue,pheno[testing,]["SEMEXT"],use = "pairwise.complete.obs")#calculate accuracy with all markers


###############
##Implementing paper method
set.seed(123)

genoPheno = cbind.data.frame(genoImputed[training,],pheno[training,]$TWT)
dim(genoPheno)
colnames(genoPheno)[13162] = "TWT"
genoPheno[1:5,c(1:5,13155:13162)]


#install.packages("FSelector")
library("FSelector")

corMatrix = cor(genoImputed)

subsetSample = sample(1187,50)

cfsTest = cfs(TWT~.,genoPheno[subsetSample,10000:13162])

a = proc.time()
cfs(TWT~.,genoPheno[,11000:13162])
b = proc.time()
print(b-a)

#50 ind 1161 markers - 1559.03  
#all ind 1161 markers ~ 5533
#all indi 2161 markers ~ 32918
