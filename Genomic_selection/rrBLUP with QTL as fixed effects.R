getwd()
setwd("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/Spring Wheat FHB resistance manuscript/Spring_wheat_marker_size_validation")
setwd("C:/Users/Administrator/Downloads")

library("rrBLUP")
phenotype = read.table("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/Spring Wheat FHB resistance manuscript/4th_June_2019/phenotype_estimated blups of HSEVR across all trials_439.txt",
             row.names = 1, header = T)

genotype=read.csv("genotype_AYT11-16_439GBSna50ABD_Bi_LDKNNimp.F01_meanImp.csv",header = T)

genotypeUnimputed = read.table("FHB-nursery spring wheat lines_439 GBS SNPna50ABD_Bi_LDKNNimp_MAF5.F01.txt",header=T,
                               row.names = 1)

yphenotype[1:5,1:2]
genotype[1:5,1:5]

dim(phenotype)
dim(genotype)

individualList=rownames(phenotype)
row.names(genotype)=individualList
genotype=as.matrix(genotype)

##############################################################
#regular GS
predictBreedingValue=function(TPgeno,TPpheno,VPgeno){
  predictedRidgeRegression=mixed.solve(TPpheno,Z=TPgeno)
  predictedValue=VPgeno%*%predictedRidgeRegression$u
  return(predictedValue)
}



testPopulation="16" #define which year as testing population
testing=individualList[grep(pattern = "^16-",individualList)]
training=individualList[-grep(pattern = "^16-",individualList)]
phenoColumn = 4 #4:blue4  8:blup4 

predicted=predictBreedingValue(genotype[training,],phenotype[training,phenoColumn],genotype[testing,])#test with all markers
correlation=cor(predicted,phenotype[testing,phenoColumn])
correlation


#################################################################
#define fixed effects
markersAsFixedFactor=c(which(colnames(genotype)=="S1A_471476097"),which(colnames(genotype)=="S6B_474184564"))

#GS with markers as fixed effect
predictBreedingValue=function(TPgeno,TPpheno,VPgeno){
  predictedRidgeRegression=mixed.solve(TPpheno,Z=TPgeno[,-markersAsFixedFactor],X=TPgeno[,markersAsFixedFactor])
  predictedValue=VPgeno[,-markersAsFixedFactor]%*%predictedRidgeRegression$u
  return(predictedValue)
}

testPopulation="16"
testing=individualList[grep(pattern = "^16-",individualList)]
training=individualList[-grep(pattern = "^16-",individualList)]
phenoColumn = 4 #4:blue4  8:blup4 

predicted=predictBreedingValue(genotype[training,],phenotype[training,phenoColumn],genotype[testing,])#test with all markers
correlation=cor(predicted,phenotype[testing,phenoColumn])
correlation

##################################################################
#Estimate identity by state similarity
#Use Demerelate
library("Demerelate")
# data("demerelpop")
# data("demerelref")
# data("demereldist")
# 
# demerelpop.sp <- split(demerelpop,demerelpop[,2])
# 
# demerelpop.sp[[1]]
# 
# testRxy = Emp.calc(demerelpop.sp[[1]],value = "Mxy", ref.pop = "NA")
# 
# genotypeUnimputed[1:5,1:5]

#Population = c(rep(1,200),rep(2,239))

DemerelateFormatGenotype = cbind.data.frame(genotypeUnimputed[,1],1,genotypeUnimputed[,2:102148])

BxyEstimator = Emp.calc(DemerelateFormatGenotype,value = "Mxy",ref.pop = 2)

mainFunctionTest = Demerelate(DemerelateFormatGenotype[,1:102148],value = "Bxy",
                              file.output = FALSE,object = TRUE,ref.pop = "NA",NA.rm = FALSE)
str(mainFunctionTest)

#####Split genotype into two values to accomodate Demerelate package
#0 as first genotype; 1 as second genotype
#E.g. 0 -> 1  0; 1 -> 0   1.
#Two genotype calls table will be created and then merged into single one.

originalGenotype = genotypeUnimputed[,2:ncol(genotypeUnimputed)]
firstGenotypeTable = apply(originalGenotype, c(1,2), 
                           function(x) if(is.na(x)) x = NA else if(x==0){x = 1}else if(x==1){x = 0})

firstGenotypeTableName = colnames(firstGenotypeTable)

firstGenotypeTableName = sapply(firstGenotypeTableName, FUN = function(x) paste(x,"0",sep = "_"))

colnames(firstGenotypeTable) = firstGenotypeTableName

secondGenotypeTable = apply(originalGenotype, c(1,2), 
                            function(x) if(is.na(x)) x = NA else if(x==0){x = 0}else if(x==1){x = 1})

secondGenotypeTableName = colnames(secondGenotypeTable)

secondGenotypeTableName = sapply(secondGenotypeTableName, FUN = function(x) paste(x,"1",sep = "_"))

colnames(secondGenotypeTable) = secondGenotypeTableName

finalGenotypeTable = cbind.data.frame(firstGenotypeTable,secondGenotypeTable)

DemerelateFormatGenotype = cbind.data.frame(genotypeUnimputed[,1],1,finalGenotypeTable)

mainFunctionTest = Demerelate(DemerelateFormatGenotype[,1:10],value = "Mxy",
                              file.output = FALSE,object = TRUE,ref.pop = "NA",NA.rm = FALSE,iteration = 500,pairs = 500)
str(mainFunctionTest)

mainFunctionTest

########
#USe snpRelate from Bioconductor

genotypeUnimputed[1:5,1:5]

genotypeTest = apply(genotypeUnimputed,c(1,2),FUN = function(x) if(is.na(x)){x=9}else if(x==1){x=2}else{x=0})

genotypeTest[1:5,1:5]

snpgdsCreateGeno("sampler.gds",genmat = as.matrix(genotypeTest),sample.id = row.names(genotypeTest),
                 snp.id = colnames(genotypeTest),snpfirstdim = F)

genofile = snpgdsOpen("sampler.gds")

ibs = snpgdsIBS(genofile,num.thread = 1)

snpgdsClose(genofile)

identityByState = ibs$ibs

row.names(identityByState) = row.names(genotypeTest)

colnames(identityByState) = row.names(genotypeTest)

identityByState[1:5,1:5]

getwd()

write.csv(identityByState,file = "identify by state -spring wheat.csv",quote = F)
