setwd("C:/Users/liu.yuan/OneDrive - North Dakota University System/Research/GenomicPredictionWithNeuralNetwork")
list.files()
phenotype=read.table("phenotype_AYT12-16_ND12-quality_BLUPs.txt",header = T,stringsAsFactors = F,row.names = 1)
phenotype[1:5,1:5]
phenotype=phenotype[sort(row.names(phenotype)),]

PCInfo=read.table("PC_info_1413_marker.txt",header = T,row.names = 1)
PCInfo[1:5,1:5]

PCEigenvalue=read.table("PC_info_Eigenvalue_1413_marker.txt",header = T,row.names = 1)
PCEigenvalue[1:20,]
plot(PCEigenvalue$proportion_of_total[1:20])

PCVector=read.table("PC_info_EigenVector_1413_marker.txt",header = T,row.names = 1)
PCVector[1:5,1:5]

geneticInfo=read.table("numberic_geno.txt",header = T,row.names = 1,na.strings = "NA")
geneticInfo[1:5,1:5]

intersect(row.names(phenotype),row.names(PCInfo)) #get individual intersection
individualList=intersect(row.names(phenotype),row.names(PCInfo))

phenotype=phenotype[individualList,]
PCInfo=PCInfo[individualList,]
geneticInfo=geneticInfo[individualList,]


#########Principle component regression
training=individualList[1:1014]
testing=individualList[1015:1184]

dataPCRtraining=cbind(as.numeric(phenotype[training,]$SedVblup),PCInfo[training,1:20])
colnames(dataPCRtraining)[1]="SedVblup"
PCLinearRegression=lm(SedVblup~.,data = dataPCRtraining)
summary(PCLinearRegression)
anova(PCLinearRegression)
pchisq(PCLinearRegression$df.residual,df.residual(PCLinearRegression),lower.tail = F)
plot(PCLinearRegression$fitted.values,PCLinearRegression$residuals,main = "Residual plot in PCR")
PCLinearRegression$coefficients
predict(PCLinearRegression,PCInfo[testing,1:20])
cor(predict(PCLinearRegression,PCInfo[testing,1:20]),as.numeric(phenotype[testing,3]),use = "pairwise.complete.obs")


##############Nerual network
#install.packages("neuralnet")
library("neuralnet")

dataPCRtraining=cbind(as.numeric(phenotype[training,]$SedVblup),PCInfo[training,1:200])

colnames(dataPCRtraining)[1]="SedVblup"
n=names(dataPCRtraining)
f=as.formula(paste("SedVblup~",paste(n[!n %in% "SedVblup"],collapse="+")))
f
nerualNetworkTraining=neuralnet(f,dataPCRtraining,hidden = c(1000,1000,1000),linear.output = T,stepmax = 1e+06)
nerualNetworkTraining
predictNeurualNetwork=compute(nerualNetworkTraining,PCInfo[testing,1:200])
predictNeurualNetwork$net.result
cor(predictNeurualNetwork$net.result,as.numeric(phenotype[testing,3]),use = "pairwise.complete.obs")


#Marker instead of PC
hundredRandomMarker=sample(13161,200)
#geneticInfo[is.na(geneticInfo)]=0.5
dataMultinomialRegression=geneticInfo[training,hundredRandomMarker]




############Multinomial model
library("nnet")

plot(dataPCRtraining$SedVblup)

hundredRandomMarker=sample(13161,200)
abline(h=10,col=2)
abline(h=-10,col=2)

geneticInfo[is.na(geneticInfo)]=0.5

dataMultinomialRegression=geneticInfo[training,hundredRandomMarker]

dataMultinomialRegression= cbind(cut(dataPCRtraining$SedVblup,c(-30,-10,10,30)),dataMultinomialRegression)

dataMultinomialRegression[1:5,1:5]
colnames(dataMultinomialRegression)[1]="SedVblup"
dataMultinomialRegression[1:5,1:5]

multinomialModel=multinom(SedVblup~.,data = dataMultinomialRegression,Hess = T)

table(dataMultinomialRegression$SedVblup)
#chisq.test(table(dataMultinomialRegression$SedVblup),table(dataMultinomialRegression$SedVblup)/1014)

summary(multinomialModel)

predict(multinomialModel,geneticInfo[testing,hundredRandomMarker],"probs")
predict(multinomialModel,geneticInfo[testing,hundredRandomMarker])

cor(predict(multinomialModel,geneticInfo[testing,hundredRandomMarker]),
    cut(as.numeric(phenotype[testing,3]),c(-30,-10,10,30)))

predictionComparison=cbind(predict(multinomialModel,geneticInfo[testing,hundredRandomMarker]),cut(as.numeric(phenotype[testing,3]),c(-30,-10,10,30)))

nrow(predictionComparison[predictionComparison[,1]==predictionComparison[,2],])/length(testing)
#############Ridge regression
library("rrBLUP")
geneticInfoImputed=A.mat(geneticInfo)
gMatrix=as.matrix(geneticInfoImputed)%*%t(as.matrix(geneticInfoImputed))
gMatrix[1:5,1:5]
max(gMatrix)
min(gMatrix)
scalingFactor=max(abs(min(gMatrix)),max(gMatrix)) #Naive way to normalize G matrix
normalizedGMatrix=gMatrix/scalingFactor

heatmap(gMatrix,keep.dendro = F)
normalizedGMatrix[testing[1],]
plot(normalizedGMatrix[testing[5],])



geneticInfoTraining=geneticInfoImputed[training,]

dataRidgeRegression=geneticInfoTraining[complete.cases(dataPCRtraining$SedVblup),]

trainingTrait=dataPCRtraining$SedVblup[complete.cases(dataPCRtraining$SedVblup)]

predictedRidgeRegression=mixed.solve(trainingTrait,Z=dataRidgeRegression)

geneticInfoImputed[testing,]%*%predictedRidgeRegression$u

cor(geneticInfoImputed[testing,]%*%predictedRidgeRegression$u,as.numeric(phenotype[testing,]$SedVblup),use = "pairwise.complete.obs")

#####rrBLUP with relatives
library("rrBLUP")

geneticInfoImputed=A.mat(geneticInfo)
gMatrix=as.matrix(geneticInfoImputed)%*%t(as.matrix(geneticInfoImputed))
gMatrix[1:5,1:5]
max(gMatrix)
min(gMatrix)
scalingFactor=max(abs(min(gMatrix)),max(gMatrix)) #Naive way to normalize G matrix
normalizedGMatrix=gMatrix/scalingFactor

heatmap(gMatrix,keep.dendro = F)
normalizedGMatrix[testing[1],]
plot(normalizedGMatrix[testing[5],])

function(TPgeno,TPpheno,VPgeno){
  
  
}

geneticInfoTraining=geneticInfoImputed[training,]

dataRidgeRegression=geneticInfoTraining[complete.cases(dataPCRtraining$SedVblup),]

trainingTrait=dataPCRtraining$SedVblup[complete.cases(dataPCRtraining$SedVblup)]

predictedRidgeRegression=mixed.solve(trainingTrait,Z=dataRidgeRegression)

geneticInfoImputed[testing,]%*%predictedRidgeRegression$u

cor(geneticInfoImputed[testing,]%*%predictedRidgeRegression$u,as.numeric(phenotype[testing,]$SedVblup),use = "pairwise.complete.obs")

