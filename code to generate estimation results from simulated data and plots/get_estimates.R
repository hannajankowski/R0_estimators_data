#####
## install the package pomp and source these files (may need to change working directory)

library(pomp) # version 2.7.1.0 used
source("computeLL.R")
source("ID.R")
source("IDEA.R")
source("plugnplay.R")
source("seqB.R")
source("WP_known.R")
source("WP_unknown.R")
source("WP.R")
source("fullBayes.R")
source("plugnplay.R")

#####
## read in and clean data (may need to change working directory)

data		<-	read.csv("Flu1_Weekly_SEIR.csv")
#data		<-	read.csv("Flu1_Weekly_SEIR.csv")
#data		<-	read.csv("Flu1_Weekly_SEAIR.csv")
#data		<-	read.csv("COVID_Weekly_SIR.csv")
#data		<-	read.csv("COVID_Weekly_SEIR.csv")
#data		<-	read.csv("COVID_Weekly_SEAIR.csv")
#data		<-	read.csv("Flu2_Weekly_SIR.csv")
#data		<-	read.csv("Flu2_Weekly_SEIR.csv")
#data		<-	read.csv("Flu2_Weekly_SEAIR.csv")
head(data)


res_Rhat		<-	matrix(0,1000,15)

for(i in 1:1000){
for(j in 1:15){

	NT					<-	as.vector(c(1,unlist(data[i,1:j])))
	res_Rhat[i,j]		<-	WP(NT=NT, method="known", mu=2/7)$Rhat	#(***)
	print(j)		
}
print(i)		
}

## write out result (may need to cange working directory)
write.table(res_Rhat, "name of your choice.csv", sep=",")


# replace (***) with 


#res_Rhat[i,j]		<-	round(WP(NT=NT, method="known", mu=2/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(WP(NT=NT, method="known", mu=5/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(WP(NT=NT, method="known", mu=8/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(WP(NT=NT, method="unknown")$Rhat,2)

#res_Rhat[i,j]		<-	round(seqB(NT=NT, mu=2/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(seqB(NT=NT, mu=5/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(seqB(NT=NT, mu=8/7)$Rhat,2)

#res_Rhat[i,j]		<-	round(ID(NT=NT, mu=2/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(ID(NT=NT, mu=5/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(ID(NT=NT, mu=8/7)$Rhat,2)

#res_Rhat[i,j]		<-	round(IDEA(NT=NT, mu=2/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(IDEA(NT=NT, mu=5/7)$Rhat,2)
#res_Rhat[i,j]		<-	round(IDEA(NT=NT, mu=8/7)$Rhat,2)

#res_Rhat[i,j]		<-	plugnplay_SIR(data=NT, sbeta=1, rbeta=2, sgamma=5, rgamma=26, N=10000)$R0[length(NT)]
#res_Rhat[i,j]		<-	plugnplay_SEIR(data=NT, sbeta=13, rbeta=11, ssigma=1, rsigma=3,sgamma=5, rgamma=11, N=10000)$R0[length(NT)]
#res_Rhat[i,j]		<-	plugnplay_SEAIR(data=NT, sbeta=26, rbeta=57, ssigma=1, rsigma=3,sgamma=5, rgamma=11, srhoA=2, rrhoA=7, N=10000)$R0[length(NT)]


#res_Rhat[i,j]		<-	fullbayes_SIR(data=NT, sbeta=1, rbeta=2, sgamma=5, rgamma=26, theta=1, N=10000)$R0[length(NT)]
#res_Rhat[i,j]		<-	fullbayes_SEIR(data=NT, sbeta=13, rbeta=11, sgamma=5, rgamma=11,ssigma=1, rsigma=3, theta=1, N=10000)$R0[length(NT)]
#res_Rhat[i,j]		<-	fullbayes_SEAIR(data=NT, sbeta=26, rbeta=57, sgamma=5, rgamma=11,ssigma=1, rsigma=3,srho=2, rrho=7, theta=1, N=10000)$R0[length(NT)]



######
## obtain estimators of serial interval (in weeks)

data		<-	read.csv("Flu1_Weekly_SEIR.csv")
#data		<-	read.csv("Flu1_Weekly_SEIR.csv")
#data		<-	read.csv("Flu1_Weekly_SEAIR.csv")
#data		<-	read.csv("COVID_Weekly_SIR.csv")
#data		<-	read.csv("COVID_Weekly_SEIR.csv")
#data		<-	read.csv("COVID_Weekly_SEAIR.csv")
#data		<-	read.csv("Flu2_Weekly_SIR.csv")
#data		<-	read.csv("Flu2_Weekly_SEIR.csv")
#data		<-	read.csv("Flu2_Weekly_SEAIR.csv")
head(data)


res_SI		<-	matrix(0,1000,15)

for(i in 1:1000){
for(j in 1:15){

	NT				<-	as.vector(c(1,unlist(data[i,1:j])))
	SD				<-	WP(NT=NT, method="unknown")$SD
	res_SI[i,j]		<-	sum(SD$supp*SD$pmf)-0.5	#(***)
	print(j)		
}
print(i)		
}

## write out result (may need to cange working directory)
write.table(res_SI, "name of your choice.csv", sep=",")


# replace (***) with 

#res_SI[i,j]		<-	plugnplay_SIR(data=NT, sbeta=1, rbeta=2, sgamma=5, rgamma=26, N=10000)$SI[length(NT)]/7
#res_SI[i,j]		<-	plugnplay_SEIR(data=NT, sbeta=13, rbeta=11, ssigma=1, rsigma=3,sgamma=5, rgamma=11, N=10000)$SI[length(NT)]/7
#res_SI[i,j]		<-	plugnplay_SEAIR(data=NT, sbeta=26, rbeta=57, ssigma=1, rsigma=3,sgamma=5, rgamma=11, srhoA=2, rrhoA=7, N=10000)$SI[length(NT)]/7

#res_SI[i,j]		<-	fullbayes_SIR(data=NT, sbeta=1, rbeta=2, sgamma=5, rgamma=26, theta=1, N=10000)$SI[length(NT)]/7
#res_SI[i,j]		<-	fullbayes_SEIR(data=NT, sbeta=13, rbeta=11, sgamma=5, rgamma=11,ssigma=1, rsigma=3, theta=1, N=10000)$SI[length(NT)]/7
#res_SI[i,j]		<-	fullbayes_SEAIR(data=NT, sbeta=26, rbeta=57, sgamma=5, rgamma=11,ssigma=1, rsigma=3,srho=2, rrho=7, theta=1, N=10000)$SI[length(NT)]/7


