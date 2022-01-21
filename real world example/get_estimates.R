#####
## install the package pomp and source these files (may need to change working directory)

library(pomp)	# version 2.7.1.0 used
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

# (**)
# CANADA
data0	<-	read.table("cases_timeseries_canada.csv", sep=",", header=TRUE)
NT_vec	<-	as.vector(unlist(data0[,4]))
times	<-	7*(0:30)
NT_vec	<-	NT_vec[times]
NT_vec	<-	c(NT_vec[1],diff(NT_vec))
data 	<- 	NT_vec[5:10]  
popltn 	<-	 38000000
data
# (**)

######
## obtain various estimators of R0
## code assumes data is as above, and return Rhat estimate for last observed week

round(WP(NT=data, method="known", mu=2/7)$Rhat,2)
round(WP(NT=data, method="known", mu=5/7)$Rhat,2)
round(WP(NT=data, method="known", mu=8/7)$Rhat,2)
round(WP(NT=data, method="unknown")$Rhat,2)

round(seqB(NT=data, mu=2/7)$Rhat,2)
round(seqB(NT=data, mu=5/7)$Rhat,2)
round(seqB(NT=data, mu=8/7)$Rhat,2)

round(ID(NT=data, mu=2/7)$Rhat,2)
round(ID(NT=data, mu=5/7)$Rhat,2)
round(ID(NT=data, mu=8/7)$Rhat,2)

round(IDEA(NT=data, mu=2/7)$Rhat,2)
round(IDEA(NT=data, mu=5/7)$Rhat,2)
round(IDEA(NT=data, mu=8/7)$Rhat,2)

res_fullBayes_SIR	<-	fullbayes_SIR(data=data, sbeta=1, rbeta=2, sgamma=5, rgamma=26, theta=1, N=popltn)
res_fullBayes_SEIR	<-	fullbayes_SEIR(data=data, sbeta=13, rbeta=11, sgamma=5, rgamma=11,ssigma=1, rsigma=3, theta=1, N=popltn)
res_fullBayes_SEAIR	<-	fullbayes_SEAIR(data=data, sbeta=26, rbeta=57, sgamma=5, rgamma=11,ssigma=1, rsigma=3,srho=2, rrho=7, theta=1, N=popltn)

round(res_fullBayes_SIR$R0[length(data)],2)
round(res_fullBayes_SEIR$R0[length(data)],2)
round(res_fullBayes_SEAIR$R0[length(data)],2)

res_plugnplay_SIR	<-	plugnplay_SIR(data=data, sbeta=1, rbeta=2, sgamma=5, rgamma=26, N=popltn)
res_plugnplay_SEIR	<-	plugnplay_SEIR(data=data, sbeta=13, rbeta=11, ssigma=1, rsigma=3,sgamma=5, rgamma=11, N=popltn)
res_plugnplay_SEAIR	<-	plugnplay_SEAIR(data=data, sbeta=26, rbeta=57, ssigma=1, rsigma=3,sgamma=5, rgamma=11, srhoA=2, rrhoA=7, N=popltn)

round(res_plugnplay_SIR$R0[length(data)],2)
round(res_plugnplay_SEIR$R0[length(data)],2)
round(res_plugnplay_SEAIR$R0[length(data)],2)


######
## obtain various estimators of serial interval (in weeks)

SD	<-	WP(NT=data, method="unknown")$SD
round(sum(SD$supp*SD$pmf)-0.5,2)

res_fullBayes_SIR	<-	fullbayes_SIR(data=data, sbeta=1, rbeta=2, sgamma=5, rgamma=26, theta=1, N=popltn)
res_fullBayes_SEIR	<-	fullbayes_SEIR(data=data, sbeta=13, rbeta=11, sgamma=5, rgamma=11,ssigma=1, rsigma=3, theta=1, N=popltn)
res_fullBayes_SEAIR	<-	fullbayes_SEAIR(data=data, sbeta=26, rbeta=57, sgamma=5, rgamma=11,ssigma=1, rsigma=3,srho=2, rrho=7, theta=1, N=popltn)

round(res_fullBayes_SIR$SI[length(data)]/7,2)
round(res_fullBayes_SEIR$SI[length(data)]/7,2)
round(res_fullBayes_SEAIR$SI[length(data)]/7,2)

res_plugnplay_SIR	<-	plugnplay_SIR(data=data, sbeta=1, rbeta=2, sgamma=5, rgamma=26, N=popltn)
res_plugnplay_SEIR	<-	plugnplay_SEIR(data=data, sbeta=13, rbeta=11, ssigma=1, rsigma=3,sgamma=5, rgamma=11, N=popltn)
res_plugnplay_SEAIR	<-	plugnplay_SEAIR(data=data, sbeta=26, rbeta=57, ssigma=1, rsigma=3,sgamma=5, rgamma=11, srhoA=2, rrhoA=7, N=popltn)

round(res_plugnplay_SIR$SI[length(data)]/7,2)
round(res_plugnplay_SEIR$SI[length(data)]/7,2)
round(res_plugnplay_SEAIR$SI[length(data)]/7,2)


#####
## how to read in and clean data for provinces (may need to change working directory)
## change code below for code between (**) and (**) above

#BC
data0		<-	read.table("cases_timeseries_prov.csv", sep=",", header=TRUE)
i_BC		<-	which(data0$province=="BC")
NT_vec_BC	<-	as.vector(unlist(data0[i_BC[1:80],4]))
times		<-	7*(1:10)
NT_vec_BC	<-	NT_vec_BC[times]
NT_vec_BC	<-	c(NT_vec_BC[1],diff(NT_vec_BC))
data		<-	NT_vec_BC[5:10]
popltn 		<-	5200000
data

#ON
data0		<-	read.table("cases_timeseries_prov.csv", sep=",", header=TRUE)
i_ON		<-	which(data0$province=="Ontario")
NT_vec_ON	<-	as.vector(unlist(data0[i_ON[1:80],4]))
times		<-	7*(1:10)
NT_vec_ON	<-	NT_vec_ON[times]
NT_vec_ON	<-	c(NT_vec_ON[1],diff(NT_vec_ON))
data		<-	NT_vec_ON[5:10]
popltn 		<-	14800000
data

#QC
data0		<-	read.table("cases_timeseries_prov.csv", sep=",", header=TRUE)
i_QC		<-	which(data0$province=="QC")
NT_vec_QC	<-	as.vector(unlist(data0[i_QC[1:80],4]))
times		<-	7*(1:10)
NT_vec_QC	<-	NT_vec_QC[times]
NT_vec_QC	<-	c(NT_vec_QC[1],diff(NT_vec_QC))
data		<-	NT_vec_QC[5:10]
popltn 		<-	8600000
data

