# OVERALL PROJECT: Feasibility of applying pharmacoepidemiologic drug-drug interaction screening 
# methods to a population of nursing home residents: An application to clopidogrel
# 
#   DESCRIPTION: 
#   Screen potential precipitant drugs for DDIs with clopidogrel resulting in 
#   major bleed and fall related injury (FRI) in a sample of NH residents using 
#   the self controlled case series (SCCS) design. 
# 
# PROGRAM: semibayes_adjustment_NHDDI_screening_github
# 
#   DESCRIPTION: 
#   This code uses the set the estimates and variances from the 
#   conditional poisson models for major bleed and (separately)
#   FRI.
#
#   To account for the multiple estimation inherent in calculating many rate ratios and 
#   CIs, we used a semi-Bayes shrinkage method. This increases the validity of effect estimates 
#   and preserves the nominal type-1 error rate. Operationally, we prespecified a variance 
#   (0.25) to assume that 95% of true rate ratios would be within an unspecified 
#   7-fold range, then shrunk outlying effect estimates toward their geometric mean.
#
# Programmer: Adam DAmico
# 
# Date: 28FEB25
# 
# Version History:


#Clear all objects
rm(list=ls())

#set directory
setwd('P:/nhddi/a5d/Aim1_Screening/Antiplatelets/Output/7_SCCS_models')

# SemiBayes function ####

# Note, we did not write this function - it was shared with us 
# by colleagues at UPenn

# This function runs Semi-Bayes adjustment towards the global mean for groups of estimates.
# For each group of estimates within which Semi-Bayes adjustment needs to be performed,
# a .csv file must be created under the working directory containing 3 columns:
# The first one with variable names (e.g. precipitant),
# the second one with ML beta estimates and
# the third one with the variances of the ML estimates.
#
# There must be one file for each group within which SB adjustment is performed.

SemiBayes <- function(vart, infile, outfile){
	for (i in 1:1){
		Data<-read.csv(infile, header=TRUE)
		bhat<-as.numeric(c(Data[[2]]))
		vhat<-as.numeric(c(Data[[3]]))
		name<-as.character(Data[[1]])
		newbeta<-vector(length=length(bhat))
		newvar<-vector(length=length(bhat))
		SB_RR<-vector(length=length(bhat))
		SB_RR_LL<-vector(length=length(bhat))
		SB_RR_UL<-vector(length=length(bhat))
		n <- length(bhat)
		w <- 1/(vhat + vart)
		pi <- sum(w*bhat)/sum(w)
		D <- bhat - pi
		varo <- sum(w*D*D)/sum(w)
		varm <- sum(w*vhat)/sum(w)
		newbeta<- w*((vart*bhat)+(vhat*pi))
		E <- w* vhat*D*sqrt(varm/vhat)
		newvar <- vhat * (1-vhat*w) + (2*E*t(E))/n
		SB_RR<-c(exp(newbeta))
		SB_RR_LL<-c(exp(newbeta-1.96*sqrt(newvar)))
		SB_RR_UL<-c(exp(newbeta+1.96*sqrt(newvar)))
		newstd<-t(sqrt(newvar))

                stderr<-sqrt(vhat)
		RR<-exp(bhat)
		RR_lower<-exp((bhat-1.96*sqrt(vhat)))
		RR_upper<-exp((bhat+1.96*sqrt(vhat)))

		out<-cbind(name,bhat,vhat,RR,RR_lower,RR_upper,newbeta,t(newvar),SB_RR,SB_RR_LL,SB_RR_UL)
		b<-as.data.frame(out)
		names(b)[c(1)] <- c("Precipitant")
		names(b)[c(2)] <- c("beta")
		names(b)[c(3)] <- c("var")
		names(b)[c(4)] <- c("RR")
		names(b)[c(5)] <- c("RR_lower")
		names(b)[c(6)] <- c("RR_upper")
		names(b)[c(7)] <- c("SB_beta")
		names(b)[c(8)] <- c("SB_var")
		names(b)[c(9)] <- c("SB_RR")
		names(b)[c(10)] <- c("SB_RR_lower")
		names(b)[c(11)] <- c("SB_RR_upper")
	
		write.table(b, file=outfile, sep=",", col.names=TRUE, row.names=FALSE, quote=F, na="NA")		

	}
}


#Runs the function SemiBayes for major bleeds and fall-related injuries		

SemiBayes(0.25, "sumstats_SAS_simplest_models_mb_02OCT24.csv", "sumstats_SAS_simplest_models_mb_28FEB25_adj_0p25.csv")
SemiBayes(0.25, "sumstats_SAS_simplest_models_fri_02OCT24.csv", "sumstats_SAS_simplest_models_fri_28FEB25_adj_0p25.csv")

