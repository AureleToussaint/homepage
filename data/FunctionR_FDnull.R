###########################################################################################################################
# R function to compare the functional volume of a community based on functional traits to a random a random pool assemblage                                                #
#																																																				                             #
#	Functional Diversity (FD) was tested against a null model																																					                             #
#	Null Model was conservative to the choosing taxonomic level                                                                                                                                                                        #
#																																																																                             #
#						Inputs : 	- traits : a matrix (S x T) with S the names of the species (lines) and T the functional traits (in colones)                                                          #
#										- comm : a matrix (C x S) with C the community	and S the presence/absence of the species in C communities                                            #
#										- taxo : a matrix of S lines with the taxonomy of each species 																										                             #
#										- boot : the number of iterations																																								                             #
# 																																																																                             #
#						Output: a list object with :	- Observed : a matrix of 2 lines and C columns with in the first line the number of species (TD)                                   #
#																									and the second line the functional volume (FD),  																	                             #
#																									both  expessed as a percentage of the total TD and FD, respectively									                             #
#																																																																                             #
#																			- Randomized : a matrix of 2 lines and C columns with the first line the MEAN and SD in the second line	             #
#																																																																                             #
#																			- SESandPval : a matrix of 2 lines and C columns with the first line the values of Standardized Effect Size               #
#																										 and the second line the values of the significant SES (alpha = 5%)               				                             #
#																																																																                             #
###########################################################################################################################


FD.null <- function (traits,comm,taxo,boot = 999, save = T){

# loading required libraries 
library(ape) 
library(geometry) 

# C = number of communities
C <- nrow(comm)

# T = number of traits
T<-dim(traits)[2] 

# 'traits' matrix needs to be scaled and centered
traitsCS <- scale(traits, center=TRUE, scale=TRUE) 

# definition of vector for results
rch <- rep(NA,C)
observed <- matrix(NA,nr=2,nc=C,dimnames=list(c("TD","FD"),rownames(comm)))			
randomize <- matrix(NA,nr=C,nc=boot,dimnames=list(rownames(comm),c(paste("boot",1:(boot),sep=""))))  

############################################################################################################ 
# loop to compute on each community the Functional volume and the number of species (Observed and randomize)						

for(i in 1:C){

# selection of species present in the community		
sppres<-colnames(comm)[which(comm[i,] ==1)]

# check if the number of species is higher than the number of trait in the C community
if(length(sppres) > (T+1)){
	
# traits of the species present only in the community 
tr<-traitsCS[sppres,]

# taxonomic and functional observed indicies

	# TD
	observed["TD",i]<-dim(traits[which(rownames(traits) %in% sppres),])[1]
	
	# FD
	 tryCatch(observed["FD",i]<-round(convhulln(tr,"FA")$vol,6), error=function(e){}) 


# randomization		
for (rand in 1:(boot)){
	
# conservative null model : the number of sppre
spec.boot<-NULL

# random pool 
for(taxonomy in levels(factor(taxo)) ) {
	lgt<-length(which(sppres %in%names(which(taxo == taxonomy)) ))
	spec.boot<-c(spec.boot,sample(names(taxo[taxo == taxonomy]),lgt))
} # end of taxonomy

# traits of the species sorting by the null model
trboot<-traitsCS[spec.boot,]

# functional randomize indicies
	
	# FD random
	 tryCatch(randomize[i,paste("boot",rand,sep="")]<-round(convhulln(trboot,"FA")$vol,6), error=function(e){}) 

} # end of boot
} # end of if length > T+1
print(paste("Community :",rownames(comm)[i],"- Done"))
} # end of C

# Compute the SES and pval values

		# SES
		SES <- apply(cbind(observed[2,], randomize),1,function(y){(y[1]-mean(y[-c(1)],na.rm=T))/sd(y[-c(1)],na.rm=T)})
		
		# pval
		pval <- apply(cbind(observed[2,], randomize),1,function(y){length(y[y[1]<y[-c(1)]][-1])/length(y[-c(1)])})
		
	random <- rbind(apply(randomize,1,mean),apply(randomize,1,sd)) ; colnames(random)<-rownames(comm) ; rownames(random)<-c("FD Mean","FD SD")
if (save == T){
result <- list(observed = observed,randomizeMeanSD = random, randomize = randomize ,SESandPval = rbind(SES,pval))	
}else{
result <- list(observed = observed,randomizeMeanSD = random, SESandPval = rbind(SES,pval))	
}
return(result)
}# end of FD.null

############################################################################################################ 
# Example
############################################################################################################ 

showexample <- T

if ( showexample == T) {
	
	taxo <- rbind(c("GenA","AA"),c("GenA","AB"),c("GenA","AC"),c("GenA","AD"),c("GenB","BA"),c("GenB","BB"),c("GenB","BC"),c("GenB","BD"),c("GenC","CA"),c("GenC","CB"),c("GenC","CC"),c("GenD","DD"),c("GenE","EA"),c("GenE","EB"),c("GenE","EC"),c("GenE","ED"),c("GenF","FA"),c("GenF","FB"),c("GenF","FC"),c("GenF","FD")) ; colnames(taxo) <- c("Genus","species"); rownames(taxo) <- taxo[,2]
	traits <- cbind(sample(rnorm(100,0.7,0.10),nrow(taxo)),sample(rnorm(100,0.1,0.25),nrow(taxo))) ; rownames(traits)<- rownames(taxo) 
	comm<-rbind (round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1))) ; colnames(comm)<-rownames(taxo) ; rownames(comm) <- paste("Comm",1:3)
	
	print(FD.null(traits,comm,taxo[,1],boot = 999,save=T))
}

############################################################################################################ 
# END
############################################################################################################ 