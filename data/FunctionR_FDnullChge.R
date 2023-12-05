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


FD.null <- function (traits,comm_histo,comm_curr,realmsp,IUCN,taxo,boot = 999, save = T){
  
  # loading required libraries 
  library(ape) 
  library(geometry) 
  
  # C = number of communities
  C <- nrow(comm_histo)
  
  # T = number of traits
  T<-dim(traits)[2] 
  
  # 'traits' matrix needs to be scaled and centered
  traitsCS <- scale(traits, center=TRUE, scale=TRUE) 
  
  # definition of vector for results
  rch <- rep(NA,C)
  nb_species_chg<-matrix(NA,nr=6,nc=C,dimnames=list(c("Histo","Curr","Intro","Translo","Exo","Extirp"),rownames(comm_histo)))
  observed <- matrix(NA,nr=2,nc=C,dimnames=list(c("TD","FD"),rownames(comm_histo)))			
  randomize <- matrix(NA,nr=C,nc=boot,dimnames=list(rownames(comm_histo),c(paste("boot",1:(boot),sep=""))))  
  
  ############################################################################################################ 
  # loop to compute on each community the Functional volume and the number of species (Observed and randomize)						
  
  for(i in 1:C){
    
    # selection of species present in the community		
    sppres_histo<-colnames(comm_histo)[which(comm_histo[i,] ==1)]
    sppres_curr<-colnames(comm_curr)[which(comm_curr[i,] ==1)]
    IUCN_nat<-IUCN[which(names(IUCN)%in%sppres_histo)]
    
    # check if the number of species is higher than the number of trait in the C community
    if(length(sppres_histo) > (T+1)){
      
      # traits of the species present only in the community 
      tr_histo<-traitsCS[sppres_histo,]
      tr_curr<-traitsCS[sppres_curr,]
      
      # taxonomic and functional observed indicies
      td_histo=dim(traits[which(rownames(traits) %in% sppres_histo),])[1]
      td_curr=dim(traits[which(rownames(traits) %in% sppres_curr),])[1]
      # TD
      observed["TD",i]<-td_curr/td_histo
      
      # introduction species
      
      nb_intro<-sppres_curr[!(sppres_curr %in% sppres_histo)]
      id_nat<-rep(NA,length(nb_intro)) ; names(id_nat)<-nb_intro
      id_nat[which(names(id_nat) %in% realmsp)]<-"translo"
      id_nat[!(names(id_nat) %in% realmsp)]<-"exo"
      nb_translo<-nb_intro[which(nb_intro %in% names(id_nat[which(id_nat == "translo")]))]
      nb_exo<-nb_intro[which(nb_intro %in% names(id_nat[which(id_nat == "exo")]))]
      
      #extirpation species
      nb_extirp<-sppres_histo[!(sppres_histo %in% sppres_curr)]
      
      nb_species_chg["Histo",i]<-td_histo
      nb_species_chg["Curr",i]<-td_curr
      nb_species_chg["Intro",i]<-length(nb_intro)
      nb_species_chg["Translo",i]<-length(nb_translo)
      nb_species_chg["Exo",i]<-length(nb_exo)
      nb_species_chg["Extirp",i]<-length(nb_extirp)
      
      # FD
      #tryCatch(fd_histo=round(convhulln(tr_histo,"FA")$vol,6), error=function(e){})
      fd_histo=round(convhulln(tr_histo,"FA")$vol,6)
      #tryCatch(fd_curr=round(convhulln(tr_curr,"FA")$vol,6), error=function(e){}) 
      fd_curr=round(convhulln(tr_curr,"FA")$vol,6)
      
      observed["FD",i]<-fd_curr/fd_histo
      # randomization		
      for (rand in 1:(boot)){
        # random pool 
        
        rd_translo<-unique(sample( names(id_nat[which(id_nat == "translo")]),length(nb_translo),replace=T))
        rd_exo<-unique(sample( names(id_nat[which(id_nat == "exo")]),length(nb_exo),replace=T))
        if(length(names(IUCN_nat[which(IUCN_nat %in% c("CR","EN","VU","NT"))])) !=0){
          rd_extirp<-unique(sample( names(IUCN_nat[which(IUCN_nat %in% c("CR","EN","VU","NT"))]),length(nb_extirp),replace=T)) 
        }else{
          rd_extirp<-unique(sample( names(IUCN_nat),length(nb_extirp),replace=T))  
        }
        
        
        # traits of the species sorting by the null model
        trboot<-tr_histo[!(rownames(tr_histo) %in% rd_extirp),]
        trboot<-rbind(trboot,traitsCS[c(rd_translo,rd_exo),])
        
        # functional randomize indicies
        td_curr_rd=dim(traits[which(rownames(traits) %in% rownames(trboot)),])[1]
        
        # FD random
        while(td_curr_rd != td_curr){
          rd_translo<-unique(sample( names(id_nat[which(id_nat == "translo")]),length(nb_translo),replace=T))
          rd_exo<-unique(sample( names(id_nat[which(id_nat == "exo")]),length(nb_exo),replace=T))
          rd_extirp<-unique(sample( names(IUCN_nat[which(IUCN_nat %in% c("CR","EN","VU","NT"))]),length(nb_extirp),replace=T))
          
          
          # traits of the species sorting by the null model
          trboot<-tr_histo[!(rownames(tr_histo) %in% rd_extirp),]
          trboot<-rbind(trboot,traitsCS[c(rd_translo,rd_exo),])
          
          # functional randomize indicies
          td_curr_rd=dim(traits[which(rownames(traits) %in% rownames(trboot)),])[1]
        }
        
        #tryCatch(fd_curr_rd=round(convhulln(trboot,"FA")$vol,6), error=function(e){}) 
        fd_curr_rd=round(convhulln(trboot,"FA")$vol,6)
        randomize[i,paste("boot",rand,sep="")]<-fd_curr_rd/fd_histo          
        
      } # end of boot
    } # end of if length > T+1
    print(paste("Community :",rownames(comm_histo)[i],"- Done"))
  } # end of C
  
  # Compute the SES and pval values
  
  # SES
  SES <- apply(cbind(observed[2,], randomize),1,function(y){(y[1]-mean(y[-c(1)],na.rm=T))/sd(y[-c(1)],na.rm=T)})
  
  # pval
  pval <- apply(cbind(observed[2,], randomize),1,function(y){length(y[y[1]<y[-c(1)]][-1])/length(y[-c(1)])})
  
  random <- rbind(apply(randomize,1,mean),apply(randomize,1,sd)) ; colnames(random)<-rownames(comm_histo) ; rownames(random)<-c("FD Mean","FD SD")
  if (save == T){
    result <- list(observed = observed,nb_species=nb_species_chg,randomizeMeanSD = random, randomize = randomize ,SESandPval = rbind(SES,pval))	
  }else{
    result <- list(observed = observed,nb_species=nb_species_chg,randomizeMeanSD = random, SESandPval = rbind(SES,pval))	
  }
  return(result)
}# end of FD.null

############################################################################################################ 
# Example
############################################################################################################ 

showexample <- T

if ( showexample == T) {
  
  taxo <- rbind(c("GenA","AA"),c("GenA","AB"),c("GenA","AC"),c("GenA","AD"),c("GenB","BA"),c("GenB","BB"),c("GenB","BC"),c("GenB","BD"),c("GenC","CA"),c("GenC","CB"),c("GenC","CC"),c("GenD","DD"),c("GenE","EA"),c("GenE","EB"),c("GenE","EC"),c("GenE","ED"),c("GenF","FA"),c("GenF","FB"),c("GenF","FC"),c("GenF","FD")) ; colnames(taxo) <- c("Genus","species"); rownames(taxo) <- taxo[,2]
  IUCN<-sample(c("CR","EN","VU","NC","LC","edn"),nrow(taxo),replace=T)
  names(IUCN)<-rownames(taxo)
  traits <- cbind(sample(rnorm(100,0.7,0.10),nrow(taxo)),sample(rnorm(100,0.1,0.25),nrow(taxo))) ; rownames(traits)<- rownames(taxo) 
  id_nat<-sample(c("exo","translo"),nrow(taxo),replace=T)
  names(id_nat)<-rownames(traits)
  comm_histo<-rbind (round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1))) ; colnames(comm_histo)<-rownames(taxo) ; rownames(comm_histo) <- paste("Comm",1:3)
  comm_curr<-rbind (round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1)),round(runif(nrow(taxo),0,1))) ; colnames(comm_curr)<-rownames(taxo) ; rownames(comm_curr) <- paste("Comm",1:3)
  print(FD.null(traits,comm_histo,comm_curr,IUCN,taxo,boot = 999, save = T))
}

############################################################################################################ 
# END
############################################################################################################ 