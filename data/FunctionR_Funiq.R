###################################################################################################################################################
#                                               R function to compute the FUNCTIONAL UNIQUENESS                                                   #
#																							              Aurele Toussaint (PhD)                                                                #
#                                                                                                                                                 #
# 'Funiq': function to compute and illustrate the proportion of the area (in 2dimensions) or volume of a taxonomic group not filled               #
#                                                                                                   by any others species of the community        #
#         For details see Toussaint et al. 2016, Scientific Report (6:22125)                                                                      #
#																																																																                  #
#						INPUTS : 	- traits : a matrix (S x T) with S the species (rows) and T the functional traits (in columns)                              #
#										  - comm : a vector (S) of the presence/absence of the species in the community                                               #
#										  - foc.sp : focal taxonomic level (F) for which functional uniqueness will be calculated 																		#
# 																																																																                #
#						OUTPUTS : - Funiq : a value of the functional uniqueness of the foc.sp in the C community                                             #
#																- expressed as a percentage of the functional volume of F species 					                                      #
#																																										                                                              #
#                     - Graphic representation on the functional space of the community and the focal taxonomic group (2 first dimensions)        #
#																			(if T > 2, only the two firsts axes will be represented, for plot the other axes,                           #
#                                                                           change the 'i' and 'j' in the function arguments )                    #
#																																																																                  #
#																			- polygon and points "blue" = functional space and species of the community without the focal species       #
#																			- polygon and points "red"  = functional space and species of the focal species              				        #
#																																															                                                    #
#   N.B. If the number of traits and/or the number of species is high, this function is time consuming                                            #
#                                                                                                                                                 #
###################################################################################################################################################
rm(list=ls())

###################################################################################################################################################
# Funiq function
Funiq <- function (traits,comm,foc.sp,i=1,j=2, 
                   graph = T, col_bck="#D6D6D699", col_frt="#CD853F99",wd_save_graph=getwd(),
                   res=150, width=900, height=900, chck_imput_data=T)
  {

  # saving name of current working directory
  current_wd<-getwd()
  
  # loading required libraries 
  library(betapart)
  
  ##########################
  # check inputs data
  if (chck_imput_data==T)
    {
    if( ncol(traits)<2 ) stop ("error: 'traits' must have at least 2 columns") 
    if( is.numeric(traits)==F ) stop ("error: traits values must be numeric") 
    if( is.character(foc.sp)==F ) stop ("error: 'foc.sp' must be a character chain")
    if( length(foc.sp)!=sum(comm[foc.sp])) stop ("error: all 'foc.sp' must be occur in 'comm'")
    if( sum(comm[foc.sp])<(ncol(traits)+1) ) stop("error: 'foc.spec' length must be higher than number of 'traits' columns")
  }
  
  ##########################
  # Focal taxonomic group (F)
  F<-foc.sp
  
  ##########################
  # others species occuring in the community (OS)
  OS<-names(comm)[!(names(comm)%in% F)]

  ##########################
  # compute functional beta diversity between F and OS
  comm.beta<-matrix(0,2,length(unique(c(F,OS))),dimnames=list( c(1,2) ,c(F,OS) ) )
  comm.beta[2, F]<-1
  comm.beta[1, OS]<-1
  tr1<-traits[match(colnames(comm.beta),rownames(traits)),]

  fctbeta<-functional.betapart.core(x=comm.beta, traits=as.matrix(tr1))

  ##########################
  # Functional volume (or area) shared (F.shared) and unshared (F.unshared) between F and OS
  F.shared<-fctbeta$shared[2,1] 
  F.unshared<-fctbeta$not.shared[2,1] 
  
  ##########################
  # Compute functional uniqueness of F species (Funiq)
  Funiq<-round((F.unshared/(F.unshared+F.shared))*100,6)
  
  ###########################################################
  # End of Funiq computation
  ###########################################################
  
  ###########################################################
  # if graphical output
  if (graph == T)
  {
    
  #########################
  # Saving plot in a jpeg files
  setwd(wd_save_graph)
  nmjpeg<-paste("Funiq_",strsplit(foc.sp," ")[[1]][[1]],".jpeg",sep="")
  jpeg(file=nmjpeg, res=res, width=width, height=height)

  Xmin<- min(traits[,i]) ; Xmax<-max(traits[,i])
  Ymin<- min(traits[,j]) ; Ymax<-max(traits[,j])
  limX=c(Xmin,Xmax); limY=c(Ymin,Ymax)
  
  par(mar=c(5,5,3,3))
  plot(limX,limY,type="n",xlab=paste("Trait ",i,sep=""),ylab=paste("Trait ",j,sep=""),xlim=limX,ylim=limY) 
  title(main=paste(strsplit(foc.sp," ")[[1]][[1]],"- Funiq =",round(Funiq,2),"%"))
  
  #########################
  # Functional space of OS (blue)
  pntsEnd.realm <- cbind(traits[OS,i],traits[OS,j]) # get traits of the species in 2 dimensions ('i' and 'j')
  con.hull.pos.realm <- chull(pntsEnd.realm) # find positions of convex hull  in the 2 dimensional plan ('i' and 'j')
  con.hull.realm <- rbind(pntsEnd.realm[con.hull.pos.realm,], pntsEnd.realm[con.hull.pos.realm[1],]) # get coordinates for convex hull
  polygon(con.hull.realm[,1], con.hull.realm[,2], col=col_bck, border="blue",angle=95) # draw polygon
  points(pntsEnd.realm[,1], pntsEnd.realm[,2],pch=16,cex=0.8,col= "blue") # plot species position
  
  #########################
  # Functional space of F (red)
  pntsEnd.realm <- cbind(traits[F,i],traits[F,j]) # get traits of the species in 2 dimensions ('i' and 'j')
  con.hull.pos.realm <- chull(pntsEnd.realm) # find positions of convex hull in the 2 dimensional plan ('i' and 'j')
  con.hull.realm <- rbind(pntsEnd.realm[con.hull.pos.realm,], pntsEnd.realm[con.hull.pos.realm[1],]) # get coordinates for convex hull
  polygon(con.hull.realm[,1], con.hull.realm[,2], col=col_frt, border="red",angle=90) # draw polygon
  points(pntsEnd.realm[,1], pntsEnd.realm[,2],pch=16,cex=0.8,col= "red") # plot species position
  
  #########################
  # Legend  
  leg.txt<-c(paste("All sp. -",strsplit(foc.sp," ")[[1]][[1]],"sp."),paste(strsplit(foc.sp," ")[[1]][[1]],"sp."))
  legend("bottomright",leg.txt, pch = 16, title = NULL,col=c("blue","red"),ncol=1,cex=1)
  
  # closing jpeg
  graphics.off()
  ################################

  }# end of graphical output 
  ###########################################################
  
  # returning to current working directory
  setwd(current_wd)
  
  #########################
  # returning resuts  
  return(Funiq=Funiq)
}  # end of Funiq function

###################################################################################################################################################
# Example
###################################################################################################################################################
  
  showexample <- T
  
  if ( showexample == T) {
    
    taxo<- rbind(cbind(rep("GenA",10),seq(1,10,1)),cbind(rep("GenB",10),seq(1,10,1)),cbind(rep("GenC",10),seq(1,10,1)),
                 cbind(rep("GenD",10),seq(1,10,1))); colnames(taxo) <- c("Genus","species"); rownames(taxo) <- paste(taxo[,1],taxo[,2])
    traits <- cbind(sample(rnorm(100,0.7,0.10),nrow(taxo)),sample(rnorm(100,0.1,0.25),nrow(taxo))) ; rownames(traits)<- rownames(taxo) 
  
    comm<-round(runif(nrow(taxo),0,1)) ; names(comm)<-rownames(taxo)

    # This example is designed to calculate the Funiq of the species belonging to the genus B ('GenB') whithin the community 'comm'
    foc.sp<-names(which(comm[paste(rep("GenB",10),seq(1,10,1))] ==1))

    print(Funiq(traits,comm,foc.sp,i=1,j=2, graph = T, col_bck="#D6D6D6", col_frt="#CD853F99"))
  }
  
###################################################################################################################################################
# END
###################################################################################################################################################