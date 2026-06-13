###################################################################################################################################################
#                                               R function to compute the FUNCTIONAL VULNERABILITY                                                #
#																							              Aurele Toussaint (PhD)                                                                #
#                                                                                                                                                 #
# 'Fvulne': function to compute and illustrate the proportion of the area (in 2dimensions) or volume loss if all vulnerable species disapeared    #
#         For details see Toussaint et al. 2016, Scientific Report (6:22125)                                                                      #
#																																																																                  #
#						Inputs : 	- traits : a matrix (S x T) with S the names of the species (lines) and T the functional traits (in colums)                 #
#										  - comm : a vector (length S) of S the presence/absence of the species in the focal C community                              #
#										  - vuln.sp : character chain of vulnerable species (e.g. IUCN, rare species)                       													#
#               N.B. vuln.sp must occur in 'comm'                                                                                                 #
# 																																																																                #
#						Outputs : - Fvuln : a value of the functional vulnerability in the C community                                                        #
#																- expressed as a percentage of the functional diversity considering all species																		#
#																																										                                                              #
#                     - Graphic representation on the functional space of the community and the proportion loss if all vulnerable species         #
#                                                                                                                            are removed          #
#																			(if T > 2, only the two firsts axes will be represented, for plot the other axes,                           #
#                                                                           change the 'i' and 'j' in the function arguments )                    #
###################################################################################################################################################
#rm(list=ls())

###################################################################################################################################################
# 'Fvulne' function
###################################################################################################################################################

Fvulne <- function (traits,comm,vuln.sp,i=1,j=2, 
                    graph = T, wd_save_graph=getwd(),
                    res=150, width=900, height=900, chck_imput_data=T,comnm=NULL,leg.pos="bottomright")
  {
  
  # saving name of current working directory
  current_wd<-getwd()
  
  # loading required libraries 
  library(betapart)
  library(geometry)
  
  ##########################
  # check inputs data
  if (chck_imput_data==T)
  {
    if( ncol(traits)<2 ) stop ("error: 'traits' must have at least 2 columns") 
    if( is.numeric(traits)==F ) stop ("error: traits values must be numeric") 
    if( length(which(is.na(traits)==T))!=0) stop ("error : NA in 'traits' matrix") 
    if( is.character(vuln.sp)==F ) stop ("error: 'foc.sp' must be a character chain")
    if( length(vuln.sp)!=sum(comm[vuln.sp])) stop ("error: all 'foc.sp' must be occur in 'comm'")
    if( sum(comm[vuln.sp])<(ncol(traits)+1) ) stop ("error: 'foc.spec' length must be higher than number of 'traits' columns")
  }
  
  ##########################
  # species present in the community
  comm.foc<-names(comm)[which(comm == 1)]
  
  ##########################
  # vulnerable species 
  VS<-vuln.sp
  
  ##########################
  # NOT vulnerable species
  OS<-comm.foc[!(comm.foc%in% VS)]
  
  #########################
  # scaling and centering of each trait according to all species values 
  traitsCS<-scale(traits, center=TRUE, scale=TRUE) 
  
  ##########################
  # calculate the functional diversity 
  tryCatch(F.all<-round(convhulln(traitsCS[c(VS,OS),],"FA")$vol,6), error=function(e){}) 
  tryCatch(F.allvuln<-round(convhulln(traitsCS[c(OS),],"FA")$vol,6), error=function(e){}) 
  
  ##########################
  # calculate the functional vulnerability
  Fvuln<-round(((F.all-F.allvuln)/F.allvuln)*100,6)
  
  ###########################################################
  # end of Fvulnerability computation
  ###########################################################
  
  ###########################################################
  # if graphical output
  if (graph == T)
  {
    ##########################
    # check graphical inputs 
    if(length(comnm) == 0){comnm="unnamed_graph"}
    
    #########################
    # saving plot in a jpeg files
    setwd(wd_save_graph)
    nmjpeg<-paste("Fvuln_",comnm,".jpeg",sep="")
    jpeg(file=nmjpeg, res=res, width=width, height=height)
    
    Xmin<- min(traits[,i]) ; Xmax<-max(traits[,i])
    Ymin<- min(traits[,j]) ; Ymax<-max(traits[,j])
    limX=c(Xmin,Xmax); limY=c(Ymin,Ymax)
  
    par(mar=c(5,5,3,3))
    plot(limX,limY,type="n",xlab=paste("Trait ",i,sep=""),ylab=paste("Trait ",j,sep=""),xlim=limX,ylim=limY) 
    title(main=paste("Nb Sp.vuln = ",length(VS)," (",round((length(VS)/(length(VS)+length(OS)))*100,2),"%) - Fvuln = ",round(Fvuln,2),"%",sep=""))
    
    #########################
    # functional space of all species (blue)
    pntsEnd.realm <- cbind(traits[comm.foc,i],traits[comm.foc,j]) # get traits of the species in 2 dimensions (i.e. 'i' and 'j')
    con.hull.pos.realm <- chull(pntsEnd.realm) # find positions of convex hull
    con.hull.realm <- rbind(pntsEnd.realm[con.hull.pos.realm,], pntsEnd.realm[con.hull.pos.realm[1],]) # get coordinates for convex hull
    polygon(con.hull.realm[,1], con.hull.realm[,2], col="#FF450095", border="black",angle=95) # draw polygon
  
    #########################
    # functional space of all species - vulnerable species (red)
    pntsEnd.realm <- cbind(traits[OS,i],traits[OS,j]) # get traits of the species in 2 dimensions (i.e. 'i' and 'j')
    con.hull.pos.realm <- chull(pntsEnd.realm) # find positions of convex hull
    con.hull.realm <- rbind(pntsEnd.realm[con.hull.pos.realm,], pntsEnd.realm[con.hull.pos.realm[1],]) # get coordinates for convex hull
    polygon(con.hull.realm[,1], con.hull.realm[,2], col="#D6D6D6", border="blue",angle=90) # draw polygon
    points(pntsEnd.realm[,1], pntsEnd.realm[,2],pch=16,cex=0.8,col= "blue") # plot species position
  
    # position of the Vulnerable species
    points(traits[VS,i],traits[VS,j],pch=16,cex=0.8,col= "red") # plot species position
    
    #########################
    # plot legend  
    leg.txt<-c("Not vulnerable sp.","Vulnerable sp.","Vulnerable funct. space")
    legend(leg.pos,leg.txt, pch = c(16,16,15), title = NULL,col=c("blue","red","#FF4500"),ncol=1,cex=0.7)
     
    # closing jpeg
    graphics.off()
    ################################
    
  }# end of graphical output 
  ###########################################################
  
  # returning to current working directory
  setwd(current_wd)
  
  #########################
  # returning resuts  
  return (Fvuln)
} # end of function

###################################################################################################################################################
# example
###################################################################################################################################################

showexample <- T

if ( showexample == T) {
  
  taxo<- rbind(cbind(rep("GenA",10),seq(1,10,1)),cbind(rep("GenB",10),seq(1,10,1)),cbind(rep("GenC",10),seq(1,10,1)),
               cbind(rep("GenD",10),seq(1,10,1)),cbind(rep("GenE",10),seq(1,10,1)),cbind(rep("GenF",10),seq(1,10,1)),
               cbind(rep("GenG",10),seq(1,10,1))); colnames(taxo) <- c("Genus","species"); rownames(taxo) <- paste(taxo[,1],taxo[,2])
         
               comm<-round(runif(nrow(taxo),0,1)) ; names(comm)<-rownames(taxo)
               comnm<-"comm_ex"
               traits <- cbind(sample(rnorm(100,0.7,0.10),nrow(taxo)),sample(rnorm(100,0.1,0.25),nrow(taxo))) ; rownames(traits)<- rownames(taxo)
               vuln.sp<-names(sample(which(comm ==1),0.3*length(which(comm ==1))))
               
               print(Fvulne(traits,comm,vuln.sp,i=1,j=2, graph = T, comnm=comnm))
}

###################################################################################################################################################
# end of example
###################################################################################################################################################