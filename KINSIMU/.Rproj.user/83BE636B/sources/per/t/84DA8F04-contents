#' @title Simulation in kinship analysis
#' @description Simulation and calculation in kinship analysis, in which genotype combinations of two groups of individual pairs with different relationships would be generated, specific likelihood ratio, as well as the CIBS score would be calculated for each pair.
#'
#' @param afmatrix name of allele frequency list, which can be loaded with "EvaluatePanel" function
#' @param ss sample size
#' @param tdelta distribution of the actual IBD or Jacquard coefficient of the individual pairs, i.e., kappa(0-2) or Delta(1-9), respectively. Which should be input in form of single row of data with c() function. It should be noted that the data should be in order of kappa0 to kappa2 or Delta1 to Delta9.
#' @param adelta3 distributions of the IBD coefficient of the outbred plaintiff's hypotheses in LR calculation, which should be a data.frame with 3 columns and x rows, where x stood for the number of such LR being calculated. The names of columns should be "k0", "k1" and "k2", and those of rows the name of LRs
#' @param adelta9 distributions of the Jacquard coefficient of the inbred plaintiff's hypotheses in LR calculation, which should be a data.frame with 9 columns and x rows, where x stood for the number of such LR being calculated. The names of columns should be "D1-D9", and those of rows the name of LRs
#' @param pedname name of the real relationship, defaults to "SimPed".
#' @param stepPI If TRUE, empirical decreasing model of STR mutation would be taken when paternity index is needed to be calculated, otherwise, mutation rate would be taken as PI if IBS=0 between an alleged PC pair.
#' @param mu mutation rate when paternity index is needed to be calculated, defualts to 0.002.
#'
#' @return a data.frame containing multiple columns: relationship, CIBS and multiple CLR for each simulation
#' @export
#'
#' @examples
#' adelta3<-data.frame(k0=c(0,0.25,0.5),k1=c(1,0.5,0.5),k2=c(0,0.25,0),row.names = c("PC","FS","HS"))
#' adelta9<-data.frame(D1=0,D2=0,D3=0,D4=0,D5=0.25,D6=0,D7=0.25,D8=0.5,D9=0,row.names = "FIMCpair")
#' data(FortytwoSTR)
#' results<-testsimulation(afmatrix = FortytwoSTR[["afmatrix"]],ss = 10000,tdelta = c(0,1,0),
#' adelta3 = adelta3, adelta9 = adelta9,pedname="PC")
#'
#'

testsimulation<-function(afmatrix,ss,tdelta,adelta3=NULL,adelta9=NULL,pedname="SimPed",stepPI=FALSE,mu=0.002){
  nl<-length(afmatrix)
  if (length(tdelta)==9 & all(tdelta[1:6]==0)) {
    tdelta=c(tdelta[9],tdelta[8],tdelta[7])
  }
  if (sum(tdelta)!=1 || any(tdelta<0) || isFALSE(length(tdelta) %in% c(3,9))) {
    stop("Error in actual coefficient distribution")
  }
  if (is.null(adelta3) & is.null(adelta9)) {
    warning("No alleged hypothesis was input, and only CIBS was output")
    results<-data.frame(Relationship=rep(pedname,ss),CIBS=rep(0,ss))
    for (i in 1:nl) {
      AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta,allelename = FALSE)
      results[,2]<-results[,2]+
        LRparas(AB = AB)$ibs
    }
  } else if (is.null(adelta9)) {
    if (any(apply(adelta3,1,sum)!=1) || any(adelta3<0)) {
      stop("Error in alleged IBD coefficient")
    }
    results<-as.data.frame(matrix(data = 0,nrow = ss,ncol = 2+nrow(adelta3)))
    results[,1]<-pedname
    colnames(results)<-c("Relationship","CIBS",paste("CLR_",row.names(adelta3),sep = ""))
    if (any(adelta3$k1==1) && isTRUE(stepPI)) {
      steppi<-TRUE
    } else {
      steppi<-FALSE
    }
    for (i in 1:nl) {
      AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta)
      para<-LRparas(AB = AB,af = afmatrix[[i]],stepwisePI = steppi,mu = mu)
      results[,2]<-results[,2]+para$ibs
      for (j in 1:nrow(adelta3)) {
        results[,j+2]<- results[,j+2]+log10(adelta3$k2[j]*para$LRid+
                                              as.numeric(adelta3$k1[j]==1)*adelta3$k1[j]*para$PImu+
                                              as.numeric(adelta3$k1[j]<1)*adelta3$k1[j]*para$PInomu+
                                              adelta3$k0[j])
      }
    }
  } else if (is.null(adelta3)) {
    if (any(apply(adelta9,1,sum)!=1) || any(adelta9<0)) {
      stop("Error in alleged Jacquard coefficient")
    }
    results<-as.data.frame(matrix(data = 0,nrow = ss,ncol = 2+nrow(adelta9)))
    results[,1]<-pedname
    colnames(results)<-c("Relationship","CIBS",paste("CLR_",row.names(adelta9),sep = ""))
    if (any(adelta9$D8==1) && isTRUE(stepPI)) {
      steppi<-TRUE
    } else {
      steppi<-FALSE
    }
    if (sum(adelta9$D1,adelta9$D2,adelta9$D3,adelta9$D4,adelta9$D5,adelta9$D6)==0) {
      for (i in 1:nl) {
        AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta)
        para<-LRparas(AB = AB,af = afmatrix[[i]],stepwisePI = steppi,mu = mu)
        results[,2]<-results[,2]+para$ibs
        for (j in 1:nrow(adelta9)) {
          results[,j+2]<-results[,j+2]+log10(adelta9$D7[j]*para$LRid+
                                               as.numeric(adelta9$D8[j]==1)*adelta9$D8[j]*para$PImu+
                                               as.numeric(adelta9$D8[j]<1)*adelta9$D8[j]*para$PInomu+
                                               adelta9$D9[j])
        }
      }
    } else{
      for (i in 1:nl) {
        AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta)
        para<-LRparas(AB = AB,af = afmatrix[[i]],stepwisePI = steppi,bred = TRUE,mu = mu)
        results[,2]<-results[,2]+para$ibs
        for (j in 1:nrow(adelta9)) {
          results[,j+2]<-results[,j+2]+log10(adelta9$D1[j]*para$FD1+
                                               adelta9$D2[j]*para$FD2+
                                               adelta9$D3[j]*para$FD3+
                                               adelta9$D4[j]*para$FD4+
                                               adelta9$D5[j]*para$FD5+
                                               adelta9$D6[j]*para$FD6+
                                               adelta9$D7[j]*para$LRid+
                                               as.numeric(adelta9$D9[j]==0)*adelta9$D8[j]*para$PImu+
                                               as.numeric(adelta9$D9[j]>0)*adelta9$D8[j]*para$PInomu+
                                               adelta9$D9[j])
        }
      }
    }
  } else {
    if (any(apply(adelta9,1,sum)!=1) || any(adelta9<0)) {
      stop("Error in alleged Jacquard coefficient")
    }
    if (any(apply(adelta3,1,sum)!=1) || any(adelta3<0)) {
      stop("Error in alleged IBD coefficient")
    }
    results<-as.data.frame(matrix(data = 0,nrow = ss,ncol = 2+nrow(adelta3)+nrow(adelta9)))
    results[,1]<-pedname
    colnames(results)<-c("Relationship","CIBS",paste("CLR_",row.names(adelta3),sep = ""),paste("CLR_",row.names(adelta9),sep = ""))
    if ((any(adelta9$D8==1) || any(adelta3$k1==1)) && isTRUE(stepPI)) {
      steppi<-TRUE
    } else {
      steppi<-FALSE
    }
    if (sum(adelta9$D1,adelta9$D2,adelta9$D3,adelta9$D4,adelta9$D5,adelta9$D6)==0) {
      for (i in 1:nl) {
        AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta)
        para<-LRparas(AB = AB,af = afmatrix[[i]],stepwisePI = steppi,mu = mu)
        results[,2]<-results[,2]+para$ibs
        for (j in 1:nrow(adelta3)) {
          results[,j+2]<- results[,j+2]+log10(adelta3$k2[j]*para$LRid+
                                                as.numeric(adelta3$k1[j]==1)*adelta3$k1[j]*para$PImu+
                                                as.numeric(adelta3$k1[j]<1)*adelta3$k1[j]*para$PInomu+
                                                adelta3$k0[j])
        }
        for (j in 1:nrow(adelta9)) {
          results[,j+2+nrow(adelta3)]<-results[,j+2+nrow(adelta3)]+log10(adelta9$D7[j]*para$LRid+
                                                                           as.numeric(adelta9$D8[j]==1)*adelta9$D8[j]*para$PImu+
                                                                           as.numeric(adelta9$D8[j]<1)*adelta9$D8[j]*para$PInomu+
                                                                           adelta9$D9[j])
        }
      }
    } else{
      for (i in 1:nl) {
        AB<-pairsimu(af = afmatrix[[i]],ss = ss,delta = tdelta)
        para<-LRparas(AB = AB,af = afmatrix[[i]],stepwisePI = steppi,bred = TRUE,mu = mu)
        results[,2]<-results[,2]+para$ibs
        for (j in 1:nrow(adelta3)) {
          results[,j+2]<- results[,j+2]+log10(adelta3$k2[j]*para$LRid+
                                                as.numeric(adelta3$k1[j]==1)*adelta3$k1[j]*para$PImu+
                                                as.numeric(adelta3$k1[j]<1)*adelta3$k1[j]*para$PInomu+
                                                adelta3$k0[j])
        }
        for (j in 1:nrow(adelta9)) {
          results[,j+2+nrow(adelta3)]<-results[,j+2+nrow(adelta3)]+log10(adelta9$D1[j]*para$FD1+
                                                                           adelta9$D2[j]*para$FD2+
                                                                           adelta9$D3[j]*para$FD3+
                                                                           adelta9$D4[j]*para$FD4+
                                                                           adelta9$D5[j]*para$FD5+
                                                                           adelta9$D6[j]*para$FD6+
                                                                           adelta9$D7[j]*para$LRid+
                                                                           as.numeric(adelta9$D9[j]==0)*adelta9$D8[j]*para$PImu+
                                                                           as.numeric(adelta9$D9[j]>0)*adelta9$D8[j]*para$PInomu+
                                                                           adelta9$D9[j])
        }
      }
    }
  }

  return(results)
}
