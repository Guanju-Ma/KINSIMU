#' @title Calculating parameters in different pairwise LR calculation
#' @description Counting and calculating parameters used in the calculation process of different pairwise LR calculation
#'
#' @param AB Genotypes of two individuals of each case, which should be data.frame with 4 columns (2 for each individual) and \code{ss} rows, where \code{ss} stand for sample size;
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param stepwisePI If TRUE, empirical decreasing model of STR mutation would be taken when paternity index is needed to be calculated, otherwise, mutation rate would be taken as PImu if IBS=0 between an alleged PC pair.
#' @param bred if TRUE, parameters used in the calculation of LR in inbred relationship would be calculated
#' @param mu mutation rate for PI calculation, with a default of 0.002
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#'
#' @return a data frame with one column and ss rows, containing log10 value of the CLR of each case
#' @export
#'

LRparas<-function(AB,af=NULL,rare=NULL,allelename=FALSE,stepwisePI=FALSE,bred=FALSE,mu=0.002){
  if (ncol(AB)!=4) {
    stop(paste("false in individual data"))
  }
  ac<-as.double(AB[,1]==AB[,3])
  ad<-as.double(AB[,1]==AB[,4])
  bc<-as.double(AB[,2]==AB[,3])
  bd<-as.double(AB[,2]==AB[,4])
  para<-data.frame(ibs=as.double(ac+ad+bc+bd>0)+as.double(ac*bd+ad*bc>0))
  if (isTRUE(allelename)) {
    if (is.null(af) || is.null(rare)) {
      if (isTRUE(stepwisePI) || isTRUE(bred)) {
        stop("Please input the frequency data of regular and rare alleles")
      }
    } else {
      pc<-af[as.character(AB[,3]),]
      pd<-af[as.character(AB[,4]),]
      pc[is.na(pc)]<-rare
      pd[is.na(pd)]<-rare
      pc<-as.numeric(pc)
      pd<-as.numeric(pd)
      para$LRid<-(ac*bd+ad*bc)/(2*pc*pd)
      para$PInomu<-((ac+bc)/pc+(ad+bd)/pd)/4 # Paternity index without considering mutation
      if (all(para$PInomu>0)) {
        para$PImu=para$PInomu
      } else if (isTRUE(stepwisePI)) {
        acm=abs(AB[,1]-AB[,3])+abs(AB[,1]-AB[,3])%%1*100000
        adm=abs(AB[,1]-AB[,4])+abs(AB[,1]-AB[,4])%%1*100000
        bcm=abs(AB[,2]-AB[,3])+abs(AB[,2]-AB[,3])%%1*100000
        bdm=abs(AB[,2]-AB[,4])+abs(AB[,2]-AB[,4])%%1*100000
        steps<-pmin(acm,adm,bcm,bdm)
        dc<-as.double(acm==steps)+as.double(bcm==steps)
        dd<-as.double(adm==steps)+as.double(bdm==steps)
        para$PImu<-as.double(steps==0)*(dc/pc+dd/pd)/4+ # ignore mutation if IBS>0
          as.double(steps>10000)*mu+ # take mutation rate as PI if there is no integer step of mutation
          as.double(steps>0 & steps<10000)*mu*10^(1-steps)*(dc/pc+dd/pd)/8 # step wise model of mutation
      } else {
        para$PImu<-para$PInomu*as.double(para$ibs>0)+mu*as.double(para$ibs==0) # Paternity index taking mutation rate as result if IBS=0
      }
      if (isTRUE(bred)) {
        ab<-as.double(AB[,1]==AB[,2])
        cd<-as.double(AB[,4]==AB[,3])
        pa<-af[as.character(AB[,1]),]
        pa[is.na(pa)]<-rare
        pa<-as.numeric(pa)
        para$LD1<-ab*ac*ad/(pa^3)
        para$LD2<-ab*cd/(pa*pc)
        para$LD3<-ab*(ac+ad)/(2*pa^2)
        para$LD4<-ab/pa
        para$LD5<-(ac+bc)*cd/(2*pc^2)
        para$LD6<-cd/pc
      }
    }
  } else {
    if (is.null(af)) {
      if(isTRUE(stepwisePI) || isTRUE(bred)) {
        stop("Please input the frequency data")
      }
    } else {
      pc<-af[AB[,3],]
      pd<-af[AB[,4],]
      para$LRid<-(ac*bd+ad*bc)/(2*pc*pd)
      para$PInomu<-((ac+bc)/pc+(ad+bd)/pd)/4 # Paternity index without considering mutation
      if (all(para$PInomu>0)) {
        para$PImu<-para$PInomu
      } else if (isTRUE(stepwisePI)) {
          an<-as.data.frame(as.numeric(row.names(af)))
          acm=abs(an[AB[,1],]-an[AB[,3],])+abs(an[AB[,1],]-an[AB[,3],])%%1*100000
          adm=abs(an[AB[,1],]-an[AB[,4],])+abs(an[AB[,1],]-an[AB[,4],])%%1*100000
          bcm=abs(an[AB[,2],]-an[AB[,3],])+abs(an[AB[,2],]-an[AB[,3],])%%1*100000
          bdm=abs(an[AB[,2],]-an[AB[,4],])+abs(an[AB[,2],]-an[AB[,4],])%%1*100000
          steps<-pmin(acm,adm,bcm,bdm)
          dc<-as.double(acm==steps)+as.double(bcm==steps)
          dd<-as.double(adm==steps)+as.double(bdm==steps)
          para$PImu<-as.double(steps==0)*para$PInomu+ # ignore mutation if IBS>0
            as.double(steps>10000)*mu+ # take mutation rate as PI if there is no integer step of mutation
            as.double(steps>0 & steps<10000)*mu*10^(1-steps)*(dc/pc+dd/pd)/8 # step wise model of mutation
      } else {
        para$PImu<-para$PInomu*as.double(para$ibs>0)+mu*as.double(para$ibs==0) # Paternity index taking mutation rate as result if IBS=0
      }
      if (isTRUE(bred)) {
        ab<-as.double(AB[,1]==AB[,2])
        cd<-as.double(AB[,4]==AB[,3])
        pa<-af[AB[,1],]
        para$FD1<-ab*ac*ad/(pa^3)
        para$FD2<-ab*cd/(pa*pc)
        para$FD3<-ab*(ac+ad)/(2*pa^2)
        para$FD4<-ab/pa
        para$FD5<-(ac+bc)*cd/(2*pc^2)
        para$FD6<-cd/pc
      }
    }
  }
return(para)
}


