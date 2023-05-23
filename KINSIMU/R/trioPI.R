#' @title Calculating LR in standard trio cases
#' @description Calculating LR in cases where 3 participants being available, a child, his/her biological mother (or father) whose parentage is confirmed, and a male (or female) who is unrelated to the confirmed parent and alleged to be specific relative of the child, usually the father. Null hypothesis, i.e., that the alleged participant is unrelated to the child, is taken as Hd.
#'
#' @param AR Genotype of the alleged relative, which should be data.frame with 2 columns and ss rows, where ss stand for sample size;
#' @param C Genotype of the child, see in \code{AR} for the data form
#' @param TP Genotype of the confirmed parent of the child, see in \code{AR} for the data form
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param muAtoC mutation rate from \code{AR} to \code{C} if AR is alleged to be the child's parent, with a default of 0.002, please note that mistakes would be introduced if the mutation rate is larger than 0.2
#' @param muTtoC mutation rate from \code{TP} to \code{C}, with a default of 0.002/3.5, please note that mistakes would be introduced if the mutation rate is larger than 0.2
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#' @param kappa1 kappa1 of the alleged relationship between AR and C, with a default of 1, meaning the AR is alleged to be the other parent of C
#'
#'
#' @details If any one of TP or AR cannot provide any of the C's alleles through integer steps, \code{muAtoC} would be output when calculating PI, and under other situations, the case with minimum mutation steps, (AR to C)+(TP to C) under Hp and (TP to C) under Hd, would be considered. Much more is required for further discussion in PI calculation involving mutations.
#'
#' @return a list of multiple parameters
#' @export
#'

trioPI<-function(AR,C,TP,af,rare=NULL,allelename=FALSE,muAtoC=0.002,muTtoC=0.002/3.5,kappa1=1){
  if (ncol(AR)!=2 || ncol(C)!=2 || ncol(TP)!=2 || ncol(AR)!=ncol(C) || ncol(AR)!=ncol(TP) || ncol(C)!=ncol(TP)) {
    stop(paste("false in individual data"))
  }
  if (kappa1>1) {
    stop(paste("kappa1 > 1"))
  }
  if (isTRUE(allelename)) {
    pc<-af[as.character(C[,1]),]
    pd<-af[as.character(C[,2]),]
    pc[is.na(pc)]<-rare
    pd[is.na(pd)]<-rare
    pc<-as.numeric(pc)
    pd<-as.numeric(pd)
  } else{
    pc<-af[C[,1],]
    pd<-af[C[,2],]
    an<-as.data.frame(as.numeric(row.names(af)))
    AR[,1]<-an[AR[,1],]
    AR[,2]<-an[AR[,2],]
    C[,1]<-an[C[,1],]
    C[,2]<-an[C[,2],]
    TP[,1]<-an[TP[,1],]
    TP[,2]<-an[TP[,2],]
  }
  dmc<-as.double(TP[,1]==C[,1])+as.double(TP[,2]==C[,1])
  dmd<-as.double(TP[,1]==C[,2])+as.double(TP[,2]==C[,2])
  dfc<-as.double(AR[,1]==C[,1])+as.double(AR[,2]==C[,1])
  dfd<-as.double(AR[,1]==C[,2])+as.double(AR[,2]==C[,2])
  if (kappa1<1) {
    if (isTRUE(all(dmc+dmd>0))) {
      results<-log10(1-kappa1+kappa1*(dfc*dmd+dfd*dmc)/(2*(dmc*pd+dmd*pc)))
    } else {
      ec<-round(abs(TP[,1]-C[,1])+abs(TP[,1]-C[,1])%%1*100000,1)
      ed<-round(abs(TP[,1]-C[,2])+abs(TP[,1]-C[,2])%%1*100000,1)
      fc<-round(abs(TP[,2]-C[,1])+abs(TP[,2]-C[,1])%%1*100000,1)
      fd<-round(abs(TP[,2]-C[,2])+abs(TP[,2]-C[,2])%%1*100000,1)
      msteps<-pmin(ec,ed,fc,fd)
      dmc<-as.double(ec==msteps)+as.double(fc==msteps)
      dmd<-as.double(ed==msteps)+as.double(fd==msteps)
      mstepc<-pmin(ec,fc)
      mstepd<-pmin(ed,fd)
      dmc2<-1+as.double(ec==fc)
      dmd2<-1+as.double(ed==fd)
      results<-log10(1-kappa1+
                       kappa1*as.double(msteps==0)*(dfc*dmd+dmc*dfd)/(2*(dmc*pd+dmd*pc))+
                       kappa1*as.double(msteps>0)*(dfc*dmd2*10^(msteps-mstepd)+dfd*dmc2*10^(msteps-mstepc))/(2*(dmc*pd+dmd*pc)))
    }
  } else {
    if (isTRUE(all(dmc*dfd+dmd*dfc>0))) {
      results<-log10((dfc*dmd+dfd*dmc)/(2*(dmc*pd+dmd*pc)))
    } else {
        ac<-round(abs(AR[,1]-C[,1])+abs(AR[,1]-C[,1])%%1*100000,1)
        ad<-round(abs(AR[,1]-C[,2])+abs(AR[,1]-C[,2])%%1*100000,1)
        bc<-round(abs(AR[,2]-C[,1])+abs(AR[,2]-C[,1])%%1*100000,1)
        bd<-round(abs(AR[,2]-C[,2])+abs(AR[,2]-C[,2])%%1*100000,1)
        ec<-round(abs(TP[,1]-C[,1])+abs(TP[,1]-C[,1])%%1*100000,1)
        ed<-round(abs(TP[,1]-C[,2])+abs(TP[,1]-C[,2])%%1*100000,1)
        fc<-round(abs(TP[,2]-C[,1])+abs(TP[,2]-C[,1])%%1*100000,1)
        fd<-round(abs(TP[,2]-C[,2])+abs(TP[,2]-C[,2])%%1*100000,1)
        f1m1<-ac+ed
        f1m2<-ac+fd
        f2m1<-bc+ed
        f2m2<-bc+fd
        m1f1<-ec+ad
        m1f2<-ec+bd
        m2f1<-fc+ad
        m2f2<-fc+bd

        ac[ac==0]<-1-log10(2/muAtoC)
        ad[ad==0]<-1-log10(2/muAtoC)
        bc[bc==0]<-1-log10(2/muAtoC)
        bd[bd==0]<-1-log10(2/muAtoC)
        ec[ec==0]<-1-log10(2/muTtoC)
        ed[ed==0]<-1-log10(2/muTtoC)
        fc[fc==0]<-1-log10(2/muTtoC)
        fd[fd==0]<-1-log10(2/muTtoC) # transforming 0 to make no mutation case suitable for mutation formula, less than 0 if mutation rate less than 0.2
        # Thus, if mutation rate is larger than 0.2, the formula would be misleading

        minall<-pmin(f1m1,f1m2,f2m1,f2m2,m1f1,m1f2,m2f1,m2f2)
        minm<-pmin(ac,ad,bc,bd)
        minf<-pmin(ec,ed,fc,fd)

        x<-as.double(f1m1==minall)*muAtoC*muTtoC*10^(2-ac-ed)/4+
          as.double(f1m2==minall)*muAtoC*muTtoC*10^(2-ac-fd)/4+
          as.double(f2m1==minall)*muAtoC*muTtoC*10^(2-bc-ed)/4+
          as.double(f2m2==minall)*muAtoC*muTtoC*10^(2-bc-fd)/4+
          as.double(m1f1==minall)*muAtoC*muTtoC*10^(2-ec-ad)/4+
          as.double(m1f2==minall)*muAtoC*muTtoC*10^(2-ec-bd)/4+
          as.double(m2f1==minall)*muAtoC*muTtoC*10^(2-fc-ad)/4+
          as.double(m2f2==minall)*muAtoC*muTtoC*10^(2-fc-bd)/4
        y<-as.double(ec==minm)*muTtoC*10^(1-ec)+
          as.double(ed==minm)*muTtoC*10^(1-ed)+
          as.double(fc==minm)*muTtoC*10^(1-fc)+
          as.double(fd==minm)*muTtoC*10^(1-fc)
        results<-log10(as.double(minall+minf+minm>10000)*muAtoC+
                         as.double(minf+minm+minall<10000)*x/y)
        }
    }
  return(results)
}


