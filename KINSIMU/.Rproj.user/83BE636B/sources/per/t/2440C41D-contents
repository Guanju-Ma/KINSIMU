#' @title LR in grandparent-child identification with the reference of a full-sibling of the alleged parent between them
#' @description LR when a child (B) is alleged to be a grand-child of another individual (GP), and an offspring(A) of the alleged grandparent participated, with or without the assistant of B's other parent (M). Hp is that, B is an offspring of A's full-sibling; and Hd that, B is unrelated to GP and A. Inbreeding is not considered.
#'
#'
#'
#' @param A Genotype data of the alleged uncle/aunt, should be data.frame with 2 columns and ss rows, where ss stand for sample size;
#' @param B Genotype data of the child, with the same form with \code{A}
#' @param GP Genotype data of the alleged grandparent of B, with the same form with \code{A}
#' @param P Genotype data of the other parent of B, which can be \code{NULL} or as that of \code{A}
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names
#' @param rare frequency of rare allele on the locus
#' @param allelename if TRUE, the input genotype data would be regarded as allelenames, otherwise, the position in the af matrix
#'
#' @return a data frame with one column and ss rows, containing log10 value of the CLR of each case
#' @details There might be no allele sharing between GP and A, or between M and B, if so, the part in LR would be output as 0, which can be further optimized in future version.
#' @export
#'

LRgpgcam<-function(A,B,GP,af,rare=NULL,allelename=FALSE,P=NULL){
  if (ncol(A)!=2 || ncol(B)!=2 || ncol(GP)!=2 || nrow(A)!=nrow(B) || nrow(A)!=nrow(GP) || nrow(B)!=nrow(GP)) {
    stop("false in individual data")
  }
  if (isTRUE(allelename)) {
    pa<-af[as.character(A[,1]),]
    pb<-af[as.character(A[,2]),]
    pc<-af[as.character(B[,1]),]
    pd<-af[as.character(B[,2]),]
    if (any(is.na(pa)) || any(is.na(pb)) || any(is.na(pc)) || any(is.na(pd))) {
      if (is.null(rare)) {
        stop("please input frequency data of rare alleles")
      }
      pa[is.na(pa)]<-rare
      pb[is.na(pb)]<-rare
      pc[is.na(pc)]<-rare
      pd[is.na(pd)]<-rare
    }
    pa<-as.numeric(pa)
    pb<-as.numeric(pb)
    pc<-as.numeric(pc)
    pd<-as.numeric(pd)
  } else {
    pa<-af[A[,1],]
    pb<-af[A[,2],]
    pc<-af[B[,1],]
    pd<-af[B[,2],]
  }
  if (isFALSE(is.null(P))) {
    if (ncol(P)!=2 || nrow(P)!=nrow(P)) {
      stop("false in indiviudal data")
    }
    dmc<-as.double(P[,1]==B[,1])+as.double(P[,2]==B[,1])
    dmd<-as.double(P[,1]==B[,2])+as.double(P[,2]==B[,2])
  } else {
    dmc<-pc
    dmd<-pd
  }
  dgpa<-as.double(P[,1]==A[,1])+as.double(P[,2]==A[,1])
  dgpb<-as.double(P[,1]==A[,2])+as.double(P[,2]==A[,2])
  dgpc<-as.double(P[,1]==B[,1])+as.double(P[,2]==B[,1])
  dgpd<-as.double(P[,1]==B[,2])+as.double(P[,2]==B[,2])
  ac<-as.double(A[,1]==B[,1])
  ad<-as.double(A[,2]==B[,1])
  bc<-as.double(A[,1]==B[,2])
  bd<-as.double(A[,2]==B[,2])

  LR1<-(pa*dgpb*ac*dmd+dgpa*pb*bc*dmd+pa*dgpb*ad*dmc+dgpa*pb*bd*dmc)/
    (4*(pa*dgpb+pb*dgpa)*(pc*dmd+pd*dmc))
  LR2<-(dmc*dgpd+dmd*dgpc)/(4*(pc*dmd+pd*dmc))
  LR1[is.na(LR1)]=0
  LR2[is.na(LR2)]=0
  results=data.frame(Log10CLR=log10(1/4+LR1+LR2))

  return(results)
}
