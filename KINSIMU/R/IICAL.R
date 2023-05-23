#' @title calculation of Incest index
#' @description Calculation of the ratio for a parent-child pair between the likelihood that the child's father is a relative of the parent to that is unrelated
#'
#' @param Parent Genotypes of individual A of each case, which should be data.frame with 2 columns and ss rows, where ss stand for sample size
#' @param Child Genotypes of individual B of each case, which should be data.frame with 2 columns and ss rows, where ss stand for sample size
#' @param af name of allele frequency matrix, a data.frame of 1 column with the allele name being row names
#' @param rare frequency of rare allele
#' @param theta kinship coefficent between the parents under Hp (that under Hd euquals to 0)
#'
#' @return a data.frame containing 3 columns: relationship, IBS and PI for each simulation
#' @export
#'

IICAL<-function(Parent,Child,af,rare,theta){
  if (ncol(Parent)!=2 || ncol(Child)!=2 || nrow(Parent)!=nrow(Child)) {
    stop(paste("false in individual data"))
  }
  colnames(Parent)=c("P","M")
  colnames(Child)=c("P","M")
  pc<-af[as.character(Child$P),]
  pd<-af[as.character(Child$M),]
  pc[is.na(pc)]<-rare
  pd[is.na(pd)]<-rare
  pc<-as.numeric(pc)
  pd<-as.numeric(pd)
  dc<-as.double(Parent$P==Child$P)+as.double(Parent$M==Child$P)
  dd<-as.double(Parent$P==Child$M)+as.double(Parent$M==Child$M)
  Ngs<-1-as.double(dc*dd==0)
  LR<-as.double(Ngs>0)*(2*theta*dc*dd/(dc*pc+dd*pd))+1-2*theta
  results<-data.frame(Ngs=Ngs,IItheta=LR)
  return(results)
}
