#' @title Pedigree simulation
#' @description Generating genotype combinations of multiple pedigrees with specific relationships on an autosomal marker.
#'
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names, which can be loaded with "EvaluatePanel" function, not necessary if Parent is not NULL
#' @param ss sample size, i.e., how many individual pairs do you want simulate
#' @param pedi a matrix in data.frame form containing the pedigree structure information
#' @param random_name the name of random individual, with a default of "RI"
#' @param muf father-child mutation rate
#' @param mum mother-child mutation rate
#' @param allelename if TRUE the output data would be the allele name, which should be the row names of af matrix, otherwise, it would be the position in that matrix
#'
#'
#' @return a list of n data.frame of 2 columns and ss rows containing genotypes of individuals, in which n equals to the row number of "pedi" frame
#' @export
#'
#' @examples
#' #simulating a first cousin pedigree
#' pedi<-pediexample
#' af<-FortytwoSTR$afmatrix[[1]]
#' pedisimu(af = af,ss = 10000,pedi = pedi)
#'
#'

pedisimu<-function(af,ss,pedi,random_name='RI',muf=0,mum=0,allelename=FALSE){
  if (ncol(pedi)!=3) {
    stop("Too many columns in pedi matrix")
  }
  colnames(pedi)<-c("Person","Father","Mother")
  if (any(duplicated(pedi$Person))) {
    stop("Duplicate personnel names")
  }
  if (length(which(pedi$Person==random_name))>0) {
    stop("Random individual should not be generated")
  }
  if (isFALSE(allelename) & (muf>0 || mum>0)) {
    stop("Allele name should be in/output if mutation rate is set as larger than 0")
  }
  np<-nrow(pedi)
  results<-as.data.frame(matrix(data = 0,nrow = ss,ncol = 2*np))
  colname<-data.frame(a=rep(pedi$Person,each=2),b=rep(c(1,2),np))
  colname$c<-paste('indi_',colname$a,"_allele",colname$b,sep = "")
  colnames(results)<-colname$c

  for (i in 1:np) {
    if (length(which(pedi$Person==pedi$Pather[i]))==0) {
      if (isTRUE(allelename)) {
        results[,2*i-1]<-sample(x = as.numeric(row.names(af)),size = ss,replace = TRUE,prob = af$Freq)
      } else {
        pop<-1:nrow(af)
        results[,2*i-1]<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
      }
      # results[,2*i-1]<-ptoc(Parent = NULL,af = af,ss = ss,mu = muf,allelename=allelename)
      if (pedi$Father[i]!=random_name) {
        print(paste("Individual",pedi$Father[i],"was not defined and set as random individual"))
      }
    } else if (which(pedi$Person==pedi$Father[i])>i) {
      stop(paste("Father of Individual",pedi$Person[i],"should be defined before him/her"))
    } else {
      f<-which(pedi$Person==pedi$Father[i])
      # results[,2*i-1]<-ptoc(Parent = as.data.frame(results[,(2*f-1):(2*f)]),mu = muf,allelename=allelename)
      RN<-sample(x = c(1,2),size = ss,replace = TRUE,prob = c(0.5,0.5))
      results[,2*i-1]<-results[,2*f-1]*as.double(RN==1)+results[,2*f]*as.double(RN==2)
      if (muf>0) {
        results[,2*i-1]=results[,2*i-1]+sample(x = c(-1,0,1),size = ss,replace = TRUE,prob = c(muf/2,1-muf,muf/2))
      }
    }
    if (length(which(pedi$Person==pedi$Mother[i]))==0) {
      if (isTRUE(allelename)) {
        results[,2*i]<-sample(x = as.numeric(row.names(af)),size = ss,replace = TRUE,prob = af$Freq)
      } else {
        pop<-1:nrow(af)
        results[,2*i]<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
      }
      # results[,2*i]<-ptoc(Parent = NULL,af = af,ss = ss,mu = mum,allelename=allelename)
      if (pedi$Mother[i]!=random_name) {
        print(paste("Individual",pedi$Mother[i],"was not defined and set as random individual"))
      }
    } else if (which(pedi$Person==pedi$Mother[i])>i) {
      stop(paste("Mother of Individual",pedi$Person[i],"should be defined before him/her"))
    } else {
      m<-which(pedi$Person==pedi$Mother[i])
      # results[,2*i]<-ptoc(Parent = as.data.frame(results[,(2*m-1):(2*m)]),mu = mum,allelename=allelename)
      RN<-sample(x = c(1,2),size = ss,replace = TRUE,prob = c(0.5,0.5))
      results[,2*i]<-results[,2*m-1]*as.double(RN==1)+results[,2*m]*as.double(RN==2)
      if (mum>0) {
        results[,2*i-1]=results[,2*i-1]+sample(x = c(-1,0,1),size = ss,replace = TRUE,prob = c(mum/2,1-mum,mum/2))
      }
    }
  }
  return(results)
}
