#' @title Parent-child inherence
#' @description Generating alleles of multiple children on a single locus with one of their parents' genotypes
#'
#' @param Parent information of Parents, with a default of NULL, standing for random male, otherwise it should be the name of a existing data.frame containing the Parents' genotypes, 2 columns × ss rows;
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names, which can be loaded with "EvaluatePanel" function, not necessary if Parent is not NULL
#' @param ss sample size, i.e., how many children do you want simulate, not necessary if Parent is not NULL
#' @param mu mutation rate, with a default of 0
#' @param allelename if TRUE the output data would be the allele name, which should be the row names of af matrix, otherwise, it would be the position in that matrix
#'
#' @return a list of two data.frame of 2 columns and ss rows containing genotypes of individual pairs with specific relationships
#' @export
#'

ptoc<-function(Parent = NULL,af = NULL,ss = NULL,mu = 0,allelename=FALSE){
  if (is.null(Parent)) {
    if (is.null(ss) || is.null(af)) {
        stop("Please input the parent data or the population data")
    }
  } else {
    if (ncol(Parent)!=2) {
      stop("False in Parent data")
      }
    colnames(Parent)=c("P","M")
    if (isFALSE(is.null(ss))) {
      if (ss!=nrow(Parent)) {
        warning("Sample size did not equal to number of Parents and was reset")
      }
    }
    ss=nrow(Parent)
  }
  if (isFALSE(allelename) & mu>0) {
    stop("Allele name should be in/output if mutation rate is set as larger than 0")
  }
  results<-data.frame(rep(0,ss))
  if (is.null(Parent)) {
    if (isTRUE(allelename)) {
      results<-sample(x = as.numeric(row.names(af)),size = ss,replace = TRUE,prob = af$Freq)
    } else {
      pop<-1:nrow(af)
      results<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
    }
  } else{
    results<-(Parent$P+Parent$M)/2+sample(x = c(-0.5,0.5),size = ss,replace = TRUE,prob = c(0.5,0.5))*(Parent$P-Parent$M)
  }
  if (mu>0) {
    results=results+sample(x = c(-1,0,1),size = ss,replace = TRUE,prob = c(mu/2,1-mu,mu/2))
  }
  return(results)
}
