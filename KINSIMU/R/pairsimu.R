#' @title Pairwise simulation
#' @description Generating genotype combinations of multiple individual pairs with specific relationships on an autosomal marker, ignoring mutaion
#'
#' @param af name of allele frequency matrix, a data.frame of 1 column containing frequencies with allele names being row names, which can be loaded with "EvaluatePanel" function, not necessary if Parent is not NULL
#' @param ss sample size, i.e., how many individual pairs do you want simulate
#' @param delta distribution of IBD or Jacquard coefficient of specific relationship, i.e., kappa(0-2) or Delta(1-9), respectively. Which should be input in form of single row of data with c() function. It should be noted that the data should be in order of kappa0 to kappa2 or Delta1 to Delta9.
#' @param allelename if TRUE, outputing the names of alleles, otherwise, the positions of them in the af matrix
#'
#' @return a data.frame containing four columns and \code{ss} rows
#' @export
#'
#' @examples
#' # Take the first STR in the 42 STR as example
#' af = FortytwoSTR$afmatrix[[1]]
#' # simulating 10,000 unrelated pairs
#' a<-pairsimu(af = af,ss = 10000,delta = c(1,0,0),allelename = FALSE)
#' # simulating 10,000 parent-child pairs
#' b<-pairsimu(af = af,ss = 10000,delta = c(0,1,0),allelename = FALSE)
#' # simulating 10,000 full-sibling pairs
#' c<-pairsimu(af = af,ss = 10000,delta = c(0.25,0.5,0.25),allelename = FALSE)
#'

pairsimu<-function(af,ss,delta,allelename=FALSE){
  if (isFALSE(all(delta>=0)) || sum(delta)!=1 || isFALSE(length(delta) %in% c(3,9))) {
    stop(paste("false in delta distribution"))
  }
  if (ncol(af)>1) {
    stop("Error in allele frequency matrix")
  }
  colnames(af)<-"Freq"
  pop<-1:nrow(af)
  if (length(delta)==9 & all(delta[1:6]==0)) {
    delta=c(delta[9],delta[8],delta[7])
  }
  results<-data.frame(A1=sample(x = pop,size = ss,replace = TRUE,prob = af$Freq))
  if (length(delta)==3) {
    results$A2<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
    if (all(delta==c(1,0,0))) {
      # unrelated individuals
      results$B1<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
      results$B2<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
    } else if (all(delta==c(0,1,0))) {
      results$B1<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)
      # results$B2<-(results$A1+results$A2)/2+sample(x = c(-0.5,0.5),size = ss,replace = TRUE,prob = c(0.5,0.5))*(results$A1-results$A2)
      RN<-sample(x = c(1,2),size = ss,replace = TRUE,prob = c(0.5,0.5))
      results$B2<-results$A1*as.double(RN==1)+results$A2*as.double(RN==2)
    } else {
      # individuals related with outbred relationships
      RN<-sample(x = 1:4,size = ss,replace = TRUE,prob = c(delta[1],delta[2]/2,delta[2]/2,delta[3]))
      T1<-as.double(RN>2)
      T2<-as.double(RN%%2==0)
      results$B1<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)*(1-T1)+results$A1*T1
      results$B2<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)*(1-T2)+results$A2*T2
    }
  } else {
    # individuals related with inbred relationships
    RN<-sample(x = 1:11,size = ss,replace = TRUE,
               prob = c(delta[1],delta[2],delta[3],delta[4],delta[5]/2,delta[5]/2,delta[6],delta[7],
                        delta[8]/2,delta[8]/2,delta[9]))
    T0<-as.double(RN>4)
    T1A<-as.double(RN %in% c(1,3,5,8,9))
    T1B<-as.double(RN==6)
    T2A<-as.double(RN %in% c(1,8,10))
    T2B<-as.double(RN %in% c(2,5,6,7))
    results$A2<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)*T0+results$A1*(1-T0)
    results$B1<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)*(1-T1A-T1B)+results$A1*T1A+results$A2*T1B
    results$B2<-sample(x = pop,size = ss,replace = TRUE,prob = af$Freq)*(1-T2A-T2B)+results$A2*T2A+results$B1*T2B
  }
  if (isTRUE(allelename)) {
    an<-as.data.frame(as.numeric(row.names(af)))
    results$A1<-an[results$A1,]
    results$A2<-an[results$A2,]
    results$B1<-an[results$B1,]
    results$B2<-an[results$B2,]
  }
  return(results)
}
