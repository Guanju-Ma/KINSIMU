#' @docType data
#' @name pediexample
#' @title Example of pedi matrix
#' @description Example of pedi, a matrix containing inherence information of a pedigree, used in function of "pedisimu"
#' @format a data.frame containg 3 columns
#'
#'

pediexample<-data.frame(Person=c("GF","GM","F1","F2","M1","M2","A","B"),
                 Father=c("RI","RI","GF","GF","RI","RI","F1","F2"),
                 Mother=c("RI","RI","GM","GM","RI","RI","M1","M2"))
