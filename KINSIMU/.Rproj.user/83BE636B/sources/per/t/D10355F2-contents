#' @title Evaluate a panel with systemic effectiveness indicators
#'
#' @param strpath name of .csv file containing af matrix, with the marker names listed in the first row
#' @param raremode mode of the calculation method of rare alleles,with a default of "ISFG", indicating the recommended method from ISFG, i.e., (X+1)/(2N+1), where X and N stood for the number of allele types detected in a survey and the sample size, respectively. Two alternative methods are given: "MAF", take the minimum allele frequency as such frequency; and "1/2N", take the minimum of possible frequency a survey can achieave.
#' @param Nind mode of sample size, with a default of "lastrow", meaning that the sample size was presented in the last row of the .csv file. An alternative method is given, i.e., input a unified sample size.
#'
#' @return list of four vectors: afmatrix, a list of allele frequency data of each locus; rare, a data.frame containing the frequency of rare allele on each locus; indicators, a data.frame containing parameters of system efficiency for each locus; panelparas, a data.frame containing system efficiency parameters for the whole panel, with the form of log10(1-paramter) to avoid the situation that the parameters being displayed as 1 because they were too close to 1
#' @export
#' @importFrom utils read.csv
#'
EvaluatePanel<-function (strpath, raremode = "ISFG", Nind = "lastrow"){
  allelefreq <- read.csv(file = strpath, header = T)
  allelefreq[is.na(allelefreq)]=0
  if (isFALSE(colnames(allelefreq)[1] %in% c("Allele", "allele"))) {
    stop("The first column should be allele names with the first cell being the letter Allele")
  }
  if (any(duplicated(colnames(allelefreq)))) {
    stop("Duplicate locus names")
  }
  nl <- ncol(allelefreq) - 1
  if (is.numeric(Nind)) {
    if (Nind < 1 || as.integer(Nind) != as.numeric(Nind)) {
      stop("Error in sample size data")
    }
    allelename <- as.numeric(allelefreq[, 1])
    allelefreq <- allelefreq[, 2:(nl + 1)]
    row.names(allelefreq) = allelename
    n_of_indi <- as.data.frame(matrix(nrow = 1, ncol = nl,
                                      data = as.vector(rep(Nind, nl))))
  } else if (Nind == "lastrow") {
    allelename <- as.numeric(allelefreq[1:(nrow(allelefreq) -
                                             1), 1])
    n_of_indi <- as.numeric(allelefreq[nrow(allelefreq),
                                       2:(nl + 1)])
    if (isFALSE(all(n_of_indi >= 1)) || isFALSE(all(as.integer(n_of_indi) ==
                                                    as.numeric(n_of_indi)))) {
      stop("Error in sample size data")
    }
    allelefreq <- allelefreq[1:(nrow(allelefreq) - 1), 2:(nl +
                                                            1)]
    row.names(allelefreq) = allelename
  } else {
    stop("False in N-of-individual mode")
  }
  if (isFALSE(all(allelefreq >= 0))) {
    stop("error: negative number exists")
  }
  if (isFALSE(all(allelefreq < 1))) {
    stop("error: there are frequencies larger than 1")
  }
  if (isFALSE(all(apply(allelefreq, 2, sum) == 1))) {
    warning(paste("allele frequencies of loci",
                colnames(allelefreq[,which(apply(allelefreq, 2, sum) != 1)]),
                "does not equal to 1, there may be some errors in the calculation results of this marker"))
  }
  indicators <- as.data.frame(matrix(nrow = 12, ncol = nl))
  colnames(indicators) <- colnames(allelefreq)
  row.names(indicators) <- c("n of alleles", "Het_unadjusted",
                             "Het_adjusted", "MAF", "PM", "DP", "PED", "PET", "PE double doubt",
                             "RGE", "RGE without mother", "PIC")
  af2 <- apply(allelefreq^2, 2, sum)
  af3 <- apply(allelefreq^3, 2, sum)
  af4 <- apply(allelefreq^4, 2, sum)
  af5 <- apply(allelefreq^5, 2, sum)
  af6 <- apply(allelefreq^6, 2, sum)
  af7 <- apply(allelefreq^7, 2, sum)
  afn0 <- allelefreq
  afn0[afn0 == 0] <- 2
  indicators[1, ] <- colSums(allelefreq != 0)
  indicators[2, ] <- 1 - af2
  indicators[3, ] <- indicators[2, ] * 2 * n_of_indi/
    (2 * n_of_indi - 1)
  indicators[4, ] <- apply(afn0, 2, min)
  indicators[5, ] <- 2 * af2^2 - af4
  indicators[6, ] <- 1 - indicators[5, ]
  indicators[7, ] <- 1 - 4 * af2 + 4 * af3 - 3 * af4 + 2 *
    af2^2
  indicators[8, ] <- 1 - 2 * af2 + af3 + 2 * af4 - 3 * af5 -
    2 * af2^2 + 3 * af2 * af3
  indicators[9, ] <- 1 + 4 * af4 - 4 * af5 - 3 * af6 - 8 *
    af2^2 + 2 * af3^2 + 8 * af2 * af3
  indicators[10, ] <- 1 - 4 * af2 + 6 * af3 - 17 * af5 + 28 *
    af6 - 15 * af7 - 4 * af2^2 + 18 * af2 * af3 - 16 * af2 *
    af4 + 5 * af2 * af5 - 12 * af3^2 + 10 * af3 * af4
  indicators[11, ] <- 1 - 8 * af2 + 16 * af3 - 26 * af4 +
    30 * af5 - 15 * af6 + 12 * af2^2 - 24 * af2 * af3 +
    8 * af2 * af4 + 6 * af3^2
  indicators[12, ] <- 1 - af2 - af2^2 + af4
  panelparas <- data.frame(parameters = c("TDP", "CPED", "CPET",
                                          "CPE double doubt", "CRGE", "CRGE without mother"),
                           values = rep(0, 6))
  for (i in 1:6) {
    panelparas$values[i] = log10(prod(1 - indicators[i +
                                                       5, ]))
  }
  rare <- as.data.frame(matrix(nrow = 1, ncol = nl))
  row.names(rare) = "rareallelefreq"
  colnames(rare) <- colnames(allelefreq)
  if (raremode == "MAF") {
    rare <- indicators[4, ]
    row.names(rare) = "rareallelefreq"
  } else if (raremode == "1/2N") {
    rare <- 0.5/n_of_indi
    row.names(rare) = "rareallelefreq"
  } else if (raremode == "ISFG") {
    rare <- (indicators[1, ] + 1)/(2 * n_of_indi + 1)
    row.names(rare) = "rareallelefreq"
  } else if (is.numeric(raremode)) {
    rare[] <- raremode
  } else {
    warning("Please set frequency data of rare allele on each locus in rare matrix, which is blank now")
  }
  afmatrix <- vector(mode = "list", length = nl)
  names(afmatrix) <- colnames(allelefreq)
  for (i in 1:nl) {
    afmatrix[[i]] <- data.frame(Freq = allelefreq[which(allelefreq[,
                                                                   i] != 0), i])
    row.names(afmatrix[[i]]) <- allelename[which(allelefreq[,
                                                            i] != 0)]
  }
  results <- list(afmatrix, as.data.frame(rare), as.data.frame(indicators),
                  as.data.frame(panelparas))
  names(results) <- c("afmatrix", "rare", "indicators", "panelparas")
  return(results)
}
