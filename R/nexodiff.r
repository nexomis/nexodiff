#' nexodiff: Nexomis package to establish RNA-seq differential
#' signatures
#'
#' This package provides different classes that represents different steps in
#' expression data analysis.
#'
#' Those steps are specified in the schematic below with the names of the
#' different class.
#'
#' * line arrow means that a class is used and embeded in another.
#' * dashed arrow means that the class is used but not embeded.
#'
#' \if{html}{\figure{schema.svg}{options: width=100% alt="Figure: schema.svg"}}
#' \if{latex}{\figure{schema.png}{options: width=6in}}
#'
#' Below are listed the different class with links to the documentation:
#' * \href{./PairwiseDesignWithAnnotation.html}{PairwiseDesignWithAnnotation}
#' * \href{./ExprData.html}{ExprData}
#'   * \href{./ExprDataTranscript.html}{ExprDataTranscript}
#'   * \href{./ExprDataGene.html}{ExprDataGene}
#' * \href{./PairwiseLFC.html}{PairwiseLFC}
#' * \href{./PairwiseDESeq2.html}{PairwiseDESeq2}
#'
#' @docType package
#' @name nexodiff
NULL
