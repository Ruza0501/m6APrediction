#' Encode DNA sequences into factor data frame
#'
#' Convert DNA strings into a data frame with each nucleotide position as a factor.
#'
#' @param dna_strings Character vector of DNA sequences (same length).
#' @return Data frame with columns `nt_pos1`, `nt_pos2`, ... (factors: "A", "T", "C", "G").
#' @examples
#' if (FALSE) {
#'   m6APrediction:::dna_encoding(c("ATGAT", "CGTAC"))
#' }
#' @import randomForest
dna_encoding <- function(dna_strings){
  nn <- nchar( dna_strings[1] )
  seq_m <- matrix( unlist( strsplit(dna_strings, "") ), ncol = nn, byrow = TRUE)
  colnames(seq_m) <- paste0("nt_pos", 1:nn)
  seq_df <- as.data.frame(seq_m)
  seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
  return(seq_df)
}

#' Predict m6A for multiple sequences
#'
#' Use a trained model to predict m6A status for multiple sequences.
#'
#' @param ml_fit Trained model (e.g., random forest).
#' @param feature_df Data frame with columns: `gc_content`, `RNA_type`, `RNA_region`, `exon_length`, `distance_to_junction`, `evolutionary_conservation`, `DNA_5mer`.
#' @param positive_threshold Probability threshold for "Positive" (default: 0.5).
#' @return Input data frame with added `predicted_m6A_prob` and `predicted_m6A_status`.
#' @examples
#'
#' library(m6APrediction)
#'
#'
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#' feature_df <- read.csv(system.file("extdata", "m6A_input_example.csv", package = "m6APrediction"))
#'
#'
#' prediction_multiple(ml_fit, feature_df, positive_threshold = 0.6)
#' @export
prediction_multiple <- function(ml_fit, feature_df, positive_threshold = 0.5){
  stopifnot(all(c("gc_content", "RNA_type", "RNA_region", "exon_length", "distance_to_junction", "evolutionary_conservation", "DNA_5mer") %in% colnames(feature_df)))
  dna_encoding <- function(dna_strings) {
    nn <- nchar(dna_strings[1])
    seq_m <- matrix(unlist(strsplit(dna_strings, "")), ncol = nn, byrow = TRUE)
    colnames(seq_m) <- paste0("nt_pos", 1:nn)
    seq_df <- as.data.frame(seq_m)
    seq_df[] <- lapply(seq_df, factor, levels = c("A", "T", "C", "G"))
    return(seq_df)
  }

  feature_df$RNA_type <- factor(feature_df$RNA_type,
                                levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  feature_df$RNA_region <- factor(feature_df$RNA_region,
                                  levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  encoded_dna <- dna_encoding(feature_df$DNA_5mer)
  feature_df_encoded <- cbind(feature_df, encoded_dna)

  predicted_probs <- predict(ml_fit, newdata = feature_df_encoded, type = "prob")[, "Positive"]

  predicted_status <- ifelse(predicted_probs > positive_threshold, "Positive", "Negative")

  feature_df$predicted_m6A_prob <- predicted_probs
  feature_df$predicted_m6A_status <- predicted_status
  return(feature_df)
}

#' Predict m6A for a single sequence
#'
#' Use a trained model to predict m6A status for a single sequence.
#'
#' @param ml_fit Trained model (e.g., random forest).
#' @param gc_content Numeric GC content value.
#' @param RNA_type RNA type (one of: "mRNA", "lincRNA", "lncRNA", "pseudogene").
#' @param RNA_region RNA region (one of: "CDS", "intron", "3'UTR", "5'UTR").
#' @param exon_length Numeric exon length.
#' @param distance_to_junction Numeric distance to junction.
#' @param evolutionary_conservation Numeric conservation score.
#' @param DNA_5mer 5-mer DNA sequence string.
#' @param positive_threshold Probability threshold for "Positive" (default: 0.5).
#' @return Named vector with `predicted_m6A_prob` and `predicted_m6A_status`.
#' @examples
#'
#' ml_fit <- readRDS(system.file("extdata", "rf_fit.rds", package = "m6APrediction"))
#'
#' prediction_single(ml_fit, gc_content = 0.5, RNA_type = "mRNA",
#'                   RNA_region = "CDS", exon_length = 12,
#'                   distance_to_junction = 50, evolutionary_conservation = 0.7,
#'                   DNA_5mer = "ATGAT", positive_threshold = 0.4)
#' @export
prediction_single <- function(ml_fit, gc_content, RNA_type, RNA_region, exon_length, distance_to_junction, evolutionary_conservation, DNA_5mer, positive_threshold = 0.5){
  single_df <- data.frame(
    gc_content = gc_content,
    RNA_type = RNA_type,
    RNA_region = RNA_region,
    exon_length = exon_length,
    distance_to_junction = distance_to_junction,
    evolutionary_conservation = evolutionary_conservation,
    DNA_5mer = DNA_5mer,
    stringsAsFactors = FALSE
  )

  single_df$RNA_type <- factor(single_df$RNA_type,
                               levels = c("mRNA", "lincRNA", "lncRNA", "pseudogene"))
  single_df$RNA_region <- factor(single_df$RNA_region,
                                 levels = c("CDS", "intron", "3'UTR", "5'UTR"))

  result_df <- prediction_multiple(ml_fit, single_df, positive_threshold)
  returned_vector <- c(
    predicted_m6A_prob = result_df$predicted_m6A_prob,
    predicted_m6A_status = result_df$predicted_m6A_status
  )

  return(returned_vector)
}
