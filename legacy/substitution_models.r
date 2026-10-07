# Nucleotide substitution models (JC, K80, K81, F81, HKY, GTR) plus the
# gamma-based per-position rate helpers used for target and non-target sites.
# Every *_sub_rate_mat() returns a 4x4 matrix of per-timepoint substitution
# probabilities whose rows are the source base and whose columns are the
# destination base, both ordered A, G, C, T (the 1, 2, 3, 4 encoding used
# throughout the simulator). Sourced by sim5_code.R; attaches its own packages.
suppressPackageStartupMessages({
  library(phangorn)
  library(ggplot2)
  library(reshape2)
  library(docstring)
  library(ggpubr)
  library(psych)  
})



# all nucleotide matrices have order AGCT

#' Build a Jukes-Cantor (JC) substitution rate matrix
#'
#' One free parameter: all 12 non-self substitutions share a single rate and
#' the base composition is ignored.
#'
#' @param overall_sub_rate Numeric scalar. Per-timepoint probability assigned
#'   to each of the 12 off-diagonal from-base/to-base pairs.
#' @return A 4x4 numeric matrix. Rows are the source base and columns the
#'   destination base, both ordered A, G, C, T. Every off-diagonal entry is
#'   `overall_sub_rate` and the diagonal is 0. Entries are per-substitution
#'   probabilities, not a normalized transition distribution, so rows do not
#'   sum to 1.
jc_sub_rate_mat <- function(overall_sub_rate){
  # Jukes-Cantor model
  # overall_sub_rate: numeric. overall probability of a substitution occurring
  # return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively 
  # JC assumes equal base frequencies and equally likely substitutions between all combinations of bases
  
  # equal rate for all substitutions
  
  # assumes equal nucleotide frequencies
  
  rate_mat <- matrix(overall_sub_rate, nrow = 4, ncol = 4)
  diag(rate_mat) <- 0
  
  return(rate_mat)
}


#' Build a Kimura two-parameter (K80) substitution rate matrix
#'
#' Two free rates: one shared transition rate for A<->G and C<->T, and one
#' shared transversion rate for the four remaining base pairs. The three
#' arguments are a redundant parametrisation of those two rates, so any two of
#' them determine the third; whichever one is left `NULL` is derived. Base
#' composition is ignored.
#'
#' @param transition_to_transversion_ratio Numeric scalar or `NULL`. Equals
#'   `transition_rate / transversion_rate`.
#' @param transition_rate Numeric scalar or `NULL`. Probability of an A<->G or
#'   C<->T substitution.
#' @param transversion_rate Numeric scalar or `NULL`. Probability of any one of
#'   the four transversion substitutions (A<->C, A<->T, G<->C, G<->T).
#' @return A 4x4 numeric matrix of substitution probabilities, rows = source
#'   base and columns = destination base, both ordered A, G, C, T. The diagonal
#'   is 0, the (A,G), (G,A), (C,T) and (T,C) entries hold `transition_rate` and
#'   the remaining eight hold `transversion_rate`; rows do not sum to 1.
#' @note At least two of the three arguments must be supplied, all supplied
#'   values must be finite and non-negative, and a zero ratio cannot be used to
#'   derive `transversion_rate`; each of these raises an error.
k80_sub_rate_mat <- function(transition_to_transversion_ratio = NULL, transition_rate = NULL,
                             transversion_rate = NULL){
  # Kimura's 2-parameter (K80) model 
  # transition_to_transversion_ratio: numeric. ratio of transition mutation rate to transversion mutation rate (transversion rate in denominator)
  # transition_rate: numeric. Baseline uniform transition rate
  # transversion_rate: numeric. Baseline uniform transversion rate
  # return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively
  # K80 assumes equal base frequencies, one transition rate, and one transversion rate. 
  # At least two of transition_to_transversion_ratio, transition_rate, and transversion_rate must be specified
  
  # different rates for transitions and transversions
  

  # set rate of tranasversions to constant 1, and use transition_to_transversion ratio
  # to specify relative substitution rates
  
  # abbreviate input param
  ttv <- transition_to_transversion_ratio

  supplied <- c(!is.null(ttv), !is.null(transition_rate), !is.null(transversion_rate))
  if(sum(supplied) < 2){
    stop('K80 requires at least two of ratio, transition_rate, and transversion_rate.')
  }
  supplied_values <- c(ttv, transition_rate, transversion_rate)
  if(any(!is.finite(supplied_values)) || any(supplied_values < 0)){
    stop('K80 parameters must be finite and non-negative.')
  }
  if(is.null(transversion_rate) && ttv == 0){
    stop('Cannot derive transversion_rate from a zero transition/transversion ratio.')
  }
  
  # if the ratio is unspecified, calculate it
  if(is.null(ttv)){
    ttv <- transition_rate / transversion_rate
  }
  # if transition rate is unspecified, calculate it
  if(is.null(transition_rate)){
    transition_rate <- ttv * transversion_rate
  }
  # if transversion rate is unspecified, calculate it
  if(is.null(transversion_rate)){
    transversion_rate <- transition_rate / ttv
  }
  
  rate_mat <- rbind(
    c(0, transition_rate, transversion_rate, transversion_rate),
    c(transition_rate, 0, transversion_rate, transversion_rate),
    c(transversion_rate, transversion_rate, 0, transition_rate),
    c(transversion_rate, transversion_rate, transition_rate, 0)
  )
  
  return(rate_mat)
}


#' Build a Kimura three-parameter (K81) substitution rate matrix
#'
#' Three free rates: one transition rate and two transversion rates, one per
#' pair of complementary transversion classes. Base composition is ignored and
#' no normalization is applied.
#'
#' @param transition_rate Numeric scalar. Probability placed on A<->G and
#'   C<->T.
#' @param transversion_rate_weakstrong_conserved Numeric scalar. Probability
#'   placed on the A<->C and G<->T entries of the returned matrix.
#' @param transversion_rate_aminoketo_conserved Numeric scalar. Probability
#'   placed on the A<->T and G<->C entries of the returned matrix.
#' @return A 4x4 numeric matrix of substitution probabilities, rows = source
#'   base and columns = destination base, both ordered A, G, C, T. The matrix
#'   is symmetric with a 0 diagonal; rows do not sum to 1.
#' @note A<->C and G<->T are the amino- and keto-conserving transversions while
#'   A<->T and G<->C are the weak/strong-conserving ones, so each transversion
#'   argument is placed on the class opposite to the one its name describes.
#'   This only matters when the two transversion rates differ.
k81_sub_rate_mat <- function(transition_rate, 
                             transversion_rate_weakstrong_conserved, 
                             transversion_rate_aminoketo_conserved){
  # Kimura's 3-parameter (K81) model
  # transition_rate: numeric. The overall probability of a transition occurring.
  # transversion_rate_weakstrong_conserved: numeric. The overall probability of a A<->T or C<->G transversion occurring
  # transversion_rate_aminoketo_conserved: numeric. The overall probability of a A<->C or T<->G transversion occurring
  # return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively
  # K81 assumes equal base frequencies, two different transition rates, and two different transversion rates
  
  # different rates for (1) transitions, 
  # (2) transversions that maintain strength pairing properties (i.e. A/T, C/G)
  # (3) transversions that maintain certain chemical structures (i.e. A/C, G/T)
  
  # assumes equal nucleotide frequencies
  
  # rename for brevity
  ti <- transition_rate
  tv_ws <- transversion_rate_weakstrong_conserved
  tv_ak <- transversion_rate_aminoketo_conserved
  
  rate_mat <- rbind(
    c(0, ti, tv_ws, tv_ak),
    c(ti, 0, tv_ak, tv_ws),
    c(tv_ws, tv_ak, 0, ti),
    c(tv_ak, tv_ws, ti, 0)
  )
  
  return(rate_mat)
}

#' Build a Felsenstein 1981 (F81) substitution rate matrix
#'
#' Five arguments: the four nucleotide fractions of the sequence plus one
#' baseline rate. Each destination base is weighted by how far its fraction
#' deviates from 0.25 (`baseline * (1 + frac_dest - 0.25)`), which makes the
#' matrix constant down each column; the 12 off-diagonal entries are then
#' rescaled so that their mean is exactly `baseline_overall_sub_rate`.
#'
#' @param frac_a Numeric scalar. Fraction of the sequence that is A.
#' @param frac_g Numeric scalar. Fraction of the sequence that is G.
#' @param frac_c Numeric scalar. Fraction of the sequence that is C.
#' @param frac_t Numeric scalar. Fraction of the sequence that is T. The four
#'   fractions are expected to sum to 1; this is not checked.
#' @param baseline_overall_sub_rate Numeric scalar. Target mean probability
#'   across the 12 non-self substitutions.
#' @return A 4x4 numeric matrix of substitution probabilities, rows = source
#'   base and columns = destination base, both ordered A, G, C, T. The diagonal
#'   is 0, the three off-diagonal entries within a column are identical (the
#'   rate depends only on the destination base), and the mean of the 12
#'   off-diagonal entries equals `baseline_overall_sub_rate`; rows do not sum
#'   to 1. A `baseline_overall_sub_rate` of 0 returns the 4x4 zero matrix.
f81_sub_rate_mat <- function(frac_a,
                             frac_g,
                             frac_c,
                             frac_t,
                             baseline_overall_sub_rate){
  # Felsenstein 1981 (F81) model 
  # frac_a: numeric. The fraction of all nucleotides that are A
  # frac_g: numeric. The fraction of all nucleotides that are G
  # frac_c: numeric. The fraction of all nucleotides that are C
  # frac_t: numeric. The fraction of all nucleotides that are T
  # baseline_overall_sub_rate numeric. The average probability of a substitution occurring at a base at a given mutation timepoint
  # To-base-specific mutation rates will be greater or less than this value but will average to it
  # return 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively
  # F81 assumes variable base frequencies and equally substitution rates proportional to these nucleotide ratios
  
  # fracs refer to what fraction of the sequence each nucleotide comprises
  # substitution rates weighted by to-base's fraction of the entire sequence. 
  
  
  if(baseline_overall_sub_rate == 0){
    return(matrix(0, nrow = 4, ncol = 4))
  }

  # general approach is to compute the deviation away from 0.25 that each base fraction is,
  # then adjust the mean overall substitution rate accordingly 
  to_a_rate <- baseline_overall_sub_rate*(1 + (frac_a - 0.25))
  to_g_rate <- baseline_overall_sub_rate*(1 + (frac_g - 0.25))
  to_c_rate <- baseline_overall_sub_rate*(1 + (frac_c - 0.25))
  to_t_rate <- baseline_overall_sub_rate*(1 + (frac_t - 0.25))
  
  rate_mat <- matrix(rep(c(to_a_rate, to_g_rate, to_c_rate, to_t_rate), 4), nrow = 4, byrow = TRUE)
  diag(rate_mat) <- 0
  
  # normalize according to the 12 permissible substitutions such that the mean of the 12 permissible 
  # substitutions is the baseline_overall_sub_rate
  prenorm_mean <- sum(rate_mat)/12
  scale_factor <- (baseline_overall_sub_rate / prenorm_mean)
  norm_mat <- rate_mat * scale_factor
  
  diag(norm_mat) <- 0
  
  return(norm_mat)
}


#' Build a Hasegawa-Kishino-Yano (HKY) substitution rate matrix
#'
#' Combines base composition with a transition/transversion distinction: every
#' off-diagonal entry is a baseline rate (the transition rate for A<->G and
#' C<->T, the transversion rate for the other four pairs) multiplied by the
#' destination weight `1 + frac_dest - 0.25`. The whole matrix is then rescaled
#' so that the mean of the 12 off-diagonal entries equals the harmonic mean of
#' the two baseline rates (`psych::harmonic.mean`).
#'
#' @param frac_a Numeric scalar. Fraction of the sequence that is A.
#' @param frac_g Numeric scalar. Fraction of the sequence that is G.
#' @param frac_c Numeric scalar. Fraction of the sequence that is C.
#' @param frac_t Numeric scalar. Fraction of the sequence that is T. The four
#'   fractions are expected to sum to 1; this is not checked.
#' @param transition_to_transversion_ratio Numeric scalar or `NULL`. Equals
#'   `baseline_transition_rate / baseline_transversion_rate`.
#' @param baseline_transition_rate Numeric scalar or `NULL`. Pre-weighting rate
#'   for A<->G and C<->T.
#' @param baseline_transversion_rate Numeric scalar or `NULL`. Pre-weighting
#'   rate for the four transversion pairs.
#' @return A 4x4 numeric matrix of substitution probabilities, rows = source
#'   base and columns = destination base, both ordered A, G, C, T, with a 0
#'   diagonal; rows do not sum to 1. Rows are named `from_a_rate`,
#'   `from_g_rate`, `from_c_rate` and `from_t_rate` (a by-product of `rbind`)
#'   while columns are unnamed. Two zero baseline rates return the 4x4 zero
#'   matrix.
#' @note As in K80, at least two of the ratio and the two baseline rates must
#'   be supplied, supplied values must be finite and non-negative, and a zero
#'   ratio cannot be used to derive `baseline_transversion_rate`; each of these
#'   raises an error.
hky_sub_rate_mat <- function(frac_a, frac_g, frac_c, frac_t, transition_to_transversion_ratio = NULL,
                             baseline_transition_rate = NULL, baseline_transversion_rate = NULL){
  # Hasegawa-Kishino-Yano (HKY) model
  # frac_a: numeric. The fraction of all nucleotides that are A
  # frac_g: numeric. The fraction of all nucleotides that are G
  # frac_c: numeric. The fraction of all nucleotides that are C
  # frac_t: numeric. The fraction of all nucleotides that are T
  # transition_to_transversion_ratio numeric. The ratio of transition mutation rate to transversion mutation rate
  # baseline_transition_rate numeric. See details. Baseline transition rate that is further modified by to-base proportions in sequence
  # baseline_transversion_rate numeric. See details. Baseline transversion rate that is further modified by to-base proportions in sequence
  # returns 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively
  # HKY assumes variable base frequencies, one transition rate, and one transversion rate
  #' Note that mean_transition_rate and mean_transversion_rate are modified such that to-base substitution rates vary according to
  #' whether the to-base is a transition or transversion and according to to-base rates
  #' At least two of transition_to_transversion_ratio, baseline_transition_rate, and baseline_transversion_rate must be provided
  #' Non-self substitution rates are scaled such that overall mean substitution rate across 12 non-self substitutions is equal to
  #' the harmonic mean of baseline_transition_rate and baseline_transversion_rate
  
  
  # does not assume equal frequencies of nucleotides in sequence
  # accounts for different transition/transverison rates by fixing transversion rates at 1
  # and multiplying appropriate rates by transition/transversion ratios
  
  
  # rename
  ttv <- transition_to_transversion_ratio

  supplied <- c(
    !is.null(ttv),
    !is.null(baseline_transition_rate),
    !is.null(baseline_transversion_rate)
  )
  if(sum(supplied) < 2){
    stop(
      paste(
        'HKY requires at least two of ratio, baseline_transition_rate,',
        'and baseline_transversion_rate.'
      )
    )
  }
  supplied_values <- c(ttv, baseline_transition_rate, baseline_transversion_rate)
  if(any(!is.finite(supplied_values)) || any(supplied_values < 0)){
    stop('HKY rate parameters must be finite and non-negative.')
  }
  if(is.null(baseline_transversion_rate) && ttv == 0){
    stop('Cannot derive baseline_transversion_rate from a zero ratio.')
  }
  
  if(is.null(ttv)){
    ttv <- baseline_transition_rate / baseline_transversion_rate
  }
  if(is.null(baseline_transition_rate)){
    baseline_transition_rate <- ttv * baseline_transversion_rate
  }
  if(is.null(baseline_transversion_rate)){
    baseline_transversion_rate <- baseline_transition_rate / ttv
  }
  if(baseline_transition_rate == 0 && baseline_transversion_rate == 0){
    return(matrix(0, nrow = 4, ncol = 4))
  }
  
  # Account for both nucleotide composition and the transition/transversion
  # rates. Each off-diagonal rate is weighted by the destination nucleotide;
  # A<->G and C<->T use baseline_transition_rate and all other pairs use
  # baseline_transversion_rate.
  destination_weights <- 1 + (c(frac_a, frac_g, frac_c, frac_t) - 0.25)
  from_a_rate <- destination_weights * c(
    0, baseline_transition_rate, baseline_transversion_rate, baseline_transversion_rate
  )
  from_g_rate <- destination_weights * c(
    baseline_transition_rate, 0, baseline_transversion_rate, baseline_transversion_rate
  )
  from_c_rate <- destination_weights * c(
    baseline_transversion_rate, baseline_transversion_rate, 0, baseline_transition_rate
  )
  from_t_rate <- destination_weights * c(
    baseline_transversion_rate, baseline_transversion_rate, baseline_transition_rate, 0
  )
  
  prenormalized_mat <- rbind(from_a_rate,
                             from_g_rate,
                             from_c_rate,
                             from_t_rate)
  
  harmonic_mean <- harmonic.mean(c(baseline_transition_rate, baseline_transversion_rate))
  
  prenorm_mean <- sum(prenormalized_mat) / 12
  scale_factor <- (harmonic_mean / prenorm_mean)
  norm_mat <- prenormalized_mat * scale_factor
  
  
  
  
  return(norm_mat)
  
 
}


#' Build a General Time Reversible (GTR) substitution rate matrix
#'
#' Ten arguments: six symmetric exchangeabilities, one per unordered base pair,
#' and the four nucleotide fractions. The entry for a from/to pair is that
#' pair's exchangeability multiplied by the fraction of the destination base,
#' so the matrix satisfies `frac_from * rate[from, to] == frac_to *
#' rate[to, from]`. No normalization is applied and no argument is validated.
#'
#' @param ag_rate Numeric scalar. Exchangeability of the A/G pair.
#' @param ac_rate Numeric scalar. Exchangeability of the A/C pair.
#' @param at_rate Numeric scalar. Exchangeability of the A/T pair.
#' @param gc_rate Numeric scalar. Exchangeability of the G/C pair.
#' @param gt_rate Numeric scalar. Exchangeability of the G/T pair.
#' @param ct_rate Numeric scalar. Exchangeability of the C/T pair.
#' @param frac_a Numeric scalar. Fraction of the sequence that is A.
#' @param frac_g Numeric scalar. Fraction of the sequence that is G.
#' @param frac_c Numeric scalar. Fraction of the sequence that is C.
#' @param frac_t Numeric scalar. Fraction of the sequence that is T.
#' @return A 4x4 numeric matrix of substitution probabilities, rows = source
#'   base and columns = destination base, both ordered A, G, C, T, with a 0
#'   diagonal. The matrix is only symmetric when all four fractions are equal,
#'   and rows do not sum to 1.
gtr_sub_rate_mat <- function(ag_rate,
                             ac_rate,
                             at_rate,
                             gc_rate,
                             gt_rate,
                             ct_rate,
                             frac_a,
                             frac_g,
                             frac_c,
                             frac_t){
  # General Time Reversible (GTR) model
  # ag_rate: numeric. The overall probability of a A<->G substitution occurring.
  # ac_rate: numeric. The overall probability of a A<->C substitution occurring.
  # at_rate: numeric. The overall probability of a A<->T substitution occurring.
  # gc_rate: numeric. The overall probability of a G<->C substitution occurring.
  # gt_rate: numeric. The overall probability of a G<->T substitution occurring.
  # ct_rate: numeric. The overall probability of a C<->T substitution occurring.
  # frac_a: numeric. The fraction of all nucleotides that are A
  # frac_g: numeric. The fraction of all nucleotides that are G
  # frac_c: numeric. The fraction of all nucleotides that are C
  # frac_t: numeric. The fraction of all nucleotides that are T
  # returns 4x4 matrix containing substitution probabilities from A, G, C, T to A, G, C, T, respectively. 
  # GTR assumes variable base frequencies and pairwise-specific substitution rates between nucleotides
  
  # substitution rates specific to dinucleotide pairs
  # allows non-uniform distribution of nucleotides in sequence
  
 
  rate_mat <- rbind(
    c(0, ag_rate*frac_g, ac_rate*frac_c, at_rate*frac_t),
    c(ag_rate*frac_a, 0, gc_rate*frac_c, gt_rate*frac_t),
    c(ac_rate*frac_a, gc_rate*frac_g, 0, ct_rate*frac_t),
    c(at_rate*frac_a, gt_rate*frac_g, ct_rate*frac_c, 0)
  )
  
  return(rate_mat)
}



#' Draw per-target edit rates from a binned gamma distribution
#'
#' Unlike the *_sub_rate_mat() models, this assigns a rate to a position rather
#' than to a base pair: it bootstraps `num_bootstrap_draws` gamma values, splits
#' them at quantiles whose widths are the relative counts of the low, medium and
#' high target classes, then samples (with replacement) one draw from the
#' matching bin for each target position. Low-rate targets therefore draw from
#' the bottom of the gamma, medium from the middle and high from the top, so
#' baseline rate and rate heterogeneity are estimated together.
#'
#' @param sequence_length Integer. Length of the barcode sequence; part of the
#'   signature but not used by the body.
#' @param h_pos Integer vector. Positions of the high-edit-rate targets.
#' @param m_pos Integer vector. Positions of the medium-edit-rate targets.
#' @param l_pos Integer vector. Positions of the low-edit-rate targets.
#' @param shape_param Numeric. Shape of the gamma the rates are drawn from.
#' @param scale_param Numeric. Scale of that gamma.
#' @param num_bootstrap_draws Numeric. Number of gamma draws to bin and sample
#'   from.
#' @return A named list of per-position edit rates: names are target positions
#'   as characters, values are numeric per-timepoint rates. An empty list when
#'   there are no targets. Positions that appear in more than one class are
#'   collapsed to their maximum rate, which also leaves the list ordered by
#'   position name as a character sort rather than in H/M/L order.
#' @section Side effects: Consumes draws from the R random number stream
#'   (`rgamma` and `sample`).
#' @note Empty classes are dropped before the quantile cutpoints are computed,
#'   so a missing class does not shift the remaining ones. If every gamma draw
#'   comes back identical the function exits early and gives every target that
#'   one value, keeping duplicate position names if any exist.
SIMPLIFY_target_site_gamma_based_sub_rates <- function(sequence_length, h_pos, m_pos, l_pos, 
                                              shape_param = 0.5, scale_param = 0.001,
                                              num_bootstrap_draws = 1000){
  # Gamma-distributed mutation rate variation by categorical mutation rate class
  # sequence_length integer. The length of the barcode sequence
  # h_pos integer. A vector of integer barcode positions corresponding to targets with High edit rates
  # m_pos integer. A vector of integer barcode positions corresponding to targets with Medium edit rates
  # l_pos integer. A vector of integer barcode positions corresponding to targets with Low edit rates
  # shape_param numeric. The shape parameter of the edit rate gamma distribution.
  # scale_param numeric. The scale parameter of the edit rate gamma distribution.
  # num_bootstrap_draws numeric. The number of bootstrap draws 
  # returns list of position:edit rate for all BE targets
  # this function differs from other nucleotide substitution models in that baseline substitution rates AND heterogeneity
  # are simultaneously estimated using a discretized gamma distribution, based on the provided positions and counts of high-, 
  # medium-, and low- edit-rate targets
  
  # draw from a gamma distribution sequence_length number of times
  rgam_vals <- rgamma(n = num_bootstrap_draws, shape = shape_param, scale = scale_param)
  
  # if we didn't encode heterogeneity, exit early
  if(length(unique(rgam_vals)) == 1){
    all_target_positions <- c(h_pos, m_pos, l_pos)
    all_rates <- rep(rgam_vals[1], length(all_target_positions))
    
    pos_er_list <- as.list(all_rates)
    names(pos_er_list) <- all_target_positions
    return(pos_er_list)
  }
  
  
  
  total_num_targets <- length(l_pos) + length(m_pos) + length(h_pos)
  
  # if there are no targets, exit early
  if(total_num_targets == 0){
    return(list())
  }

  
  # account for potentially missing target classes
  class_lengths <- c('L' = length(l_pos),
                     'M' = length(m_pos),
                     'H' = length(h_pos))
  nonempty_classes <- names(class_lengths)[class_lengths > 0]
  prob_breakpoints <- c(0, cumsum(class_lengths[nonempty_classes]) / total_num_targets)
  
  
  # set quantile cutpoints at the levels corresponding to the relative numbers of HML targets
  cuts_quant <- quantile(rgam_vals, probs = prob_breakpoints, names = FALSE)
  new_cuts <- cut(rgam_vals, breaks = cuts_quant, labels = nonempty_classes, include.lowest = TRUE)
  
  val_df <- data.frame('gam' = rgam_vals,
                       'bin' = new_cuts)
  
  # sample high/medium/low edit rates from the rgam values 
  all_rates <- c()
  for(erc in c('H', 'M', 'L')){
    if(erc %in% nonempty_classes){
      all_rates <- append(all_rates, sample(val_df$gam[val_df$bin == erc], size = as.integer(class_lengths[erc]), replace = TRUE))
    }
  }

  # create list of pos:edit_rate or pos:edit_rate_class list
  all_target_positions <- c(h_pos, m_pos, l_pos)

  pos_er_list <- as.list(all_rates)
  names(pos_er_list) <- all_target_positions
  
  # if there are multiple rates per position, retain the maximum:
  filt_pos_er_list <- tapply(unlist(pos_er_list),
                             names(unlist(pos_er_list)),
                             max)
  pos_er_list <- as.list(filt_pos_er_list)
 
  
  return(pos_er_list)
  

}


#' Choose the positions that are forced to be invariant
#'
#' Selects `round(frac_invariant * length(eligible_invariant_sites))` positions
#' without replacement; callers zero the substitution probability at these
#' positions. Restricting the eligible set is how target sites are kept out of
#' the invariant pool.
#'
#' @param eligible_invariant_sites Integer vector of positions allowed to
#'   become invariant.
#' @param frac_invariant Numeric in [0, 1]. Fraction of those positions to make
#'   invariant.
#' @return An integer vector of the sampled positions; empty when the rounded
#'   count is 0.
#' @section Side effects: Consumes draws from the R random number stream
#'   (`sample`).
#' @note `sample()` is called on `eligible_invariant_sites` directly, so a
#'   single-element eligible set of value n is interpreted as `1:n`.
nontarget_get_invariant_inds <- function(eligible_invariant_sites, 
                                         frac_invariant){
  # Add invariant sites to background mutational processes
  # Force the substitution probability to zero at a specified fraction of non-target genomic sites
  # position_er_list: list. List with names == position numbers, values == edit rate at respective position
  # eligible_invariant_sites: integer. Vector of integers corresponding to the positions that are permitted to be invariant.
  # Possible use case would be preventing target sites from being becoming invariant.
  # frac_invariant: numeric. Numeric value indicating the fraction of non-target sites that are not permitted to mutate. 
  # returnsindices of barcode or mt that will be forced to zero (ie are invariant)
  
  # determine the number of genomic positions that will have forced-zero mutation rate
  num_invariant <- round(frac_invariant * length(eligible_invariant_sites))
  
  # randomly sample the inds to get the positions numbers of invariant sites
  invariant_inds <- sample(eligible_invariant_sites, size = num_invariant, 
                           replace = FALSE)
  
  return(invariant_inds)
  

}


#' Scale per-position rates by discretized gamma heterogeneity
#'
#' Bootstraps 10000 gamma draws, optionally reduces them to `num_discrete_bins`
#' equal-probability bins summarized by their mean or median, and multiplies
#' every rate in `position_er_list` by one scaling factor sampled from those
#' values. The sampling is done independently for each element inside a
#' position, so a position holding several rates (for example separate
#' transition and transversion rates) gets a separate factor per rate.
#'
#' @param position_er_list List whose names are position numbers and whose
#'   elements are the rate (or list of rates) at that position.
#' @param shape_param Numeric. Shape of the heterogeneity gamma; `0` returns
#'   `position_er_list` unchanged.
#' @param scale_param Numeric. Scale of that gamma; the default `1/shape_param`
#'   makes its expectation 1, so rates are preserved in expectation. The draws
#'   themselves are not re-normalized.
#' @param num_discrete_bins Integer >= 0. Number of equal-area bins the gamma is
#'   collapsed into; `0` samples scaling factors from the raw draws instead.
#' @param bin_agg_metric Character, `'mean'` or `'median'`. Statistic used to
#'   summarize each bin.
#' @return A copy of `position_er_list` with the same names and nesting, every
#'   rate multiplied by its own scaling factor.
#' @section Side effects: Consumes draws from the R random number stream
#'   (`rgamma` and `sample`).
#' @note Returns the input untouched for an empty list or `shape_param == 0`
#'   before any validation runs; otherwise a non-finite or non-positive gamma
#'   parameter, a `bin_agg_metric` other than mean/median, and a
#'   `num_discrete_bins` that is not one non-negative integer each raise an
#'   error.
nontarget_scale_gamma_heterogeneity <- function(position_er_list, shape_param = 0.5, scale_param = 1/shape_param,
                                                num_discrete_bins = 4, bin_agg_metric = 'mean'){
  # Scale substitution rates with stochastic gamma-distribution-based heterogeneity
  # Following specification of a baseline substitution probability matrix, scale substitution probability values 
  # by multiplying by gamma distribution draw. 
  # position_er_list list. List with names == position numbers, values == edit rate at respective position
  # shape_param numeric. The shape parameter of the edit rate gamma distribution.
  # scale_param numeric. The scale parameter of the edit rate gamma distribution; defaults to 1/shape_param so expected value is 1.
  # num_discrete_bins integer. The number of equal-area bins into which the gamma distribution should be divided. 
  # Increasing num_discrete_bins increases the resolution of the gamma distribution, thus increasing the number of possible scaling 
  # factors by which substitution rates can be multiplied. If num_discrete_bins == 0, then the gamma distribution is not discretized. 
  # bin_agg_metric character. Either 'mean' or 'median.' For each discretized bin, if applicable, summarize those values in 
  #' the respective bin by finding either the mean or median of the values falling in that bin. 
  # returns list of position:edit rate for all BE targets
  # This function generates a bootstrapped gamma distribution as specified by inputted shape and scale parameters. 
  # If specified, this distribution is divided into num_discrete_bins equally-weighted bins, each of which is characterized
  # by the bin_agg_metric function. Each edit rate in position_er_list is scaled by multiplying the originally-specified 
  # substitution rate by a random draw from these aggregated metrics (or a random draw from the entire distribution, if 
  # num_discrete_bins == 0)
  
  if(length(position_er_list) == 0 || shape_param == 0){
    return(position_er_list)
  }
  if(!is.finite(shape_param) || shape_param < 0 ||
     !is.finite(scale_param) || scale_param <= 0){
    stop('shape_param and scale_param must define a positive finite gamma distribution.')
  }
  if(!(bin_agg_metric %in% c('mean', 'median'))){
    stop("bin_agg_metric must be either 'mean' or 'median'.")
  }
  if(length(num_discrete_bins) != 1 || is.na(num_discrete_bins) ||
     num_discrete_bins < 0 || num_discrete_bins %% 1 != 0){
    stop('num_discrete_bins must be one non-negative integer.')
  }

  rgam_vals <- rgamma(n = 10000, shape = shape_param, scale = scale_param)
  
  if(length(unique(rgam_vals)) == 1){
    # if we're encoding zero heterogeneity, break early
    for(i in seq_along(position_er_list)){
      for(j in seq_along(position_er_list[[i]])){
        position_er_list[[i]][[j]] <- position_er_list[[i]][[j]] * rgam_vals[1]
      }
    }
    
    return(position_er_list)
  }
  
  
  if(num_discrete_bins != 0){
    # generate equally-spaced probability-breakpoints that will be used to 
    # find equal-area-under-the-curve breakpoints for gamma distribution
    prob_breakpoints <- seq(0, 1, length.out = num_discrete_bins + 1)
    
    # find numeric quantiles according to these breakpoints
    quantiles <- quantile(rgam_vals, probs = prob_breakpoints)
    
    # cut the data into labeled bins
    cuts <- cut(rgam_vals, breaks = quantiles, labels = seq(1, num_discrete_bins), include.lowest = TRUE)
    
    # empty vector to which either mean or median values of each bin will be added
    hetero_scales <- numeric()
    
    # cuts will be processed in ascending order from 1 to number of cuts
    for(breaknum in sort(unique(cuts))){
      if(bin_agg_metric == 'mean'){
        stat <- mean(rgam_vals[which(cuts == breaknum)])
      }
      else if(bin_agg_metric == 'median'){
        stat <- median(rgam_vals[which(cuts == breaknum)])
      }
      
      # append stat to growing vector
      hetero_scales <- c(hetero_scales, stat)

    }
    
  }
  
  # if the user chooses not to discretize the gamma distribution, the scaling factors will 
  # just be draws from the gamma distribution (again, normalized such that mean == 1)
  else if(num_discrete_bins == 0){
    hetero_scales <- rgam_vals
  }
  
  
  # now for each position, we randomly choose the stat from one of these bins as the scaling factor
  # for the substitution rate
  # when we have transversions, there are two possible rates. so we add noise to each of them with SEPARATE scaling factors
  # this is the relevance of the nested loop
  # in other cases, only one scaling factor per position will be necessary 
  # print(as.numeric(position_er_list))
  for(i in seq_along(position_er_list)){
    for(j in seq_along(position_er_list[[i]])){
      scaling_factor <- sample(hetero_scales, size = 1)[1]
      position_er_list[[i]][[j]] <- position_er_list[[i]][[j]] * scaling_factor
    }
  }
  
  
  return(position_er_list)
  
  
}
