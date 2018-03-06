df <- data.table::fread("~/AndersenLab/Github_Repos/cegwas2/inst/extdata/test_phenotype.tsv")

# extract strain, isotype dataframe from database
generate_isotype_lookup <- function(species = "ce") {

  if( species == "ce" ) {
    isotype_lookup <- dplyr::collect(get_db("strain")) %>%
      dplyr::mutate(strain_names = ifelse(!is.na(previous_names),
                                          paste(strain, previous_names, sep="|"),
                                          strain)) %>%
      tidyr::separate_rows(strain_names, sep = "\\|") %>%
      dplyr::select(strain = strain_names, isotype)
  }

  return( isotype_lookup )
}

# strains2resolve is a vector of strain names
# isotype_lookup is a dataframe of strain, isotype
resolve_isotypes <- function( strains2resolve, isotype_lookup = generate_isotype_lookup() ){

  resolved_strains <- data.frame( strain = as.character(strains2resolve) ) %>%
    dplyr::left_join( ., isotype_lookup, by = "strain") %>%
    dplyr::select( -strain ) %>%
    dplyr::rename( strain = isotype ) %>%
    dplyr::pull( strain )

  if( sum(is.na(resolved_strains)) > 0 ){
    unresolved_strains <- sum(is.na(resolved_strains))
    message(glue::glue("~ ~ ~ WARNING ~ ~ ~
                       \n{unresolved_strains} strains were not resolved.
                       \n~ ~ ~ WARNING ~ ~ ~"))
  }

  return( resolved_strains )
}

# data = strain, trait, phenotype
BAMF_prune <- function( data, remove_outliers = TRUE ){

  categorize1 <- function(data) {
    with(data,
         ( (sixhs >= 1 & ( (s6h + s5h + s4h ) / numst) <= .05))
         | ((sixls >= 1 & ( (s6l + s5l + s4l) / numst) <= .05))
    )
  }

  # If the 5 innermost bins are discontinuous by more than a 1 bin gap, the
  # observation is in the fifth bin (between 7 and 10x IQR outside the
  # distribution), and the four outermost bins make up less than 5% of the
  # population, mark the observation an outlier

  categorize2 <- function(data) {
    with(data,
         ( (fivehs >= 1 & ( (s6h + s5h + s4h + s3h) / numst) <= .05))
         | ((fivels >= 1 & ( (s6l + s5l + s4l + s3l) / numst) <= .05))
    )
  }

  # If the 4 innermost bins are discontinuous by more than a 1 bin gap, the
  # observation is in the fourth bin (between 5 and 7x IQR outside the
  # distribution), and the four outermost bins make up less than 5% of the
  # population, mark the observation an outlier

  categorize3 <- function(data) {
    with(data,
         ( (fourhs >= 1 & (s5h + s4h + s3h + s2h) / numst <= .05))
             | ((fourls >= 1 & (s5l + s4l + s3l + s2l) / numst <= .05))
    )
  }

  napheno <- data[is.na(data$phenotype), ] %>%
    dplyr::mutate(bamfoutlier1 = NA, bamfoutlier2 = NA, bamfoutlier3 = NA)

  datawithoutliers <- data %>%
    # Filter out all of the wash and/or empty wells
    dplyr::filter(!is.na(strain)) %>%
    # Group by trait, the, calculate the first and third
    # quartiles for each of the traits
    dplyr::group_by(trait) %>%
    dplyr::summarise(iqr = IQR(phenotype, na.rm = TRUE),
                     q1 = quantile(phenotype, probs = .25, na.rm = TRUE),
                     q3 = quantile(phenotype, probs = .75,
                                   na.rm = TRUE)) %>%
    # Add a column for the boundaries of each of the bins
    dplyr::mutate(cut1h = q3 + (iqr * 2),
                  cut1l =q1 - (iqr * 2),
                  cut2h = q3 + (iqr * 3),
                  cut2l =q1 - (iqr * 3),
                  cut3h = q3 + (iqr * 4),
                  cut3l =q1 - (iqr * 4),
                  cut4h = q3 + (iqr * 5),
                  cut4l =q1 - (iqr * 5),
                  cut5l = q1 - (iqr * 7),
                  cut5h = q3 + (iqr * 7),
                  cut6l = q1 - (iqr * 10),
                  cut6h = q3 + (iqr * 10)) %>%
    # Join the bin boundaries back to the original data frame
    dplyr::left_join(data, ., by=c("trait")) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    # Add columns tallying the total number of points in each of the bins
    dplyr::mutate(onehs = ifelse( cut2h > phenotype & phenotype >= cut1h,
                                  1, 0),
                  onels = ifelse( cut2l < phenotype & phenotype <= cut1l,
                                  1, 0),
                  twohs = ifelse( cut3h > phenotype & phenotype >= cut2h,
                                  1, 0),
                  twols = ifelse( cut3l < phenotype & phenotype <= cut2l,
                                  1, 0),
                  threehs = ifelse(cut4h > phenotype & phenotype >= cut3h,
                                   1, 0),
                  threels = ifelse(cut4l < phenotype & phenotype <= cut3l,
                                   1, 0),
                  fourhs = ifelse(cut5h > phenotype &  phenotype >= cut4h,
                                  1, 0),
                  fourls = ifelse(cut5l < phenotype &  phenotype <= cut4l,
                                  1, 0),
                  fivehs = ifelse(cut6h > phenotype & phenotype >= cut5h,
                                  1, 0),
                  fivels = ifelse(cut6l < phenotype & phenotype <= cut5l,
                                  1, 0),
                  sixhs = ifelse(phenotype >= cut6h, 1, 0),
                  sixls = ifelse(phenotype <= cut6l, 1, 0)) %>%
    # Group on condition and trait, then sum the total number of data points
    # in each of the IQR multiple bins
    dplyr::group_by(trait) %>%
    dplyr::mutate(s1h = sum(onehs, na.rm = TRUE),
                  s2h = sum(twohs, na.rm = TRUE),
                  s3h = sum(threehs, na.rm = TRUE),
                  s4h = sum(fourhs, na.rm = TRUE),
                  s5h = sum(fivehs, na.rm = TRUE),
                  s1l = sum(onels, na.rm = TRUE),
                  s2l = sum(twols, na.rm = TRUE),
                  s3l = sum(threels, na.rm = TRUE),
                  s4l = sum(fourls, na.rm = TRUE),
                  s5l = sum(fivels, na.rm = TRUE),
                  s6h = sum(sixhs, na.rm = TRUE),
                  s6l = sum(sixls, na.rm = TRUE))%>%
    # Group on condition and trait, then check to see if the number of
    # points in each bin is more than 5% of the total number of data points
    dplyr::group_by(trait) %>%
    dplyr::mutate(p1h = ifelse(sum(onehs, na.rm = TRUE) / n() >= .05,1,0),
                  p2h = ifelse(sum(twohs, na.rm = TRUE) / n() >= .05,1,0),
                  p3h = ifelse(sum(threehs, na.rm = TRUE) / n() >= .05,1,0),
                  p4h = ifelse(sum(fourhs, na.rm = TRUE) / n() >= .05,1,0),
                  p5h = ifelse(sum(fivehs, na.rm = TRUE) / n() >= .05,1,0),
                  p6h = ifelse(sum(sixhs, na.rm = TRUE) / n() >= .05,1,0),
                  p1l = ifelse(sum(onels, na.rm = TRUE) / n() >= .05,1,0),
                  p2l = ifelse(sum(twols, na.rm = TRUE) / n() >= .05,1,0),
                  p3l = ifelse(sum(threels, na.rm = TRUE) / n() >= .05,1,0),
                  p4l = ifelse(sum(fourls, na.rm = TRUE) / n() >= .05,1,0),
                  p5l = ifelse(sum(fivels, na.rm = TRUE) / n() >= .05,1,0),
                  p6l = ifelse(sum(sixls,
                                   na.rm = TRUE) / n() >= .05,1,0)) %>%
    # Count the number of observations in each condition/trait combination
    dplyr::mutate(numst = n()) %>%
    # Group on trait, then filter out NAs in any of the added
    # columns
    dplyr::group_by(trait) %>%
    dplyr::filter(!is.na(trait), !is.na(phenotype), !is.na(iqr), !is.na(q1),
                  !is.na(q3), !is.na(cut1h), !is.na(cut1l), !is.na(cut2h),
                  !is.na(cut2l), !is.na(cut3h), !is.na(cut3l),
                  !is.na(cut4h), !is.na(cut4l), !is.na(cut5l),
                  !is.na(cut5h), !is.na(cut6l), !is.na(cut6h),
                  !is.na(onehs), !is.na(onels), !is.na(twohs),
                  !is.na(twols), !is.na(threehs), !is.na(threels),
                  !is.na(fourhs), !is.na(fourls), !is.na(fivehs),
                  !is.na(fivels), !is.na(sixhs), !is.na(sixls),
                  !is.na(s1h), !is.na(s2h), !is.na(s3h), !is.na(s4h),
                  !is.na(s5h), !is.na(s1l), !is.na(s2l), !is.na(s3l),
                  !is.na(s4l), !is.na(s5l), !is.na(s6h), !is.na(s6l),
                  !is.na(p1h), !is.na(p2h), !is.na(p3h), !is.na(p4h),
                  !is.na(p5h), !is.na(p6h), !is.na(p1l), !is.na(p2l),
                  !is.na(p3l), !is.na(p4l), !is.na(p5l), !is.na(p6l),
                  !is.na(numst)) %>%
    # Add three columns stating whether the observation is an outlier
    # based the three outlier detection functions below
    dplyr::ungroup() %>%
    dplyr::mutate(cuts = categorize1(.),
                  cuts1 = categorize2(.),
                  cuts2 = categorize3(.))

  if( remove_outliers == T){
    outliers_removed <- datawithoutliers %>%
      dplyr::rename(bamfoutlier1 = cuts,
                    bamfoutlier2 = cuts1,
                    bamfoutlier3 = cuts2)%>%
      dplyr::filter(!bamfoutlier1 & !bamfoutlier2 & !bamfoutlier3)%>%
      dplyr::select(trait, strain, phenotype) %>%
      tidyr::spread(trait, phenotype)

    return(outliers_removed)
  } else {
      with_outliers <- datawithoutliers %>%
          dplyr::mutate(bamfoutlier1 = dplyr::case_when( cuts == TRUE ~ trait,
                                                         cuts == FALSE ~ "FALSE"),
                        bamfoutlier2 = dplyr::case_when( cuts1 == TRUE ~ trait,
                                                         cuts1 == FALSE ~ "FALSE" ),
                        bamfoutlier3 = dplyr::case_when( cuts2 == TRUE ~ trait,
                                                         cuts2 == FALSE ~ "FALSE") ) %>%
          dplyr::select(trait, strain, phenotype, bamfoutlier1, bamfoutlier2, bamfoutlier3) %>%
          tidyr::spread(trait, phenotype)

    return(with_outliers)
  }
}


#' Process phenotypes for mapping
#'
#' \code{process_phenotypes} takes input raw phenotype data and converts strain names to
#' isotype names, summarizes replicate data, and eliminates outliers.
#'
#' @param df a dataframe with a strain column and columns for each trait for processing.
#' @param summarize_replicates summarization method, currently limited to "mean" or "median"
#' @param prune_method method for eliminating outliers, currently limited to "BAMF"
#' @param remove_outliers boolean to specify if outliers should be eliminated within the function.
#' If FALSE, additional columns will be output specifying if the strain phenotype is an outier.
#' @return Output is a dataframe with
#' @seealso \link{generate_kinship} \link{generate_mapping}
#' @export
#'
process_phenotypes <- function(df,
                               summarize_replicates = "mean",
                               prune_method = "BAMF",
                               remove_outliers = TRUE){

  if( sum(grepl(colnames(df)[1], "Strain", ignore.case = T)) == 0 ) {
    message(glue::glue("~ ~ ~ WARNING ~ ~ ~
                       \nCheck input data format, strain should be the first column.
                       \n~ ~ ~ WARNING ~ ~ ~"))
  }

  # ~ ~ ~ # resolve strain isotypes # ~ ~ ~ #
  # get strain isotypes
  strain_isotypes_db <- generate_isotype_lookup()
  # identify strains that were phenotyped, but are not part of an isotype
  non_isotype_strains <- dplyr::filter(df,
                                       !(strain %in% strain_isotypes_db$strain),
                                       !(strain %in% strain_isotypes_db$isotype))
  # remove any strains identified to not fall into an isotype
  if( nrow(non_isotype_strains) > 0 ) {

    strains_to_remove <- unique(non_isotype_strains$strain)

    message(glue::glue("~ ~ ~ WARNING ~ ~ ~
                       \nRemoving strain(s) {strains_to_remove} because they do not fall into a defined isotype.
                       \n~ ~ ~ WARNING ~ ~ ~"))

    df_non_isotypes_removed <- dplyr::filter( df, !(strain %in% strains_to_remove) )
  } else {
    df_non_isotypes_removed <- df
  }

  # resolve isotypes
  df_isotypes_resolved <- df_non_isotypes_removed %>%
    dplyr::mutate( isotype = resolve_isotypes( strain, strain_isotypes_db ) ) %>%
    tidyr::gather( trait, phenotype, -strain, -isotype ) %>%
    dplyr::filter( !is.na(phenotype) )

  # ~ ~ ~ # summarize replicate data # ~ ~ ~ #

  df_replicates_summarized <- df_isotypes_resolved %>%
    dplyr::group_by( isotype, trait ) %>%
    { if (summarize_replicates == "mean") dplyr::summarise(., phenotype = mean( phenotype, na.rm = T ) )
      else if (summarize_replicates == "median") dplyr::summarise(., phenotype = median( phenotype, na.rm = T ) )
      else  message(glue::glue("~ ~ ~ WARNING ~ ~ ~
                       \nPlease choose mean or median as options for summarizeing replicate data.
                                  \n~ ~ ~ WARNING ~ ~ ~")) } %>%
    dplyr::rename(strain = isotype)

  # included for testing to make sure BAMF_prune is removing something
  # df_replicates_summarized_with_outs <- dplyr::mutate(df_replicates_summarized, new_pheno = ifelse(strain == "N2", 1000, phenotype))%>%
  #   dplyr::select(-phenotype)%>%
  #   dplyr::select(strain, trait, phenotype = new_pheno)

  # ~ ~ ~ # Perform outlier removal # ~ ~ ~ #
  if( prune_method == "BAMF" ){
    if(remove_outliers == TRUE ){
      df_BAMF <- BAMF_prune(df_replicates_summarized, remove_outliers = T) %>%
          tidyr::spread(trait, phenotype)
    } else {
      df_BAMF <- BAMF_prune(df_replicates_summarized, remove_outliers = F) %>%
          tidyr::spread(trait, phenotype)
      df_BAMF <- BAMF_prune(data = df_replicates_summarized_with_outs, remove_outliers = F) %>%
          tidyr::spread(trait, phenotype)
    }
  }
}






