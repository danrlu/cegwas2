df <- data.table::fread("~/AndersenLab/Github_Repos/cegwas2/inst/extdata/test_phenotype.tsv")

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


process_phenotypes <- function(df){
  
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
    
    df_non_isotypes_removed <- dplyr::filter(df, !(strain %in% strains_to_remove) )
  }
  
  # resolve isotypes 
  df_isotypes_resolved <- df_non_isotypes_removed %>%
    dplyr::mutate(isotype = resolve_isotypes( strain, strain_isotypes_db )) %>%
    tidyr::gather(trait, phenotype, -strain, -isotype) %>%
    dplyr::filter(!is.na(phenotype))
  
}









