ceGWAS <- R6Class("ceGWAS",
                  # use cegwas2 defaults, unless specified by user
                  private = list(
                      snvs = cegwas2::snp_set,
                      K = cegwas2::kinship_matrix
                  ),

                  public = list(
                      phenotype = NULL,
                      processed_phenotype = NULL,
                      mapping = NULL,
                      peak_intervals = NULL,
                      fine_mapping = NULL,

                      initialize = function(phenotype = NA,
                                            outlier_method = NA,
                                            remove_outliers = FALSE) {
                          self$phenotype = phenotype
                          self$processed_phenotype = process_phenotypes()
                      },

                      set_markers = function(genotype_matrix) {
                          private$snvs <- genotype_matrix
                      },

                      set_kinship = function() {
                          private$K <- generate_kinship(private$snvs)
                      },

                      # always perform mappings on objects processed phenotypes
                      # combine gwas_mapping and process_mapping from cegwas
                      perform_mapping = function(significance_threshold = "Bonferroni",
                                                 confidence_interval = "someMethod") {
                          self$mapping = gwas_mapping( phenotype = self$processed_phenotype,
                                                       genotype = private$snvs,
                                                       kinship = private$K,
                                                       threshold = significance_threshold,
                                                       CI = confidence_interval )

                          self$peak_intervals <- extract_mapping_peaks(self$mapping)
                      }

                      perform_fine_mapping = function(variant_severity = "ALL") {
                          self$fine_mapping = interval_fine_mapping(processed_phenotype,
                                                                    peak_intervals,
                                                                    variant_severity = variant_severity)
                      }
                  )
)
