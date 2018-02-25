#' Get a dataframe from the CeNDR database
#'
#' \code{get_db} downloads and caches the CeNDR database and fetches requested tables.
#'
#' The database has five tables:
#'
#' \itemize{
#'      \item \strong{homologs} - A table of homologs and orthologs.
#'      \item \strong{metadata} - data about the data in this database.
#'      \item \strong{strain} - A replicate of the WI Strain Database. Contains information on their names, isotypes, isolation location, and more.
#'      \item \strong{wormbase_gene} - Database of the gff-based gene models.
#'      \item \strong{wormbase_gene_summary} - A summary of the genes in \emph{C. elegans}.
#' }
#'
#' @param table The table to request. [\strong{Default:} \code{"wormbase_gene"}]
#' @param renew Force download of the latest database. Then return the requested table. [\strong{Default:} \code{FALSE}]
#' @examples get_db()
#' @return Dataframe with the requested table.
#' @export


get_db <- function(table="wormbase_gene", renew=FALSE) {
    # Function for fetching the variant database
    base::dir.create("~/.cegwas", showWarnings=F)
    file_path <- "~/.cegwas/cegwas.db" # nolint

    table_list <-  c("homologs",
                     "metadata",
                     "wormbase_gene",
                     "wormbase_gene_summary",
                     "strain")

    assertthat::assert_that(table %in% table_list,
                            msg = glue::glue("{table} not an option. Must be one of:\n\n{paste0(table_list, collapse='\n')}\n"))

    if (file.info(file_path)$size < 128 | is.na(file.info(file_path)$size == 0) | renew) {
        message(paste0("Downloading Gene Database to ", file_path))
        url <- "https://storage.googleapis.com/elegansvariation.org/db/_latest.db"
        utils::download.file(url, file_path)
    }
    dplyr::tbl(dplyr::src_sqlite(file_path), table)
}
