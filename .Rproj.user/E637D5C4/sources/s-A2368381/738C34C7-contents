#' is.available.genome
#'
#' @description A function that allows to check whether genome sequence and transcript model data can be downloaded
#'
#' @return A boolean value
#'
#' @noRd
is.available.genome <- function (db = "refseq", organism, details = FALSE) 
{
  if (is.element(db, c("refseq", "genbank"))) {
    if (file.exists(file.path(tempdir(), paste0("AssemblyFilesAllKingdoms_", 
                                                db, ".txt")))) {
      suppressWarnings(AssemblyFilesAllKingdoms <- readr::read_tsv(file.path(tempdir(), 
                                                                             paste0("AssemblyFilesAllKingdoms_", db, ".txt")), 
                                                                   col_names = c("assembly_accession", "bioproject", 
                                                                                 "biosample", "wgs_master", "refseq_category", 
                                                                                 "taxid", "species_taxid", "organism_name", 
                                                                                 "infraspecific_name", "isolate", "version_status", 
                                                                                 "assembly_level", "release_type", "genome_rep", 
                                                                                 "seq_rel_date", "asm_name", "submitter", "gbrs_paired_asm", 
                                                                                 "paired_asm_comp", "ftp_path", "excluded_from_refseq"), 
                                                                   comment = "#", col_types = readr::cols(assembly_accession = readr::col_character(), 
                                                                                                          bioproject = readr::col_character(), biosample = readr::col_character(), 
                                                                                                          wgs_master = readr::col_character(), refseq_category = readr::col_character(), 
                                                                                                          taxid = readr::col_integer(), species_taxid = readr::col_integer(), 
                                                                                                          organism_name = readr::col_character(), infraspecific_name = readr::col_character(), 
                                                                                                          isolate = readr::col_character(), version_status = readr::col_character(), 
                                                                                                          assembly_level = readr::col_character(), release_type = readr::col_character(), 
                                                                                                          genome_rep = readr::col_character(), seq_rel_date = readr::col_date(), 
                                                                                                          asm_name = readr::col_character(), submitter = readr::col_character(), 
                                                                                                          gbrs_paired_asm = readr::col_character(), 
                                                                                                          paired_asm_comp = readr::col_character(), 
                                                                                                          ftp_path = readr::col_character(), excluded_from_refseq = readr::col_character())))
    }
    else {
      kgdoms <- biomartr::getKingdoms(db = db)(db = db)
      storeAssemblyFiles <- vector("list", length(kgdoms))
      for (i in seq_along(kgdoms)) {
        storeAssemblyFiles[i] <- list(getSummaryFile(db = db, 
                                                     kingdom = kgdoms[i]))
      }
      AssemblyFilesAllKingdoms <- dplyr::bind_rows(storeAssemblyFiles)
      readr::write_tsv(AssemblyFilesAllKingdoms, file.path(tempdir(), 
                                                           paste0("AssemblyFilesAllKingdoms_", db, ".txt")))
    }
    organism_name <- assembly_accession <- taxid <- NULL
    orgs <- stringr::str_replace_all(AssemblyFilesAllKingdoms$organism_name, 
                                     "\\(", "")
    orgs <- stringr::str_replace_all(orgs, "\\)", "")
    AssemblyFilesAllKingdoms <- dplyr::mutate(AssemblyFilesAllKingdoms, 
                                              organism_name = orgs)
    organism <- stringr::str_replace_all(organism, "\\(", 
                                         "")
    organism <- stringr::str_replace_all(organism, "\\)", 
                                         "")
    if (!biomartr:::is.taxid(organism)) {
      FoundOrganism <- dplyr::filter(AssemblyFilesAllKingdoms, 
                                     stringr::str_detect(organism_name, organism) | 
                                       assembly_accession == organism)
    }
    else {
      FoundOrganism <- dplyr::filter(AssemblyFilesAllKingdoms, 
                                     taxid == as.integer(organism))
    }
    if (nrow(FoundOrganism) == 0) {
      message("Unfortunatey, no entry for '", organism, 
              "' was found in the '", db, "' database. ", 
              "Please consider specifying ", paste0("'db = ", 
                                                    dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                     "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
              "' to check whether '", organism, "' is available in these databases.")
      return(FALSE)
    }
    if (nrow(FoundOrganism) > 0) {
      if (!details) {
        if (all(FoundOrganism$refseq_category == "na")) {
          message("Only a non-reference genome assembly is available for '", 
                  organism, "'.", " Please make sure to specify the argument 'reference = FALSE' when running any get*() function.")
        }
        else {
          message("A reference or representative genome assembly is available for '", 
                  organism, "'.")
        }
        if (nrow(FoundOrganism) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
      if (details) {
        if (all(FoundOrganism$refseq_category == "na")) {
          message("Only a non-reference genome assembly is available for '", 
                  organism, "'.", " Please make sure to specify the argument 'reference = FALSE' when running any get*() function.")
        }
        return(FoundOrganism)
      }
    }
  }
  if (db == "ensembl") {
    name <- accession <- accession <- assembly <- taxon_id <- NULL
    new.organism <- stringr::str_replace_all(organism, " ", 
                                             "_")
    ensembl.available.organisms <- biomartr:::get.ensembl.info()
    ensembl.available.organisms <- dplyr::filter(ensembl.available.organisms, 
                                                 !is.na(assembly))
    if (!biomartr:::is.taxid(organism)) {
      selected.organism <- dplyr::filter(ensembl.available.organisms, 
                                         stringr::str_detect(name, stringr::str_to_lower(new.organism)) | 
                                           accession == organism, !is.na(assembly))
    }
    else {
      selected.organism <- dplyr::filter(ensembl.available.organisms, 
                                         taxon_id == as.integer(organism), !is.na(assembly))
    }
    if (!details) {
      if (nrow(selected.organism) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      }
      if (nrow(selected.organism) > 0) {
        message("A reference or representative genome assembly is available for '", 
                organism, "'.")
        if (nrow(selected.organism) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
    }
    if (details) 
      return(selected.organism)
  }
  if (db == "ensemblgenomes") {
    new.organism <- stringr::str_replace_all(organism, " ", 
                                             "_")
    name <- accession <- assembly <- taxon_id <- NULL
    selected.organism <- get.ensemblgenome.info(new.organism)
    if (!details) {
      if (length(selected.organism) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      } else {
        return(TRUE)
      }
    }
    if (details) 
      return(selected.organism)
  }
  if (db == "uniprot") {
    if (biomartr:::is.taxid(organism)) {
      unipreot_rest_url <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&taxid=", 
                                  as.integer(organism))
      rest_status_test <- curl::curl_fetch_memory(unipreot_rest_url)
      if (rest_status_test$status_code != 200) {
        message("Something went wrong when trying to access the API 'https://www.ebi.ac.uk/proteins/api/proteomes'.", 
                " Sometimes the internet connection isn't stable and re-running the function might help. Otherwise, could there be an issue with the firewall? ", 
                "Is it possbile to access the homepage 'https://www.ebi.ac.uk/' through your browser?")
      }
      uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url))
    }
    else {
      organism_new <- stringr::str_replace_all(organism, 
                                               " ", "%20")
      unipreot_rest_url_name <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&name=", 
                                       organism_new)
      rest_status_test_name <- curl::curl_fetch_memory(unipreot_rest_url_name)
      unipreot_rest_url_upid <- paste0("https://www.ebi.ac.uk/proteins/api/proteomes?offset=0&size=-1&upid=", 
                                       organism)
      rest_status_test_upid <- curl::curl_fetch_memory(unipreot_rest_url_upid)
      if ((rest_status_test_upid$status_code != 200) & 
          (rest_status_test_name$status_code != 200)) {
        message("Something went wrong when trying to access the API 'https://www.ebi.ac.uk/proteins/api/proteomes'.", 
                " Sometimes the internet connection isn't stable and re-running the function might help. Otherwise, could there be an issue with the firewall? ", 
                "Is it possbile to access the homepage 'https://www.ebi.ac.uk/' through your browser?")
      }
      if (rest_status_test_name$status_code == 200) {
        uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url_name))
      }
      if (rest_status_test_upid$status_code == 200) {
        uniprot_species_info <- tibble::as_tibble(jsonlite::fromJSON(unipreot_rest_url_upid))
      }
    }
    if (!details) {
      if (nrow(uniprot_species_info) == 0) {
        message("Unfortunatey, no entry for '", organism, 
                "' was found in the '", db, "' database. ", 
                "Please consider specifying ", paste0("'db = ", 
                                                      dplyr::setdiff(c("refseq", "genbank", "ensembl", 
                                                                       "ensemblgenomes", "uniprot"), db), collapse = "' or "), 
                "' to check whether '", organism, "' is available in these databases.")
        return(FALSE)
      }
      if (nrow(uniprot_species_info) > 0) {
        message("A reference or representative genome assembly is available for '", 
                organism, "'.")
        if (nrow(uniprot_species_info) > 1) {
          message("More than one entry was found for '", 
                  organism, "'.", " Please consider to run the function 'is.available.genome()' and specify 'is.available.genome(organism = ", 
                  organism, ", db = ", db, ", details = TRUE)'.", 
                  " This will allow you to select the 'assembly_accession' identifier that can then be ", 
                  "specified in all get*() functions.")
        }
        return(TRUE)
      }
    }
  }
  if (details) 
    return(uniprot_species_info)
}

#' getsequence
#'
#' @description A function that allows download the genome sequence of an organism
#'
#' @return A character string giving the location of the file encoding the genome sequence downloaded
#'
#' @noRd
getsequence <- function (organism, release = NULL, type = "dna", id.type = "toplevel"){
  if (!is.element(type, c("dna", "cds", "pep", "ncrna"))) 
    stop("Please a 'type' argument supported by this function: \n                 'dna', 'cds', 'pep', 'ncrna'.")
  name <- NULL
  if (!suppressMessages(is.available.genome(organism = organism, 
                                            db = "ensemblgenomes", details = FALSE))) {
    warning("Unfortunately organism '", organism, "' is not available at ENSEMBLGENOMES. ", 
            "Please check whether or not the organism name is typed correctly or try db = 'ensembl'.", 
            " Thus, download of this species has been omitted. ", 
            call. = FALSE)
    return(FALSE)
  }
  else {
    taxon_id <- assembly <- accession <- NULL
    new.organism <- stringr::str_to_lower(stringr::str_replace_all(organism, 
                                                                   " ", "_"))
    ensembl_summary <- suppressMessages(as.data.frame(unlist(is.available.genome(organism = organism, 
                                                            db = "ensemblgenomes", details = TRUE))))
    if (nrow(ensembl_summary) == 0) {
      message("Unfortunately, organism '", organism, "' does not exist in this database. Could it be that the organism name is misspelled? Thus, download has been omitted.")
      return(FALSE)
    }
    new.organism <- paste0(stringr::str_to_upper(stringr::str_sub(ensembl_summary$name[1], 
                                                                  1, 1)), stringr::str_sub(ensembl_summary$name[1], 
                                                                                           2, nchar(ensembl_summary$name[1])))
  }
  get.org.info <- ensembl_summary
  rest_url <- paste0("http://rest.ensembl.org/info/assembly/", 
                     ensembl_summary["name",], "?content-type=application/json")
  rest_api_status <- biomartr:::test_url_status(url = rest_url, organism = organism)
  if (is.logical(rest_api_status)) {
    return(FALSE)
  }
  else {
      release_api <- jsonlite::fromJSON("http://rest.ensembl.org/info/eg_version?content-type=application/json")
      if (!is.null(release)) {
        if (!is.element(release, seq_len(as.integer(release_api)))) 
          stop("Please provide a release number that is supported by ENSEMBLGENOMES.", 
               call. = FALSE)
      }
      if (is.null(release)) 
        core_path <- "http://ftp.ensemblgenomes.org/pub/current/"
      if (!is.null(release)) 
        core_path <- paste0("http://ftp.ensemblgenomes.org/pub/current", 
                            release, "/")
      ensembl.qry <- paste0(core_path, stringr::str_to_lower(stringr::str_replace(get.org.info$division[1], 
                                                                                  "Ensembl", "")), "plants/fasta/",
                            stringr::str_to_lower(ensembl_summary["name",]), 
                            "/", type, "/", paste0(ensembl_summary["url_name",], ".", rest_api_status$default_coord_system_version, 
                                                   ".", "dna_rm", ifelse(id.type == "none", "", "."), 
                                                   ifelse(id.type == "none", "", id.type), ".fa.gz"))

    if (file.exists(file.path(path, paste0(new.organism, 
                                           ".", rest_api_status$default_coord_system_version, 
                                           ".", type, ifelse(id.type == "none", "", "."), ifelse(id.type == 
                                                                                                 "none", "", id.type), ".fa.gz")))) {
      message("File ", file.path(path, paste0(new.organism, 
                                              ".", rest_api_status$default_coord_system_version, 
                                              ".", type, ifelse(id.type == "none", "", "."), 
                                              ifelse(id.type == "none", "", id.type), ".fa.gz")), 
              " exists already. Thus, download has been skipped.")
    }
    else {
      path <- getwd()
      file.fasta <- file.path(path, 
                              paste0(ensembl_summary["name",], ".", rest_api_status$default_coord_system_version, 
                              ".", type, ifelse(id.type == "none", "", "."), 
                                          ifelse(id.type == "none", "", id.type), ".fa.gz"))
      utils::download.file(url = ensembl.qry, destfile = file.fasta)
    }
    return(file.fasta)
  }
}

