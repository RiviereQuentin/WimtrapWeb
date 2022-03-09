#' carepat 
#'
#' @description A function that makes use of prebuilt models
#'
#' @return Data table describing the predicted transcription factor binding sites
#'
#' @noRd
carepat <- function(organism = c("Arabidopsis thaliana", "Solanum lycopersicum"),
                    condition = c("seedlings", "flowers", "roots", "roots_non_hairs",
                                  "seed_coats", "seedlings_dark7d", "seedlings_dark7dLD24h",
                                  "seedlings_dark7dlight3h", "seedlings_dark7dlight30min", "seedlings_heatshock",
                                  "ripening_fruits", "immature_fruits"),
                    TFnames = NULL,
                    pfm = NULL,
                    show_annotations = FALSE,
                    score_threshold = 0.5){

  #Download the required data and models from the github
  #repository RiviereQuentin/caretrap if necessary
  package.dir <- getwd()
  dir.models <- paste0(package.dir, "/carepat-main/models/")
  if (organism == "Solanum lycopersicum"){
    dir.data <- paste0(package.dir, "/carepat-main/data/Solanum_lycopersicum/")
    if (condition == "ripening_fruits"){
      genomic_data <- c("CE_sl.bed", "CDS_sl.bed", "Intron_sl.bed", "X5UTR_sl.bed",
                               "X3UTR_sl.bed", "ripeningfruits/DHS_sl_ripening.bed",
                               "ripeningfruits/DGF_sl_ripening.bed", "ripeningfruits/H3K27me3_sl_ripening.bed",
                               paste0("ripeningfruits/Methylome_sl_ripening_part", seq(1,14), ".bed"))
      names(genomic_data) <- c("phastcons", "CDS", "Intron", "X5UTR", "X3UTR", "DHS", "DGF", "H3K27me3", paste0("Methylome_", seq(1,14)))
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      TFBSmodel <- paste0(dir.models, "TFBSmodel_sl_ripening_full.RData")
      load(TFBSmodel)
      TFBSmodel <- TFBSmodel.sl.ripening_full
    } else if (condition == "immature_fruits"){
      genomic_data <- c("CE_sl.bed", "CDS_sl.bed", "Intron_sl.bed", "X5UTR_sl.bed",
                               "X3UTR_sl.bed", "immaturefruits/DHS_sl_immature.bed",
                               "immaturefruits/H3K27me3_sl_immature.bed",
                               paste0("ripeningfruits/Methylome_sl_ripening_part", seq(1,14), ".bed"))
      names(genomic_data) <- c("phastcons", "CDS", "Intron", "X5UTR", "X3UTR", "DHS", "H3K27me3", paste0("Methylome_", seq(1,14)))
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      TFBSmodel <- paste0(dir.models, "TFBSmodel_sl_ripening_reduced.RData")
      load(TFBSmodel)
      TFBSmodel <- TFBSmodel.sl.ripening_reduced
    }

    imported_genomic_data <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                                        tts = paste0(dir.data, "TTS_sl.bed"),
                                                        tss = paste0(dir.data, "TSS_sl.bed"))
    imported_genomic_data$Methylome_1 <- imported_genomic_data[grep(pattern = "Methylome", names(imported_genomic_data))]
    imported_genomic_data$Methylome_1 <- lapply(imported_genomic_data$Methylome_1, function(x){x <- data.table::as.data.table(x); colnames(x)[ncol(x)] <- "CpG"; return(x)})
    imported_genomic_data$Methylome_1 <- do.call(rbind, imported_genomic_data$Methylome_1)
    imported_genomic_data$Methylome_1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data$Methylome_1)
    names(imported_genomic_data)[which(names(imported_genomic_data) == "Methylome_1")] <- "Cme"
    imported_genomic_data <- imported_genomic_data[!(seq(1, length(imported_genomic_data)) %in% grep(pattern = "Methylome", names(imported_genomic_data)))]
    genome_sequence <- Biostrings::readDNAStringSet(paste0(dir.data, c("genome_sl_chr0.fa.gz",
                                                                       "genome_sl_chr1_part1.fa.gz",
                                                                       "genome_sl_chr1_part2.fa.gz",
                                                                        paste0("genome_sl_chr", seq(2,12),".fa.gz"))))
    genome_sequence[2] <- Biostrings::xscat(genome_sequence[2], genome_sequence[3])
    genome_sequence <- genome_sequence[c(1,2,seq(4,14))]
    PWMs.file <- paste0(dir.data, "PWMs_sl.meme")

  } else if (organism == "Arabidopsis thaliana"){
    dir.data <- paste0(package.dir, "/carepat-main/data/Arabidopsis_thaliana/")
    if (condition == "seedlings"){
      genomic_data <- c(DHS = "seedlings/DHS_athal_seedlings.bed",
                        DGF = "seedlings/DGF_athal_seedlings.bed",
                        Dloops = "seedlings/Dloops_athal_seedlings.bed",
                        H2AZ = "seedlings/H2AZ_athal_seedlings.bed",
                        H2BuB = "seedlings/H2BUB_athal_seedlings.bed",
                        H3K4me1 = "seedlings/H3K4me1_athal_seedlings.bed",
                        H3K4me2 = "seedlings/H3K4me2_athal_seedlings.bed",
                        H3K9me2 = "seedlings/H3K9me2_athal_seedlings.bed",
                        H3K14ac = "seedlings/H3K14ac_athal_seedlings.bed",
                        H3K18ac = "seedlings/H3K18ac_athal_seedlings.bed",
                        H3K27ac = "seedlings/H3K27ac_athal_seedlings.bed",
                        H3K27me1 = "seedlings/H3K27me1_athal_seedlings.bed",
                        H3K56ac = "seedlings/H3K56ac_athal_seedlings.bed",
                        H4K5ac = "seedlings/H4K5ac_athal_seedlings.bed",
                        H4K8ac = "seedlings/H4K8ac_athal_seedlings.bed",
                        H4K12ac = "seedlings/H4K12ac_athal_seedlings.bed",
                        H4K16ac = "seedlings/H4K16ac_athal_seedlings.bed",
                        Methylome = "seedlings/Methylome_athal_seedlings.bed",
                        phast1 = "phastcons_athal_part1.bed",
                        phast2 = "phastcons_athal_part2.bed",
                        CDS = "CDS_athal.bed",
                        Intron = "Intron_athal.bed",
                        X3UTR = "X3UTR_athal.bed",
                        X5UTR = "X5UTR_athal.bed",
                        CNS = "CNS_athal.bed")
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      load(file = paste0(dir.models, "TFBSmodel_athal_seedlings_full.RData"))
      TFBSmodel <- TFBSmodel.athal.seedlings_full
    } else if (condition == "flowers"){
      genomic_data <- c(DHS = "flowers/DHS_athal_flowers.bed",
                        DGF = "flowers/DGF_athal_flowers.bed",
                        Methylome = "flowers/Methylome_athal_flowers.bed",
                        phast1 = "phastcons_athal_part1.bed",
                        phast2 = "phastcons_athal_part2.bed",
                        CDS = "CDS_athal.bed",
                        Intron = "Intron_athal.bed",
                        X3UTR = "X3UTR_athal.bed",
                        X5UTR = "X5UTR_athal.bed",
                        CNS = "CNS_athal.bed"
                      )
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp
      load(file = paste0(dir.models, "TFBSmodel_athal_flowers.RData"))
      TFBSmodel <- TFBSmodel.athal.flowers
    } else if (condition %in% c("roots", "roots_non_hairs", "seed_coats", "seedlings_dark7d", "seedlings_dark7dLD24h",
                                "seedlings_dark7dlight3h", "seedlings_dark7dlight30min", "seedlings_heatshock")){
      load(file = paste0(dir.models, "TFBSmodel_athal_seedlings_reduced.RData"))
      TFBSmodel <- TFBSmodel.athal.seedlings_reduced
      genomic_data <- c(phast1 = "phastcons_athal_part1.bed",
                                         phast2 = "phastcons_athal_part2.bed",
                                         CDS = "CDS_athal.bed",
                                         Intron = "Intron_athal.bed",
                                         X3UTR = "X3UTR_athal.bed",
                                         X5UTR = "X5UTR_athal.bed",
                                         CNS = "CNS_athal.bed")
      if (condition == "roots"){
        genomic_data <- c(genomic_data,
                          c(DGF = "roots/DGF_athal_roots.bed",
                                   DHS = "roots/DHS_athal_roots.bed"))
      } else if (condition == "roots_non_hairs"){
        genomic_data <- c(genomic_data,
                          c(DGF = "roots/DGF_athal_roots_non_hairs.bed",
                                   DHS = "roots/DHS_athal_roots_non_hairs.bed"))
      } else if (condition == "seed_coats"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seed_coats/DGF_athal_seed_coats.bed",
                                   DHS = "seed_coats/DHS_athaliana_seed_coats.bed"))
      } else if (condition == "seedlings_dark7d"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7d.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7d.bed"))
      } else if (condition == "seedlings_dark7dLD24h"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7dLD24h.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dLD24h.bed"))
      } else if (condition == "seedlings_dark7dlight3h"){
        genomic_data <- c(genomic_data,
                          c(DGF = "seedlings/DGF_athal_seedlings_dark7dlight3h.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dlight3h.bed"))
      } else if (condition == "seedlings_dark7dlight30min"){
        genomic_data <- c(genomic_data,
                         c(DGF = "seedlings/DGF_athal_seedlings_dark7dlight30min.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_dark7dlight30min.bed"))
      } else if (condition == "seedlings_heatshock"){
        genomic_data <- c(genomic_data,
                         c(DGF = "seedlings/DGF_athal_seedlings_heatshock.bed",
                                   DHS = "seedlings/DHS_athal_seedlings_heatshock.bed"))
      }
      tmp <- paste0(dir.data, genomic_data)
      names(tmp) <- names(genomic_data)
      genomic_data <- tmp

    }

    imported_genomic_data <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                               tts = paste0(dir.data, "TTS_athal.bed"),
                                               tss = paste0(dir.data, "TSS_athal.bed"))

    imported_genomic_data$phast1 <- imported_genomic_data[grep(pattern = "phast", names(imported_genomic_data))]
    imported_genomic_data$phast1 <- lapply(imported_genomic_data$phast1,
                                           function(x){x <- as.data.frame(x);
                                           colnames(x)[ncol(x)] <- "phastcons";
                                           return(x)})
    imported_genomic_data$phast1 <- do.call(rbind, imported_genomic_data$phast1)
    imported_genomic_data$phast1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data$phast1, keep.extra.columns = TRUE)
    imported_genomic_data <- imported_genomic_data[!(names(imported_genomic_data)=="phast2")]
    names(imported_genomic_data)[which(names(imported_genomic_data)=="phast1")] <- "phastcons"
    genome_sequence <- Biostrings::readDNAStringSet(paste0(dir.data, paste0("genome_athal_chr", seq(1,5), ".fa.gzip")))
    PWMs.file <- paste0(dir.data, "PWMs_athal.meme")

  }

  if(length(pfm) == 1){
    PWMs.file <- pfm
  }
  TFBSdata <- Wimtrap::getTFBSdata(pfm = PWMs.file,
                          TFnames = TFnames,
                          genome_sequence = genome_sequence,
                          imported_genomic_data = imported_genomic_data)
  TFBSpredictions <- Wimtrap::predictTFBS(TFBSmodel = TFBSmodel,
                                 TFBSdata = TFBSdata,
                                 studiedTFs = TFnames,
                                 score_threshold = score_threshold,
                                 show_annotations = show_annotations)
  return(TFBSpredictions)
}