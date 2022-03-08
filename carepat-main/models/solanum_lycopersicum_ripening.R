genomic_data <- paste0("/home/quentin/Documents/carepat/data/Solanum_lycopersicum/",
                       c("CE_sl.bed", "CDS_sl.bed", "Intron_sl.bed", "X5UTR_sl.bed",
                         "X3UTR_sl.bed", "ripeningfruits/DHS_sl_ripening.bed",
                         "ripeningfruits/DGF_sl_ripening.bed", "ripeningfruits/H3K27me3_sl_ripening.bed",
                         paste0("ripeningfruits/Methylome_sl_ripening_part", seq(1,14), ".bed")))
names(genomic_data) <- c("phastcons", "CDS", "Intron", "X5UTR", "X3UTR", "DHS", "DGF", "H3K27me3", paste0("Methylome_", seq(1,14)))

imported_genomic_data.sl.ripening <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                                                tts = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/TTS_sl.bed",
                                                                tss = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/TSS_sl.bed")
imported_genomic_data.sl.ripening$Methylome_1 <- imported_genomic_data.sl.ripening[grep(pattern = "Methylome", names(imported_genomic_data.sl.ripening))]
imported_genomic_data.sl.ripening$Methylome_1 <- lapply(imported_genomic_data.sl.ripening$Methylome_1, function(x){x <- data.table::as.data.table(x); colnames(x)[ncol(x)] <- "CpG"; return(x)})
imported_genomic_data.sl.ripening$Methylome_1 <- do.call(rbind, imported_genomic_data.sl.ripening$Methylome_1)
imported_genomic_data.sl.ripening$Methylome_1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data.sl.ripening$Methylome_1)
names(imported_genomic_data.sl.ripening)[which(names(imported_genomic_data.sl.ripening) == "Methylome_1")] <- "Cme"
imported_genomic_data.sl.ripening <- imported_genomic_data.sl.ripening[!(seq(1, length(imported_genomic_data.sl.ripening)) %in% grep(pattern = "Methylome", names(imported_genomic_data.sl.ripening)))]

genome_sequence <- Biostrings::readDNAStringSet(c("/home/quentin/Documents/carepat/data/Solanum_lycopersicum/genome_sl_chr0.fa.gz",
                                                  "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/genome_sl_chr1_part1.fa.gz",
                                                  "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/genome_sl_chr1_part2.fa.gz",
                                                  paste0("/home/quentin/Documents/carepat/data/Solanum_lycopersicum/genome_sl_chr", seq(2,12),".fa.gz")))
genome_sequence[2] <- Biostrings::xscat(genome_sequence[2], genome_sequence[3])
genome_sequence <- genome_sequence[c(1,2,seq(4,14))]
TFBSdata.sl.ripening <- Wimtrap::getTFBSdata(pfm = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/PFMs_sl.pfm",
                                             TFnames = c("RIN", "EIN3", "TAGL1", "FUL1", "FUL2"),
                                             genome_sequence = genome_sequence,
                                             imported_genomic_data = imported_genomic_data.sl.ripening)

TFBSmodel.sl.ripening_full <- Wimtrap::buildTFBSmodel(TFBSdata = TFBSdata.sl.ripening,
                                                 ChIPpeaks = c(RIN = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/RIN_sl_ripening.bed",
                                                               EIN3 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/EIN3_sl_ripening.bed",
                                                               TAGL1 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/TAGL1_sl_ripening.bed",
                                                               FUL1 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/FUL1_sl_ripening.bed",
                                                               FUL2 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/FUL2_sl_ripening.bed"))
save(TFBSmodel.sl.ripening_full, file="/home/quentin/Documents/carepat/models/TFBSmodel_sl_ripening_full.RData")

TFBSdata.sl.ripening <- c(FUL2 = "FUL2_47_18_annotations.tsv",
FUL1 = "FUL1_46_44_annotations.tsv",
TAGL1 = "TAGL1_46_11_annotations.tsv",
EIN3 = "EIN3_45_34_annotations.tsv",
RIN = "RIN_44_50_annotations.tsv")

TFBSdata.sl.ripening.reduced <- paste0("reduced_", TFBSdata.sl.ripening)
names(TFBSdata.sl.ripening.reduced) <- names(TFBSdata.sl.ripening)
for (TF in names(TFBSdata.sl.ripening)){
  data <- data.table::fread(TFBSdata.sl.ripening[TF])
  data <- data[,-grep("DGF", colnames(data)), with = FALSE]
  data.table::fwrite(data, paste0("reduced_", TFBSdata.sl.ripening[TF]))
}
TFBSmodel.sl.ripening_reduced <- Wimtrap::buildTFBSmodel(TFBSdata = TFBSdata.sl.ripening.reduced,
                                                 ChIPpeaks = c(RIN = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/RIN_sl_ripening.bed",
                                                               EIN3 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/EIN3_sl_ripening.bed",
                                                               TAGL1 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/TAGL1_sl_ripening.bed",
                                                               FUL1 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/FUL1_sl_ripening.bed",
                                                               FUL2 = "/home/quentin/Documents/carepat/data/Solanum_lycopersicum/ripeningfruits/FUL2_sl_ripening.bed"))
save(TFBSmodel.sl.ripening_reduced, file="/home/quentin/Documents/carepat/models/TFBSmodel_sl_ripening_reduced.RData")
