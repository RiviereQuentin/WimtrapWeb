genomic_data <- c(DHS = "data/Arabidopsis_thaliana/seedlings/DHS_athal_seedlings.bed",
                  DGF = "data/Arabidopsis_thaliana/seedlings/DGF_athal_seedlings.bed",
                  Dloops = "data/Arabidopsis_thaliana/seedlings/Dloops_athal_seedlings.bed",
                  H2AZ = "data/Arabidopsis_thaliana/seedlings/H2AZ_athal_seedlings.bed",
                  H2BuB = "data/Arabidopsis_thaliana/seedlings/H2BUB_athal_seedlings.bed",
                  H3K4me1 = "data/Arabidopsis_thaliana/seedlings/H3K4me1_athal_seedlings.bed",
                  H3K4me2 = "data/Arabidopsis_thaliana/seedlings/H3K4me2_athal_seedlings.bed",
                  H3K9me2 = "data/Arabidopsis_thaliana/seedlings/H3K9me2_athal_seedlings.bed",
                  H3K14ac = "data/Arabidopsis_thaliana/seedlings/H3K14ac_athal_seedlings.bed",
                  H3K18ac = "data/Arabidopsis_thaliana/seedlings/H3K18ac_athal_seedlings.bed",
                  H3K27ac = "data/Arabidopsis_thaliana/seedlings/H3K27ac_athal_seedlings.bed",
                  H3K27me1 = "data/Arabidopsis_thaliana/seedlings/H3K27me1_athal_seedlings.bed",
                  H3K56ac = "data/Arabidopsis_thaliana/seedlings/H3K56ac_athal_seedlings.bed",
                  H4K5ac = "data/Arabidopsis_thaliana/seedlings/H4K5ac_athal_seedlings.bed",
                  H4K8ac = "data/Arabidopsis_thaliana/seedlings/H4K8ac_athal_seedlings.bed",
                  H4K12ac = "data/Arabidopsis_thaliana/seedlings/H4K12ac_athal_seedlings.bed",
                  H4K16ac = "data/Arabidopsis_thaliana/seedlings/H4K16ac_athal_seedlings.bed",
                  Methylome = "data/Arabidopsis_thaliana/seedlings/Methylome_athal_seedlings.bed",
                  phast1 = "data/Arabidopsis_thaliana/phastcons_athal_part1.bed",
                  phast2 = "data/Arabidopsis_thaliana/phastcons_athal_part2.bed",
                  CDS = "data/Arabidopsis_thaliana/CDS_athal.bed",
                  Intron = "data/Arabidopsis_thaliana/Intron_athal.bed",
                  X3UTR = "data/Arabidopsis_thaliana/X3UTR_athal.bed",
                  X5UTR = "data/Arabidopsis_thaliana/X5UTR_athal.bed",
                  CNS = "data/Arabidopsis_thaliana/CNS_athal.bed")


imported_genomic_data.arabidopsis.seedlings <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                                                          tts = "data/Arabidopsis_thaliana/TTS_athal.bed",
                                                                          tss = "data/Arabidopsis_thaliana/TSS_athal.bed")

imported_genomic_data.arabidopsis.seedlings$phast1 <- imported_genomic_data.arabidopsis.seedlings[grep(pattern = "phast", names(imported_genomic_data.arabidopsis.seedlings))]
imported_genomic_data.arabidopsis.seedlings$phast1 <- lapply(imported_genomic_data.arabidopsis.seedlings$phast1,
                                                             function(x){x <- as.data.frame(x);
                                                                         colnames(x)[ncol(x)] <- "phastcons";
                                                                         return(x)})
imported_genomic_data.arabidopsis.seedlings$phast1 <- do.call(rbind, imported_genomic_data.arabidopsis.seedlings$phast1)
imported_genomic_data.arabidopsis.seedlings$phast1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data.arabidopsis.seedlings$phast1, keep.extra.columns = TRUE)
imported_genomic_data.arabidopsis.seedlings <- imported_genomic_data.arabidopsis.seedlings[!(names(imported_genomic_data.arabidopsis.seedlings)=="phast2")]
names(imported_genomic_data.arabidopsis.seedlings)[which(names(imported_genomic_data.arabidopsis.seedlings)=="phast1")] <- "phastcons"

ChIPpeaks <- c(PIF3 ="data/Arabidopsis_thaliana/seedlings/PIF3_athal_seedlings.bed",
               ABF3 = "data/Arabidopsis_thaliana/seedlings/ABF3_athal_seedlings.narrowPeak",
               PIF5 = "data/Arabidopsis_thaliana/seedlings/PIF5_athal_seedlings.narrowPeak",
               NAC52 = "data/Arabidopsis_thaliana/seedlings/SGS1_athal_seedlings.narrowPeak",
               NAC50 = "data/Arabidopsis_thaliana/seedlings/NAC50_athal_seedlings.narrowPeak",
               CCA1 = "data/Arabidopsis_thaliana/seedlings/CCA1_athal_seedlings.narrowPeak",
               HB7 = "data/Arabidopsis_thaliana/seedlings/ATHB7_athal_seedlings.narrowPeak",
               PIF4 = "data/Arabidopsis_thaliana/seedlings/PIF4_athal_seedlings.narrowPeak",
               FBH3 = "data/Arabidopsis_thaliana/seedlings/FBH3_athal_seedlings.narrowPeak",
               PRR5 = "data/Arabidopsis_thaliana/seedlings/PRR5_athal_seedlings.bed",
               TOC1 = "data/Arabidopsis_thaliana/seedlings/TOC1_athal_seedlings.bed",
               PRR7 = "data/Arabidopsis_thaliana/seedlings/PRR7_athal_seedlings.bed",
               ABF1="data/Arabidopsis_thaliana/seedlings/ABF1_athal_seedlings.narrowPeak",
               GBF2="data/Arabidopsis_thaliana/seedlings/GBF2_athal_seedlings.narrowPeak",
               GBF3="data/Arabidopsis_thaliana/seedlings/GBF3_athal_seedlings.narrowPeak",
               HAT22="data/Arabidopsis_thaliana/seedlings/HAT22_athal_seedlings.narrowPeak",
               HB5="data/Arabidopsis_thaliana/seedlings/HB5_athal_seedlings.narrowPeak",
               HB6="data/Arabidopsis_thaliana/seedlings/HB6_athal_seedlings.narrowPeak",
               HBI1="data/Arabidopsis_thaliana/seedlings/HBI1_athal_seedlings.narrowPeak",
               LHY="data/Arabidopsis_thaliana/seedlings/LHY_athal_seedlings.narrowPeak",
               MYB3="data/Arabidopsis_thaliana/seedlings/MYB3_athal_seedlings.narrowPeak",
               MYB44="data/Arabidopsis_thaliana/seedlings/MYB44_athal_seedlings.narrowPeak",
               PIF1="data/Arabidopsis_thaliana/seedlings/PIF1_athal_seedlings.narrowPeak",
               WRKY18="data/Arabidopsis_thaliana/seedlings/WRKY18_athal_seedlings.narrowPeak",
               WRKY33="data/Arabidopsis_thaliana/seedlings/WRKY33_athal_seedlings.narrowPeak",
               WRKY40="data/Arabidopsis_thaliana/seedlings/WRKY40_athal_seedlings.narrowPeak",
               RD26="data/Arabidopsis_thaliana/seedlings/RD26_athal_seedlings.narrowPeak",
               IBH1 = "data/Arabidopsis_thaliana/seedlings/IBH1_athal_seedlings.narrowPeak",
               ABF4 = "data/Arabidopsis_thaliana/seedlings/ABF4_athal_seedlings.narrowPeak",
               REVOLUTA = "data/Arabidopsis_thaliana/seedlings/REVOLUTA_athal_seedlings.narrowPeak"
)
TFBSdata.athal.seedlings <- Wimtrap::getTFBSdata(pfm = "data/Arabidopsis_thaliana/PFMs_athal.pfm",
                                             TFnames = names(ChIPpeaks),
                                             genome_sequence = paste0("data/Arabidopsis_thaliana/genome_athal_chr", seq(1,5), ".fa.gzip"),
                                             imported_genomic_data = imported_genomic_data.arabidopsis.seedlings)
rm(imported_genomic_data.arabidopsis.seedlings)

TFBSmodel.athal.seedlings_full <- Wimtrap:::buildTFBSmodel(TFBSdata.athal.seedlings,
                                                      ChIPpeaks = ChIPpeaks,
                                                      model_assessment = FALSE)

save(TFBSmodel.athal.seedlings_full, file = "models/TFBSmodel_athal_seedlings_full.RData")

gd <- c("seqnames", "start", "end", "width", "strand", "matchScore", "matchLogPval", "ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
        "X3UTR", "Downstream", "ClosestTTS", "DistToClosestTTS", "ClosestTSS", "DistToClosestTSS", "Matches_400bp", "Matches_1000bp",
        paste0(rep(c("phastcons", "DGF", "DHS", "CNS"), each = 3), rep(c("_20bp", "_400bp", "_1000bp"), 4)))
for (i in seq_along(TFBSdata.athal.seedlings)){
  TFdata <- data.table::fread(TFBSdata.athal.seedlings[i])
  TFdata <- TFdata[,colnames(TFdata) %in% gd, with= FALSE]
  data.table::fwrite(TFdata, file = paste0("reduced_", TFBSdata.athal.seedlings[i]))
}

TFBSdata.athal.seedlings_reduced <- paste0("reduced_", TFBSdata.athal.seedlings)
names(TFBSdata.athal.seedlings_reduced) <- names(TFBSdata.athal.seedlings)
TFBSmodel.athal.seedlings_reduced <- Wimtrap:::buildTFBSmodel(TFBSdata.athal.seedlings_reduced,
                                                           ChIPpeaks = ChIPpeaks,
                                                           model_assessment = FALSE)
save(TFBSmodel.athal.seedlings_reduced, file = "models/TFBSmodel_athal_seedlings_reduced.RData")

gd <- c("seqnames", "start", "end", "width", "strand", "matchScore", "matchLogPval", "ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
        "X3UTR", "Downstream", "ClosestTTS", "DistToClosestTTS", "ClosestTSS", "DistToClosestTSS", "Matches_400bp", "Matches_1000bp",
        paste0(rep(c("DGF", "DHS"), each = 3), rep(c("_20bp", "_400bp", "_1000bp"), 2)))
for (i in seq_along(TFBSdata.athal.seedlings)){
  TFdata <- data.table::fread(TFBSdata.athal.seedlings[i])
  TFdata <- TFdata[,colnames(TFdata) %in% gd, with= FALSE]
  data.table::fwrite(TFdata, file = paste0("reduced_", TFBSdata.athal.seedlings[i]))
}

TFBSdata.athal.seedlings_reduced <- paste0("reduced", TFBSdata.athal.seedlings)
TFBSmodel.athal.seedlings_full <- Wimtrap:::buildTFBSmodel(TFBSdata.athal.seedlings_reduced,
                                                           ChIPpeaks = ChIPpeaks[names(TFBSdata.athal.seedlings_reduced)],
                                                           model_assessment = FALSE)
save(TFBSmodel.athal.seedlings_reduced, file = "models/TFBSmodel_general_DHSDGF.RData")

gd <- c("seqnames", "start", "end", "width", "strand", "matchScore", "matchLogPval", "ProximalPromoter", "Promoter", "X5UTR", "CDS", "Intron",
        "X3UTR", "Downstream", "ClosestTTS", "DistToClosestTTS", "ClosestTSS", "DistToClosestTSS", "Matches_400bp", "Matches_1000bp",
        paste0(rep(c("DHS"), each = 3), rep(c("_20bp", "_400bp", "_1000bp"), 2)))
for (i in seq_along(TFBSdata.athal.seedlings)){
  TFdata <- data.table::fread(TFBSdata.athal.seedlings[i])
  TFdata <- TFdata[,colnames(TFdata) %in% gd, with= FALSE]
  data.table::fwrite(TFdata, file = paste0("reduced_", TFBSdata.athal.seedlings[i]))
}

TFBSdata.athal.seedlings_reduced <- paste0("reduced", TFBSdata.athal.seedlings)
TFBSmodel.athal.seedlings_full <- Wimtrap:::buildTFBSmodel(TFBSdata.athal.seedlings_reduced,
                                                           ChIPpeaks = ChIPpeaks[names(TFBSdata.athal.seedlings_reduced)],
                                                           model_assessment = FALSE)
save(TFBSmodel.athal.seedlings_reduced, file = "models/TFBSmodel_general_DHS.RData")
