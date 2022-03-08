genomic_data <- c(DHS = "data/Arabidopsis_thaliana/flowers/DHS_athal_flowers.bed",
                  DGF = "data/Arabidopsis_thaliana/flowers/DGF_athal_flowers.bed",
                  Methylome = "data/Arabidopsis_thaliana/flowers/Methylome_athal_flowers.bed",
                  phast1 = "data/Arabidopsis_thaliana/phastcons_athal_part1.bed",
                  phast2 = "data/Arabidopsis_thaliana/phastcons_athal_part2.bed",
                  CDS = "data/Arabidopsis_thaliana/CDS_athal.bed",
                  Intron = "data/Arabidopsis_thaliana/Intron_athal.bed",
                  X3UTR = "data/Arabidopsis_thaliana/X3UTR_athal.bed",
                  X5UTR = "data/Arabidopsis_thaliana/X5UTR_athal.bed",
                  CNS = "data/Arabidopsis_thaliana/CNS_athal.bed"
                  )


imported_genomic_data.arabidopsis.flowers <- Wimtrap::importGenomicData(genomic_data = genomic_data,
                                                                          tts = "data/Arabidopsis_thaliana/TTS_athal.bed",
                                                                          tss = "data/Arabidopsis_thaliana/TSS_athal.bed")

imported_genomic_data.arabidopsis.flowers$phast1 <- imported_genomic_data.arabidopsis.flowers[grep(pattern = "phast", names(imported_genomic_data.arabidopsis.flowers))]
imported_genomic_data.arabidopsis.flowers$phast1 <- lapply(imported_genomic_data.arabidopsis.flowers$phast1,
                                                             function(x){x <- as.data.frame(x);
                                                             colnames(x)[ncol(x)] <- "phastcons";
                                                             return(x)})
imported_genomic_data.arabidopsis.flowers$phast1 <- do.call(rbind, imported_genomic_data.arabidopsis.flowers$phast1)
imported_genomic_data.arabidopsis.flowers$phast1 <- GenomicRanges::makeGRangesFromDataFrame(imported_genomic_data.arabidopsis.flowers$phast1, keep.extra.columns = TRUE)
imported_genomic_data.arabidopsis.flowers <- imported_genomic_data.arabidopsis.flowers[!(names(imported_genomic_data.arabidopsis.flowers)=="phast2")]
names(imported_genomic_data.arabidopsis.flowers)[which(names(imported_genomic_data.arabidopsis.flowers)=="phast1")] <- "phastcons"
ChIPpeaks <- c(PI = "data/Arabidopsis_thaliana/flowers/PI_athal_flowers.bed",
               AP3 = "data/Arabidopsis_thaliana/flowers/AP3_athal_flowers.bed",
               AG = "data/Arabidopsis_thaliana/flowers/AG_athal_flowers.bed")
TFBSdata.athal.flowers <- Wimtrap::getTFBSdata(pfm = "data/Arabidopsis_thaliana/PFMs_athal.pfm",
                                               TFnames = names(ChIPpeaks),
                                               genome_sequence = paste0("data/Arabidopsis_thaliana/genome_athal_chr", seq(1,5), ".fa.gzip"),
                                               imported_genomic_data = imported_genomic_data.arabidopsis.flowers)
TFBSmodel.athal.flowers <- Wimtrap:::buildTFBSmodel(TFBSdata.athal.flowers,
                                                    ChIPpeaks = ChIPpeaks,
                                                    model_assessment = TRUE)
save(TFBSmodel.athal.flowers, file = "models/TFBSmodel_athal_flowers.RData")

