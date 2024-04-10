#_______________________________________________________#
#                                                       #
#     ~*~                                      ~*~      #
#        code for analysis of sRNA-seq of               #
#   different breeds of pigs, supplementary file for    #
#                 Asratian et. al. 2024                 #
#              Code written by Signe Skog               #
#                                                       #
#     ~*~                                      ~*~      #
#_______________________________________________________#

library(seqpac)
library(tidyverse)
library(patchwork)


# Preparation of files ----------------------------------


merge_lanes(in_path = "D:/2023/pig_sperm_2023/fastq",
            out_path ="D:/2023/pig_sperm_2023/merged",
            threads = 6)

make_trim(input="D:/2023/pig_sperm_2023/merged",
          output="D:/2023/pig_sperm_2023/trimmed",
          adapt_3 = "AGATCGGAAGAGCACACGTCTGAACTCCA",
          threads=1)

counts<-make_counts(input="D:/2023/pig_sperm_2023/trimmed",
                    threads=1)

pheno <- data.frame(sample_ID=colnames(counts$counts), 
                    id=colnames(counts$counts), 
                    row.names = colnames(counts$counts), 
                    breed=rep(c("Duroc", "Hampshire", "Landrace", "Yorkshire"), each=6),
                    lines=rep(c("Father", "Mother"), each=12))

pac <- make_PAC(counts=counts, pheno=pheno)
save(pac, file="pac_pigs_master.Rdata")
load("pac_pigs_master.Rdata")

pac
pac_f <- PAC_filter(pac, size=c(16,75), threshold=5, coverage=10, stat=TRUE)
pac_f <- PAC_norm(pac_f)
pac_h <- PAC_filter(pac_f, threshold = 1, coverage = 5, norm="cpm", stat=TRUE)


## Mapping ---------------------------------------------


### Genomes --------------------------------------------

map_reanno(pac_h, 
           output_path="D:/OUT",
           ref_paths=list(duroc="D:/Genomes/sus/ensembl_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa",
                          landrace="D:/Genomes/sus/landrace/GCA_001700215.1_Landrace_pig_v1_genomic.fa",
                          hampshire="D:/Genomes/sus/hamp/GCA_001700165.1_Hampshire_pig_v1_genomic.fa",
                          yorkshire="D:/Genomes/sus/largewhite/GCA_001700135.1_Large_White_v1_genomic.fa"),
           mismatches=3)

reanno_genome <- make_reanno(reanno_path="D:/OUT", 
                             PAC=pac_h, 
                             mis_fasta_check = TRUE)

pac_h <- add_reanno(reanno_genome,
                    mismatches = 3, merge_pac=pac_h)

head(pac_h@Anno)


### Biotypes -------------------------------------------


map_reanno(pac_h, 
           ref_paths=list(ens="D:/Genomes/sus/ncRNA_ensembl_2023/Sus_scrofa.Sscrofa11.1.ncrna.fa",
                          mattRNA="D:/Genomes/sus/GtRNAdb_2023/susScr11-mature-tRNAs.fa",
                          piRNA="D:/Genomes/sus/piRBase_2023/ssc.v3.0.fa",
                          pc="D:/Genomes/sus/pc/Sus_scrofa.Sscrofa11.1.cds.rm_df.fa"), 
           output_path="D:/Genomes/out",
           type="internal", 
           mismatches=3,  
           import="biotype", 
           threads=1)

reanno_biotype <- make_reanno(reanno_path="D:/Genomes/out",
                              PAC=pac_h, 
                              mis_fasta_check = TRUE)

bio_search <- list(ens=c("snRNA", "snoRNA","Mt_rRNA", "rRNA", "miRNA", "tRNA",
                         "Mt_tRNA","scaRNA", "lncRNA", "misc_RNA", "Y_RNA", 
                         "ribozyme","vault_RNA", "miscRNA"),
                   mattRNA=c("tRNA"),
                   piRNA=c("piR"),
                   pc=c("ENS"))

pac_h <- add_reanno(reanno_biotype, bio_search=bio_search, 
                    type="biotype", bio_perfect=TRUE, 
                    mismatches = 3, merge_pac=pac_h)

pac_h<-simplify_reanno(input=pac_h, 
                       hierarchy = list(mt_rRNA="ens_Mt_rRNA",
                                        rRNA="ens_rRNA",
                                        mt_tRNA="ens_Mt_tRNA",
                                        tRNA="mattRNA|ens_tRNA",
                                        miRNA="ens_miRNA",
                                        piRNA="piRNA",
                                        pc="pc|pc_ENS"),
                       mismatches = 3,
                       merge_pac=TRUE,
                       bio_name = "Biotypes_3mm")

pac_h<-simplify_reanno(input=pac_h, 
                       hierarchy = list(mt_rRNA="ens_Mt_rRNA",
                                        rRNA="ens_rRNA",
                                        mt_tRNA="ens_Mt_tRNA",
                                        tRNA="mattRNA|ens_tRNA",
                                        miRNA="ens_miRNA",
                                        piRNA="piRNA",
                                        pc="pc|pc_ENS"),
                       mismatches = 0,
                       merge_pac=TRUE,
                       bio_name = "Biotypes_0mm")

head(pac_h@Anno)
table(pac_h@Anno$Biotypes_0mm)
unique((pac_h@Anno$Biotypes_0mm))


### Mitochondria --------------------------------------


map_reanno(pac_h, 
           output_path="D:/OUT",
           ref_paths=list(mito="D:/Genomes/sus/NC_000845.1_mito.fa"),
           mismatches=3)

reanno_genome <- make_reanno(reanno_path="D:/OUT", 
                             PAC=pac_h, 
                             mis_fasta_check = TRUE)

pac_h <- add_reanno(reanno_genome,
                    mismatches = 3, merge_pac=pac_h)

head(pac_h@Anno)


#readjust the biotypes

table(pac_h@Anno$Biotypes_0mm)
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_0mm %in% "piRNA"),]$Biotypes_0mm <- "mt_piRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_0mm %in% "rRNA"),]$Biotypes_0mm <- "mt_rRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_0mm %in% "tRNA"),]$Biotypes_0mm <- "mt_tRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_0mm %in% "pc"),]$Biotypes_0mm <- "mt_pc"
table(pac_h@Anno$Biotypes_0mm)


table(pac_h@Anno$Biotypes_3mm)
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_3mm %in% "piRNA"),]$Biotypes_3mm <- "mt_piRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_3mm %in% "rRNA"),]$Biotypes_3mm <- "mt_rRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_3mm %in% "tRNA"),]$Biotypes_3mm <- "mt_tRNA"
pac_h@Anno[(pac_h@Anno$Any_genome2 %in% "Hit" & 
              pac_h@Anno$Biotypes_3mm %in% "pc"),]$Biotypes_3mm <- "mt_pc"
table(pac_h@Anno$Biotypes_3mm)
unique(pac_h@Anno$Biotypes_3mm)


### Save and load progress -----------------------------

save(pac_h, file="pac_pig_anno.Rdata")
load("C:/Users/sigsk47/OneDrive - Linköpings universitet/Documents/2023/PigBreeds_sRNA/pac_pig_anno.Rdata")
