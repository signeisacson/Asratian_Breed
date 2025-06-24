
#_______________________________________________________#
#                                                       #
#     ~*~                                      ~*~      #
#        code for analysis of sRNA-seq of               #
#   different breeds of pigs, supplementary file for    #
#                 Asratian et. al. 2025                 #
#              Code written by Signe Skog               #
#                                                       #
#     ~*~                                      ~*~      #
#_______________________________________________________#

library(seqpac)
library(tidyverse)
library(patchwork)
library(circlize)
library(ggplot2)


# Preparation of files ----------------------------------


merge_lanes(in_path = "D:/2023/pig_sperm_2023/fastq",
            out_path ="D:/2023/pig_sperm_2023/merged",
            threads = 6)

make_trim(input="D:/2023/pig_sperm_2023/merged",
          output="D:/2023/pig_sperm_2023/trimmed",
          adapt_3 = "AGATCGGAAGAGCACACGTCTGAACTCCA",
          threads=6)

counts<-make_counts(input="D:/2023/pig_sperm_2023/trimmed",
                    threads=6)

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

rm(pac)
rm(pac_f)
gc()

## Mapping ---------------------------------------------


### Genomes --------------------------------------------

map_reanno(pac_h, 
           output_path="D:/OUT",
           ref_paths=list(duroc="D:/Genomes/sus/ensembl_genome/Sus_scrofa.Sscrofa11.1.dna.toplevel.fa",
                          landrace="D:/Genomes/sus/landrace/GCA_001700215.1_Landrace_pig_v1_genomic.fa",
                          hampshire="D:/Genomes/sus/hamp/GCA_001700165.1_Hampshire_pig_v1_genomic.fa",
                          yorkshire="D:/Genomes/sus/largewhite/GCA_001700135.1_Large_White_v1_genomic.fa"),
           mismatches=3, threads = 6)

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
                          pc="D:/Genomes/sus/pc/Sus_scrofa.Sscrofa11.1.cds.all.fa"), 
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
load("E:/Pig_BreedData/PigBreeds_sRNA/pac_pig_anno.Rdata")


## Preprocessing ---------------------------------------


pac_h <- PAC_summary(pac_h, norm="cpm", pheno_target = list("breed"))
pac_h <- PAC_summary(pac_h, norm="cpm", pheno_target = list("lines"))
pac_h <- PAC_summary(pac_h, norm="cpm", pheno_target = list("breed"), type = "log2FC")
pac_h <- PAC_summary(pac_h, norm="cpm", pheno_target = list("lines"), type="log2FC")



# Figures ----------------------------------------------

col_5 <- c("#2355DC","#1DB4D7","#FFB100","#BD484B")

## Figure 1 --------------------------------------------



#Fig 1 c - pie chart of biotypes in 3 mm


pie <- PAC_pie(pac_h, summary="all", 
               anno_target = list("Biotypes_3mm",
                                  c("miRNA", "mt_pc","mt_piRNA", "mt_rRNA", "mt_tRNA",
                                    "no_anno", "other", "pc", "piRNA", "rRNA", "tRNA")),
               colors=c("#8121d0","#dc00a3","#ff1e70",
                        "#ff5e64","#00c4a2", "#7dca44","#d4bf00","#7d7d7d", "#ff952d",
                        "#d3d3da","#6b6b6e"),
                labels="all")
pie

###Fig S1  --------------------------------------------

### Fig S1 a

pie <- PAC_pie(pac_h, summary="all", 
               anno_target = list("Biotypes_0mm",
                                  c("miRNA", "mt_pc","mt_piRNA", "mt_rRNA", "mt_tRNA",
                                    "no_anno", "other", "pc", "piRNA", "rRNA", "tRNA")),
               colors=c("#8121d0","#dc00a3","#ff1e70",
                        "#ff5e64","#00c4a2", "#7dca44","#d4bf00","#7d7d7d", "#ff952d",
                        "#d3d3da","#6b6b6e"),
               labels="all")
pie



## Figure 2 --------------------------------------------


library(factoextra)
library(FactoMineR)


### Fig 2 a

tst.pca <- PCA(t(pac_h@norm$cpm), graph = FALSE)
fviz_pca_ind(tst.pca, label = "ind",
             repel = TRUE, 
             pointsize = 3, 
             habillage = as.factor(pac_h@Pheno$breed),
             palette = col_5, 
             addEllipses = TRUE, 
             ellipse.level=0.6)+ 
  scale_shape_manual(values=c(15,15,16,16))


### Fig 2 b

pc1 <- fviz_contrib(tst.pca, choice ="var", axes = 1)
ggplot(pc1$data, aes(x=contrib))+
  geom_histogram(fill="grey")+
  theme_bw()+
  labs(y= "Number of unique sRNA ( x 103)", x = "contribution (%)")


pc1<-pc1$data[order(pc1$data$contrib, decreasing = TRUE),]
contr_10000 <- pc1[1:10000,]
pac_c10000 <- PAC_filter(pac_h, anno_target = rownames(contr_10000))
sb_10000<-PAC_stackbar(pac_c10000, anno_target = list("Biotypes_3mm"),
                       summary="all", color =c("#d3d3da","#6b6b6e",
                                               "#ff1e70","#ff952d","#dc00a3",
                                               "#d4bf00", "#ff5b60", "#00c4a2",
                                               "#8121d0"))
sb_10000
clipr::write_clip(contr_10000)

### Fig 2 d

pc2 <- fviz_contrib(tst.pca, choice ="var", axes = 2)
ggplot(pc2$data, aes(x=contrib))+
  geom_histogram(fill="grey")+
  theme_bw()+
  labs(y= "Number of unique sRNA ( x 103)", x = "contribution (%)")


pc2<-pc2$data[order(pc2$data$contrib, decreasing = TRUE),]
contr_10000 <- pc2[1:10000,]
pac_c10000 <- PAC_filter(pac_h, anno_target = rownames(contr_10000))
sb_10000<-PAC_stackbar(pac_c10000, anno_target = list("Biotypes_3mm"),
                       summary="all", color=c("#6b6b6e", "#d6d6dd",
                                              "#7dca44", "#00c4a2", 
                                              "#ff1e70", "#d4bf00",
                                              "#8121d0", "#ff952d",
                                              "#dc00a3", "#ff5b60",
                                               "#ff6348"))
sb_10000


###Fig S2  --------------------------------------------

### Fig S2 a

dq <- PAC_deseq(pac_h, model=~breed, pheno_target = list("breed", c("Yorkshire", "Landrace")))
res_m<-dq$result
res_m$lbl <- ""
res_m$FC<-res_m[,3]
res_m$diffexpressed <- " "
res_m$diffexpressed <- ifelse(res_m$padj<0.1 & res_m$log2FC_breed_Yorkshire_vs_Landrace > 1, "UP", "")
res_m$diffexpressed <- ifelse((res_m$padj<0.1 & res_m$log2FC_breed_Yorkshire_vs_Landrace < (-1)), "DOWN", res_m$diffexpressed)
res_m[res_m$diffexpressed %in% "UP",]$lbl <- res_m[res_m$diffexpressed %in% "UP",]$Biotypes_3mm
res_m[res_m$diffexpressed %in% "DOWN",]$lbl <- res_m[res_m$diffexpressed %in% "DOWN",]$Biotypes_3mm
table(res_m$lbl)
res_m$lbl <- factor(res_m$lbl, levels=c("miRNA", "mt_pc","mt_piRNA", "mt_rRNA", "mt_tRNA",
                                    "no_anno", "other", "pc", "piRNA", "rRNA", "tRNA"))

ggplot(res_m, aes(FC, -log10(padj), color=lbl)) + 
  geom_point(size=2.5) +
  scale_color_manual(values=c("#ff952d", "#ff6348", 
                              "#d4bf00", "#7dca44", "#00c4a2", 
                              "#6b6b6e", "#d6d6dd", "#ff5b60", 
                              "#ff1e70", "#dc00a3", "#8121d0"))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_vline(xintercept = c(-1,1), linetype="dashed")+
  theme(text=element_text(size=15), legend.position="bottom")+
  labs(y= "-log10 adjusted p value", x = "log2FC Sire versus Dame lines")

###Figure 3  --------------------------------------------

### Fig 3 a

dsq <- PAC_deseq(pac_h, model=~lines)
res<-dsq$result
res$lbl <- ""
res$FC<-res[,3]
res$diffexpressed <- " "
res[res[,3]>1 & res$padj<0.1,]$diffexpressed <- "UP"
res[res[,3]<(-1) & res$padj<0.1,]$diffexpressed <- "DOWN"
res[res$diffexpressed %in% "UP",]$lbl <- res[res$diffexpressed %in% "UP",]$Biotypes_3mm
res[res$diffexpressed %in% "DOWN",]$lbl <- res[res$diffexpressed %in% "DOWN",]$Biotypes_3mm
table(res$lbl)
res$lbl <- factor(res$lbl, levels=c("miRNA", "mt_pc","mt_piRNA", "mt_rRNA", "mt_tRNA",
                                    "no_anno", "other", "pc", "piRNA", "rRNA", "tRNA"))

ggplot(res, aes(FC, -log10(padj), color=lbl)) + 
  geom_point(size=2.5) +
  scale_color_manual(values=c("#ff952d", "#ff6348", 
                              "#d4bf00", "#7dca44", "#00c4a2", 
                              "#6b6b6e", "#d6d6dd", "#ff5b60", 
                              "#ff1e70", "#dc00a3", "#8121d0"))+
  theme_classic()+
  geom_hline(yintercept = 1, linetype="dashed")+
  geom_vline(xintercept = c(-1,1), linetype="dashed")+
  theme(text=element_text(size=15), legend.position="bottom")+
  labs(y= "-log10 adjusted p value", x = "log2FC Sire versus Dame lines")


### Fig 3 c

result_df<-readxl::read_excel("SupTab_Asratian2025.xlsx", sheet="Table S7")

result_df <- result_df[!result_df$family %in% "unk",]
result_df$show <- ""
result_df[result_df$family %in% "tRNA",]$show <- "tRNA-derived"
result_df[result_df$family %in% "SINE/tRNA",]$show <- "SINE/tRNA"
result_df[result_df$family %in% "LINE/L1",]$show <- "LINE/L1"

ggplot(result_df, aes(x=Sum.cpm.Sire.vs.Dam, y=Sum.cpm.York.vs.Land, 
                      color=show, size=Total.no.of.piRNAs))+
  geom_point()+
  theme_classic()+
  scale_color_manual(values=c("#a8a8a8", "#2c5aa0", "#192a55", "#986159"))+
  labs(y= "Total CPM Landrace vs Yorkshire", x = "Total CPM Sire vs dam lines")+
  theme(text=element_text(size=15), legend.position="bottom")


###Fig S3  --------------------------------------------


pac_mt <- PAC_filter(pac_h, anno_target = list("Any_genome2", "Hit"))
map_m <- PAC_mapper(pac_mt, 
                    ref="D:/Genomes/sus/mitochr_pigs/tRNA_mito.fa",
                    mismatches = 3)

cplt <- PAC_covplot(pac_mt, map_m, summary_target = list("cpmMeans_breed"), 
                    colors = c("#2765b6", "#2dc0e2", "#f8b720", "#b74e55"))

cplt

###Figure 4  --------------------------------------------


### Fig 4 a

pac_mt

#manually adding mito info from ncbi (https://www.ncbi.nlm.nih.gov/nuccore/NC_024511.2)
df_chr = data.frame(
  name  = c("D-loop","tRNA-Phe",  "12S rRNA",    "tRNA-Val","16s rRNA", 
            "tRNA-Leu1", "ND1", "tRNA-Ile", "tRNA-Gln", "tRNA-Met", 
            "ND2", "tRNA-Trp",  "tRNA-Ala", "tRNA-Asn", "tRNA-Cys", 
            "tRNA-Tyr", "COX1", "tRNA-Ser1", "tRNA-Asp", "COX2", 
            "tRNA-Lys", "ATP8", "ATP6", "COX3", "tRNA-Gly",
            "ND3", "tRNA-Arg", "ND4L", "ND4", "tRNA-His", 
            "tRNA-Ser2", "tRNA-Leu2",  "ND5", "ND6", "tRNA-Glu", 
            "CYTB", "tRNA-Thr", "tRNA-Pro"),
  start = c(1,1296,1366,2327,2395,
            3965,4042, 4997,5063,5137,
            5207,6249,6323,6392,6499, 
            6565, 6631, 8179, 8255,  8323, 
            9011, 9079, 9240, 9920, 10704,
            10773, 11120, 11189, 11479, 12857,           
            12926, 12985,13055, 14859, 15387,         
            15460, 16600, 16667),
  end   = c(1295,1365,2325,2394,3964,
            4039,4996, 5065,5135,5206,
            6248,6316,6390,6466, 6564, 
            6629, 8175, 8247, 8322, 9010, 
            9077, 9282, 9920, 10703, 10772, 
            11118, 11188, 11485, 12856, 12925, 
            12984, 13054,14875, 15386,15455, 
            16599, 16667, 16731))
df_chr$name<-factor(df_chr$name, ordered=TRUE, levels = df_chr$name)

h <- rep(c("D-loop","tRNA-Phe",  "12S rRNA", "",
           "tRNA-Val","16s rRNA",  "tRNA-Leu1", "",
           "ND1", "tRNA-Ile", "tRNA-Gln", "", 
           "tRNA-Met", "ND2", "tRNA-Trp",  "",
           "tRNA-Ala", "",
           "tRNA-Asn", "",
           "tRNA-Cys", 
           "tRNA-Tyr","", 
           "COX1", "",
           "tRNA-Ser1","", 
           "tRNA-Asp", "COX2", 
           "tRNA-Lys", "",
           "ATP8", "ATP6", "COX3", 
           "tRNA-Gly",
           "ND3", "",
           "tRNA-Arg", "ND4L", "ND4", "tRNA-His", 
           "tRNA-Ser2", "tRNA-Leu2",  "ND5", "ND6", "tRNA-Glu", "", 
           "CYTB", "tRNA-Thr", "tRNA-Pro"),
         c(1295,70,960,1,
           68,1570, 75, 2,
           955, 69,70,1, 
           70,1042,68, 6,   
           68, 1, 
           75, 32, 
           66, 
           65, 1, 
           1545,3, 
           69,7, 
           68, 688,67,1, 
           204,637,784,69,346,1, 
           69,297,1371,69,59,70,1821,511,69, 4,
           1140,68,64))


covpi<-PAC_mapper(pac_mt, 
                  ref="E:/Pig_BreedData/PigBreeds_sRNA/mito_suscr.fa",
                  mismatches = 3)
cov <- PAC_covplot(pac_mt, covpi, summary_target = list("cpmMeans_breed"),
                   xseq=TRUE, style="solid")
df<-cov$KX094894.1$data

dim(df)
length(h)

df2 <- as.data.frame(pivot_wider(df, names_from="Group", values_from="Coverage"))
df2$chr<-h
df2$end_column <- df2$Position
df2$start_column <- (as.numeric(df2$end_column) - 1)
df22 <- df2[,c(6,8,7,1,2,3,4,5)]
df23<-df22
cols <- c("Duroc", "Hampshire", "Landrace", "Yorkshire")

df23[cols]<-log(df23[cols])
df23$Duroc<-sub("-Inf", "0", df23$Duroc)
df23$Duroc<-as.numeric(df23$Duroc)
df23$Hampshire<-sub("-Inf", "0", df23$Hampshire)
df23$Hampshire<-as.numeric(df23$Hampshire)
df23$Landrace<-sub("-Inf", "0", df23$Landrace)
df23$Landrace<-as.numeric(df23$Landrace)
df23$Yorkshire<-sub("-Inf", "0", df23$Yorkshire)
df23$Yorkshire<-as.numeric(df23$Yorkshire)

dfmean <- df23[,c(1:3)]
dfmean$mean <- rowMeans(df23[,5:8])

#PREPARE GRAPH
circos.clear()
circos.par("track.height"=0.28)
circos.genomicInitialize(df_chr,  
                         tickLabelsStartFromZero=FALSE)

circos.track(ylim = c(0, 1),
             bg.col = c("#b2b2b2",
                        "#f4dd91", "#81B185","#f4dd91", "#81B185",
                        "#f4dd91","#b1cecf","#f4dd91","#f4dd91","#f4dd91",
                        "#b1cecf","#f4dd91","#f4dd91","#f4dd91","#f4dd91","#f4dd91",
                        "#b1cecf","#f4dd91","#f4dd91",
                        "#b1cecf","#f4dd91",
                        "#b1cecf","#b1cecf","#b1cecf","#f4dd91",
                        "#b1cecf","#f4dd91",
                        "#b1cecf","#b1cecf","#f4dd91","#f4dd91","#f4dd91",
                        "#b1cecf","#b1cecf","#f4dd91",
                        "#b1cecf","#f4dd91","#f4dd91"),
             bg.border = NA, track.height = 0.05)

ylimm<-c(max(df23[,c(5,6,7,8)])+0.5)
#ylimm<-c(0, ylimm)

dfmean$chr<-as.character(dfmean$chr)
dfmean$start_column <- as.integer(dfmean$start_column)
dfmean$end_column <- as.integer(dfmean$end_column)
dfmean$mean <- as.numeric(dfmean$mean)

circos.genomicTrack(dfmean, 
                    panel.fun = function(region, value, ...) {
                      circos.genomicLines(region, value, area=TRUE,
                                          type="l", col=c("#8121d0"), lwd=0.4)
                    })


### Fig 4 b

cplt$`tRNA-Ser_NC_000845.1:12806-12864`

### Fig 4 c

cplt$`tRNA-Arg_NC_000845.1:11000-11068`

### Fig 4 d

map <- PAC_mapper(pac_mt, 
                   ref="D:/Genomes/sus/ncRNA_ensembl_2023/Sus_scrofa.Sscrofa11.1.ncrna_rRNA.fa",
                   mismatches = 3)
cpt <- PAC_covplot(pac_mt, map, summary_target = list("cpmMeans_breed"),
                   colors = c("#2765b6", "#2dc0e2", "#f8b720", "#b74e55"))

cpt$`ENSSSCT00000019658.3_ncrna_primary_assembly:Sscrofa11.1:MT:2274:3844:1_gene:ENSSSCG00000018063.3_gene_biotype:Mt_rRNA_transcript_biotype:Mt_rRNA`


### Fig 4 g

cpt$`ENSSSCT00000019656.3_ncrna_primary_assembly:Sscrofa11.1:MT:1246:2205:1_gene:ENSSSCG00000018061.3_gene_biotype:Mt_rRNA_transcript_biotype:Mt_rRNA`



### Figure 5 -------------------------------------------


df_circ <- readxl::read_excel("SupTab_Asratian2025.xlsx", sheet="Table S8", skip=3)

chordDiagram(df_circ, annotationTrack =  c("grid"))
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  xlim = get.cell.meta.data("xlim")
  ylim = get.cell.meta.data("ylim")
  sector.name = get.cell.meta.data("sector.index")
  circos.text(mean(xlim), ylim[1] + .1, sector.name, facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5))
}, bg.border = NA)


