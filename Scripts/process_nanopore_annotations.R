
# process and filter annotations generated from Nanopore

# This script requires the minimap2.intersectWS271.gff, Intron_minimap2.gff, WS271 geneIDs, WS271 annotations and *.SJ.out.tab files of STAR alignment
# c_elegans.PRJNA13758.WS271.geneIDs.txt
# caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf
# Intron_caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf
# STAR_DRR003389.SJ.out.tab, STAR_DRR024071.SJ.out.tab, STAR_DRR024072.SJ.out.tab, STAR_SRR5202807.SJ.out.tab, STAR_SRR5202808.SJ.out.tab, STAR_SRR5202811.SJ.out.tab, STAR_SRR5202812.SJ.out.tab

# packages
library(dplyr)
library(reshape2)
library(roperators)
options(scipen = 100)

# This script requires the WS271 geneIDs, WS271 annotations and .SJ.out.tab files of STAR alignment
# "c_elegans.PRJNA13758.WS271.geneIDs.txt"
# "caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf"
# "Intron_caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf"
# "STAR_DRR003389.SJ.out.tab","STAR_DRR024071.SJ.out.tab","STAR_DRR024072.SJ.out.tab","STAR_SRR5202807.SJ.out.tab","STAR_SRR5202808.SJ.out.tab", "STAR_SRR5202811.SJ.out.tab","STAR_SRR5202812.SJ.out.tab"

# inputs: "minimap2.intersectWS271.gff" and "Intron_minimap2.gff"
args <- commandArgs(trailingOnly = TRUE)






# process WS271-intersected gff annotations --------------------------------------------------------

# label transcripts with most-overlapping genes
t <- read.table(args[1], sep = "\t")
e <- t[t$V3 == "exon",]

# transcripts and gene V9 = Nanopore transcript id, V18 = Wormbase gene id
e$V9 <- gsub("transcript_id |;", "", e$V9)
e$V18 <-  gsub("gene_id |;.*$", "", e$V18)

# count overlap exons
e <- e %>% group_by(V9, V18) %>% mutate(overlapping_exons = n())
e <- e %>% distinct(V9, V18, overlapping_exons)
e[e$V18 == ".", "overlapping_exons"] <- 0

# most-overlapping
e <- e %>% group_by(V9) %>% mutate(most_overlapping = ifelse(overlapping_exons == max(overlapping_exons), "yes", "no"))
e <- e[e$most_overlapping == "yes",]
e <- e %>% group_by(V9) %>% slice(which.max(nchar(V18)))

# Not-annotated transcripts were labeled with Pinfish read IDs
g <- t[t$V3 == "mRNA", ]
g <- cbind(g$V9, colsplit(g$V9, "; transcript_id ", names = c("gene","transcript")))
g$gene <- gsub("gene_id ", "", g$gene)
g$transcript <- gsub(";", "", g$transcript)
g <- g %>% distinct(gene, transcript)
eg <- left_join(e, g, by=c("V9" = "transcript"))
g <- eg[,c("V18","gene")] %>% group_by(gene) %>% slice(which.max(nchar(V18)))
eg <- left_join(eg, g, by=c("gene"))
eg$gene <- ifelse(eg$V18.y != ".", eg$V18.y, eg$gene)
eg <- eg[,c("V9","gene")]
colnames(eg) <- c("transcript","wb")

# label with CGC names
n <- read.table("c_elegans.PRJNA13758.WS271.geneIDs.txt", sep = ",")
colnames(n) <- c("species","wb","cgc","seq","live","type")
n$cgc <- as.character(n$cgc)
n$seq <- as.character(n$seq)
n$gene <- ifelse(n$cgc != "", n$cgc, n$seq)
n <- n %>% distinct(wb, gene)
gene_name <- left_join(eg, n, by="wb")
gene_name$gene <- ifelse(!is.na(gene_name$gene), gene_name$gene, gene_name$wb)

# merge
n <- read.table(args[1], sep = "\t")
n <- n[n$V3 == "mRNA",]
n$transcript <- gsub(".*transcript_id |;", "", n$V9)
n <- left_join(n, gene_name, by="transcript")






# filter 5' truncated transcripts ---------------------------

# WS271 data (5'-end annotations)
g <- read.table("caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf", sep = "\t") 
g <- g[g$V3 == "transcript",]
g$V1 <- gsub("MtDNA", "M", g$V1)
g$V1 <- paste("chr", g$V1, sep = "")
g$rnaStart <- ifelse(g$V7 == "+", g$V4, g$V5)
g$left <- g$rnaStart - 100
g$right <- g$rnaStart + 100
g$wb <- gsub("^gene_id |;.*$", "", g$V9)
g <- g[,c("V1","V7","V9","left","right","wb")]
colnames(g) <- c("chr","strand","transcript_id","left","right","wb")

# whether the 5'end is within +- 100 nt of annotated 5'ends
d <- left_join(n, g, by="wb") %>% filter(ifelse(V7 == "+", left <= V4 & V4 <= right, left <= V5 & V5 <= right))
d$annotatedStart <- "yes"
d <- d[,c("transcript","annotatedStart")]

# add novel transcripts
m <- data.frame(transcript = n[n$wb == n$gene, "transcript"], annotatedStart = "yes")
d <- rbind(d, m)
d <- d %>% distinct(transcript, annotatedStart)

# filtered transcripts
n <- read.table(args[1], sep = "\t")
n$transcript <- gsub("^.*transcript_id |;.*$", "", n$V9)
d <- left_join(n, d, by="transcript")
d <- d[!is.na(d$annotatedStart),]
filtered <- d[,1:9]






# prepare suppoprted junctions ----------------------------------------------

# Illumina junctions
junc_files <- list.files(pattern="*.SJ.out.tab", recursive = TRUE)
d <- data.frame()
for (f in junc_files) {
  
  t <- read.table(f, sep = "\t")
  d <- rbind(d, t)
  
}
colnames(d) <- c("chr",
                 "intronStart",
                 "intronEnd",
                 "strand",     # 0: undefined, 1: +, 2: -
                 "GTAG",       # 0: non-canonical, 1: GT/AG, 2: CT/AC, 3: GC/AG, 4: CT/GC, 5: AT/AC, 6: GT/AT
                 "annotated",
                 "number of uniquely mapping reads",
                 "number of multi-mapping reads",
                 "maximum spliced alignment overhang")
d <- d[,c("chr","intronStart","intronEnd")]

# merge Illumina junctions and WS271 junctions
w <- read.table("Intron_caenorhabditis_elegans.PRJNA13758.WBPS14.canonical_geneset.gtf", sep = "\t")
w <- w[,c("V1", "V4", "V5")]
w$V1 <- gsub("MtDNA", "M", w$V1)
w$V1 <- paste("chr", w$V1, sep = "")
colnames(w) <- c("chr", "intronStart", "intronEnd")
d <- rbind(d, w)
d <- d %>% distinct(chr,intronStart,intronEnd)
d$supported <- "yes"
supported <- d






# filter out transcripts that containing not-supported junctions --------------

# introns from collapsed transcripts
i <- read.table(args[2], sep = "\t")
i$transcript <- gsub("^.*transcript_id |; intron_number.*$", "", i$V9)
i <- i[c("V1","V4","V5","transcript")]
colnames(i) <- c("chr","intronStart","intronEnd","transcript")

# list of not-supported transcripts
d <- supported
m <- left_join(i, d, by=c("chr","intronStart","intronEnd"))
m <- m[is.na(m$supported),]
m$supported <- "no"
m <- m %>% distinct(transcript, supported)

# filter not-supported transcripts
n <- filtered
n$transcript <- gsub("^.*transcript_id |;$", "", n$V9)
n <- left_join(n, m, by="transcript")
n <- n[is.na(n$supported), 1:9]
outname <- strsplit(args[1], "\\.")[[1]][1]
write.table(n, paste(outname, ".filtered.gff", sep = ""), sep = "\t", quote = F, col.names = F, row.names = F)





