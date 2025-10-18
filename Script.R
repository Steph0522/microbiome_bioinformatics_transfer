
# Parte 1 -----------------------------------------------------------------


# Ajuste de rutas y directorio de trabajo
# ruta que contiene los archivos fastq
ruta_fastq <- "C:/Users/shere/Documents/microbiome_bioinformatics_transfer/Datos/seqs_transfer/"

# Crear carpeta para guardar resultados y ajusta directorio de trabajo
ruta_resultados <-"C:/Users/shere/Documents/microbiome_bioinformatics_transfer/resultados"
setwd(ruta_resultados)

fnFs <- sort(list.files(ruta_fastq, pattern = "_1.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(ruta_fastq, pattern = "_2.fastq.gz", full.names = TRUE))


head(fnFs)


head(fnRs)

# Parte 2 -----------------------------------------------------------------

setwd("../")

# install.packages("tidyverse")
# 
# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")
# 
# BiocManager::install("phyloseq")
# BiocManager::install("Biostrings")
# BiocManager::install("dada2")
# BiocManager::install("ShortRead")
# 


library(dada2)
library(Biostrings)
library(ShortRead)


# Primers usados en PCR
FWD <- "CCTACGGGNGGCWGCAG"
REV <- "GACTACHVGGGTATCTAATCC"

allOrients <- function(primer) {
  require(Biostrings)
  dna <- DNAString(primer) 
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convertir de nuevo a vector de caracteres
}

# Determinar todas las orientaciones posibles de los primers
FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)
FWD.orients


# #esta parte toma aprox: 9 min
# 
# Crear un directorio nuevo con filtrado
fnFs.filtN <- file.path(ruta_resultados, "filtN", basename(fnFs))
fnRs.filtN <- file.path(ruta_resultados, "filtN", basename(fnRs))

filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, multithread = TRUE)


primerHits <- function(primer, fn) {
# Cuenta el numero de lecturas en donde se encontraron primers
 nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
 return(sum(nhits > 0))
}
 
# Crear tabla para la deteccion de los primers
primeres_detected <- rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[1]]),
                          FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[1]]),
                          REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[1]]),
                          REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[1]]))

primeres_detected


#write.csv(primeres_detected, "primers_detected.csv")


## py -m pip install cutadapt

## py -m cutadapt --version

# #poner ruta aquí
# cutadapt <- "C:/Users/shere/AppData/Local/Programs/Python/Python314/Scripts/cutadapt.exe" # en mi pc está aquí
# 
# system2(cutadapt, args = "--version")


# #correr cutadapt
path.cut <- file.path(ruta_fastq, "cutadapt")
if(!dir.exists(path.cut)) dir.create(path.cut)
fnFs.cut <- file.path(path.cut, basename(fnFs))
fnRs.cut <- file.path(path.cut, basename(fnRs))
 
FWD.RC <- dada2:::rc(FWD)
REV.RC <- dada2:::rc(REV)
R1.flags <- paste("-g", FWD, "-a", REV.RC)
R2.flags <- paste("-G", REV, "-A", FWD.RC)
 
# Correr Cutadapt
for (i in seq_along(fnFs)) {
 system2(
   cutadapt,
   args = c(
     R1.flags,
     R2.flags,
     "-n",2,
     "-o", fnFs.cut[i],
     "-p",fnRs.cut[i],
     fnFs.filtN[i],
     fnRs.filtN[i],
     "--minimum-length=1"
   )
 )
}


primers_detected_despues = rbind(
 FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[1]]),
 FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[1]]),
 REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[1]]),
 REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[1]])
)
 
 
#write.csv(primers_detected_despues ,"primers_detected_despues.csv")


# Parte 3 -----------------------------------------------------------------
ruta_fastq <- "C:/Users/shere/Documents/microbiome_bioinformatics_transfer/Datos/seqs_transfer/"
path.cut <- file.path(ruta_fastq, "cutadapt")


# Ruta
cutFs <- sort(list.files(path.cut, pattern = "_1.fastq.gz", full.names = TRUE))
cutRs <- sort(list.files(path.cut, pattern = "_2.fastq.gz", full.names = TRUE))

# Extraer nombres
get.sample.name <- function(fname) strsplit(basename(fname), "_")[[1]][1]
sample.names <- unname(sapply(cutFs, get.sample.name))
head(sample.names)


vector =c("SRR24359981", "SRR24359982", "SRR24359983", "SRR24359984",
          "SRR24359985", "SRR24359986")
print(vector)


plotQualityProfile(cutFs[1:2])
plotQualityProfile(cutRs[1:2])

#duración: 10 min aprox
filtFs <- file.path(ruta_resultados, "filtered2", basename(fnFs))
filtRs <- file.path(ruta_resultados, "filtered2", basename(fnRs))
filtFs

names(filtFs) <- sample.names
names(filtRs) <- sample.names

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(280,200),maxN=0, maxEE=c(2,2),
                     truncQ=2, rm.phix=TRUE,compress=TRUE, multithread=TRUE)

out
#write.csv(out, "out.csv")


#duración = 11 min aprox.
errF <- learnErrors(filtFs, multithread = TRUE)
errR <- learnErrors(filtRs, multithread = TRUE)

png("error_model.png",
    units = "in",
    height = 7,
    width = 10,
    res = 300)
plotErrors(errF)
dev.off()

#saveRDS(errR, "errR.RDS")


#duración = 2 min aprox
derepFs <- derepFastq(filtFs, verbose = TRUE)
derepRs <- derepFastq(filtRs, verbose = TRUE)
#saveRDS(derepFs, "derepFs.RDS")


#duración = 4 min aprox
dadaFs <- dada(derepFs, err = errF, multithread = TRUE)
dadaRs <- dada(derepRs, err = errR, multithread = TRUE)
#saveRDS(dadaRs, "dadaRs.RDS")


mergers <- mergePairs(dadaFs, derepFs, dadaRs, derepRs, verbose=TRUE)
#saveRDS(mergers, "mergers.RDS")


mergers = readRDS("mergers.RDS")

#tabla
seqtab <- makeSequenceTable(mergers)
dim(seqtab)


#remoción quimeras
seqtab_nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = T,
                                    verbose = T)
dim(seqtab_nochim)


getN <- function(x)sum(getUniques(x))
stats <- cbind(
  out,
  sapply(dadaFs, getN),
  sapply(dadaRs, getN),
  sapply(mergers, getN),
  rowSums(seqtab_nochim)
)

colnames(stats) <- c("input",
                     "filtered",
                     "denoisedF",
                     "denoisedR",
                     "merged",
                     "nonchim")

stats
write.csv(stats, "denoising-stats.csv")


# Parte 4 -----------------------------------------------------------------


#duración : aprox 7 min

classifier <- "greengenes2_trainset.fa.gz"
taxa <- assignTaxonomy(seqtab_nochim, classifier, multithread = TRUE, tryRC = TRUE)

# Visualizar lo que se genero despues de la asignacion
taxa_print <- taxa
rownames(taxa_print) <- NULL

head(taxa_print)
dim(taxa_print)

#write.csv(taxa, "taxa.csv")

# Exportar objetos generados durante el preprocesamiento
# save( errF, dadaFs, dadaRs,seqtab_nochim, taxa,
#       file = "data.RData")
# write.csv(taxa, "taxonomy.csv")
# write.csv(seqtab_nochim, "table.csv")
# write.csv(stats, "stats.csv")

# Parte 5 -----------------------------------------------------------------

library(phyloseq)
library(Biostrings)
library(tidyverse)
theme_set(theme_bw())

samples.out <- rownames(seqtab_nochim)
samdf <-read.delim("Datos/seqs_transfer/sra-metadata.tsv") %>% 
  dplyr::rename(Origen="Host.Life.Stage..sample.", Tipo="Host.Tissue.Sampled..sample.")
seqtab.nochim <- seqtab_nochim[match(samdf$ID, rownames(seqtab_nochim)),]
rownames(samdf) <- samdf$ID

ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps


dna <- Biostrings::DNAStringSet(taxa_names(ps))
names(dna) <- taxa_names(ps)
ps <- merge_phyloseq(ps, dna)
taxa_names(ps) <- paste0("ASV", seq(ntaxa(ps)))
ps

ps.prop <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))

my_phyla <- tax_glom(ps.prop ,taxrank = "Phylum")
my_bar_phyla <- plot_bar(my_phyla, fill="Phylum") + 
  facet_grid(~Origen+Tipo, scales="free_x")+
  ggtitle("Abundancia relativa de Filos") 
my_bar_phyla

my_Class <- tax_glom(ps.prop ,taxrank = "Class")
my_bar_Class <- plot_bar(my_Class, fill="Class") + 
  facet_grid(~Origen+Tipo, scales="free_x")+
  ggtitle("Abundancia relativa de clases") 
my_bar_Class


#madre
madre <- subset_samples(ps, Origen=="adult")
madre

#embrion
embrion <- subset_samples(ps, Origen=="embryo")
embrion

plot_richness(madre, x="Tipo", measures=c("Observed","Shannon", "Simpson"), color="Tipo")

plot_richness(madre, x="Tipo", measures=c("Observed","Shannon", "Simpson"), color="Tipo")+
  geom_boxplot()

plot_richness(embrion, x="Tipo", measures=c("Observed","Shannon", "Simpson"), color="Tipo")+
  geom_boxplot()


#ps_rare <- rarefy_even_depth(ps_fil, sample.size = min(sample_sums(ps_fil)), rngseed = 19)
#ps_rare

# Transformar la data  proporciones
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
plot_ordination(ps.prop, ord.nmds.bray, color="Tipo", shape="Origen",title="Bray NMDS")

ord.nmds.jac <- ordinate(ps.prop, method="PCoA", distance="jaccard")
plot_ordination(ps, ord.nmds.jac, color="Tipo", shape="Origen",title="Bray NMDS")

# install.packages("devtools")
# library(devtools)
# install_github("Steph0522/MicroBioMeta")
# library(MicroBioMeta)

taxonomys <- taxa %>%
  as.data.frame() %>%
  tidyr::unite("taxonomy", Kingdom:Species, sep = ";")

tabla = merge_feature_taxonomy(t(seqtab.nochim), taxonomys)

abundance_barplot(
  table = tabla,
  metadata = samdf %>% dplyr::rename("SAMPLEID"="ID"),
  level = "phylum",
  top_n_groups = 10,
  x_col = "Tipo",
  facet_col = "Origen"
)+theme(axis.text.x = element_text(angle = 50, hjust = 0.95))

abundance_barplot(
  table = tabla,
  metadata = samdf %>% dplyr::rename("SAMPLEID"="ID"),
  level = "genus",
  top_n_groups = 20,
  x_col = "Tipo",
  facet_col = "Origen",add_remained = TRUE,
)+theme(axis.text.x = element_text(angle = 50, hjust = 0.95))


abundance_heatmap_plot(
  table = tabla,
  metadata = samdf,
  condition1 = "Origen",
  condition2 = "Tipo",
  top_n = 20,
  cluster = TRUE,
  show_column_names = FALSE
)


alpha_hill_plot(
  table = tabla,
  metadata = samdf,
  x_col = "Tipo",
  fill_col = "Tipo",
  stat = "anova"
) + theme(axis.text.x = element_blank(),
          axis.ticks.x = element_blank())


beta_div_plot(
  table = tabla,
  metadata = samdf,
  group_col = "Tipo",
  shape_col = "Origen",
  ordination = "PCoA",
  distance = "aitchison"
)



