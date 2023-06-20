#Hector Gutierrez
#Camila Rodriguez

#--------------Parte 1 -------------------------------------------

library(Biostrings)
library(ade4)
library(seqinr)
library(adegenet)
library(ape)


setwd("C:/Users/alang/Documents/Tec/Tareas/Analisis de Biologia Computacional")


corona_virus <- c("NC_045512", "OP435368", "OQ918256", "BS007312", "OQ913932", "OP848485", "ON291271", "MT994849",
                  "OK096766", "MW466791")  


virus_sequences <- read.GenBank(corona_virus)
virus_sequences


write.dna(virus_sequences,file = "coronavirus_seqs.fasta",format = "fasta")

#5
virus_seq_not_align <- readDNAStringSet("coronavirus_seqs.fasta", format = "fasta")
class(virus_seq_not_align)

virus_seq_not_align


#------------Parte 2 - 3-------------------------------------------------


library(DECIPHER)


# Alineamiento de las primeras 150 posiciones
virus_seq_not_align_150 <- virus_seq_not_align[,1:150]
virus_seq_not_align_150 <- OrientNucleotides(virus_seq_not_align_150)
virus_seq_align_150 <- AlignSeqs(virus_seq_not_align_150)

# Alineamiento de los nucleótidos 500 al 650
virus_seq_not_align_500_650 <- virus_seq_not_align[,500:650]
virus_seq_not_align_500_650 <- OrientNucleotides(virus_seq_not_align_500_650)
virus_seq_align_500_650 <- AlignSeqs(virus_seq_not_align_500_650)


BrowseSeqs(virus_seq_align_150)
BrowseSeqs(virus_seq_align_500_650)

#------------ Parte 4 - 5---------------------------------------------------------------------------------------------------------------



writeXStringSet(virus_seq_align_150, file = "coronavirus_seq_align_150.fasta")
writeXStringSet(virus_seq_align_500_650, file = "coronavirus_seq_align_500_650.fasta")


virus_aligned_150 <- read.alignment("coronavirus_seq_align_150.fasta", format = "fasta")
virus_aligned_500_650 <- read.alignment("coronavirus_seq_align_500_650.fasta", format = "fasta")


matriz_distancia_150 <- dist.alignment(virus_aligned_150, matrix = "similarity")
matriz_distancia_500_650 <- dist.alignment(virus_aligned_500_650, matrix = "similarity")
as.data.frame(as.matrix(matriz_distancia_150))
as.data.frame(as.matrix(matriz_distancia_500_650))

tablas_grises_150 <- as.data.frame(as.matrix(matriz_distancia_150))
tablas_grises_500_650 <- as.data.frame(as.matrix(matriz_distancia_500_650))

table.paint(tablas_grises_150, cleg = 0, clabel.row = .5, clabel.col = .5)
table.paint(tablas_grises_500_650, cleg = 0, clabel.row = .5, clabel.col = .5)

#--------------Parte 6 - 7 -------------------------------------------------------------

library(phytools)
library(maps)
library(viridis)
library(viridisLite)
library(ggtree)
library(ggplot2)

virus_tree <- nj(matriz_distancia_150)
virus_tree2 <- nj(matriz_distancia_500_650)

virus_colors <- c("red", "blue", "#2E8B57", "purple","orange", "#008B8B",
                  "#8B795E", "#CD6090", "brown", "black")

virus_tree <- ladderize(virus_tree)




#Titulo
plot(virus_tree, main = "Arbol Filogenetico del virus SARS-COV2", tip.color=virus_colors)

# Asignar nombres y ubicaciones a cada virus
virus_id <- c("NC_045512", "OP435368", "OQ918256", "BS007312", "OQ913932", "OP848485", "ON291271", "MT994849", "OK096766", "MW466791")
virus_date <- c("China", "Finlandia", "India", "Japón", "USA", "Australia", "Francia", "Irán", "Alemania", "Corea del Sur")


tip_dates <- data.frame(tips=virus_tree$tip.label, date = virus_date)


legend("bottomright", legend = paste(tip_dates$tips, " - (", tip_dates$date, ")", sep = ""),
       pch = 18, col = c("red", "blue", "#2E8B57", "purple", "orange", "#008B8B", "#8B795E", "#CD6090", "brown",  "black" ), pt.bg = "white", title = "Codigo de Accesion - Ubicacion")



                                                 
                                    
