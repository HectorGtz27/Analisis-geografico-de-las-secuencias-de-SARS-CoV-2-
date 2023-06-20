#Descargamos los paquetes que contienen las funciones y herramientas necesarias 
#para trabajar con los datos biológicos

library(viridisLite)
library(ape)
library(ade4)
library(seqinr)
library(adegenet)
library(Biostrings)
library(DECIPHER)
library(ggtree)
library(ggplot2)
library(tidyr)

#Utilizamos el getwd() para usar el directorio de trabajo actual
#y el setwd() le da a R una dirección a un directorio diferente 
#para que ahi mismo se almacene y guarde la lectura y escritura de archivos.

getwd()
setwd("C:/Users/alang/Documents/Tec/Tareas/Analisis de Biologia Computacional")

#EJERCICIO 1 
#Obtener los 10 genomas de diferentes países 

corona_virus <- c("NC_045512", "OP435368", "OQ918256", "BS007312", "OQ913932", "OP848485", "ON291271", "MT994849",
                  "OK096766", "MW466791")        


# Leer archivo GenBank
virus_sequences <- read.GenBank(corona_virus)
write.dna(virus_sequences,file = "virus_coronavirus",format = "fasta")


#EJERCICIO 2
#Longitud de las secuencias de cada variante 

enumerar <- c(1:10)


longitud <- sapply(virus_sequences, length)
gen_long <- data.frame(Numero_Genoma = enumerar, Secuencia_ID = corona_virus, Longitud = longitud)

print(gen_long,row.names = FALSE)


#EJERCICIO 3
#Gráfico de longitudes 

grafica_longitud <- ggplot(gen_long, aes(x = Secuencia_ID, y = Longitud, fill = Secuencia_ID)) +
  geom_bar(stat = "identity") +
  labs(x = "Virus", y = "Longitud", title = "Longitud de Sars-Cov-2 de diferentes países") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 8),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

grafica_longitud


#EJERCICIO 5
#Comosición de nucleotidos de cada genoma 

virus_sequences_character <- c(as.character(virus_sequences))
nucleotidos_gen <- sapply(virus_sequences_character,count,1)

nucleotidos_gen 


#EJERCICIO 6 
#Gráfica nucleótidos de cada genoma 

nucleotidos <- as.data.frame(nucleotidos_gen)

nucleotidos

#-------NC_045512------------------------------------------

# Definir colores
colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = NC_045512, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma NC_045512") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))


# Imprimir gráfico
print(grafica_nucleotidos)

#------OP435368---------------------------------------------------------------

# Definir colores
colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos2 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = OP435368, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma OP435368") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos2)

#--------OQ918256-------------------------------------------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos3 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = OQ918256, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma OQ918256") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos3)


#----------BS007312-----------------------------------------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos4 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = BS007312, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma BS007312") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos4)


#---------OQ913932-----------------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos5 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = OQ913932, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma OQ913932") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos5)


#---------OP848485----------------------------------
colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos6 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = OP848485, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma OP848485") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos6)

#---------ON291271--------------------------------
colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos7 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = ON291271, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma ON291271") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos7)


#----------MT994849----------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos8 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = MT994849, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma MT994849") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))

# Imprimir gráfico
print(grafica_nucleotidos8)

#----------OK096766--------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos9 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = OK096766, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma OK096766") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))


# Imprimir gráfico
print(grafica_nucleotidos9)

#---------MW466791-----------------------------

colores <- c("skyblue", "steelblue", "navy", "darkgreen")

# Crear gráfico de barras
grafica_nucleotidos10 <- ggplot(nucleotidos, aes(x = rownames(nucleotidos), y = MW466791, fill = rownames(nucleotidos))) +
  geom_bar(stat = "identity") +
  labs(x = "Base nitrogenada", y = "Cantidad", title = "Composición de nucleótidos del genoma MW466791") +
  scale_x_discrete(labels = c("A", "C", "G", "T")) +
  scale_fill_manual(values = colores) +
  theme(axis.text.x = element_text(face = "bold"),
        axis.text.y = element_text(face = "bold"),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 13))


# Imprimir gráfico
print(grafica_nucleotidos10)

#EJERCICIO 8
#Calcula el %GC de cada variante

NC_045512_GC <- GC(virus_sequences_character[[1]])*100
OP435368_GC <- GC(virus_sequences_character[[2]])*100
OQ918256_GC <- GC(virus_sequences_character[[3]])*100
BS007312_GC <- GC(virus_sequences_character[[4]])*100
OQ913932_GC <- GC(virus_sequences_character[[5]])*100
OP848485_GC <- GC(virus_sequences_character[[6]])*100
ON291271_GC <- GC(virus_sequences_character[[7]])*100
MT994849_GC <- GC(virus_sequences_character[[8]])*100
OK096766_GC <- GC(virus_sequences_character[[9]])*100
MW466791_GC <- GC(virus_sequences_character[[10]])*100

gc_porcentaje <- c(NC_045512_GC, OP435368_GC, OQ918256_GC, BS007312_GC, OQ913932_GC, OP848485_GC, ON291271_GC, MT994849_GC, OK096766_GC, MW466791_GC)
genomas <- c("NC_045512", "OP435368", "OQ918256", "BS007312", "OQ913932", "OP848485", "ON291271", "MT994849", "OK096766", "MW466791")
gc_final <- data.frame(Genoma = genomas, GC_Porcentaje = gc_porcentaje)

gc_final 


#EJERCICIO 9    

ggplot(gc_final, aes(x = Genoma, y = GC_Porcentaje, fill = Genoma)) +
  geom_bar(stat = "identity") +
  labs(x = "Genoma", y = "Porcentaje de GC", title = "Porcentaje de GC de cada genoma de SARS-CoV-2") +
  theme(axis.text.x = element_text(size = 8),
        axis.text.y = element_text(size = 9),
        axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold", size = 14))


#EJERCICIO 11

virus_nombre <- rep("SARS-coV-2",10)
pais <- c("China","Finlandia","India","Japón","Estdos Unidos","Australia","Francia","Iran","Alemania","Corea del Sur")
gc_genomas <- sapply(virus_sequences_character, function(x) GC(x) * 100)

dataFrameFinal <- data.frame(Virus = virus_nombre, ID = corona_virus,País_de_Origen = pais, Longitud = longitud, Porcentaje_GC = gc_genomas)

print(dataFrameFinal, row.names = FALSE)











