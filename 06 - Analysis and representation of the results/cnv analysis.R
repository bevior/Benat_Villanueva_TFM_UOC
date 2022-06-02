# Introducción
#Este documento sirve para hacer varios cálculos y representaciones y análisis sobre los datos de los CNV predichos por GATK y DECoNT y sus respectivas intersecciones.
#El primero paso sera establecer el directorio de trabajo al mismo directorio donde se encuentra este script y la carpeta con los datos.
#Esto ha de hacerse manualmente.
#Para leer los datos adecuadamente, guardaremos el directorio donde se encuentran.
directorio_CNV<- "Datos/datos_reordenados/"
directorio_intersectos<- "Datos/intersecciones/"
#Además si aún no existe, crearemos un directorio para los resultados.
if (file.exists("Resultados")) {
  cat("The folder already exists")
} else {
  dir.create("Resultados")
}

#Creamos una variable que contenga todos los nombres de las muestras. Nos será útil en el futuro.
muestras<- read.table(file=paste(directorio_CNV, "DECoNT_all/", "DECoNT_cnv_chr_ordered.bed", sep=""), sep="\t")
muestras<- sort(unique(muestras[,6])) #Eliminamos duplicados, ordenamos alfabéticamente y guardamos los nombres de las muestras.

#Creamos una función para abrir los archivos en formato BED de manera más rápida.
open.bed<- function(X) {
  table<- read.table(file=X, sep="\t")
  return(table)
}

#Creamos una función para calcular los tamaños de los CNVs automáticamente.
get.sizes<- function(data_frame){
  size<- data_frame$end - data_frame$start
  return(size)
}

#Los archivos no tienen encabezado, por lo que nosotros mismos le crearemos uno depndiendo si es una archivo BED de GATK o de DECoNT.
header_gatk<- c("chrom", "start", "end", "ID_0", "ID_1", "sample_name", "CN", "GT")
header_decont<- c("chrom", "start", "end", "ID_0", "ID_1", "sample_name", "GT_XHMM", "GT")

#Además también crearemos una función que partiendo de un BED de GATK o DECoNT, directamente nos formatee las columnas de genotipo de manera adecuada.
factorize.GT<- function(vector) {
  vector<- factor(vector, levels= c("DEL", "DUP"))
  return(vector)
}

#Y otra para formatear la columna de las muestras como factor.
factorize.sample<- function(vector) {
  vector<- factor(vector, levels= muestras)
  return(vector)
}







# 1.- Pulido DECoNT.
#En este apartado analizaremos el pulido echo por DECoNT sobre las predicciones de XHMM.
#Establecemos la ruta del archivo.
file<- "Datos/DECoNT_XHMM_polished_cnvs.txt" #PC PERSONAL
decont.data<- read.table(file, header=TRUE, sep="\t")

#Hacemos una exploración rápida de los datos leídos.
head(decont.data)
class(decont.data)
str(decont.data)

## ¿Cuantos CNVs ha corregido como NO-CALL, DEL o DUP?
nrow(decont.data) #Numero de CNVs predichos por XHMM originalmente.
corrections<- decont.data[decont.data$XHMM.Prediction!=decont.data$DECoNT.Polished.Prediction,]; nrow(corrections) #Numero de correcciones hechas por DECoNT.
nrow(decont.data) - nrow(corrections) #Número de CNVs que no han sido corregidos.

#Presentamos el desglose de CNV corregidos a modo de tabla.
datos_tabla_01<- as.matrix(table(corrections$DECoNT.Polished.Prediction)) #Correcciones desglosadas.
datos_tabla_01<- t(datos_tabla_01)
datos_tabla_01<- cbind(datos_tabla_01, sum(datos_tabla_01))
rownames(datos_tabla_01)<- "CNV corrected to"
colnames(datos_tabla_01)[4]<- "Total"
datos_tabla_01



# 2.- GATK vs DECoNT, todos las muestras con todos sus CNVs.
#En este apartado analizaremos todos los CNV predichos por GATK y DECoNT. Sin consolidarlos ni separarlos por individuos.

#Abrimos el archivo de GATK.
file<- paste(directorio_CNV, "GATK_all/GATK_cnv.bed", sep="")
gatk.cnv<- open.bed(file) 
#Le ponemos el encabezado.
names(gatk.cnv)<- header_gatk
#Factorizamos el genotipo predicho.
gatk.cnv$GT<- factorize.GT(gatk.cnv$GT)
#Factorizamos la columna "sample_name".
gatk.cnv$sample_name<- factorize.sample(gatk.cnv$sample_name)
str(gatk.cnv) #Comprobamos como tenemos los datos.

#Abrimos el archivo de DECoNT.
file<- paste(directorio_CNV, "DECoNT_all/DECoNT_cnv.bed", sep="")
decont.cnv<- open.bed(file) 
#Le ponemos el encabezado.
names(decont.cnv)<- header_decont
#Factorizamos el genotipo predicho.
decont.cnv$GT<- factorize.GT(decont.cnv$GT)
#Facotirzamos el genotipo predicho por XHMM
decont.cnv$GT_XHMM<- factorize.GT(decont.cnv$GT_XHMM)
#Factorizamos la columna "sample_name"
decont.cnv$sample_name<- factorize.sample(decont.cnv$sample_name)
str(decont.cnv) #Comprobamos como tenemos los datos.


## ¿Cuantos CNV se han detectado en total en cada programa? Tabla
library (viridis) #Para la paleta de colores.
datos_tabla_02<- c(nrow(gatk.cnv), nrow(decont.cnv))
names(datos_tabla_02)<- c("GATK", "XHMM + DECoNT")
colnames(datos_tabla_02)
barplot(datos_tabla_02, col=plasma(2), main="CNVs detected by each program in total", ylab="N of CNVs")
datos_tabla_02


## ¿Cuantos DEL y DUP ha detectado cada programa? Tabla
#Crearemos una matriz con los porcentajes de DEL y DUP para cada programa.
tabla_gatk<- table(gatk.cnv$GT)
tabla_decont<- table(decont.cnv$GT)
datos_tabla_03<- cbind(tabla_gatk,tabla_decont)
colnames(datos_tabla_03)<- c("GATK", "DECoNT")
#Añadimos los totales al final de las lineas.
datos_tabla_03<-rbind(datos_tabla_03, datos_tabla_02)
rownames(datos_tabla_03)[3]<- "Total"
datos_tabla_03

#Ahora cambiamos el formato de los datos a un data frame largo para representarlo mejor.
datos_figura_01<- rbind(data.frame(prop.table(table(gatk.cnv$GT)), group="GATK"), 
                        data.frame(prop.table(table(decont.cnv$GT)), group="XHMM + DECoNT"))
names(datos_figura_01)<- c("GT", "CNVs", "group")
#Factorizamos las variables.
datos_figura_01$GT<- factorize.GT(datos_figura_01$GT)
datos_figura_01$group<- factor(datos_figura_01$group, levels= c("GATK", "XHMM + DECoNT"))

#Dinujamos el gráfico.
library(ragg) #Para un entorno gráfico con acceso directo a todas las fuentes de Windows y con opciones de antialiasing (suavizado). Mas info en: https://ragg.r-lib.org/
library(ggplot2)
library(viridis)
library(hrbrthemes)
#windowsFonts("Arial Narrow"= windowsFont("TT Arial")) #Esto sería para importar esta fuente en concreto solo durante la sesión.
#windowsFonts("Roboto Condensed"= windowsFont("TT Roboto Condensed")) #Esto sería para importar esta fuente en concreto solo durante la sesión.
#hrbrthemes::import_roboto_condensed() #Esto sería para lo mismo pero usando el propio paquete.

figura_01<- ggplot(data=datos_figura_01, aes(x=group, y=CNVs*100, fill=GT)) + 
  geom_bar(position="stack", stat="identity", colour= "black", width=0.6) +
  geom_text(aes(label=round(CNVs*100, digits=2)), position=position_stack(vjust=0.5), size=5) + 
  theme_ipsum() + 
  scale_fill_manual(values=rainbow(2)) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_text(size=10, hjust=0.5),
        axis.text.x=element_text(size= 10, angle= 0)) + 
  ggtitle(label="Genotype prediction % for each programm") + 
  xlab(label="") + 
  ylab(label="%") + 
  labs(fill= "Genotype")
figura_01



## ¿Se han detectado proporciones parecidas de DEL y DUP en cada individuo al usar distintos programas? ESTO VA CON LA PRIMERA DE LAS TABLAS
#Crearemos un data.frame ancho que cuente cuantos CNV, cuantos, DEL y cuantos DUP, ha detectado cada programa en cada muestra.
datos_figura_02<- data.frame(sample_name= muestras, 
                             gatk_cnv=as.numeric(table(gatk.cnv$sample_name)),
                             decont_cnv=as.numeric(table(decont.cnv$sample_name)))
datos_figura_02$gatk_del<- as.matrix(table(gatk.cnv$sample_name, gatk.cnv$GT))[,1]
datos_figura_02$gatk_dup<- as.matrix(table(gatk.cnv$sample_name, gatk.cnv$GT))[,2]
datos_figura_02$decont_del<- as.matrix(table(decont.cnv$sample_name, decont.cnv$GT))[,1]
datos_figura_02$decont_dup<- as.matrix(table(decont.cnv$sample_name, decont.cnv$GT))[,2]
head(datos_figura_02)
tail(datos_figura_02)
#Dibujamos un diagrama de barras en el que mostraremos las diferencias entre % de DEL y DUP para cada individuo.
#Tomaremos como total el número de CNVs detectados en cada individuo.
#Todos estos cálculos se realizarán justo antes de dibujar y no se guardarán como variables.
library(cowplot)

#Diferencias en DEL
figura_02a<- ggplot(data=datos_figura_02, aes(x=sample_name, y=((gatk_del/gatk_cnv) - (decont_del/decont_cnv)) * 100)) + 
  geom_bar(position="dodge", stat="identity", fill="#FF0000", color="black", size=0.2) + 
  ylim(-60,60) + 
  #scale_fill_viridis(discrete=TRUE) + 
  theme_ipsum() + 
  theme(legend.position="none", 
        plot.title=element_text(size=12, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_text(size=10, hjust=0.5),
        axis.text.x=element_blank()) + 
  ggtitle(label="Differencial of DEL % in each sample",
          subtitle="(GATK DELs %) - (XHMM + DECoNT DELs %)") + 
  xlab(label="Samples") + 
  ylab(label="Differential of DEL %")

#Diferencias en DUP
figura_02b<- ggplot(data=datos_figura_02, aes(x=sample_name, y=((gatk_dup/gatk_cnv) - (decont_dup/decont_cnv)) * 100)) + 
  geom_bar(position="dodge", stat="identity", fill="#00FFFF", color="black", size=0.2) +
  ylim(-60,60) + 
  #scale_fill_viridis(discrete=TRUE) + 
  theme_ipsum() + 
  theme(legend.position="none", 
        plot.title=element_text(size=12, hjust=0.5),
        plot.subtitle=element_text(size=10, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_text(size=10, hjust=0.5),
        axis.text.x=element_blank()) + 
  ggtitle(label="Differencial of DUP % in each sample",
          subtitle="(GATK DUPs %) - (XHMM + DECoNT DUPs %)") + 
  xlab(label="Samples") + 
  ylab(label="Differential of DUP %")

#Conclusión: En porcentaje por cada individuo generalmente se detectan más DEL con GATK y más DUP con DECoNT.
figura_02<- plot_grid(figura_02a,figura_02b, labels="AUTO")
figura_02
#Contamos cuantas veces ha aparecido cada numero de copias.
datos_tabla_04<- table(gatk.cnv$CN) #Resulta que con diferencia, GATK está detectando DELs con 0 copias.


## ¿Que tamaños tienen los CNV predichos por cada programa (DEL y DUP)?
#Primero tenemos que calcular los tamaños y añadirlos a los data.frame originales.
gatk.cnv$size<- get.sizes(gatk.cnv)
decont.cnv$size<- get.sizes(decont.cnv)
#Podemos ver que en GATK se encuentran CNV de tamaño máximo y medio mucho mayores que en DECoNT.
summary(gatk.cnv$size)
summary(decont.cnv$size)

#Hay pocas longitudes con valores demasiado altos, no sé pueden representar bien.
plot(density(gatk.cnv$size))
#De modo que haremos una transformación logarítmica (base 10), para representar mejor los datos.
#Pero no se lo aplicaremos a los datos antes de guardarlos, si no justo antes de dibujar el gráfico.
#Creamos el data.frame largo para la representación.
#Los saltos de linea sirven tan solo para que al representarlos aparezcan separados.
datos_figura_03<- rbind(data.frame(size=gatk.cnv$size, GT=gatk.cnv$GT, group=paste("GATK \n", gatk.cnv$GT, sep=" ")),
                        data.frame(size=decont.cnv$size, GT=decont.cnv$GT, group=paste("XHMM + DECoNT \n", decont.cnv$GT, sep=" ")))
head(datos_figura_03)
tail(datos_figura_03)
str(datos_figura_03)

#Ahora realizaremos un boxplot y una diagrama de densidad de las longitudes de los CNVs de cada programa, pero distinguiendo entre DEL o DUP.
library(ggplot2)
library(latex2exp)
#library(dplyr)
library(hrbrthemes)
library(viridis)
#Establecemos el grupo como factor y determinamos su orden.
datos_figura_03$GT<- factor(datos_figura_03$GT, levels= c("DEL", "DUP"))
datos_figura_03$group<- factor(datos_figura_03$group, levels= c("GATK \n DEL", "GATK \n DUP", "XHMM + DECoNT \n DEL", "XHMM + DECoNT \n DUP"))

figura_03<- ggplot(data= datos_figura_03, aes(x=group, y= log10(size), fill=GT)) +
  geom_violin(width=1, size=0.7) + 
  geom_boxplot(width=0.1, color="black", alpha=0.2) + 
  stat_summary(fun=mean, geom="point", shape=21, size=2, stroke= 1, color="black", fill="grey") + 
  #scale_fill_viridis(discrete=TRUE) + 
  theme_ipsum() + 
  theme(legend.position="none",
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_text()) + 
  ggtitle("Distributions of CNV length by genotype and program") + 
  xlab(label="Program and CNV type") + 
  ylab(TeX("$\\log_{10}(bp)$"))
figura_03
#Referencia en: https://r-graph-gallery.com/violin_and_boxplot_ggplot2.html
#El "size" predeterminado del violín es 0.5. Ajustarlo para que quede bien en el documento final.


## ¿Cuantos CNVs exclusivos y en común se han detectado?
#Cargaremos los datos de las intersecciones para responder estas preguntas.

#Exclusivos de GATK
file<- paste(directorio_intersectos, "all_samples/", "GATK_cnv_exclusive.bed", sep="")
gatk.exclusive<- open.bed(file)
#Le ponemos el encabezado.
names(gatk.exclusive)<- header_gatk
#Factorizamos el genotipo predicho.
gatk.exclusive$GT<- factorize.GT(gatk.exclusive$GT)
#Factorizamos la columna "sample_name".
gatk.exclusive$sample_name<- factorize.sample(gatk.exclusive$sample_name)
#Añadimos la columna de la longitud de cada CNV.
gatk.exclusive$size<- get.sizes(gatk.exclusive)
str(gatk.exclusive) #Comprobamos como tenemos los datos.

#Exclusivos de DECoNT
file<- paste(directorio_intersectos, "all_samples/", "DECoNT_cnv_exclusive.bed", sep="")
decont.exclusive<- open.bed(file)
#Le ponemos el encabezado.
names(decont.exclusive)<- header_decont
#Factorizamos el genotipo predicho.
decont.exclusive$GT<- factorize.GT(decont.exclusive$GT)
#Factorizamos la columna "sample_name".
decont.exclusive$sample_name<- factorize.sample(decont.exclusive$sample_name)
#Añadimos la columna de la longitud de cada CNV.
decont.exclusive$size<- get.sizes(decont.exclusive)
str(decont.exclusive) #Comprobamos como tenemos los datos.

#En común por ambos.
file<- paste(directorio_intersectos, "all_samples/", "cnv_common.bed", sep="")
cnv.common<- open.bed(file)
#Le ponemos el encabezado.
names(cnv.common)<- c("chrom", "start","end", "sample_name", "GT", "merged_cnv")
#Esta vez no factorizamos el genotipo predicho, ya que hay muchos concatenados en una sola linea.
#Factorizamos la columna "sample_name".
cnv.common$sample_name<- factorize.sample(cnv.common$sample_name)
#Añadimos la columna de la longitud de cada CNV.
cnv.common$size<- get.sizes(cnv.common)
str(cnv.common) #Comprobamos como tenemos los datos.

#Usaremos los datos de la tabla 1 que contiene el número total de CNVs.
datos_tabla_02
datos_tabla_05<- matrix(data=c(datos_tabla_02[1] - nrow(gatk.exclusive), datos_tabla_02[2] - nrow(decont.exclusive),
                                nrow(gatk.exclusive), nrow(decont.exclusive),
                               datos_tabla_02),
                         ncol=2, byrow=TRUE)
colnames(datos_tabla_05)<- c("GATK", "DECoNT")
rownames(datos_tabla_05)<- c("Overlapping CNV", "Exclusive CNV", "Total CNV")
datos_tabla_05


## ¿Una vez consolidados los CNVs en común de ambos programas, cuantos han quedado?
datos_tabla_06<- nrow(cnv.common)
names(datos_tabla_06)<- c("Sum of the total of CNVs consolidated by saple")
datos_tabla_06 #Estos van a ser los que usemos con X-CNV.


## ¿Qué tamaño tienen los CNVs detectados por ambos programas? Teniendo en cuenta si son exclusivos, en común o en común consolidados.
#Definiremos una función para encontrar los CNVs en común de ambos programas, partiendo de un data.frame con todos los CNVs y de otro solo con los CNVs exclusivos.
find.common<- function(all, exclusives){
  ind_exclusives<- c() #Aquí guardaremos los números de línea.
  for(i in 1:nrow(exclusives)){
    ind_exclusives[i]<- which(all$ID_0 == exclusives$ID_0[i])
  }
  common<- all[-ind_exclusives,]
  return(common)
}

#Crearemos un data.frame largo con todos los datos de las longitudes. Incluiremos otra columna especificando el tipo de CNV (exlcusivo etc.) y otra que los agrupe por el programa.
datos_figura_04<- rbind(data.frame(size=gatk.exclusive$size, type="exclusive", group="GATK \n exclusive"),
                        data.frame(size=decont.exclusive$size, type="exclusive", group="XHMM + DECoNT \n exclusive"),
                        data.frame(size=find.common(gatk.cnv, gatk.exclusive)$size, type="overlapping", group="GATK \n overlapping"),
                        data.frame(size=find.common(decont.cnv, decont.exclusive)$size, type="overlapping", group="XHMM + DECoNT \n overlapping"),
                        data.frame(size=cnv.common$size, type="overlapping consolidated", group="Overlapping \n consolidated"))
#Factorizamos las variables del tipo y del grupo.
datos_figura_04$type<- factor(datos_figura_04$type, levels=c("exclusive", "overlapping", "overlapping consolidated"))
datos_figura_04$group<- factor(datos_figura_04$group, levels=c("GATK \n exclusive", "GATK \n overlapping", "Overlapping \n consolidated", "XHMM + DECoNT \n overlapping", "XHMM + DECoNT \n exclusive"))

figura_04<- ggplot(data= datos_figura_04, aes(x=group, y=log10(size), fill=type)) +
  geom_violin(width=1, size=0.7) + 
  geom_boxplot(width=0.1, color="grey", alpha=0.2) + 
  stat_summary(fun=mean, geom="point", shape=21, size=2, stroke=1, color="black", fill="grey") + 
  scale_fill_viridis(discrete=TRUE) +
  theme_ipsum() + 
  theme(legend.position="none", 
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10)) + 
  ggtitle("Distributions of CNV lengths by exclusivity and program") + 
  xlab(label="Program CNV and type of CNV") + 
  ylab(TeX("$\\log_{10}(bp)$"))
figura_04
#De aquí se puede reforzar al conclusión de que los CNV de GATK son más pequeños y a la hora de consolidarlos simplemente se unen a los grandes de DECoNT.
#Los GATK common quizá tiendan a a ser un poco más grandes que los GATK exclusive, pero poco.


## ¿De los CNVs en común consolidados cuantos coinciden en el genotipo?
#Para contestar esta pregunta tenemos que añadir información al archivo con todos los CNVs en común y consolidados.

#Añadimos una columna que establezca si ha habido o no acuerdo entre los dos programas.
cnv.common$consensus<- as.numeric(
  lapply(X= lapply(X= strsplit(x=cnv.common$GT, split=","), FUN= unique), FUN= length))
#Añadimos otra columna con el genotipo acordado por ambos programas.
cnv.common$consensus_GT<- factor(lapply(X= strsplit(x=cnv.common$GT, split=","), FUN= unique),
                                 levels=c("DEL", "DUP"))

#Factorizamos la columna del consenso.
cnv.common$consensus<- factor(cnv.common$consensus, levels=c(1,2), labels=c("Agreement", "Disagreement"))
#Comprobamos los resultados
head(cnv.common)

#Creamos una matriz que nos resuma cuantos CNVs en común y consolidados hay en desacuerdo o no.
datos_tabla_07<- as.matrix(summary(cnv.common$consensus_GT))
datos_tabla_07<- rbind(datos_tabla_07, sum(datos_tabla_07))
rownames(datos_tabla_07)<- c("DEL", "DUP", "Disagreement", "total")
colnames(datos_tabla_07)<- "CNVs in common consolidated"
datos_tabla_07

#Ahora realizaremos una representación a modo de diagrama de barras apiladas.
#Creamos el data.frame largo con los datos.
datos_figura_05<- data.frame(CNVs=summary(cnv.common$consensus_GT))
datos_figura_05$consensus<- c("Agreement", "Agreement", "Disagreement")
datos_figura_05$GT<- c("DEL", "DUP", "ANY")
#Factorizamos variables.
datos_figura_05$consensus<- factor(datos_figura_05$consensus, levels=c("Agreement", "Disagreement"))
datos_figura_05$GT<- factor(datos_figura_05$GT, levels=c("DEL", "DUP", "ANY"), labels=c("DELs", "DUPs", "Discordants"))
str(datos_figura_05)
datos_figura_05

figura_05<- ggplot(data=datos_figura_05, aes(x=consensus, y=CNVs, fill=GT)) + 
  geom_bar(position="stack", stat="identity", color="black", width=0.6) +
  geom_text(aes(label=CNVs), position=position_stack(vjust=0.5), size=5) + 
  ylim(0, sum(datos_figura_05$CNVs)) + 
  theme_ipsum() + 
  scale_fill_manual(values= c("#FF0000", "#00FFFF", "#808B96")) +
  theme(plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size= 10, angle= 0)) + 
  ggtitle(label="Genotypes of common consolidated CNVs") + 
  xlab(label="") + 
  ylab(label="N of CNV") + 
  labs(fill= "Genotypes")
figura_05



# 3.- GTAK vs DECoNT, CNV consolidados por exoma.
#En esta apartado analizaremos las regiones de CNV tras haber consolidado los CNVs entre todas las muestras para cada uno de los dos programas.


## ¿Cuantos CNV consolidados por exoma tiene cada programa y de ellos cuantos son exclusivos de cada programa o comunes entre programas?
#Reutilizaremos la matriz de la figura 4 para crear una tabla similar.
datos_tabla_08<- datos_tabla_05

#Vamos abriendo los archivos para contar el número de lineas (CNVs) que tienen.
file<- paste(directorio_CNV, "GATK_by_exome/GATK_cnv_exome.bed", sep="")
datos_tabla_08[3,1]<- nrow(open.bed(file))
file<- paste(directorio_CNV, "DECoNT_by_exome/DECoNT_cnv_exome.bed", sep="")
datos_tabla_08[3,2]<- nrow(open.bed(file))
file<- paste(directorio_intersectos, "by_exome/GATK_cnv_exome_exclusive.bed", sep="")
datos_tabla_08[2,1]<- nrow(open.bed(file))
file<- paste(directorio_intersectos, "by_exome/DECoNT_cnv_exome_exclusive.bed", sep="")
datos_tabla_08[2,2]<- nrow(open.bed(file))
datos_tabla_08[1,1]<- datos_tabla_08[3,1] - datos_tabla_08[2,1]
datos_tabla_08[1,2]<- datos_tabla_08[3,2] - datos_tabla_08[2,2]
datos_tabla_08

#Una vez consolidados los CNV exomicos en común obtenemos el siguiente número de CNVs únicos por exoma detectados por ambos programas.
file<- paste(directorio_intersectos, "by_exome/GATKandDECoNT_cnv_exome_common.bed", sep="")
datos_tabla_09<- nrow(open.bed(file))
names(datos_tabla_09)<- "CNVs by exome in common and consolidated"
datos_tabla_09


## Visualización de los intervalos originales del Manifiesto, de los intervalos por exoma de GATK, intervalos por exoma de DECoNT e intervalos consolidados por exoma en común.
#Ver esto: http://compgenomr.github.io/book/visualizing-and-summarizing-genomic-intervals.html
library(GenomicRanges)
library(karyoploteR)

#Abrimos el archivo del manifiesto.
file<- "Datos/Manifiesto_2_no_chrXYM.bed"
manifiesto<- open.bed(file)
head(manifiesto)

#Abrimos los CNV consolidados por exoma de GATK.
file<- paste(directorio_CNV, "GATK_by_exome/GATK_cnv_exome.bed", sep="")
gatk.exome<- open.bed(file)[,1:3]
head(gatk.exome)

#Abrimos los CNV consolidados por exoma de DECoNT.
file<- paste(directorio_CNV, "DECoNT_by_exome/DECoNT_cnv_exome.bed", sep="")
decont.exome<- open.bed(file)[,1:3]
head(decont.exome)

#Abrimos los CNV consolidados por exoma en común para ambos programas y consolidados también entre programas.
file<- paste(directorio_intersectos, "by_exome/GATKandDECoNT_cnv_exome_common.bed", sep="")
exome.common<- open.bed(file)[,1:3]
head(exome.common)


#Cariograma del manifiesto.
kp_m<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal", main="Kariogram of targeted exomes")
kpPlotRegions(kp_m, data=toGRanges(manifiesto), col= "black")
#Cariograma de las regiones CNV de GATK.
kp_g<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal", main="Kariogram of GATK CNVs")
kpPlotRegions(kp_g, data=toGRanges(gatk.exome), col= "#BC1DBF")
#Cariograma de las regiones CNV de DECoNT.
kp_d<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal", main="Kariogram of XHMM + DECoNT CNVs")
kpPlotRegions(kp_d, data=toGRanges(decont.exome), col= "#FFC300")
#Cariograma de las regiones CNV en común consolidadas para ambos programas.
kp_c<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal", main="Kariogram of consolidated CNVs")
kpPlotRegions(kp_c, data=toGRanges(exome.common), col= "#2C9A2B")
#NOTA: los ideogramas se descargan de internet y se pueden quitar si establacemos ideogram.plotter=NULL


## ¿Tienen estos CNVs consolidados por exoma el mismo tamaño que los no consolidados?
#Abrimos los datos de los CNVs exclusivos por exoma.

#Exoma exclusivos de GATK.
file<- paste(directorio_intersectos, "by_exome/GATK_cnv_exome_exclusive.bed", sep="")
gatk.exome_exclusive<- open.bed(file)[,1:3]
head(gatk.exome_exclusive)

#Exoma exclsuivos de DECoNT.
file<- paste(directorio_intersectos, "by_exome/DECoNT_cnv_exome_exclusive.bed", sep="")
decont.exome_exclusive<- open.bed(file)[,1:3]
head(decont.exome_exclusive)

#Les ponemos los encabezados a todos los archivos.
names(exome.common)=names(gatk.exome_exclusive)=names(decont.exome_exclusive)=names(gatk.exome)=names(decont.exome)=header_gatk[1:3]
#Crearemos IDs para todos los CNVs, para poder identificarlos y buscarlos. Lo haremos mediante la siguiente función.
get.ID_0<- function(data_frame){
  ID<- paste(data_frame$chrom, data_frame$start, data_frame$end, sep="_")
  return(ID)
}
gatk.exome$ID_0<- get.ID_0(gatk.exome)
gatk.exome_exclusive$ID_0<- get.ID_0(gatk.exome_exclusive)
decont.exome$ID_0<- get.ID_0(decont.exome)
decont.exome_exclusive$ID_0<- get.ID_0(decont.exome_exclusive)
exome.common$ID_0<- get.ID_0(exome.common)
#Calculamos los tamaños de los CNVs para los 5 archivos de exoma.
gatk.exome$size<- get.sizes(gatk.exome)
gatk.exome_exclusive$size<- get.sizes(gatk.exome_exclusive)
decont.exome$size<- get.sizes(decont.exome)
decont.exome_exclusive$size<- get.sizes(decont.exome_exclusive)
exome.common$size<- get.sizes(exome.common)

#Crearemos un data.frame largo con todos los datos de las longitudes. Incluiremos otra columna especificando el tipo de CNV (exlcusivo etc.) y otra que los agrupe por el programa.
datos_figura_07<- rbind(data.frame(size=gatk.exome_exclusive$size, type="exclusive", group="GATK \n exclusive"),
                        data.frame(size=decont.exome_exclusive$size, type="exclusive", group="XHMM + DECoNT \n exclusive"),
                        data.frame(size=find.common(gatk.exome, gatk.exome_exclusive)$size, type="overlapping", group="GATK \n overlapping"),
                        data.frame(size=find.common(decont.exome, decont.exome_exclusive)$size, type="overlapping", group="XHMM + DECoNT \n overlapping"),
                        data.frame(size=exome.common$size, type="overlapping consolidated", group="Overlapping \n consolidated"))
#Factorizamos las variables del tipo y del grupo.
datos_figura_07$type<- factor(datos_figura_07$type, levels=c("exclusive", "overlapping", "overlapping consolidated"))
datos_figura_07$group<- factor(datos_figura_07$group, levels=c("GATK \n exclusive", "GATK \n overlapping", "Overlapping \n consolidated", "XHMM + DECoNT \n overlapping", "XHMM + DECoNT \n exclusive"))

figura_07<- ggplot(data= datos_figura_07, aes(x=group, y=log10(size), fill=type)) +
  geom_violin(width=1, size=0.7) + 
  geom_boxplot(width=0.1, color="grey", alpha=0.2) + 
  stat_summary(fun=mean, geom="point", shape=21, size=2, stroke=1, color="black", fill="grey") + 
  scale_fill_viridis(discrete=TRUE) +
  theme_ipsum() + 
  theme(legend.position="none", 
        plot.title=element_text(size=14, hjust=0.5, vjust=5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_text(size=10)) + 
  ggtitle("Distributions of consolidated CNV lengths, by exclusivity and program") + 
  xlab(label="Program and type of CNV") + 
  ylab(TeX("$\\log_{10}(bp)$"))
figura_07
#De aquí podemos comentar que parece que la consolidación por exoma hace juntarse varios CNV pequeños y formar CNV, mas grandes.
#La mayoría de CNV en común son grandes ahora. Y como exlcusivos quedan los más pequeños de cáda exoma de cada programa.
#Los consolidados por exoma, en común y consolidados entre programas (amarillo), parece que tienen una distribución mut parecida a la de los en común y consolidados de antes con todos los individuos).



# 4.- Exportar las figuras y tablas.
#Fuente: https://ggplot2.tidyverse.org/reference/ggsave.html

#A la hora de exportar a documentos PDF o imágenes conviene seleccionar un dispositivo gráfico.
#En MacOS el dispositivo predeterminado soporta antialiasing, para unos bordes suavizados, pero en Windows no.
#IMPORTANTE SEGUIR ESTOS PASOS A MANO PARA QUE LAS FUENTES SE USEN CORRECTAMENTE Y LOS BORDES ESTEN SUAVIZADOS
#Dentro de RStudio ir a Tools --> Global Options --> General --> Graphics
#Seleccionar "AGG" como Backend y "Subpixel" como Antialiasing. El backend "Cairo" también soporta antialiasing pero creo que no da acceso a todas las fuentes de Windows como "AGG".
#A continuación seguimos los siguientes pasos para guardar las imágenes en formato PNG.

#Establecemos el directorio de salida.
directorio_out<- "Resultados/" #Definimos el direcotorio de salida.
#También exportaremos tablas a formato Excel, por lo que necesitaremos el siguiente paquete.
library(xlsx)

write.xlsx(datos_tabla_01, 
           file=paste(directorio_out, "tabla_01.xlsx", sep=""))
write.xlsx(datos_tabla_02,
           file=paste(directorio_out, "tabla_02.xlsx", sep=""))
write.xlsx(datos_tabla_03,
           file=paste(directorio_out, "tabla_03.xlsx", sep=""))
write.xlsx(datos_tabla_04,
           file=paste(directorio_out, "tabla_04.xlsx", sep=""))
write.xlsx(datos_tabla_05,
           file=paste(directorio_out, "tabla_05.xlsx", sep=""))
write.xlsx(datos_tabla_06,
           file=paste(directorio_out, "tabla_06.xlsx", sep=""))
write.xlsx(datos_tabla_07,
           file=paste(directorio_out, "tabla_07.xlsx", sep=""))
write.xlsx(datos_tabla_08,
           file=paste(directorio_out, "tabla_08.xlsx", sep=""))
write.xlsx(datos_tabla_09,
           file=paste(directorio_out, "tabla_09.xlsx", sep=""))


#Notese que le estamos especificando el tamaño de la imagen en milímetros para que ocupe todo el horizontal de un Din A4.
#Además le ponemos alta densidad de pixeles, adecuada para impresión.
ggsave(figura_01, 
       filename= paste("figura_01", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
ggsave(figura_02, 
       filename= paste("figura_02", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
ggsave(figura_03, 
       filename= paste("figura_03", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
ggsave(figura_04, 
       filename= paste("figura_04", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
ggsave(figura_05, 
       filename= paste("figura_05", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
#Los siguientes son los cariogramas. Están hechos con R base y no consigo editar cosas como lo margenes.
#Además no consigo que no se solapen con cowplot así que los guardare independientemente para luego juntarlos manualmente.
#Notar que los estoy rederizando igual de altos que siempre pero a un 75% de anchura respecto a los anteriores.
#Cariograma del manifiesto.
kp_m<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal")
kpPlotRegions(kp_m, data=toGRanges(manifiesto), col= "black")
KP<- recordPlot()
#Cariograma de las regiones CNV de GATK.
kp_g<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal")
kpPlotRegions(kp_g, data=toGRanges(gatk.exome), col= "#BC1DBF")
KG<- recordPlot()
#Cariograma de las regiones CNV de DECoNT.
kp_d<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal")
kpPlotRegions(kp_d, data=toGRanges(decont.exome), col= "#FFC300")
KD<- recordPlot()
#Cariograma de las regiones CNV en común consolidadas para ambos programas.
kp_c<- plotKaryotype(genome="hg19", plot.type=1, chromosomes="autosomal")
kpPlotRegions(kp_c, data=toGRanges(exome.common), col= "#2C9A2B")
KC<- recordPlot()
ggsave(plot_grid(KP, labels="A"), 
       filename= paste("figura_06A", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210*0.75, height=118, units= "mm")
ggsave(plot_grid(KG, labels="B"), 
       filename= paste("figura_06B", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210*0.75, height=118, units= "mm")
ggsave(plot_grid(KD, labels="C"), 
       filename= paste("figura_06C", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210*0.75, height=118, units= "mm")
ggsave(plot_grid(KC, labels="D"), 
       filename= paste("figura_06D", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210*0.75, height=118, units= "mm")
ggsave(figura_07, 
       filename= paste("figura_07", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")

sessionInfo()