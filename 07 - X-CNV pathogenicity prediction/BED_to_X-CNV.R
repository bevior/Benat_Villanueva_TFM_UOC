######################################################################################
#Paso 1: Reformatear los datos.


#Abrimos los CNVs de todos los individuos pero que sean en común entre ambos programas y que hayan sido consolidados entre ambos programas.
#file<- "../08 - Analizar CNV e intersecciones/Datos/intersecciones/all_samples/cnv_common.bed" #PC PERSONAL
file<- "../intersecciones/all_samples/cnv_common.bed" #SERVIDOR
cnv.common<- read.table(file, sep="\t",
                        colClasses=c("character", "integer", "integer", "character", "character", "integer"))
#Le ponemos el encabezado.
names(cnv.common)<- c("chrom", "start","end", "sample_name", "GT", "merged_cnv")
#Guardamos una variable con el nombre de las muestras.
muestras<- sort(unique(cnv.common$sample_name))
#Factorizamos la columna "sample_name".
cnv.common$sample_name<- factor(cnv.common$sample_name, levels= muestras)

#Añadimos una columna que establezca si ha habido o no acuerdo entre los dos programas.
cnv.common$consensus<- as.numeric(
  lapply(X= lapply(X= strsplit(x=cnv.common$GT, split=","), FUN= unique), FUN= length))
#Añadimos otra columna con el genotipo acordado por ambos programas. Lo pondremos de manera que X-CNV lo reconozca.
cnv.common$consensus_GT<- factor(lapply(X= strsplit(x=cnv.common$GT, split=","), FUN= unique),
                                 levels=c("DEL", "DUP"), labels=c("loss", "gain"))
#Factorizamos la columna del consenso.
cnv.common$consensus<- factor(cnv.common$consensus, levels=c(1,2), labels=c("agreement", "disagreement"))
#Quitamos el prefijo de "chr" al numero de chromosoma.
cnv.common$chrom<- gsub("chr", "", cnv.common$chrom)
#Comprobamos los resultados
head(cnv.common)

#Ahora seleccionaremos solo las filas en las que haya consenso.
cnv.common_agree<- cnv.common[cnv.common$consensus=="agreement",]

#Por último crearemos la carpeta de salida y un bucle que nos vaya guardado los CNVs de cada individuo en un archivo separado.
#SERVIDOR
for(i in muestras){
  #Notese que estamos seleccionando solo las filas de esa muestra y solo sus columnas con información sobre la ubicación y el genotipo.
  output<- cnv.common_agree[cnv.common_agree$sample_name==i,][,c(1,2,3,8)]
  file_name<- paste("Paso1_data_format/",i, "_cnv_common.bed", sep="")
  write.table(output, file=file_name, quote= FALSE, sep="\t", row.names= FALSE, col.names= FALSE, eol="\n")
}



#PC PERSONAL
if (file.exists("Resultados")) {
  cat("The folder already exists")
} else {
  dir.create("Resultados")
}

for(i in muestras){
  #Notese que estamos seleccionando solo las filas de esa muestra y solo sus columnas con información sobre la ubicación y el genotipo.
  output<- cnv.common_agree[cnv.common_agree$sample_name==i,][,c(1,2,3,8)]
  file_name<- paste("Resultados/", i, "_cnv_common.bed", sep="")
  write.table(output, file=file_name, quote= FALSE, sep="\t", row.names= FALSE, col.names= FALSE, eol="\n")
}

sessionInfo()