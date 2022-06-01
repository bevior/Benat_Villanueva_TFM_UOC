# Paso 4: Pasar los resultados de GATK a un formato compatible con BED.
#Primero necesitamos la lista de los nombres de cada muestra. Accederemos al archivo DECoNT_cnv.bed para obtenerlos.
#Guardamos laruta al archivo.
archivo<- "Paso1_DECoNT_to_BED/DECoNT_cnv.bed" #SERVIDOR
#archivo<- "Results_from_DECoNT/DECoNT_cnv.bed" #PC PERSONAL
DECoNT<- read.table(archivo, header=FALSE, sep="\t")
names(DECoNT)<- c("chrom", "start", "end", "ID_0", "ID_1", "sample_name", "XHMM", "DECoNT") #Tendremos que renombrar las columnas manualmente.
str(DECoNT) #Comprobamos que la estructura del archivo es la correcta.
#Los nombres de las muestras están en la columna sample_name

DECoNT$sample_name
muestras<- sort(unique(DECoNT$sample_name))
muestras #Ya tenemos los nombres de las muestras para iterar sobre ellos.
for (i in muestras) {
  print(i)
}

#Guardamos el directorio donde se encuentran los archivos.
directorio<-"/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/laboratorio/GATK_gCNV/Paso5_PostprocessGermlineCNVCalls/" #SERVIDOR
#directorio<- "../04 - CNV calling resultados/Paso5_PostprocessGermlineCNVCalls/" #PC PERSONAL

#Abrimos un paquete para poder leer archivos VCF más fácilmente.
library(vcfR)

#Ahora abriremos en bucle los archivos VCF de GATK y extraeremos la información necesaria para crear un archivo BED similar a los que ya tenemos para los datos de DECoNT.
#Recordemos que los VCF están en base 1 tanto el inicio como el final. Tendremos que pasar el inicio a base 0.
#Debemos tener en cuenta lo siguiente a la hora de leer la predicción para los CNVs:
# 1- El genotipo se encuentra bajo el nombre GT
# 2- El numero de copias estimado se encuentra bajo el nombre CN
#Usaremos directamente el CN para trabajar.

for (i in muestras) {
  archivo<- paste(i, ".gCNV_segments.vcf", sep="")
  vcf<- read.vcfR(file=paste(directorio, archivo, sep="")) #Leemos el archivo
  #Iremos guardando las variables por separado.
  chrom<- vcf@fix[,1]
  start<- as.numeric(vcf@fix[,2]) #Hay que pasarlo a numérico.
  start<- start-1 #Lo pasamos a base 0.
  end<- INFO2df(vcf)[,4] #Se encuentra en un ligar especial y ya en formato numérico.
  ID_0<- paste("GATK", i, chrom, start, end, sep="_")
  ID_1<- paste("GATK", i, chrom, start+1, end, sep="_")
  sample_name<- rep(i, times=length(start)) #Repetimos el nombre de la muestra el numero de veces que haga falta, hasta rellenar todas las lineas.
  GATK<- extract.gt(x=vcf, element="CN", as.numeric=TRUE, convertNA= FALSE) #Información como el CN o el GT se extrae de esta manera.
  rownames(GATK)<- NULL #Para borrar los nombres de columnas.
  GATK_GT<- ifelse(GATK<2, "DEL", ifelse(GATK==2, "DIP", ifelse(GATK>2, "DUP", "ERROR"))) #Le pedimos que nos clasifique como deleccion, diploide o duplicación según el número de copias en cada predicción.
  gatk.bed<- data.frame(chrom, start, end, ID_0, ID_1, sample_name, GATK, GATK_GT)
  names(gatk.bed)[7:8]<- c("GATK", "GATK_GT") #Renombramos algunas columnas.
  gatk.bed<- gatk.bed[gatk.bed$GATK_GT!="DIP",] #Eliminaremos las filas que correspondan a CNV diploides.
  head(gatk.bed) #Comprobaciones finales
  tail(gatk.bed)
  
  ruta.output<- paste("Paso4_GATK_to_BED/", i, "_GATK_cnv.bed", sep="") #SERVIDOR
  #ruta.output<- paste("some_output/Paso4_GATK_to_BED/", i, "_GATK_cnv.bed", sep="") #PC PERSONAL
  
  #Guardamos los archivos en formato BED, delimitado por tabulaciones. En esta ocasión no guardaremos los nombres de columna.
  write.table(gatk.bed, file=ruta.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\n")
  print(paste("Se ha procesado la muestra", i, sep=" "))
}







## Este va a ser un ejemplo de lo que hay que hacer con cada archivo. Pero solo lo haremos con uno.
i<- "AFA016"
archivo<- paste(i, ".gCNV_segments.vcf", sep="")
vcf<- read.vcfR(file=paste(directorio, archivo, sep="")) #Leemos el archivo
#Iremos guardando las variables por separado.
chrom<- vcf@fix[,1]
start<- as.numeric(vcf@fix[,2]) #Hay que pasarlo a numérico.
start<- start-1 #Lo pasamos a base 0.
end<- INFO2df(vcf)[,4] #Se encuentra en un ligar especial y ya en formato numérico.
ID_0<- paste("GATK", i, chrom, start, end, sep="_")
ID_1<- paste("GATK", i, chrom, start+1, end, sep="_")
sample_name<- rep(i, times=length(start)) #Repetimos el nombre de la muestra el numero de veces que haga falta, hasta rellenar todas las lineas.
GATK<- extract.gt(x=vcf, element="CN", as.numeric=TRUE, convertNA= FALSE) #Información como el CN o el GT se extrae de esta manera.
rownames(GATK)<- NULL #Para borrar los nombres de columnas.
GATK_GT<- ifelse(GATK<2, "DEL", ifelse(GATK==2, "DIP", ifelse(GATK>2, "DUP", "ERROR"))) #Le pedimos que nos clasifique como deleccion, diploide o duplicación según el número de copias en cada predicción.
gatk.bed<- data.frame(chrom, start, end, ID_0, ID_1, sample_name, GATK, GATK_GT)
names(gatk.bed)[7:8]<- c("GATK", "GATK_GT") #Renombramos algunas columnas.
gatk.bed<- gatk.bed[gatk.bed$GATK_GT!="DIP",] #Eliminaremos las filas que correspondan a CNV diploides.
head(gatk.bed) #Comprobaciones finales
tail(gatk.bed)

ruta.output<- paste("Paso4_GATK_to_BED/", i, "_GATK_cnv.bed", sep="") #SERVIDOR
#ruta.output<- paste("some_output/Paso4_GATK_to_BED/", i, "_GATK_cnv.bed", sep="") #PC PERSONAL

#Guardamos los archivos en formato BED, delimitado por tabulaciones. En esta ocasión no guardaremos los nobres de clumna.
write.table(gatk.bed, file=ruta.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\n")
sessionInfo()