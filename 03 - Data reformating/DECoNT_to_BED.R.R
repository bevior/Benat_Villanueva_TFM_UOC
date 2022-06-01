# Paso 1: Pasar los resultados de DECoNT a un formato compatible con BED.
#Primero importamos los datos.
#Definimos la ruta de antemano.
ruta.input<- "/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/laboratorio/XHMM_gCNV/Paso10_DECoNT_polishing/DECoNT_XHMM_polished_cnvs.txt" #SERVIDOR
#ruta.input<- "DECoNT_XHMM_polished_cnvs.txt" #PC PERSONAL
decont.data<- read.table(ruta.input, header=TRUE, sep="\t")

#Hacemos una exploración rápida de los datos leídos.
head(decont.data)
class(decont.data)
str(decont.data)

#https://bedtools.readthedocs.io/en/latest/content/overview.html#:~:text=The%20VCF%20format%20uses%201,to%20work%20correctly%20with%20bedtools.
#Las coordenadas se corresponden con las del formato VCF, por lo que tanto el inicio como el final son en base a 1.
#El estándar BED requiere que el inicio sea en base a 0 y el final en base a 1, por lo que tendremos que cambiar el primero.
#De esta manera para caluclar el tamaño de un intervalo solo hay que hacer end - start.
start<- decont.data$CNV.Start.Index - 1
end<- decont.data$CNV.End.Index
head(cbind(start, end))
table(end > start) #Comprobamos que el final siempre es mayor que el inicio.

#Creamos una variable extra y la llamaremos "ID_0". Estará formada por la cromosoma, el inicio y el final. Todo en base 0.
ID_0<- paste("DECoNT", decont.data$Sample.Name, decont.data$Chromosome, start, end, sep="_")
head(ID_0)
tail(ID_0)

#Creamos otra variable y la llamaremos "ID_1". Estará formada por la cromosoma, el inicio y el final de las coordenadas originales. Base 1.
ID_1<- paste("DECoNT", decont.data$Sample.Name, decont.data$Chromosome, decont.data$CNV.Start.Index, decont.data$CNV.End.Index, sep="_")
head(ID_1)
tail(ID_1)

#Ahora creamos un nuevo data.frame en el que juntaremos todos estos datos en un orden compatible con el formato BED.
decont.bed<- data.frame(chrom=decont.data$Chromosome, start, end, ID_0, ID_1, decont.data[,c(1,5,6)])
names(decont.bed)[6:8]<- c("sample_name", "XHMM", "DECoNT") #Cambiamos ligeramente algunos nombres de columna.
head(decont.bed)
tail(decont.bed) #Comprobamos que todo este correcto.

#Filtraremos los CNV que sean "NO-CALL" seleccionando las lineas que no coincidan con esa categoría de CNV.
decont.bed<- decont.bed[decont.bed$DECoNT!="NO-CALL",]

#Ahora establecemos la ruta de salida.
ruta.output<- "Paso1_DECoNT_to_BED/DECoNT_cnv.bed" #SERVIDOR
#ruta.output<- "Results_from_DECoNT/DECoNT_cnv.bed" #PC PERSONAL
write.table(decont.bed, file=ruta.output, sep="\t", row.names=FALSE, col.names=FALSE, quote=FALSE, eol="\n") #Aquí hay que tener cuidado que al parecer en RStudio sigue guardandolo con final de linea "\r\n".
#Se puede reconvertir fácilmente con Notepad++. En GNU/Linux no creo que tengamos ese problema.
sessionInfo()