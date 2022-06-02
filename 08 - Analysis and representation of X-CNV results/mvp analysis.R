# Introducción
#En este script leeremos los datos de X-CNV, clasificaremos los MVP según las categorías del artículo y mostraremos cuantos de cada categoría tiene cada muestra.



# 1.- Preparar los datos.


## Leer los datos.
#Los datos están en formato CSV, divididos por individuo y tienen 53 columnas.
#Crearemos un bucle que itere sobre los nombres de las muestras para abrir los archivos, extraer solo las columnas necesarias y añadirle el nombre de muestra al mismo tiempo.

#Primero obtenemos los nombres de las 100 muestras.
file<- "cnv_common.bed"
muestras<- read.table(file, sep="\t")[,4]
muestras<- sort(unique(muestras))

#Ahora creamos el bucle.
directorio<- "X-CNV_phenotiping/Paso2_phenotyping/"
contenedor<- list() #Creamos una lista vacía para almacenar todos los archivos.
cnv_common_pt<- data.frame() #Luego guardaremos los datos de todas las muestras en este data.frame.
phenotypes<- c("pathogenic", "likely pathogenic", "uncertain", "likely bening", "bening") #Guardaremos los nombres de los fenotipos en esta variable.

for(i in 1:length(muestras)){
  file<- paste(directorio, muestras[i], "_cnv_common.output.csv", sep="")
  contenedor[[i]]<- read.csv(file)[,c(1,2,3,4,53)]
  contenedor[[i]]<- cbind(contenedor[[i]][,1:4], muestras[i], contenedor[[i]][,5])
  names(contenedor[[i]])<- c("chrom", "start", "end", "GT", "sample_name", "MVP")
  contenedor[[i]]$chrom<- paste("chr", contenedor[[i]]$chrom, sep="")
  contenedor[[i]]$GT<- factor(contenedor[[i]]$GT, levels=c("loss", "gain"), labels=c("DEL", "DUP"))
  contenedor[[i]]$sample_name<- factor(contenedor[[i]]$sample_name, levels= muestras)
  contenedor[[i]]$PT<- ifelse(contenedor[[i]]$MVP>0.76, phenotypes[1],
                              ifelse(contenedor[[i]]$MVP>0.46, phenotypes[2],
                                     ifelse(contenedor[[i]]$MVP>0.16, phenotypes[3],
                                            ifelse(contenedor[[i]]$MVP>0.14, phenotypes[4],
                                                   ifelse(contenedor[[i]]$MVP<=0.14, phenotypes[5], "ERROR")))))
  contenedor[[i]]$PT<- factor(contenedor[[i]]$PT, levels= phenotypes)
  cnv_common_pt<- rbind(cnv_common_pt, contenedor[[i]])
}

#Comprobamos que el data.frame resultante tiene la estructura adecuada.
head(cnv_common_pt)
tail(cnv_common_pt)
str(cnv_common_pt)



# 2.- Analizar los datos.


## Tablas.
#Cotaremos cuantos fenotipos hay de cada tipo.
datos_tabla_09<- t(as.matrix(table(cnv_common_pt$PT)))
rownames(datos_tabla_09)<- "CNVs"
datos_tabla_09<- cbind(datos_tabla_09, sum(datos_tabla_09))
colnames(datos_tabla_09)[6]<- "Total"
datos_tabla_09


## ¿Cuantos CNV de cada tipo tiene cada individuo?.
#Creamos un data.frame largo que contenga los datos a representar. Contaremos cuantos CNVs de cada tipo tiene cada individuo.
datos_figura_08<- data.frame(sample_name=cnv_common_pt$sample_name,
                        MVP=cnv_common_pt$MVP,
                        PT=cnv_common_pt$PT)

datos_figura_08<- data.frame(table(datos_figura_08$sample_name, datos_figura_08$PT))
names(datos_figura_08)<- c("sample_name", "PT", "CNV")
head(datos_figura_08)
#Creremos una paleta de colores personalizada.
colores<- c("#FF0000", "#FF6C00", "#87918B", "#3E9F75", "#1E5BFF")
#Un barplot para enseñar cuantos CNV de cada tipo tiene cada individuo.
library(ggplot2)
library(ragg)
library(hrbrthemes)
library(viridis)

figura_08<- ggplot(data= datos_figura_08, aes(x=sample_name, y=CNV, fill=PT)) +
  geom_bar(position="stack", stat="identity", color="black") + 
  #ylim(0, max(30)) + 
  theme_ipsum() + 
  scale_fill_manual(values=colores) + 
  theme(legend.position=c(0.75,0.75),
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_text(size=10, hjust=0.5),
        axis.text.x=element_blank()) + 
  ggtitle(label="Patogenicity of consolidated CNVs in agreement") + 
  xlab(label="Samples") + 
  ylab(label="N of CNVs") + 
  labs(fill="Patogenicity class")
figura_08


## ¿Que tipos CNVs tienen los individuos con algún CNV patogenico?
#Seleccionaremos solo las muestras con CNV del tipo patogenico.
samples<- cnv_common_pt[cnv_common_pt$PT=="pathogenic", "sample_name"]
#Seleccionamos todas las lineas en las que hay un CNV para esas muestras.
lines<- c() #Aquí guardaremos los números de linea.
for(i in 1:nrow(cnv_common_pt)){
  if(is.element(cnv_common_pt$sample_name[i], samples) == TRUE){
    lines[i]<- i
  }
}
lines<- na.omit(lines) #Eliminamos los valores NA
datos_figura_09<- cnv_common_pt[lines,]
#También crearemos un indice que nos indique que numero de CNV es cada CNV (para su individuo en concreto).
ind_max<- as.numeric(table(datos_figura_09$sample_name)) #Esta variable nos indicará el indice máximo de los CNVs de cada muestra.
ind_max<- ind_max[ind_max!=0] #Seleccionamos solo los valores que no son 0
names(ind_max)<- samples #Les ponemos nombres de muestras.
#Creamos un bucle que cree secuencias ordenadas solo hasta ese número máximo.
indices<- list() #Los guardaremos en esta lista vacía.
for(i in 1:length(ind_max)){
  indices[[i]]<- seq(1, ind_max[i])
}
datos_figura_09$cnv_index<- unlist(indices) #Ahora los sacamos de la lista y lo guardamos como una nueva columna.
head(datos_figura_09)
#Codificaremos los nombres de las muestras para no mostrarlos en el artículo.
#Creamos un vector con los nuevos nombres, refactorizamos y sobreescribimos los datos.
coded_names<- paste("Sample", 1:100, sep=" ")
datos_figura_09$sample_name<- factor(datos_figura_09$sample_name, levels= muestras, labels= coded_names)

#Crearemos un scatter plot para los individuos que tienen CNV patogénicos.

figura_09<- ggplot(data= datos_figura_09, aes(x=cnv_index, y=MVP, fill=PT)) +
  geom_point(shape=21, size=2, color="black") + 
  ylim(0,1) + 
  facet_wrap(~sample_name, scales= "free_x") + 
  theme_ipsum() + 
  scale_fill_manual(values=colores) + 
  theme(legend.position=c(0.9,0.1),
        #legend.justification=c("bottom"),
        plot.title=element_text(size=14, hjust=0.5),
        axis.title.y=element_text(size=10, hjust=0.5),
        axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.text.y=element_text(size=5),
        strip.text.x=element_text(size=10, hjust=0.5),
        legend.text=element_text(size=12)) + 
  ggtitle(label="Samples with pathogenic CNVs and their MVP values") + 
  xlab(label="Per sample CNV index") + 
  ylab(label="MVP score") + 
  labs(fill="Patogenicity class")
figura_09
#De aquí concluimos que hay solo 10 individuos con CNV patogenicos. La mayoría son o indecisos o benignos.
#Pero algunos si que tienen bastantes "likely pathogenic" también.



## ¿En que regiones caen esos CNV patogenicos?
datos_tabla_10<- cnv_common_pt[cnv_common_pt$PT=="pathogenic",]
#Codificamos los nombres.
datos_tabla_10$sample_name<- factor(datos_tabla_10$sample_name, levels= muestras, labels= coded_names)
datos_tabla_10
#De aquí podriamos explicar que son todo deleciones. Y en los cromosomas 5, 6, 9, 11, 12 y 17.



# 3.- Exportar las figuras y tablas.
if (file.exists("Resultados")) {
  cat("The folder already exists")
} else {
  dir.create("Resultados")
}
directorio_out<- "Resultados/" #Definimos el direcotorio de salida.

#Exportamos las tablas.
library(xlsx)

write.xlsx(datos_tabla_09, 
           file=paste(directorio_out, "tabla_09.xlsx", sep=""))
write.xlsx(datos_tabla_10, 
           file=paste(directorio_out, "tabla_10.xlsx", sep=""))

#Exportamos las figuras.
ggsave(figura_08, 
       filename= paste("figura_08", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
ggsave(figura_09, 
       filename= paste("figura_09", ".png", sep=""), 
       path= directorio_out,
       device="png", bg="white", dpi=300, width=210, height=118, units= "mm")
sessionInfo()