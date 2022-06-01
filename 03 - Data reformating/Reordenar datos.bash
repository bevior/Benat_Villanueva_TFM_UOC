##########################################################################################################
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir reordenar_datos
cd reordenar_datos
#Una vez obtenido los CNVs a través de los dos métodos, es hora de ordenarlos de manera que podamos hacer fácilmente las comparaciones necesarias.
#Algunos pasos se llevarán acabo usando R. Crearemos un entrono de Conda para instalarlo.




##########################################################################################################
#Paso 0: Intalación de un entorno R mediante Conda.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01
mkdir R_instalation
cd R_instalation

conda create -n R_benat01 r-essentials r-base #Creamos el entorno y directamente le instalamos R y sus programas básicos.
source activate R_benat01 #Activamos el entorno "source" es necesario solo en sistemas operativos (SO) GNU/Linux. En Windows o MacOS es solo "activate".
conda install -c bioconda r-vcfr #También instalaremos un paquete de R para poder leer archivos VCF. Directamente desde bioconda. SI FALLA HABRÁ QUE INSTALARLO MANUALMENTE MAS ADELANTE.
#Para salir del entorno lo mejor es hacer un "source activate" sin especificar entorno. Entonces se sale y se vuelve al entorno base.
conda env list #Para ver todos los entornos.
conda list #Para ver todo lo instalado en este entorno.
popd

#En Conda se pueden instalar cosas desde distintos canales 0 "chanels". Al igual que en R puedes instalar desde el CRAN o desde Bioconductor por ejemplo.
#Unos canales famosos sob Bioconda o Conda Forge. Se puede especificar de donde descargar un paquete poniendo "-c" de "--channel". O se pueden activar los canales de manera predeterminada para el entorno Conda con estos comandos
#conda config --add channels defaults
#conda config --add channels bioconda
#conda config --add channels conda-forge
#Según los de Bioconda (https://bioconda.github.io/user/install.html#set-up-channels) es importante activarlos en ese orden. PERO DE MOMENTO NO LO VAMOS A HACER, ya que no espero instalar paquetes comlicados en este entorno.

#En principio ya se ha instalado correctamente. Para activarlo:
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/R_instalation
source activate R_benat01
popd



##########################################################################################################
#Paso 1: Pasar los resultados de DECoNT a un formato compatible con BED.

#Lo realizaremos en R, por lo que primero hay que abrir el entorno.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/R_instalation
source activate R_benat01
popd

mkdir Paso1_DECoNT_to_BED

#Los pasos a seguir están en el archivo "DECoNT_to_BED_Linux.R" y serían los correspondientes al paso 1 con el mismo nombre que este.
#El archivo DECoNT_to_BED.R es exactamente igual solo que los separadores de línea son "\r\n" (Windows).
#IMPORTANTE: En los archivos se establecen directorios de entrada y de salida tanto para mi PC Windows como para el servidor GNU/lINUX. Hay que revisar el archivo manualmente y dejar las rutas adecuadas o editarlas completamente.

#Debemos tener en cuenta que este archivo ya no contiene los CNV que han sido corregidos a "NO-CALL". Solo tendremos en cuenta "DEL" y "DUP"

#Abrimos el archivo para comprobar que los resultados están bien.
head Paso1_DECoNT_to_BED/DECoNT_cnv.bed
tail Paso1_DECoNT_to_BED/DECoNT_cnv.bed
gawk '{print $8;}' Paso1_DECoNT_to_BED/DECoNT_cnv.bed | sort | uniq #Comporbamos que solo quedan DEL y DUP.

#Además reordenaremos los CNVs por cromosoma y ubicación, por que algunos programas pueden requerir que los datos estén ya ordenados.
sort -k1,1 -k2,2n Paso1_DECoNT_to_BED/DECoNT_cnv.bed > Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed
head Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed #Comprobamos que el earchivo se ha guardado bien.
tail Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed



##########################################################################################################
#Paso 2: Consolidar CNV de DECoNT para todo el exoma.
#Usaremos Bedtools para ello: https://bedtools.readthedocs.io/en/latest/
#En el servidor está instalada la versión 2.26. No proporciono los comandos para instalarlo en este documento, pero noc reo que sea muy complicado.
#Este tutorial exlica como hacer lo que queremos a partir del paso 10 (bedtools merge): https://sandbox.bio/tutorials?id=bedtools-intro&step=8

mkdir Paso2_consolidate_decont_cnv

#Vamos a usar bedtools merge para juntar los CNVs que se solapen entre los individuos y cntar como que son los mismos.
#Le pedimos que nos junte los intervalos del input que se crucen en una base por lo menos. Luego que nos cuente la columna 1 por cada fila que ha sido "mergeada" (para saber cuantos intervalos se han juntado en uno solo).
#Por último le pedimos que nos pegue la columnas 4 y 5 (que son las que contienen el ID) por cada intervalo juntado, para saber que muestras se solapaban y han sido juntadas en ese intervalo.
#También le pediremos que nos imprima el estado del CNV según XHMM y DECoNT, que corresponden a la columnas 7 y 8.
#Los campos estarán delimitados por comas, pero eso se puede cambiar con el argumento "-delim ";"" o "-delim "|"" por ejemplo.
bedtools merge --help #Para ver las opcines disponibles.
#El comando sería algo así como:
bedtools merge -i Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed -c 1,4,5,7,8 -o count,collapse,collapse,collapse,collapse | head
#bedtools merge -i Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed -d 1000 -c 1,4,5,8 -o count,collapse,collapse,collapse #Si queremos juntarlos aunque se alejen 1000 bases los intervalos.
#Si todo está bien lo guardaremos como un arhcivo informativo. 
bedtools merge -i Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed -c 1,4,5,7,8 -o count,collapse,collapse,collapse,collapse > Paso2_consolidate_decont_cnv/DECoNT_cnv_exome.bed

#Vamos a contar cuantos CNV individuales teníamos antes (1718):
wc -l Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed
#Y cuantos tenemos ahora (320):
wc -l Paso2_consolidate_decont_cnv/DECoNT_cnv_exome.bed

#Ahora tenemos cada CNV de DECoNT agrupado en un intervalo nuevo. Sin duplicados.



##########################################################################################################
#Paso 3: Separar los CNV de DECoNT por individuos.

mkdir Paso3_separate_decont_cnv

#Vamos a iterar sobre los nombres de las muestras. Primero tenemos que obtenerlas. 
#Abrimos el archivo sin header y ya ordenado. Imprimimos la columna que tiene el nombre de la muestra, ordenamos y eliminamos duplicados. Deberían ser 100 individuos distintos.
gawk '{print $6;}' Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed | sort | uniq | wc -l
#Guardaremos estos nombres de individuos como una varieble.
samples=`gawk '{print $6;}' Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed | sort | uniq`
echo $samples | wc
echo $samples #Debería de imprimirnos todas bien.
#Podemos iterar sobre esta variable cada vez que queramos usar los nmobres de las muestras.

#Ahora usaremos un bucle que itera por cada nombre de las muestras para escoger las lineas que coincidan solo con la muestra correspondiente.
#ATENCION cuidado con las muestras CAM-GFF-10, CAM-GFF-104 y CAM-GFF-107. Empiezan igual y cuando hacemos un "grep CAM-GFF-10" nos pilla las tres muestras!!
#Notese el "grep -w" para buscar un exact match con el nombre de la muestra.
for i in $samples; do grep -w $i Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed > Paso3_separate_decont_cnv/${i}_DECoNT_cnv.bed; done
#Los archivos se llamarán igual que el originale pero con el prefijo de su muestra.

#Podemos ver cuantos CNVs tiene cada individuo leyendo el numero de líneas de cáda archivo.
wc -l Paso3_separate_decont_cnv/*_DECoNT_cnv.bed
wc -l Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed #El archivo original tiene 2052 CNVs.
#también podemos echarle un ojo al archivo problemático en concreto para asegurarnos de que todo ha salido bien.
cat Paso3_separate_decont_cnv/CAM-GFF-10_DECoNT_cnv.bed



##########################################################################################################
#Paso 4: Pasar los resultados de GATK a un formato compatible con BED.

#Lo realizaremos en R, por lo que primero hay que abrir el entorno.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/R_instalation
source activate R_benat01
popd

mkdir Paso4_GATK_to_BED

#Una vez abierto R, si "library(vcfR)" no funciona lo instalaremos desde el propio R.
#install.packages("vcfR") #Seleccionamos el espejo del CRAN español si es que nos preguntan. Es posible que tarde unos minutos, seguramente tendrá muhcas dependencias.
#Los pasos a seguir están en el archivo "GATK_to_BED_Linux.R" y serían los correspondientes al paso 4 con el mismo nombre que este.
#El archivo GATK_to_BED.R es exactamente igual solo que los separadores de línea son "\r\n" (Windows).
#IMPORTANTE: En los archivos se establecen directorios de entrada y de salida tanto para mi PC Windows como para el servidor GNU/lINUX. Hay que revisar el archivo manualmente y dejar las rutas adecuadas o editarlas completamente.

#Comprobamos el numero de archivos generados y su estrcutura.
ls Paso4_GATK_to_BED
ls Paso4_GATK_to_BED | wc -l
cat Paso4_GATK_to_BED/AFA016_GATK_cnv.bed
cat Paso4_GATK_to_BED/WAZ118_GATK_cnv.bed
#Comprobamos que no quedan CNVs "DIP"
grep DIP Paso4_GATK_to_BED/*_GATK_cnv.bed

#Para terminar vamos a comprobar cuantos CNV teniamos antes contando los "DIP".
directorio="/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/laboratorio/GATK_gCNV/Paso5_PostprocessGermlineCNVCalls/"
grep -v ^# $directorio*.gCNV_segments.vcf | wc -l #Escojemos las lineas que no empiecen por "#" (las del encabezado) y las contamos para todos los archivos. (90021).
wc -l Paso4_GATK_to_BED/*_GATK_cnv.bed #Ahora se nos han quedado en 52691.




##########################################################################################################
#Paso 5: Juntar todos los CNV de todas las muestras en un solo archivo.

mkdir Paso5_fuse_gatk_samples

#Vamos a juntar otra vez todos los individuos en un solo archivo BED. Además también haremos una variante del archivo en el que los CNVs estén ordenados por chromosoma y no por individuos.
#Para ello volveremos a iterar sobre los nombres de las muestras.
samples=`gawk '{print $6;}' Paso1_DECoNT_to_BED/DECoNT_cnv_chr_ordered.bed | sort | uniq`
for i in $samples; do cat Paso4_GATK_to_BED/${i}_GATK_cnv.bed >> Paso5_fuse_gatk_samples/GATK_cnv.bed; done

#Comprobamos que tiene la estrcutura adecuada
head Paso5_fuse_gatk_samples/GATK_cnv.bed
tail Paso5_fuse_gatk_samples/GATK_cnv.bed
#Y el número de líneas adecuado
wc -l Paso4_GATK_to_BED/*_GATK_cnv.bed
wc -l Paso5_fuse_gatk_samples/GATK_cnv.bed


#Ahora haremos la copia del archivo pero con los CNV ordenados por cromosoma y luego posición de inicio.
sort -k1,1 -k2,2n Paso5_fuse_gatk_samples/GATK_cnv.bed > Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed
#Hacemos las mismas comprobaciones que antes.
head Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed
tail Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed
wc -l Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed



##########################################################################################################
#Paso 6: Consolidar CNV de GATK para todo el exoma.

mkdir Paso6_consolidate_gatk_cnv

#Volveremos a usar bedtools para juntar los intervalos solapantes del archivo que contiene todos los CNV de GATK ordenador por chromosoma.
#Como el archivos de GATK tiene las columnas ligeramente distintas a el de DECoNT, pedimos las adecuadas de modo que nos de:
#Cuantos CNVs se han juntado en uno solo y los ID_0, ID_1, predicción de número de copias y predicción de phenotipo de cada uno de los CNV que han sido juntados.
bedtools merge -i Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed -c 1,4,5,7,8 -o count,collapse,collapse,collapse,collapse | head
bedtools merge -i Paso5_fuse_gatk_samples/GATK_cnv_chr_ordered.bed -c 1,4,5,7,8 -o count,collapse,collapse,collapse,collapse > Paso6_consolidate_gatk_cnv/GATK_cnv_exome.bed
#De esta manera si se han detectado CNVs en intervalos solapantes entre varios individuos, se contarán solo como uno. Obteniendo los CNVs consolidados para el exoma entero.

#GATK detecta muchos CNV para cada individuo, no es posible verlos bien. Comprobaremos solo la primeras columnas.
gawk '{print $1,$2,$3,$4;}' Paso6_consolidate_gatk_cnv/GATK_cnv_exome.bed
wc -l Paso6_consolidate_gatk_cnv/GATK_cnv_exome.bed
#Podemos ver que la mayoría de CNVs se solapaban entre ellos, por que solo quedan 44 regiones de CNV tras la consolidación.
