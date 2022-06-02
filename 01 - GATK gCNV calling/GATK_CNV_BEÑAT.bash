#Este va a ser el script queusaré para hacer el gCNV calling con GATK.
#Referencia: https://gatk.broadinstitute.org/hc/en-us/articles/360035531152--How-to-Call-common-and-rare-germline-copy-number-variants#1 

#El primer paso sería instalar GATK, nos vamos a la web oficial (https://gatk.broadinstitute.org/hc/en-us) y desde ahí a su git-hub (https://github.com/broadinstitute/gatk/releases). Podemos instalarlo desde Docker, o si ya tenemos todas las dependencias, descargar el .zip con wget.
#Lo descomprimimos y creamos el entorno conda correspondiente.
#Que no se nos olvide activar el entorno antes de empezar a trabajar.

#Comandos para la instalación de GTAK v4.2.5.0. Ha tardado cerca de 2h en instalar el entorno.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01
mkdir GATK_instalation
cd GATK_instalation
wget https://github.com/broadinstitute/gatk/releases/download/4.2.5.0/gatk-4.2.5.0.zip
unzip gatk-4.2.5.0.zip
cd gatk-4.2.5.0
./gatk --help #comprobamos que funciona
#Ahora instalamos el entorno Python mediante conda
conda env create -f gatkcondaenv.yml
#conda activate gatk #para activar el entorno, pero si no funciona usar el siguiente comando
source activate gatk
conda activate #para volver al entorno básico
conda env list #Para ver los entornos disponibles.
alias GATKv42='/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk' #guardamos la ruta como un alias para invocarlo más fácilmente
GATKv42 --help #Para comprobar que el alias también funciona
GATKv42´--list #Ver la lista de todas las herramientas

#Una vez instalado hay que ir al directorio donde esta instalado GATK, por que el entorno conda también está ahí, pero oculto, y activarlo.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0
source activate gatk
#también se podría hacer esto pero no estoy seguro de que funciona bien.
source activate /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/benat01/.conda/envs/gatk

#Ahora ya podemos empezar a trabajar
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir GATK_gCNV
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/GATK_gCNV
alias GATKv42='/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk'

#Submuestreamos el archivo bed que contiene las coordenadas de los exomas. Quitaremos los cromosomas sexuales y el ADN mitoconcdrial.
grep -v 'chrX\|chrY\|chrM' /MSA01/EXOMES/intermedios/Manifiesto_2.bed > Manifiesto_2_no_chrXYM.bed

#Me gustaría guardar ciertos directorios como variables pero no me sale.
almacen=`/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams`
alias almacen='/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams'
#Guardamos en otra variable el genoma de referencia
alias hg19='/MSA01/REF/hg19/ucsc.hg19.fasta'



###############################################################################################################
#Paso 1: Collect raw counts data with PreprocessIntervals and CollectReadCounts

#Separar el genoma de referencia en partes mas o menos iguales (binning)
#PreprocessIntervals pads exome targets and bins WGS intervals.

mkdir Paso1_PreprocessIntervlas

GATKv42 PreprocessIntervals \
-R /MSA01/REF/hg19/ucsc.hg19.fasta \
-L Manifiesto_2_no_chrXYM.bed \
--bin-length 0 \
--padding 250 \
-imr OVERLAPPING_ONLY \
-O Paso1_PreprocessIntervlas/targets_no_chrXYM.interval_list

#-R se refiere al genoma de referencia. Tenemos que especificar el directorio.
#-L coordenadas de los exomas secuenciados (cuando los datos son WES).
#--bin-length especifica el tamaño de cada bin. Esto solo se usa cuando los datos son WGS. Desactivamos el binning poniendolo en 0
#padding especifica el numero de bases extra que se van a seleccionar, tanto a derecha como a la izquierda de cada intervalo. Podemos ponerlo entre 100 y 250 para WES.
#-imr OVERLAPPING_ONLY especifica que solo junte las partes del genoma que se solapan (o algo así). Se especifica en todos los comandos de este protocolo.
#-O espeficica la el archivo de salida. En este caso la salida en un archivo .interval_list que contiene intervalos del estilo Picard y sirve para separar el genoma de referencia original en distintas partes (binning).

#Contar cuantas secuencias se solapan en cada uno de esos intervalos en los que hemos dividido el genoma de referencia.
#CollectReadCounts tabulates the raw integer counts of reads overlapping an interval.

mkdir Paso1_CollectReadCounts
mkdir historiales

for i in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do echo $i; done #Así obtenemos todos los nombres de archivos. Ojo que tienen la ruta comleta.
for i in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do output=${i#/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/}; echo $output; done #Así creamos otra variable que tiene solo el nombre del archivo, sin la ruta.

#Ejemplo de los 100 a la vez con nohup. Es especial hay que meterlo todo dentro de las '' por que contiene un loop.
nohup bash -c 'for input in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do output=${input#/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/};
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk CollectReadCounts \
-L Paso1_PreprocessIntervlas/targets_no_chrXYM.interval_list \
-R /MSA01/REF/hg19/ucsc.hg19.fasta \
-imr OVERLAPPING_ONLY \
-I $input \
-O Paso1_CollectReadCounts/${output%.bam}.hdf5;
done' > historiales/paso1_collectreadcounts.txt &

cat historiales/paso1_collectreadcounts.txt
tail -f historiales/paso1_collectreadcounts.txt #Para el seguimiento



############################################################################################################
#Paso 2: (Optional) Annotate intervals with features and subset regions of interest with FilterIntervals

#Este paso trata de anotar cierta información sobre nuestros datos y luego filtrar las secuencias que tengan valores extremos en esas anotaciones. De manera predeterminada solo se calcula el contenido en GC, pero es recomentable también adjuntar un archivo de "mappability" y "segmental duplication content".
#De momento nosotros no tenemos los dos archivos adicionales así que anotaremos solo en contenido en GC. Para obtener los archivos de mappability abría que usar el UCSCTable browser, pero solo le podemos pedir 1.000 regiones cada vez y solo nos puede devolver 10.000 bases cada vez, por lo que es un rollo: https://www.biostars.org/p/181014/
#POR SUERTE en el articulo de referencia del protocolo GATK que estoy usando aparece un enlace que me puede ayudar: https://bismap.hoffmanlab.org/ Aquí tienen los archivos de mappability para varios genomas y usando distintos k-mers.
#Subtraeremos el cromosoma 4:
for i in *umap_hg19.bed; do grep 'track\|chr4' $i> chr4${i#all}; done #No he quetado los cromosomas X ni Y, imagino que con que no estén en el interval_list será suficiente.
#El segmental duplication content en cambio ha sido fácil de descargarlo desde el UCSC Table browser. he escogido el grupo "Repeats" y el track "Segmental Dups". Lo he descargado para todo el genoma y también solo para el cromosoma 4 usando la seccion de filtros.
#Tengo una errata en el nombre del archivo pequeño, lo corrijo cambiandole el nombre:
mv map_and_segment_files/chr4_segmetns_hg19.bed map_and_segment_files/chr4_segments_hg19.bed
#Estan todos en el siguiente directorio: /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/map_and_segment_files

#Anotamos nuestra lista de intervalos en base al genoma de referencia y quizá datos de mappability y segment duplication.

mkdir Paso2_AnnotateIntervals

#Este lo voy a hacer sin archivos adicionales, solo con el contenido en GC.

GATKv42 AnnotateIntervals \
-L Paso1_PreprocessIntervlas/targets_no_chrXYM.interval_list \
-R /MSA01/REF/hg19/ucsc.hg19.fasta \
-imr OVERLAPPING_ONLY \
-O Paso2_AnnotateIntervals/targets_no_chrXYM_annotated.tsv



#Filtramos nuestros datos con los puntos de corte predeterminados por GATK. Hay que tener en cuenta que los inputs son los archivos .tsv con las lecturas que se solapan en cada uno de los intervalos seleccionados. Segunda aprte del Paso 1.

mkdir Paso2_FilterIntervals
ls Paso1_CollectReadCounts | gawk '{print "-I ","Paso1_CollectReadCounts/",$1," ", OFS="", ORS=""}' #Con este comando podemos listar adecuadamente la ruta hacia cada archivo TSV o HDF5, con el "-I" por delante y los espacios adecuados entre los campos.
input_hdf5=`ls Paso1_CollectReadCounts | gawk '{print "-I ","Paso1_CollectReadCounts/",$1," ", OFS="", ORS=""}'`
echo $input_hdf5 #Debería de mostrarnos todas las muestras en el formato adecuado.

GATKv42 FilterIntervals \
-L Paso1_PreprocessIntervlas/targets_no_chrXYM.interval_list \
--annotated-intervals Paso2_AnnotateIntervals/targets_no_chrXYM_annotated.tsv \
$input_hdf5 \
-imr OVERLAPPING_ONLY \
-O Paso2_FilterIntervals/targets_gc_filtered.interval_list
 


########################################################################################################
#Paso 3: Call autosomal and allosomal contig ploidy with DetermineGermlineContigPloidy
#En este paso se cuentan cuantan als ploidias (cuantos cromosomas hay de cada tipo) en nuestro caso no es muy importante por que no usaremos cromosomas sexuales.
#Además este paso crea entrena un modelo Bayesiano con información sobre los contigs, que explique la variabilidad entre las muestras. Luego esto será usado en el paso 4.
#Hay que inlcuir una tabla siguiendo el formato espeficicado aqui: https://gatk.broadinstitute.org/hc/en-us/articles/360051304711-DetermineGermlineContigPloidy
#Para los cromosomas diploides poner lo de la fila 1 o 2. Para los cromosomas aploides como el X o el Y poner lo de sus respectivas filas.

#Ejemplo para el cromosoma 4:
CONTIG_NAME	PLOIDY_PRIOR_0	PLOIDY_PRIOR_1	PLOIDY_PRIOR_2	PLOIDY_PRIOR_3
chr4	0.01	0.01	0.97	0.01

#Ejemplo para todos los cromosomas no sexuales:
CONTIG_NAME	PLOIDY_PRIOR_0	PLOIDY_PRIOR_1	PLOIDY_PRIOR_2	PLOIDY_PRIOR_3
chr1	0.01	0.01	0.97	0.01
chr2	0.01	0.01	0.97	0.01
chr3	0.01	0.01	0.97	0.01
chr4	0.01	0.01	0.97	0.01
chr5	0.01	0.01	0.97	0.01
chr6	0.01	0.01	0.97	0.01
chr7	0.01	0.01	0.97	0.01
chr8	0.01	0.01	0.97	0.01
chr9	0.01	0.01	0.97	0.01
chr10	0.01	0.01	0.97	0.01
chr11	0.01	0.01	0.97	0.01
chr12	0.01	0.01	0.97	0.01
chr13	0.01	0.01	0.97	0.01
chr14	0.01	0.01	0.97	0.01
chr15	0.01	0.01	0.97	0.01
chr16	0.01	0.01	0.97	0.01
chr17	0.01	0.01	0.97	0.01
chr18	0.01	0.01	0.97	0.01
chr19	0.01	0.01	0.97	0.01
chr20	0.01	0.01	0.97	0.01
chr21	0.01	0.01	0.97	0.01
chr22	0.01	0.01	0.97	0.01

#Podemos usar cat > /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/ploidy_priors/all_ploidy_priors.tsv y a continuación copiar el contenido para escribirlo en el archivo. Cerramos el archivo con Control + D.
#Lo mismo con cat > /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/ploidy_priors/chr4_ploidy_priors.tsv

#Estan guardados en el siguiente directorio: /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/ploidy_priors
#Hay uno extra con el sufijo _2 que es copiado del original de Neskuts. Por si resulta que este editor no pone bien las tabulaciones.

mkdir Paso3_DetermineGermlineContigPloidy

cat /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/ploidy_priors/all_ploidy_priors.tsv #Nos aseguramos de que el archivo está bien.

#Tenemos que tener en cuenta que dentro de nohup no hay ninguna variable definida, por lo que tendremos que volver a definir la variable input. Como he tenido problemas con los '_', los he cambiado por "_" (para nohup). 
nohup bash -c "input_hdf5=`ls Paso1_CollectReadCounts | gawk '{print "-I ","Paso1_CollectReadCounts/",$1," ", OFS="", ORS=""}'`;
echo $input;
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk DetermineGermlineContigPloidy \
-L Paso2_FilterIntervals/targets_gc_filtered.interval_list \
-imr OVERLAPPING_ONLY \
$input_hdf5 \
--contig-ploidy-priors /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/ploidy_priors/all_ploidy_priors.tsv \
--output Paso3_DetermineGermlineContigPloidy \
--output-prefix ploidy \
--verbosity DEBUG" > historiales/paso3_determinegermlinecontigploidy.txt &

cat historiales/paso3_determinegermlinecontigploidy.txt
tail -f historiales/paso3_determinegermlinecontigploidy.txt #Para el seguimiento

#En este paso el modelo (creo que es Byaesiano) ha tardado 20 epocas en converger. En el segundo intento también ha tardado 20 epocas, peor lo ha hecho en 9 minutos!



############################################################################################################################
#Paso 4: Call copy number variants with GermlineCNVCaller
#Este es el paso gordo en el que se descubren los CNV. Es el paso que más tiempo, RAM y CPU consume.

#Ejemplo sin hacer ningún shard. Las 100 muestras a ala vez.

mkdir Paso4_GermlineCNVCaller

nohup bash -c "input_hdf5=`ls Paso1_CollectReadCounts | gawk '{print "-I ","Paso1_CollectReadCounts/",$1," ", OFS="", ORS=""}'`;
echo $input;
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk GermlineCNVCaller \
--run-mode COHORT \
-L Paso2_FilterIntervals/targets_gc_filtered.interval_list \
$input_hdf5 \
--contig-ploidy-calls Paso3_DetermineGermlineContigPloidy/ploidy-calls \
--annotated-intervals Paso2_AnnotateIntervals/targets_no_chrXYM_annotated.tsv \
-imr OVERLAPPING_ONLY \
--output Paso4_GermlineCNVCaller \
--output-prefix all \
--verbosity DEBUG" > historiales/paso4_germlinecnvcaller.txt &

cat historiales/paso4_germlinecnvcaller.txt
tail -f historiales/paso4_germlinecnvcaller.txt #Para el seguimiento



#El modelo ha convergido en la epoca 10 y ha tardado 22.37 miutos. Hacerlo todo junto podría tardar 8 dias y 15 horas. Tardaría unas 10h en hacer el exoma completo de una persona. Para calcular tiempos: https://rechneronline.de/add-time/multiplication-division.php
#Runtime total memory: 2.127.560.704 unos 2GB
#Con todas las muestras entrena primero un modelo "denoising warm-up". Ha tardado unas 16h, 1 sola época dividida en 5000 partes.
#Lo siguiente es entrenar otro modelo de "denoising main" que ha tardado unos pocos minutos. Seguido se ha peusto con "sampling main", también unos pocos minutos. Luego ha emepzado con "CNV calling", parece que ha dividio la epoca en 10 partes.
#Lo máximo que he visto usar han sido 75GB de RAM.
#Al final ha tardado 4 días y medio en terminarlo todo. Ha hecho 50 epocas y el modelo no ha convergido, pero creo que se pueden usar los resultados igualmente.
-L #Debería de ser la lista de intervalos sobre los que vamos a trabajar, por lo que usaré los intervalos filtrados. Otra opción sería los intervalos originales "preproccesed".
--contig-ploidy-calls #Tiene que apuntar a la carpeta de ploidy calls generada en el paso 3.
--annotated-intervals #Tienen que ser los intervalos anotados mediante GC, mappability, segmented duplicatione etc. Del paso 2.



#######################################################################################################################
#Paso 5: Call copy number segments and consolidate sample results with PostprocessGermlineCNVCalls.
#En este paso exportamos los CNV  y obtenemos dos tipos de ficheros. Hay que hacer el proceso para cáda muestra, por lo que habrá que usar bucles.

#Omitimos el argumento "--allosomal-contig" por que hemos excluido los cromosomas X e Y del análisis.
Paso4_GermlineCNVCaller/all-calls/SAMPLE_0/sample_name.txt #Aquí podemos encotrar el nombre de cáda muestra.
/MSA01/REF/hg19/ucsc.hg19.dict #Es el diccionario que vamos a usar. Lo usó Neskuts y no es lo mismo que la secuencia entera en FASTA.

mkdir Paso5_PostprocessGermlineCNVCalls

#El bucle entero metido en nohup, con la ruta completa hacia GATK y toda la creación de variables dentro de nohup, que si no, no las reconoce el subproceso.
nohup bash -c 'sample_n=`ls Paso4_GermlineCNVCaller/all-calls | grep SAMPLE_*`;
echo $sample_n;
for line in $sample_n;
do n=${line#SAMPLE_};
sample=`cat Paso4_GermlineCNVCaller/all-calls/SAMPLE_$n/sample_name.txt`;
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk PostprocessGermlineCNVCalls \
--model-shard-path Paso4_GermlineCNVCaller/all-model \
--calls-shard-path Paso4_GermlineCNVCaller/all-calls \
--contig-ploidy-calls Paso3_DetermineGermlineContigPloidy/ploidy-calls \
--sample-index $n \
--output-genotyped-intervals Paso5_PostprocessGermlineCNVCalls/$sample.gCNV_intervals.vcf \
--output-genotyped-segments Paso5_PostprocessGermlineCNVCalls/$sample.gCNV_segments.vcf \
--output-denoised-copy-ratios Paso5_PostprocessGermlineCNVCalls/$sample.gCNV_denoised_copy_ratios.tsv \
--sequence-dictionary /MSA01/REF/hg19/ucsc.hg19.dict;
done' > historiales/paso5_postprocessgermlinecnvcalls.txt &

cat historiales/paso5_postprocessgermlinecnvcalls.txt
tail -f historiales/paso5_postprocessgermlinecnvcalls.txt #Para el seguimiento

#Funciona!! Creo que me ha tardado unos 2 o 3 minutos por muestra.
#No se muy bien por que pero en este caso solo me funciona si meto todo el bucle en nohup con '_'.
for i in Paso5_PostprocessGermlineCNVCalls/*.tsv; do wc -l $i; done #pequeña comprobación de que todos los archivos tienen el mismo número de lineas.
for i in Paso5_PostprocessGermlineCNVCalls/*gCNV_intervals.vcf; do wc -l $i; done
for i in Paso5_PostprocessGermlineCNVCalls/*gCNV_segments.vcf; do wc -l $i; done #Estos no coinciden, no se por que pero creo que es normal.



###################################################################################################################
