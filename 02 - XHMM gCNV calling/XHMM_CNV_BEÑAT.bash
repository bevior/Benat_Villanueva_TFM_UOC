############################################################################################################
#En este script haremos un CNV calling con el programa XHMM (eXome Hidden Markov Models) y después puliremos los resultados con DECoNT (Deep Exome Copy Number Tunner).



###############################################################################################################
#Paso 0: Instalación de XHMM.
#A fecha de Abril de 2022, XHMM se puede instalar desde Conda (Bioconda).
#Creamos un nuevo directorio para todo esto, por si a Conda le da por instalar el entorno en el directorio de trabajo. Que ha sido el caso. Creo que el Miniconda 3 del servidor está específicamente configurado apra hacerlo, por que normalemtne no debería de funcionar así.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/XHMM_instalation

#Seguiré estas dos webs oficiales de referencia:
#https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#specifying-location
#https://docs.anaconda.com/anaconda/user-guide/tasks/using-r-language/
conda create -n r-xhmm r-essentials r-base
source activate r-xhmm #Activamos el entorno "source" es necesario solo en sistemas operativos (SO) GNU/Linux. En Windows o MacOS es solo "activate".
#Para salir del entorno lo mejor es hacer un "source activate" sin especificar entorno. Entonces se sale y se vuelve al entorno base.
conda env list #Para ver todos los entornos.
conda list #Para ver todo lo instalado en este entorno.
#Con esto instalamos R (no sé que versión escogerá) y otros paquetes esenciales para R. R era solo necesario para hacer unos gráficos de los resultados de XHMM, pero lo instalamos igualmente.

#En Conda se pueden instalar cosas desde distintos canales 0 "chanels". Al igual que en R puedes instalar desde el CRAN o desde Bioconductor por ejemplo.
#Unos canales famosos sob Bioconda o Conda Forge. Se puede especificar de donde descargar un paquete poniendo "-c" de "--channel". O se pueden activar los canales de manera predeterminada para el entorno Conda con estos comandos
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#Según los de Bioconda (https://bioconda.github.io/user/install.html#set-up-channels) es importante activarlos en ese orden.

#Ahora que tenemos todos los canales activados, instalamos XHMM:
conda install xhmm
#Si no funcionase se podría instalar XHMM especificando que hay que usar el canal de Bioconda (https://anaconda.org/bioconda/xhmm):
#conda install -c bioconda xhmm

#En principio ya se ha instalado correctamente. Para activarlo:
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/XHMM_instalation
source activate r-xhmm



################################################################################################
#Protocolo a seguir: https://statgen.bitbucket.io/xhmm/tutorial.html
#Creamos la carpeta en la que trabajaremos.
mkdir /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV

#Copiaremos los intervalos que definen los exomas secuenciados. En esta ocasion cogeremos todos a excepción de los cromosomas sexuales y el ADN mitocondrial.
grep -v 'chrX\|chrY\|chrM' /MSA01/EXOMES-CAMEROON/intermedios/Manifiesto_cameroon2.bed > Manifiesto_cameroon2_no_chrXYM.bed

#También haremos una carpeta para los historiales.
mkdir historiales



#############################################################################################
#Paso 1: Calcular las "depth of coverage" con GATK y juntarlos todos en un archivo con forma de matriz.
#Tutorial oficial de la herramienta: https://gatk.broadinstitute.org/hc/en-us/articles/4418062816667-DepthOfCoverage-BETA-
#Creamos el alias para GATK, activamos el entorno Conda, comprobamos que funciona y volvemos:
GATKv42='/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk'
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0
source activate gatk
popd
$GATKv42 --help
$GATKv42 --list
#Ya estamos de vuelta en el directorio de trabajo.

#Creamos la carpeta de salida.
mkdir Paso1_DepthOfCoverage
#Creamos una carpeta para albergar los scripts de parallel
mkdir parallel_scripts

for i in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do echo $i; done #Así obtenemos todos los nombres de archivos. Ojo que tienen la ruta comleta.
for input in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do output=${input#/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/}; echo $output; done #Así creamos otra variable que tiene solo el nombre del archivo, sin la ruta.

#Haremos un bucle que escriba todos los comandos necesarios individualmente en un archivo.
for input in /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/*.bam; do output=${input#/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/};
echo "/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk DepthOfCoverage \
-I $input \
-L Manifiesto_cameroon2_no_chrXYM.bed \
-R /MSA01/REF/hg19/ucsc.hg19.fasta \
-imr OVERLAPPING_ONLY \
-pt sample --verbosity INFO --omit-depth-output-at-each-base --omit-locus-table \
--min-base-quality 0 --start 1 --stop 5000 --nBins 200 \
--include-ref-n-sites \
--count-type COUNT_READS \
--output-format TABLE \
-O Paso1_DepthOfCoverage/${output%.bam}.DATA" >> parallel_scripts/depthofcoverage_parallel.txt;
done
cat parallel_scripts/depthofcoverage_parallel.txt
wc -l parallel_scripts/depthofcoverage_parallel.txt #Nos aseguramos de que el archivo tiene el número de lineas correcto, por si hemos escrito mas de una vez sobre el y tenemos los comandos duplicados.

#Finalmente corremos los comandos de ese archivo con Parallel (Ojo que cada proceso va a ocupar 10GB de RAM. Usando un solo hilo cada muestra tarda entre 60 y 90 miutos en procesarse):
nohup bash -c "parallel --jobs 32 < parallel_scripts/depthofcoverage_parallel.txt" > historiales/paso1_depthofcoverage_parallel.txt &

cat historiales/paso1_depthofcoverage_parallel.txt
tail -f historiales/paso1_depthofcoverage_parallel.txt
#El input (-I) es el archivo .BAM. El intervalo (-L) son la lista de intervalos. La referencia (-R) es el genoma de referencia. Y el output (-O) es el archivo de salida.
#Ojo que de manera predeterminada el output está en format CSV y XHMM lo lee en formato TSV (Tab Separated Values).
#Este proceso crea 4 archivos por cada individuo. Para el siguiente paso necesitaremos los que terminana en .sample_interval_summary.
#Para las 100 muestras, me ha tardado entr 8 y 9 horas.

#En el siguiente paso juntaremos todas las muestras (en filas) con sus respectivos datos de read depth para cada ubicación (columnas) en un archivo grande con forma de matriz.
#Creamos el siguiente directorio de salida.
mkdir Paso1_CombinedGATKRD
#Ahora vamos a juntar todos los datos de RD para cada intervalo y todos los individuos en un solo archivo especial con forma de matriz. Pero primero tenemos que abrir el entorno Conda adecuado.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/XHMM_instalation
source activate r-xhmm
popd
xhmm --help #Debería de mostrarse si hemos activado el entorno correctamente.

#Tenemos que seleccionar varios archivos a la vez, por lo que intentaremos guardarlos como una variable.
ls Paso1_DepthOfCoverage | grep .sample_interval_summary$ | gawk '{print "--GATKdepths ","Paso1_DepthOfCoverage/",$1," ", OFS="", ORS=""}' #Con este comando podemos listar adecuadamente la ruta hacia cada archivo, con el "-I" por delante y los espacios adecuados entre los campos.
input_DATA=`ls Paso1_DepthOfCoverage | grep .sample_interval_summary$ | gawk '{print "--GATKdepths ","Paso1_DepthOfCoverage/",$1," ", OFS="", ORS=""}'`
echo $input_DATA #Debería de mostrarnos todas las muestras en el formato adecuado.

xhmm --mergeGATKdepths -o Paso1_CombinedGATKRD/DATA.RD.txt \
$input_DATA
#Este proceso dura poco más de 5 mintos.



######################################################################################################################
#Paso 2: Obtener información de GC y zonas"repeated-masked" para filtrar los datos mas a delante (Opcional).
#Debemos volver activar el entorno Conda de GATK
GATKv42='/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0/gatk' #Para llamarlo más fácilmente cuando no estemos dentro de un nohup.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/GATK_instalation/gatk-4.2.5.0
source activate gatk
popd

#Ese paso debería de hacerse con la herramienta GCContentByInterval, pero esto solo estaba en GATK3. En su lugar lo voy a intentar hacer con AnnotateIntervals en GATK4.
#A modo de intervalos a analizar no vamos a meter el .interval_list, si no el Manifiesto_cameroon2_no_chrXYM.bed, que contiene las coordenadas de todos los exomas secuenciados.
mkdir Paso2_AnnotateIntervals

$GATKv42 AnnotateIntervals \
-L Manifiesto_cameroon2_no_chrXYM.bed \
-R /MSA01/REF/hg19/ucsc.hg19.fasta \
-imr OVERLAPPING_ONLY \
-O Paso2_AnnotateIntervals/DATA.locus_GC.txt

#Ahora seleccionamos los intervalos que contengan un contenido de GC extremo (inferior a 0.1 o superior a 0.9).
grep ^chr Paso2_AnnotateIntervals/DATA.locus_GC.txt | gawk 'BEGIN {FS="\t"} {if ($4 < 0.1 || $4 > 0.9) print $0}' > Paso2_AnnotateIntervals/extreme_gc_targets.txt
cat Paso2_AnnotateIntervals/extreme_gc_targets.txt #Para echarle un ojo a la lista de exones excluidos por tener una proporción GC desorbitada.
#Con esto hemos creado un subset de lo anterior solo con los lugares en los que se condsidera que hay cantidades extremas de GC. En nuestro caso particular resulta que ningún exón tienen un contenido de GC extremo, por lo que no se filtrarán.
#Si hubiese zonas que filtrar no estoy seguro de si el formato adecuado será del tipo "chr1 120 130" o "chr1:120:130". Quizá XHMM no lo detecte bien y haya que cambiar el formato manualmente para que lo lea.
#He tenido que cambiar el comando por que le original no se ajustaaba a los datos proporcionados por GATK 4.2. El original era así:
#cat Paso2_AnnotateIntervals/DATA.locus_GC.txt | awk '{if ($2 < 0.1 || $2 > 0.9) print $1}' > Paso2_AnnotateIntervals/extreme_gc_targets.txt

#El siguiente paso opcional sería usando PLINK/Sec para calcular la proporción de "repeated-masked bases" en cada intervalo. Pero no lo haremos por que no disponemos del programa y por que probablemente los comandos no se apliquen a la versión actual del programa.
#Plink/Seq no es lo mismo que Plink. Lo ha desarrollado Hardvard y a fecha de 2022 está disponible aquí: https://zzz.bwh.harvard.edu/plinkseq/



#Ahora filtraremos en base a la información calculada anteriormente.
mkdir Paso2_Filter_GC_mean-center
#Necesitamos usar XHMM para esto.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/XHMM_instalation
source activate r-xhmm
popd

xhmm --matrix -r Paso1_CombinedGATKRD/DATA.RD.txt --centerData --centerType target \
-o Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt \
--outputExcludedTargets Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt.filtered_targets.txt \
--outputExcludedSamples Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeTargets Paso2_AnnotateIntervals/extreme_gc_targets.txt \
--minTargetSize 10 --maxTargetSize 10000 \
--minMeanTargetRD 10 --maxMeanTargetRD 500 \
--minMeanSampleRD 25 --maxMeanSampleRD 200 \
--maxSdSampleRD 150

#Los "outputExcludedTargets" son archivos que van a ser necesarios para pasos posteriores. En total nos ha generado los siguientes 3 arcivos:
#DATA.filtered_centered.RD.txt,  DATA.filtered_centered.RD.txt.filtered_samples.txt (este está vacío por que todos los GC content estaban dentro del rango requerido) y  DATA.filtered_centered.RD.txt.filtered_targets.txt



######################################################################################################################################
#Paso 3: PCA sobre los datos centrados
mkdir Paso3_PCA

xhmm --PCA -r Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt --PCAfiles Paso3_PCA/DATA.RD_PCA

#DATA.RD_PCA es el prefijo de los 3 output de este comando. Los outputs en sí son estos:
#DATA.RD_PCA.PC_LOADINGS.txt (contiene los pesos),  DATA.RD_PCA.PC_SD.txt (contiene las desviaciones estandar de cáda componente principal) y  DATA.RD_PCA.PC.txt (quizá tenga los nuevos valores proyetados sobre los componentes principales)
#Este comando crea un componente principal (CP en castellano, PC en inglés) por cáda muestra. No debería de tardar demasiado.



##########################################################################################################################################
#Paso 4: Normalizar los datos centrados usando la información del PCA.
mkdir Paso4_normalize

xhmm --normalize -r Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt --PCAfiles Paso3_PCA/DATA.RD_PCA \
--normalizeOutput Paso4_normalize/DATA.PCA_normalized.txt \
--PCnormalizeMethod PVE_mean --PVE_mean_factor 0.7

#En --PCAfiles solo indicamos el prefijo de los archivos y el solo ya encuentra los 3.
#En este paso se han quitado los CP que tienen mucha varianza siguiendo los ajustes predeterminados. El resultado son estos 2 archivos:
#DATA.PCA_normalized.txt (quizá contenga los datos proyectados en los PCA que han quedado tras la normalización) y  DATA.PCA_normalized.txt.num_removed_PC.txt (nos dice que CP han sido quitados por tener demasiada varianza).



#############################################################################################################################################
#Paso 5: Filtrar y normalizar mediante z-score los datos del PCA normalizado (sin los CP que mas varianza tenían).
mkdir Paso5_Filter_z-score_normalize

xhmm --matrix -r Paso4_normalize/DATA.PCA_normalized.txt --centerData --centerType sample --zScoreData \
-o Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt \
--outputExcludedTargets Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--outputExcludedSamples Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
--maxSdTargetRD 30

#"--outputExcludedTargets/Samples" nos indican los intervalos o los individuos que han sido eliminados en este paso.
#Este comando nos genera 3 archivos:
#DATA.PCA_normalized.filtered.sample_zscores.RD.txt (este archivo tiene los principales datos), DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt (zonas filtradas mediante z-score) y DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt (individuos filtrados mediante z-score).



#################################################################################################################################################
#Paso 6: Filtramos los datos de "read-depth" originales para que sean los mismos que los ya normalizados y filtrados.
mkdir Paso6_Filter_original_RD_data

xhmm --matrix -r Paso1_CombinedGATKRD/DATA.RD.txt \
--excludeTargets Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt.filtered_targets.txt \
--excludeTargets Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_targets.txt \
--excludeSamples Paso2_Filter_GC_mean-center/DATA.filtered_centered.RD.txt.filtered_samples.txt \
--excludeSamples Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt.filtered_samples.txt \
-o Paso6_Filter_original_RD_data/DATA.same_filtered.RD.txt

#Esto solo nos produce un archivo llamado DATA.same_filtered.RD.txt.



#####################################################################################################################################################
#Paso 7: Descubrir los CNV en los datos normalizados.
mkdir Paso7_CNV_Discover

#Antes de ejecutar el comando hay que tener a mano el archivo "params.txt" que está delimitado por tabulaciones ya va a definir los parametros de la búsqueda de CNVs.
#Aquí uso los ajustes predeterminados recomendados por los autores de XHMM.
#Lo podemos encontrar en el directorio del entorno Conda o crearlo nosotros mismos con el siguiente comando:
cat > params.txt
#Seguido del siguiente texto:
1e-08	6	70	-3	1	0	1	3	1
#Y a continuación pulsar Control + D.

nohup bash -c "xhmm --discover -p params.txt \
-r Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R Paso6_Filter_original_RD_data/DATA.same_filtered.RD.txt \
-c Paso7_CNV_Discover/DATA.xcnv -a Paso7_CNV_Discover/DATA.aux_xcnv -s Paso7_CNV_Discover/DATA" > historiales/paso7_cnv_discover.txt &

cat historiales/paso7_cnv_discover.txt
tail -f historiales/paso7_cnv_discover.txt #Para monitorizarlo.
#Aquí es donde se usa el HMM.
#A pesar de que lo he puesto en segundo plano, el proceso llevará poco menos de media hora.
#En este paso obtenemos los siguientes 5 archivos:
#DATA.aux_xcnv (parece quetiene información sobre los CNVs antes de que se descartasen algunos de ellos)
#DATA.posteriors.DEL.txt (info solo sobre las deleciones) 
#DATA.posteriors.DIP.txt (info solo sobre los diploides, donde no hay CNVs)
#DATA.posteriors.DUP.txt (info solo sobre las duplicaciones)
#DATA.xcnv (tiene los CNVs oficiales)



#########################################################################################################################################################
#Paso 8: Genotipar (conseguir los VCF) los CNV detectados.
mkdir Paso8_Genotype

nohup bash -c "xhmm --genotype -p params.txt \
-r Paso5_Filter_z-score_normalize/DATA.PCA_normalized.filtered.sample_zscores.RD.txt -R Paso6_Filter_original_RD_data/DATA.same_filtered.RD.txt \
-g Paso7_CNV_Discover/DATA.xcnv -F /MSA01/REF/hg19/ucsc.hg19.fasta \
-v Paso8_Genotype/DATA.vcf" > historiales/paso8_genotype.txt &

cat historiales/paso8_genotype.txt
tail -f historiales/paso8_genotype.txt #Para monitorizarlo
#Esto nos pasa los datos de CNV presentes en el archivo DATA.xcnv a el archivo DATA.vcf
#Este archivo ya es un VCF normal (v4.1 eso si) con los CNVs detectados para cada uno de los individuos.



#################################################################################################################
#Paso 0: Instalación de DECoNT.
#A partir de aquí detallaré el proceso de pulimineto de los resultados de XHMM mediante DECoNT, una red neuronal recurrente (RNN) para no solo filtrar, sino corregir los resultados de CNV calling.
#El articulo de refereancia es: https://www.biorxiv.org/content/10.1101/2020.05.09.086082v1
#El repositorio de GitHub es: https://github.com/ciceklab/DECoNT

#Por suerte, esta red está hecha con Python, TensorFlow y Keras, por lo que es fácilmente distribuible mediante un entorno Conda.
#Vamos a salir al directoro superior para instalar el entorno Conda.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01
mkdir DECoNT_instalation
cd DECoNT_instalation
#Descargamos el .zip de todo el repositorio de GitHub y lo descomprimimos.
wget https://github.com/ciceklab/DECoNT/archive/refs/heads/master.zip
unzip master.zip
#Creamos le entono, esto suele tardar un rato.
conda env create -f DECoNT-master/DECoNT_linux.yml
source activate DECoNT_linux
conda env list
conda list
#Comprobaciones de que se ha instalado correctamente.
popd



#Una vez instalado lo activamos con los siguientes comandos:
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/DECoNT_instalation
source activate DECoNT_linux
popd



################################################################################################
#Paso 0: Instalación de Sambamba.
#DECoNT requiere que las "read depth" estén en un formato específico. Es complicado adaptar lo que ya tenemos de DepthOfCoverage de GATK, no sé si sería del todo posible. or otra parte es cierto que usando la herramienta CollectReadCounts de GATK sería más fácil de adaptar.
#De todas maneras los más sencillo es volver a calcular los read depths usando Sambamba, por que DECoNT reconoce el formato automáticametne. Usaremos la versión 0.8.2 del GitHub: https://github.com/biod/sambamba
#Haremos una carpeta para instalar ahí Sambamba y creamos el entorno Conda.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01
mkdir Sambamba_instalation
cd Sambamba_instalation
conda create -n sambamba
source activate sambamba 
#Establecemos los canales a través de los cuales instalaremos descargaremos las cosas.
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge
#Según los de Bioconda (https://bioconda.github.io/user/install.html#set-up-channels) es importante activarlos en ese orden.
#Puede que no sea necesario, pero especificaremos que los queremos descargar desde Bioconda.
conda install -c bioconda sambamba #La instalación dura poco.
sambamba --help #Comprobación de que se ha instalado.
conda env list #Para ver todos los entornos.
conda list

#En principio ya se ha instalado correctamente. Para activarlo:
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/Sambamba_instalation
source activate sambamba
popd



########################################################################################################
#Seguiremos creando carpetas numeradas con la misma numeración de antes, como si todo fuese un solo proceso.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV
mkdir Paso9_Sambamba_RD

#Haremos un bucle para calcular las RD con sambamba para cada muestra una vez más. Usaremos el sigueinte comando de ejemplo que nos dan los de DECoNT en su GitHub:
sambamba depth window -w 1000 HG00096.wes.bam > /home/user/sambamba_read_depths/HG00096.wes.bam_read_depths.txt
#No iteraremos sobre los normbres de los archivos .BAM, si no sobre los nombres de las muestras.
#Los nombres a secas están en el archivo Paso3_PCA/DATA.RD_PCA.PC_LOADINGS.txt como nombres de columna. Leeremos el archivo y borraremos el primer campo, que no nos interesa.
head -1 Paso3_PCA/DATA.RD_PCA.PC_LOADINGS.txt | cut --complement -f1 #Con esto los tenemo separados por tabulaciones si no me equivoco.
nombres=`head -1 Paso3_PCA/DATA.RD_PCA.PC_LOADINGS.txt | cut --complement -f1` #Los guardamos a modo de variables.
for name in $nombres; do echo $name; done | sort | uniq | wc -l #Para comprobar que no tenemos ninguna repetida y que las tenemos todas.

#Al igual que antes, generaremos un archivos con todos los comandos necesarios para que luego Parallel lo lea:
for name in $nombres; do file=`ls /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/bams/$name[!0-9]*.bam`; #ATENCION! ESTE SISTEMA DE NOMBRAR LOS ARCHIVOS PROBABLEMENTE SOLO ES ADECUADO PARA ESTOS DATOS EN PARTICULAR!!
echo "sambamba depth window -w 1000 $file > Paso9_Sambamba_RD/$name.read_depths.txt" >> parallel_scripts/sambamba_parallel.txt; 
done
cat parallel_scripts/sambamba_parallel.txt #Leemos el archivo para asegurarnos de que se ha escrito correctamente.
wc -l parallel_scripts/sambamba_parallel.txt #Es muy importante que el erchivo tenga el número correcto de lineas (comandos) ya que se ha generado a base de añadir, no de sobreescribir. Si ejecutamos el bucle anterior por error otra vez, nos encontraremos con el doble de líneas de las deseadas.
grep ^sambamba parallel_scripts/sambamba_parallel.txt | sort | uniq | wc -l #Además hacemos otra comprobación para asegurarnos de que todas las lineas empiezan bien.
#Ejecutamos Sambamba con nohup, de fondo y en paralelo:
nohup bash -c "parallel --jobs 64 < parallel_scripts/sambamba_parallel.txt" > historiales/paso9_sambamba_rd_parallel.txt &
#Los procesos de Sambamba no ocupan mucha RAM, por lo que le vamos a meter los 64 hilos.

cat historiales/paso9_sambamba_rd_parallel.txt
tail -f historiales/paso9_sambamba_rd_parallel.txt #Al usar Parallel, no se registra nada en este archivo, pero bueno.
#Ahora en principio ya tenemos los archivos con los read-depth creados, el siguiente paso es aplicar DECoNT.



###################################################################################################################################
#Paso 10: Pulir los resultados de XHMM con DECoNT.
#Seguiremos las instrucciones de el repositorio de GitHub: https://github.com/ciceklab/DECoNT

#Abrimos el entorno Conda
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/DECoNT_instalation
source activate DECoNT_linux
popd

#Creamos el directorio de salida.
mkdir Paso10_DECoNT_polishing
#Necesitaremos los archivos del GitHub en el propio directorio de trabajo para que todo esto funcione, los volveremos a descompri, pero en nuestra carpeta para tener una copia mas a mano.
unzip /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/DECoNT_instalation/master.zip -d Paso10_DECoNT_polishing
#Además tendremos que meternos en un directorio del GitHub antes de ejecutarlo, de los contrario no detectará los archivos necesarios. Esto nos obliga a usar rutas absolutas (¿o relativas, pero más complejas?) para los inputs y outputs.
pushd Paso10_DECoNT_polishing/DECoNT-master/DECoNT/scripts
#Tendremos que usar el archivo .xcnv del Paso 8 como input, NO el .VCF del Paso 9.
#Además definiremos los inputs y outputs como variables, para que sean más fáciles de cambiar en un futuro. Pero ojo que lo tenemos que hacer dentro del nohup.
nohup bash -c 'input="/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV/Paso7_CNV_Discover/DATA.xcnv";
output_dir="/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV/Paso10_DECoNT_polishing";
sample_RD="/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV/Paso9_Sambamba_RD";
python DECoNT_polish.py -m pretrained -cn XHMM -i $input -o $output_dir -s $sample_RD/' > /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/XHMM_gCNV/historiales/paso10_decont_polish.txt &
#Es importante ponerle el "/" al final del directorio del los sample_RD, si no no lo detecta (lo coge literalmente como un string o algo).
#Con 100 muestras, este proceso me ha tardado poco más de 3 horas.

popd #Volvemos al directorio de trabajo habitual.
cat historiales/paso10_decont_polish.txt
tail -f historiales/paso10_decont_polish.txt #Para ver el progreso.

#Una vez terminado le echamos un ojo a los cambios que ha hecho:
cat Paso10_DECoNT_polishing/DECoNT_XHMM_polished_cnvs.txt #El archivo entero.
cat Paso10_DECoNT_polishing/DECoNT_XHMM_polished_cnvs.txt | gawk '{if ($5!=$6) print $0}' #Solo las lineas que ha corregido.
grep -v DECoNT Paso10_DECoNT_polishing/DECoNT_XHMM_polished_cnvs.txt | gawk '{if ($5!=$6) print $0}' | wc -l #El número total de cambios. Hemos cogido todas las lineas menos la primera.
grep -v DECoNT Paso10_DECoNT_polishing/DECoNT_XHMM_polished_cnvs.txt | wc -l #El número original de CNVs.

#Como referencia hay 2153-1 lineas en DATA.xncv
#Hay 1475 lineas en DATA.VCF
#En los datos corregidos hay 2052 lineas. Hay 100 lineas que han desaparecido respecto al DATA.xcnv, por que no seleccona los CNV que has sido catalogados como diploides.
#Se supone que hay 495 cambios.



###################################################################################################################################