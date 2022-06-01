##################################################################################################################
#En este script instalaremos y usaremos X-CNV para obtener posibles consecuencias funcionales y patogénicas de los CNV detectados.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir X-CNV_phenotiping



######################################################################################################################
#Paso 0: INstalación de X-CNV.
#Artículo de referecia: https://genomemedicine.biomedcentral.com/articles/10.1186/s13073-021-00945-4#Sec20
#Repositorio en GitHub: https://github.com/kbvstmd/XCNV

#Empezaremos abriendo un entorno de Conda vacío e instalando R, los paquetes básicos y los paquetes data.table y xgboost. Estos últimos se llaman r-data.table y r-xgboost en Conda.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01
mkdir X-CNV_instalation
cd X-CNV_instalation

conda create -n X-CNV 
source activate X-CNV #Activamos el entorno "source" es necesario solo en sistemas operativos (SO) GNU/Linux. En Windows o MacOS es solo "activate".
#Para salir del entorno lo mejor es hacer un "source activate" sin especificar entorno. Entonces se sale y se vuelve al entorno base.

#En Conda se pueden instalar cosas desde distintos canales 0 "chanels". Al igual que en R puedes instalar desde el CRAN o desde Bioconductor por ejemplo.
#Unos canales famosos sob Bioconda o Conda Forge. Se puede especificar de donde descargar un paquete poniendo "-c" de "--channel". O se pueden activar los canales de manera predeterminada para el entorno Conda con estos comandos
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

conda install r-essentials r-base r-data.table r-xgboost #Instalamos los paquetes necesarios.
conda env list #Para ver todos los entornos.
conda list #Para ver todo lo instalado en este entorno.


#Ahora descargaremos el .ZIP del repositorio y lo descomprimiremos.
wget https://github.com/kbvstmd/XCNV/archive/refs/heads/master.zip
unzip master.zip
#Siguiendo las instrucciones del GitHub, tendremos que ejecutar un script para completar la intalación.
cd XCNV-master
nohup bash -c "sh Install.sh" > ../historial_instalacion.txt & #Esto no debería de cambiar nada del PATH, todo se hace en el directorio de trabajo actual.
#Como son muchos satos a descargar, lo hacemos con un nohup.la instalación tarda varias horas.

cat ../historial_instalacion.txt
tail -f ../historial_instalacion.txt

#Una vez este todo instalado podemos nos tendremos el programa disponible en el siguiente directorio:
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/X-CNV_instalation/XCNV-master/bin/XCNV

popd



####################################################################################################################
#Paso X: Usar X-CNV.
#Usarlo es muy sencillo. Primero tendremos que activar el entorno Conda.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/X-CNV_instalation
source activate X-CNV
popd

#El fichero tiene que tener la siguiente estructura separada por tabulaciones:
cromosoma	inicio	final	gain/loss
1	50	90	loss
1	100	200	loss
1	150	250	gain


#Luego llamamos al programa y le damos el fichero separado por tabuacines que tiene que anotar. Tarda unos pocos minutos.
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/X-CNV_instalation/XCNV-master/bin/XCNV ejemplo.bed




###############################################################################################################
#Paso 1: Refromatear los datos.

mkdir Paso1_data_format

#A X-CNV podemos pasarle intervalos desordenados o solapantes, pero no conservará información más allá de las 3 primeras columnas. Perderemos la información de la muestra si las hacemos todas a la vez.
#Así que tenemos que separar los datos por muestra.
#Además tendremos que porcesar los datos de modo que usemos solo los CNVs en común y consolidados por ambos programas que estén de acuerdo en el genotipo.
#Para hacer todo esto abriremos R y seguiremos los pasos del archivo "BED_to_X-CNV_Linux.R". 
#El archivo "BED_to_X-CNV.R" es exactamente igual pero con los delimitadores de Windows.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/R_instalation
source activate R_benat01
popd
R
#Una vez terminado salimos de R
q()

#Le echamos un ojo a los resultados.
cat Paso1_data_format/CAM-GFF-10_cnv_common.bed
#Comprobamos cuantos CNVs en común, consolidados y en acuerdo tiene cada individuo.
wc -l Paso1_data_format/*.bed



################################################################################################################
#Paso 2: Phenotipado mediante X-CNV.
#Ahora haremos un bucle que itere por cada muestra, coja su archivo correspondiente y lo pase por X-CNV.

mkdir Paso2_phenotyping

#Abrimos el entorno conda.
pushd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/X-CNV_instalation
source activate X-CNV
popd

#Crearemos una variable que guarde los nombres de las muestras para iterar sobre ellas.
samples=`gawk '{ print $4;}' ../intersecciones/all_samples/cnv_common.bed | sort | uniq`
echo $samples

#X-CNV guarda los ouputs en el mismo directorio que los inputs, así que los copiaremos todos a la carpeta de salida.
cp Paso1_data_format/*_cnv_common.bed Paso2_phenotyping
#Ahora ejecutaremos X-CNV en bucle para todas las muestras. No lo paralelizaremos por que da problemas al igual que al usar nohup pero peor. No solo falla puede dejar la unidad entera en modo solo lectura. Se puede solucionar reiniciando el servidor.
for i in $samples; do file="${i}_cnv_common.bed";
/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/X-CNV_instalation/XCNV-master/bin/XCNV Paso2_phenotyping/$file;
echo "Se ha procesado la muestra $i";
done
#El proceso entero puede durar poco menos de hora y media.
#Ahora borraremos los archivos de input y los que son igual que el input pero ordenados.
rm Paso2_phenotyping/*_cnv_common.bed
rm Paso2_phenotyping/*_cnv_common.sort.bed
#Comprobamos que tenemos un por cada muestra. Y que tienen el mismo número de lineas que antes (942). Debemos tener en cuenta que X-CNV agrega un encabezado a los resultados por lo que tendremos una linea más por cada muestra.
ls Paso2_phenotyping
ls Paso2_phenotyping | wc -l
cat Paso2_phenotyping/*_cnv_common.output.csv | wc -l
cat Paso1_data_format/*_cnv_common.bed | wc -l


#################################################################################################################
#Paso 3: Concatenar los datos.
#A continuación concatenaremos todos los archivos en uno solo para leerlo más fácilmente con R.

mkdir Paso3_concatenate

cat Paso2_phenotyping/*_cnv_common.output.csv > Paso3_concatenate/cnv_common_mvp.csv
#Comprobamos el número de lineas del archvivo (943 + 100).
wc -l Paso3_concatenate/cnv_common_mvp.csv



#################################################################################################################