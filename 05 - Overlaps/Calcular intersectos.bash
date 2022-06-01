##########################################################################################################
#En este docuemnto calcularemos los CNVs solapantes (o intersectos) entre los dos programas. Esto nos ayudará a saber cosas como cuantos CNVs hay exclusivos de cáda programa y cuantos en común.



##########################################################################################################
#Paso 1: Descubirir los CNVs solapantes en cada uno de los individuos.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir intersecciones
cd intersecciones

#Vamos a crear un bucle que itere sobre los nombres de las muestras. Para ello usaremos el archivo con todos los CNV de DECoNT, ya que es  fácil sacar los nombres de ahí.
#Si el archivo está en otro lugar, cambiar la ruta en el comando de abajo.
directorio="../datos_reordenados/DECoNT_all"
samples=`gawk '{print $6;}' $directorio/DECoNT_cnv_chr_ordered.bed | sort | uniq`
#El objetivo será que por cada individuo, ejecute "bedtools intersect" (de tres maneras distintas) entre los CNV detectados por GATK y los detectados por DECoNT.

mkdir by_sample #Aquí guardaremos todos los CNV solapantes por muestra.

#Usaremos "bedtools intersect -v -a -b" para conseguir los CNV exclusivos de cada programa. Tendremos que hacerlo dos veces, ya que bedtools solo nos devuelve los CNVs del archivo A.
directorio="../datos_reordenados"
for i in $samples; do GATK=`ls $directorio/GATK_by_samples/${i}_GATK_cnv.bed`; DECoNT=`ls $directorio/DECoNT_by_samples/${i}_DECoNT_cnv.bed`;
bedtools intersect -v -a $GATK -b $DECoNT > by_sample/${i}_GATK_cnv_exclusive.bed;
bedtools intersect -v -a $DECoNT -b $GATK > by_sample/${i}_DECoNT_cnv_exclusive.bed;
done
#Ojo que bash no entiende "$i" como variable si va seguido de un "_". Pero se puede poner ${i} y listo.
#Notese que hemos guardado para cada individuo dos directorios uno con los CNVs solapantes de GATK y otro para los de DECoNT.
#Por cáda muestra ahora mismo tenemos dos archivos uno con los CNVs exclusivos de GATK y otro con los CNVs exclusivos de DECoNT.

#Ahora vamos a seleccionar los CNVs que han sido detectado por ambos programas. Para seleccionarlos volveremos a usar "bedtools intersect -v -a -b", pero siendo A todos los CNVs de un programa y B los CNVs exclusivos de ese programa.
#De esta manera estamos seleccionando los CNVs no exclusivos (los que han sido detectados por ambos programas).
#Una vez seleccionados los CNVs en común para ambos programas, los juntaremos en un solo archivo y los consolidaremos usando "bedtools merge -i". De esta manera si hay CNVs solapantes, se juntan en uno solo teniendo en cuneta la base mas downstream y la base mas upsetram.
#Al usar la función "merge" bedtools solo nos devuelve las coordenadas de los intervalos. Tendremos que pedirle explicitamente que nos devuelva también los valores de las demás columnas.
#En principio le pediremos que nos devuleva el nombre de cada muestra. Este se encuentra en la columna 6 y le diremos que no nos muestre duplicados.
#también añadiremos el genotipo, que corresponde a la columna 8 en los archivos BED de ambos programas y crearemos una nueva columna que nos cuente el número de CNVs que han sido consolidados en uno solo.
#Debemos tener en cuenta que en este paso estamos perdiendo el ID, de cada uno de los CNVs anes de ser consolidados, se podría incluir en el output, pero de momento no lo voy a hacer por simploficar las cosas.
directorio="../datos_reordenados"
for i in $samples; do GATK=`ls $directorio/GATK_by_samples/${i}_GATK_cnv.bed`; DECoNT=`ls $directorio/DECoNT_by_samples/${i}_DECoNT_cnv.bed`;
GATK_ex=`ls by_sample/${i}_GATK_cnv_exclusive.bed`; DECoNT_ex=`ls by_sample/${i}_DECoNT_cnv_exclusive.bed`;
bedtools intersect -v -a $GATK -b $GATK_ex > GATK_common.tmp; bedtools intersect -v -a $DECoNT -b $DECoNT_ex > DECoNT_common.tmp;
cat GATK_common.tmp DECoNT_common.tmp | sort -k1,1 -k2,2n | bedtools merge -c 6,8,8 -o distinct,collapse,count -i stdin > by_sample/${i}_cnv_common.bed;
done

#Hechamos un vistazo a los archivos generados. 
ls by_sample | wc -l #Debería de haber 3 archivos por cáda muestra.
#Vemos la estructura de los archivos para la primera muestra.
head by_sample/AFA016_GATK_cnv_exclusive.bed
head by_sample/AFA016_DECoNT_cnv_exclusive.bed
head by_sample/AFA016_cnv_common.bed



#Ahora haremos lo mismo con los datos de los CNV consolidados para todo el exoma.
mkdir by_exome
#Creamos el directorio y guardamos las ubicaciones de los archivos de input.
GATK="../datos_reordenados/GATK_by_exome/GATK_cnv_exome.bed"
DECoNT="../datos_reordenados/DECoNT_by_exome/DECoNT_cnv_exome.bed"
#Guardamos los CNVs exclusivos de cada programa
bedtools intersect -v -a $GATK -b $DECoNT > by_exome/GATK_cnv_exome_exclusive.bed
bedtools intersect -v -a $DECoNT -b $GATK > by_exome/DECoNT_cnv_exome_exclusive.bed
#Seleccionamos los CNVs en común de ambos programas, los juntamos, los ordenamos y los consolidamos.
#Le pediremos una columna extra donde nos indique cuantos CNVs se han juntado en uno solo.
GATK_ex="by_exome/GATK_cnv_exome_exclusive.bed"
DECoNT_ex="by_exome/DECoNT_cnv_exome_exclusive.bed"
bedtools intersect -v -a $GATK -b $GATK_ex > GATK_common.tmp
bedtools intersect -v -a $DECoNT -b $DECoNT_ex > DECoNT_common.tmp
cat GATK_common.tmp DECoNT_common.tmp | sort -k1,1 -k2,2n | bedtools merge -c 1 -o count -i stdin > by_exome/GATKandDECoNT_cnv_exome_common.bed
rm *.tmp

#Comprobamos los archivos generados.
ls by_exome
head by_exome/GATK_cnv_exome_exclusive.bed
head by_exome/DECoNT_cnv_exome_exclusive.bed
head by_exome/GATKandDECoNT_cnv_exome_common.bed


#Para terminar juntaremos los datos de los intersectos por individuos en un solo archivo, de manera que luego sea más fácil de abrir con R.
mkdir all_samples

cat by_sample/*_GATK_cnv_exclusive.bed > all_samples/GATK_cnv_exclusive.bed
cat by_sample/*_DECoNT_cnv_exclusive.bed > all_samples/DECoNT_cnv_exclusive.bed
cat by_sample/*_cnv_common.bed > all_samples/cnv_common.bed

#Comprobamos que estáb bien los archivos.
head all_samples/GATK_cnv_exclusive.bed
head all_samples/DECoNT_cnv_exclusive.bed
head all_samples/cnv_common.bed



###########################################################################################################
