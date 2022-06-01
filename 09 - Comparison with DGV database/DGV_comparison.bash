##################################################################################################################
#En este documento compararemos nuestros CNVs consolidados entre individuos y programas con CNVs de la base de datos DGV.



####################################################################################################################
#Paso 1: Descargar los datos del Table Browser.

#Descargaremos los datos sobre todos los CNVs de DGV (DGV Support) y los DGV Gold estándar. Aunque en principio solo usaremos los primeros datos.
#Enlace a DGV Support desde el Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1357046747_IE91GaZWIlRjW7IEFBiimlEhyH7g&clade=mammal&org=Human&db=hg19&hgta_group=varRep&hgta_track=dgvPlus&hgta_table=dgvSupporting&hgta_regionType=genome&position=chr4%3A1-191%2C154%2C276&hgta_outputType=primaryTable&hgta_outFileName=
#Enlace a DGV Gold Standard en el Table Browser: https://genome.ucsc.edu/cgi-bin/hgTables?hgsid=1357046747_IE91GaZWIlRjW7IEFBiimlEhyH7g&clade=mammal&org=Human&db=hg19&hgta_group=varRep&hgta_track=dgvPlus&hgta_table=dgvGold&hgta_regionType=genome&position=chr4%3A1-191%2C154%2C276&hgta_outputType=primaryTable&hgta_outFileName=
#Ponemos nombre a los archivos y los descargamos a un directorio en nuestro PC con Windows 10.

#En el servidor creamos el directorio de trabajo.
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir DGV_comparison
cd DGV_comparison



##################################################################################################################
#Paso 2: Subir los archivos al servidor.

mkdir Paso2_DGV

#Ahora usamos scp o sftp para subirlo a nuestro servidor, por razones de seguridad no voy a poner el código para hacer este paso.



##################################################################################################################
#Paso 3: Comparar nuestros datos con los del DGV.

mkdir Paso3_intersecting_cnvs

CNV="../intersecciones/by_exome/GATKandDECoNT_cnv_exome_common.bed" #Definiremos nuestro archivo de CNVs como una variable.
#Notar que estamos usando el archivo que contiene 283 CNVs que han sido consolidados entre todos los individuos primero y luego entre programas. En las figuras esto corresponde al cariograma de bandas verdes.



#Usaremos "bedtools intersect -v" para seleccionar los CNVs que no estén registrados en la base de datos. Pero debemos quitarle el encabezado al archivo del DGV primero.
#Para ello usaremos un "grep -v", y lo pasamos por un pipe a bedtools. Bedtools coje el input del pipe si le establecemos "stdin" como input.
grep -v \#chrom Paso2_DGV/DGV_Supp_Var_dgvsupporting.tsv | bedtools intersect -v -a $CNV -b stdin > Paso3_intersecting_cnvs/analisis_unique_CNVs.bed
#Debemos tener en cuenta que estos procesos ocupan mucha RAM, este en particular 9,4GB.
#Solo obtenemos un CNV "novel" que no se solapa con ninguno de la base de datos: chr6    159057595       159066027
cat Paso3_intersecting_cnvs/analisis_unique_CNVs.bed

#Ahora usamos "bedtools intersect -wa -wb -a -b" para crear un archivo que nos indique cuales de nuestros CNVs se solapan con cuales de la base de datos del DGV.
grep -v \#chrom Paso2_DGV/DGV_Supp_Var_dgvsupporting.tsv | bedtools intersect -wa -wb -a $CNV -b stdin > Paso3_intersecting_cnvs/overlapping_CNVs_DGV.bed
#El resultado es un archivo con una linea por cada uno de los CNVs solapantes en la que esta el CNV en cuentión seguido de todas las entradas del DGV con la que se solapa.
#Comprobamos los datos.
head Paso3_intersecting_cnvs/overlapping_CNVs_DGV.bed

#Por último crearemos seleccioanremos esas entradas del DGV que se solapan con nuestros CNVs y las guardaremos en un archivo aparte usando "bedtools intersect -wa -a -b".
grep -v \#chrom Paso2_DGV/DGV_Supp_Var_dgvsupporting.tsv | bedtools intersect -wa -a stdin -b $CNV > Paso3_intersecting_cnvs/overlapping_DGV.bed
#Notese que esta vez hemos puesto el DGV como archivo A, para que nos devuelvan solo sus entradas solapantes, una en cada linea.
#Comprobamos los datos.
head Paso3_intersecting_cnvs/overlapping_DGV.bed
wc -l Paso3_intersecting_cnvs/overlapping_DGV.bed #Nuestros 283 CNVs se solapan con 822.674 CNVs de DGV.


#No lo mostraré pero también he comprobado si nuestros CNVs se solapan con los DGV Gold, la mayoría si lo hacen, solo hay uno 20 que no se solapan.
#De todas maneras lo más interesante ese ese CNV del cromosoma 6 que no se solapa con ninguno del DGV ni DGV Gold. Futuro trabajo sería comprobar si se ha detectado en otras bases de datos como el dbVar o en algún otro estudio.



##################################################################################################################
#Paso 4: Avisar de que en esta carpeta no hay paso 1.

#Avisaremos de que en esta carpeta no hay Paso 1, para evitar confusiones.
#Escribiremos un archivo que lo explique.
cat > no_hay_paso1.txt
#Y dentro el siguiente texto:
El paso 1 sería descargar los datos a un PC personal. Por lo que no hay carpeta del paso 1 en este directorio.

#Y cerramos el archivo con Control + D.



##################################################################################################################