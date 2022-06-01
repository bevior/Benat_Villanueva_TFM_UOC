##########################################################################################################
cd /mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque
mkdir datos_reordenados
cd datos_reordenados

#Es posible que el proceso anterior de reordenar los datos haya sido un poco lioso. Así que vamos a copiar los resltados a este directorio y poner nombres decarpeta más entendibles.
#Notesé que todos los datos están en un formato propio pero compatible con el estándar BED.



##########################################################################################################
#Paso 1: Establecer el directorio de origen de los datos.

directorio="/mnt/d9fa213d-9e64-475f-8e07-bb8a8721253b/benat01/parque/reordenar_datos"
#Cambiar este directorio según se necesite.



##########################################################################################################
#Paso 2: Copiar los datos y renombrar las carpetas.

#Todos los CNV de DECoNT en un solo archivo. Ordenados por muestra y chromosoma o solo por chromosoma.
cp -R $directorio/Paso1_DECoNT_to_BED DECoNT_all
ls DECoNT_all

#CNV de DECoNT por cada individuo.
cp -R $directorio/Paso3_separate_decont_cnv DECoNT_by_samples
ls DECoNT_by_samples

#CNV de DECoNT consilodados para todo el exoma.
cp -R $directorio/Paso2_consolidate_decont_cnv DECoNT_by_exome
ls DECoNT_by_exome

#Todos los CNV de GATK en un solo archivo. Ordenados por muestra y chromosoma o solo por chromosoma.
cp -R $directorio/Paso5_fuse_gatk_samples GATK_all
ls GATK_all

#CNV de GATK por cada individuo.
cp -R $directorio/Paso4_GATK_to_BED GATK_by_samples
ls GATK_by_samples

#CNV de GATK consolidados para todo el exoma.
cp -R $directorio/Paso6_consolidate_gatk_cnv GATK_by_exome
ls GATK_by_exome


