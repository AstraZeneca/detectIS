#!/usr/bin/bash

#Script designed to run detectIS using 
#a host genome and a viral genome as references 

singimg="../utils/detectIS.simg"

gref="./Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna.toplevel_Scaffold0.fa"
extragenome1="./CelTag.fa"

reads=( ./Sim_R1.fq.gz ./Sim_R2.fq.gz ) 

mkdir -p Res
name="./Res/SimRead"
paramscpu=4


echo $(date)
echo "Starting minimap2 mapping on viral ref"

#Mapping Reads onto the plasmid and extracting the mapped ones
/usr/local/singularity/bin/singularity exec $singimg bash -c "minimap2 -x sr -c ${extragenome1} ${reads[0]} -t ${paramscpu} > ${name}_mate1_vir.paf"
/usr/local/singularity/bin/singularity exec $singimg bash -c "minimap2 -x sr -c ${extragenome1} ${reads[1]} -t ${paramscpu} > ${name}_mate2_vir.paf"
/usr/local/singularity/bin/singularity exec $singimg bash -c "parallel -k cut -f1 ::: ${name}_mate1_vir.paf ${name}_mate2_vir.paf >  ${name}_vir.lst"
/usr/local/singularity/bin/singularity exec $singimg bash -c "sort --parallel ${paramscpu} ${name}_vir.lst | uniq > ${name}_vir_mapping.lst"
/usr/local/singularity/bin/singularity exec $singimg parallel --link seqtk subseq {1} ${name}_vir_mapping.lst ">" {2} ::: ${reads[0]} ${reads[1]} ::: ${name}_vir_R1.fq  ${name}_vir_R2.fq

echo $(date)
echo "minimap2 mapping on plasmid ref completed, mapping reads on host genome"
#Mapping Reads onto the host genome and run detectIS
/usr/local/singularity/bin/singularity exec $singimg bash -c "minimap2 -x sr -c ${gref} ${name}_vir_R1.fq -t ${paramscpu} > ${name}_mate1_gnm.paf"
/usr/local/singularity/bin/singularity exec $singimg bash -c "minimap2 -x sr -c ${gref} ${name}_vir_R2.fq -t ${paramscpu} > ${name}_mate2_gnm.paf"

echo $(date)
echo "minimap2 mapping complete, runnig detectIS"
/usr/local/singularity/bin/singularity exec $singimg bash -c "perl ../Workflows/bin/detectIS.pl -h1 ${name}_mate1_gnm.paf -h2 ${name}_mate2_gnm.paf -v1 ${name}_mate1_vir.paf -v2 ${name}_mate2_vir.paf -o ${name}"

echo $(date)
echo "detectIS complete"

ln -s ../Workflows/bin/detectIS.png detectIS.png
../Workflows/bin/CreatePDFandHTML.sh Res/
rm detectIS.png


