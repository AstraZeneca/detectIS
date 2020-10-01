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
echo "Starting minimap2 mapping on plasmid ref"

#Mapping Reads to the plasmid reference and extracting the mapped ones
singularity exec $singimg bash -c "minimap2 -x sr -c ${extragenome1} ${reads[0]} -t ${paramscpu} > ${name}_mate1_vir.paf"
singularity exec $singimg bash -c "minimap2 -x sr -c ${extragenome1} ${reads[1]} -t ${paramscpu} > ${name}_mate2_vir.paf"
singularity exec $singimg bash -c "cut -f1 ${name}_mate1_vir.paf > ${name}_vir.lst"
singularity exec $singimg bash -c "cut -f1 ${name}_mate2_vir.paf >> ${name}_vir.lst"
singularity exec $singimg bash -c "less -S ${name}_vir.lst | sort | uniq > ${name}_vir_mapping.lst"
singularity exec $singimg bash -c "seqtk subseq ${reads[0]} ${name}_vir_mapping.lst | gzip -vc > ${name}_vir_R1.fq.gz"
singularity exec $singimg bash -c "seqtk subseq ${reads[1]} ${name}_vir_mapping.lst | gzip -vc > ${name}_vir_R2.fq.gz"

echo $(date)
echo "minimap2 mapping on plasmid ref completed, mapping reads on host genome"
#Mapping Reads onto the host genome and run detectIS
singularity exec $singimg bash -c "minimap2 -x sr -c ${gref} ${name}_vir_R1.fq.gz -t ${paramscpu} > ${name}_mate1_gnm.paf"
singularity exec $singimg bash -c "minimap2 -x sr -c ${gref} ${name}_vir_R2.fq.gz -t ${paramscpu} > ${name}_mate2_gnm.paf"

echo $(date)
echo "minimap2 mapping complete, runnig detectIS"
singularity exec $singimg bash -c "PERL_BADLANG=0 LANG='FOO' perl ../Workflows/bin/detectIS.pl -h1 ${name}_mate1_gnm.paf -h2 ${name}_mate2_gnm.paf -v1 ${name}_mate1_vir.paf -v2 ${name}_mate2_vir.paf -o ${name}"

echo $(date)
echo "detectIS complete"

ln -s ../Workflows/bin/detectIS.png detectIS.png
../Workflows/bin/CreatePDFandHTML.sh Res/
rm detectIS.png
