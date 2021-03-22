#!usr/bin/bash

#Script designed to simulate different plasmid integrations in a genomic reference sequence 

fastadir="../TestDataset/SimIntegration/simFASTA/"
mkdir -p ${fastadir}

perl SimulateIntegration.pl ../TestDataset/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna.toplevel_Scaffold0.fa ../TestDataset/CelTag.fa 0.51 > ${fastadir}Scaffold0_0.5C_IS.fa 
perl SimulateIntegration.pl ../TestDataset/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna.toplevel_Scaffold0.fa ../TestDataset/CelTag.fa 1.03 > ${fastadir}Scaffold0_1C_IS.fa 
perl SimulateIntegration.pl ../TestDataset/Cricetulus_griseus_chok1gshd.CHOK1GS_HDv1.dna.toplevel_Scaffold0.fa ../TestDataset/CelTag.fa 2.05 > ${fastadir}Scaffold0_2C_IS.fa 
