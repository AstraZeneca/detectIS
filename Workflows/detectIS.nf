#!/usr/bin/env nextflow

/*
===============================================================
 detectIS-nextflow
===============================================================
Pipeline to identify Integration Sites (IS) in paired-end
RNA-seq/DNA-seq experiments.
 
NextFlow main file 
Jan 2020.


#### Homepage / Documentation
https://bitbucket.astrazeneca.net/users/kdvb570/repos/detectis/browse
#### Author
Luigi Grassi  <luigi.grassi@medimmune.com>
---------------------------------------------------------------
*/


def helpMessage() {
    log.info """
 
=================================================================================
detectIS pipeline	 
https://bitbucket.astrazeneca.net/users/kdvb570/repos/detectis/browse
=================================================================================

Usage:
    
The typical command for running the pipeline is as follows:

nextflow run  detectIS-nextflow -c 20190515_downsampleRNAseq_detectISTEST.conf 

Mandatory arguments: -c Configuration file with all the parameters used in the analysis
    """.stripIndent()
}


// Show help message
params.help = false
if (params.help){
    helpMessage()
    exit 0
}

if (params.project_name==null){
    helpMessage()
    exit 0
}

Channel
    .fromFilePairs(params.reads)                                             
    .ifEmpty {exit 1, error "Cannot find any reads matching: ${params.reads}" }  
    .into { reads_minimap1 }


// STEP 1 - MAPPING WITH MINIMAP2 READS AS SINGLE END ON PLSM/VIRUS

vir_mini = Channel.fromPath(params.vir_seq)
        .ifEmpty { exit 1, "Viral/plasmid genome reference fasta file not found: please check ${params.vir_seq}" }

process minimap_vir {
        errorStrategy 'retry'
        maxRetries 3
        scratch params.scratch

        cpus = params.cpu.minimap
        memory = 20.GB
        tag "$name"
        publishDir "${params.outdir}/MinimapPAF", mode: 'copy' ,
	saveAs: {filename ->
        		if (filename.indexOf("_vir.paf") > 0) "$filename"
                        else if (filename.indexOf("_mapping.lst") > 0) "$filename"
                        else null
        }
	
        input:
        set val(name), file(reads) from reads_minimap1
        file pref from vir_mini.collect()

        output:
        file "*_R1.fq.gz" into fastq_toremap1
        file "*_R2.fq.gz" into fastq_toremap2
        file "*_mate1_vir.paf" into viral1_paf
        file "*_mate2_vir.paf" into viral2_paf

        script:

        """
		minimap2 -x sr -c ${pref} ${reads[0]} -t ${params.cpu.minimap} > ${name}_mate1_vir.paf
		minimap2 -x sr -c ${pref} ${reads[1]} -t ${params.cpu.minimap} > ${name}_mate2_vir.paf
		cut -f1 ${name}_mate1_vir.paf > ${name}_vir.lst
		cut -f1 ${name}_mate2_vir.paf >> ${name}_vir.lst
		less -S ${name}_vir.lst | sort | uniq > ${name}_vir_mapping.lst
		seqtk subseq ${reads[0]} ${name}_vir_mapping.lst | gzip -vc > ${name}_vir_R1.fq.gz
		seqtk subseq ${reads[1]} ${name}_vir_mapping.lst | gzip -vc > ${name}_vir_R2.fq.gz
	 """
}

// STEP 2 - MAPPING ALL READS WITH ANY VIRAL OVERLAP ONTO THE HOST GENOME

host_mini = Channel.fromPath(params.host_seq)
        .ifEmpty { exit 1, "Host genome reference fasta file not found: please check ${params.host_seq}" }



process minimap_genome_and_detectIS {
        errorStrategy 'retry'
        maxRetries 3
        scratch params.scratch

        cpus = params.cpu.minimap
        memory = 20.GB
        tag "$sample"
        publishDir "${params.outdir}", mode: 'copy' ,
	saveAs: {filename ->
        		if (filename.indexOf("_gnm.paf") > 0) "MinimapPAF/$filename"
                        else "detectIS/$filename"
        }

        input:
        file read1 from fastq_toremap1
        file read2 from fastq_toremap2
        file gref from host_mini.collect()
	file vpaf1 from viral1_paf
        file vpaf2 from viral2_paf	
	
        output:
	file "*_mate1_gnm.paf" into genom1_paf
        file "*_mate2_gnm.paf" into genom2_paf
	file("*.md")
        file("*.txt")

        script:	
	sample = read1.toString() - ~/(_vir)?(_R1)?(\.fq)?(\.gz)?$/

        """
                minimap2 -x sr -c ${gref} ${read1} -t ${params.cpu.minimap} > ${sample}_mate1_gnm.paf
                minimap2 -x sr -c ${gref} ${read2} -t ${params.cpu.minimap} > ${sample}_mate2_gnm.paf
		perl $baseDir/bin/detectIS.pl -h1 ${sample}_mate1_gnm.paf -h2 ${sample}_mate2_gnm.paf -v1 ${vpaf1} -v2 ${vpaf2} -o ${sample}
	"""
}

