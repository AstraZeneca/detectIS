![Maturity level-Prototype](https://img.shields.io/badge/Maturity%20Level-Prototype-red)

# detectIS

DetectIS is a pipeline specifically designed to detect exogenous DNA integration sites using DNA or RNA paired-end sequencing data.
The workflow manager [nextflow](https://www.nextflow.io/) is used with a configuration file and a Singularity image 


## Getting Started

In order to run the workflow, the user has to create a configuration file, specifying:

	 	a)fasta file with the reference host genome;
		b)fasta file with the reference exogenous sequence;
		c)the directory containing the raw data, in FASTQ format
		d)the output directory. 
The analysis can be executed locally or in an HPC environment, in the latter scenario the user has also to specify the cluster executor. 


### Prerequisites

The detectIS software requirements are:
	- [Singularity](https://www.sylabs.io/docs/) V2.6 or higher.
	- [Nextflow](https://www.nextflow.io/), the workflow has been developed and tested with version 0.32.0.4897 


### Creating a Singularity container

A Singularity container with all the necessary software is required to run the pipeline.
The image can be created by using the recipe (file: "detectIS.rec" contained in the  "utils" directory). Superuser privileges are necessary to generate a Singularity container with the command:

```
sudo singularity build detectIS.simg detectIS.rec
```

N.B. superuser privileges are necessary only to create the container but no to use it. This means you can create the container in your local pc/workstation and copy it to the system where you run analyses (e.g. your hpc or cluster). 

Alternatively, If you have problems in generating a Singularity container from the recipe you can download the image from [Singularity Hub](https://singularity-hub.org/)  


### Runnig the workflow

If you have installed Singularity, Nextflow, and [configured the Singularity](https://www.sylabs.io/guides/2.6/user-guide/faq.html?highlight=disk%20access#how-are-external-file-systems-and-paths-handled-in-a-singularity-container) granting the image access to the disk partitions to read and write you can run any workflow.

```
nextflow run Workflows/detectIS.nf -c detectIS_TestDataset.conf -with-report detectIS_TestDataset_nextflow_report.html
```

In the example Workflows/detectIS.nf is the workflow for the detectIS analysis and detectIS_TestDataset.conf is the configuration file with all the information needed for that given project. In the configuration file are specified input and output file directories, references (fasta) directories, and cluster specific parameters. 


### Test data sets

In the directory "TestDataset" are contained paired-end reads and reference files to run a detectIS analysis.
The dataset simulates the integration of a plasmid in the genome of Chinese hamster ovary cell line (CHOK1) .

The analysis can be executed using the bash script "Run_detectIS.sh", also contained in the directory "TestDataset" or using nexflow:

```
nextflow run Workflows/detectIS.nf -c detectIS_TestDataset.conf -with-report detectIS_TestDataset_nextflow_report.html
```

The configuration file and the bash script can be either used as a template for other analyses. 


## Deployment

Please notice that Singularity containers can be [kernel-dependent](https://www.sylabs.io/guides/2.6/user-guide/faq.html?highlight=disk%20access#are-singularity-containers-kernel-dependent), this implies that the image recipies contained in this project will not necessarily produce an image able to run on your HPC system. If none of the available images is compatible with your system you might need to modify the recipe using an OS with compatible kernel, please raise an issue if this is the case and you need support for it.

## Citation

If you use detectIS in your research, please cite our latest [publication](https://doi.org/10.1093/bioinformatics/btab366).
