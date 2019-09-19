workflow runDx {

	call CallableLoci
	call FilterCallableLoci {
	    input:
	        filtered_bed = CallableLoci.filtered_bed
	}
    call Dx {
        input:
            filtered_cds = FilterCallableLoci.filtered_cds
    }
	
	output {
	    File TMB = Dx.mutation_burden
	    File mutSig = Dx.signatures
	}
	
	meta {
		author: "Sehyun Oh"
        email: "shbrief@gmail.com"
        description: "This workflow extracts copy number and mutation metrics (e.g. TMB and mutSig). CallableLoci task generates a BED file with callable regions. In this workflow, mutation burden calculation is restricted to coding sequences."
    }
}

task CallableLoci {
    File tumor_bam
    File tumor_bai
    File ref_fasta
    File ref_fai
    File ref_dict
    Int minDepth   # 30
    String gatk_docker   # broadinstitute/gatk3:3.8-1
    String SAMPLEID = basename(tumor_bam, ".bam")
    
	command <<<
		java -jar -Xmx4g /usr/GenomeAnalysisTK.jar \
        -T CallableLoci \
        -R ${ref_fasta} \
        -I:tumor ${tumor_bam} \
        --summary ${SAMPLEID}_table.txt \
        -o ${SAMPLEID}_callable_status.bed \
        --minDepth ${minDepth}
        
        grep CALLABLE ${SAMPLEID}_callable_status.bed > \
        ${SAMPLEID}_callable_status_filtered.bed
	>>>
	
	output {
        File filtered_bed = "${SAMPLEID}_callable_status_filtered.bed"
	}
	
	runtime {
		docker: gatk_docker
		cpu : 4
		memory: "8 GB"
	}
}

task FilterCallableLoci {
    File filtered_bed
    String genome
    
	String fname = basename(filtered_bed)
	String SAMPLEID = sub(fname, "_callable_status_filtered.bed", "")

	command <<<
		Rscript /usr/local/lib/R/site-library/PureCN/extdata/FilterCallableLoci.R \
        --genome ${genome} \
        --infile ${filtered_bed} \
        --outfile ${SAMPLEID}_callable_status_filtered_cds.bed 
	>>>
	
	output {
        File filtered_cds = "${SAMPLEID}_callable_status_filtered_cds.bed"
	}
	
	runtime {
		docker: "quay.io/shbrief/pcn_docker"
		cpu : 4
		memory: "8 GB"
	}
}

task Dx {
    File filtered_cds
    File simpleRepeats
    File resRDS
    
	String fname = basename(filtered_cds)
	String SAMPLEID = sub(fname, "_callable_status_filtered_cds.bed", "")

	command <<<
	    R -e 'BiocManager::install("TxDb.Hsapiens.UCSC.hg19.knownGene")'

		Rscript /usr/local/lib/R/site-library/PureCN/extdata/Dx.R \
        --out ${SAMPLEID} \
        --rds ${resRDS} \
        --callable ${filtered_cds} \
        --exclude ${simpleRepeats}
        --signatures
	>>>
	
	output {
		File mutation_burden = "${SAMPLEID}_mutation_burden.csv"
		File signatures = "${SAMPLEID}_signatures.csv"
	}
	
	runtime {
		docker: "quay.io/shbrief/pcn_docker"
		cpu : 4
		memory: "8 GB"
	}
}
