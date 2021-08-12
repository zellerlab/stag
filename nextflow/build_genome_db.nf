#!/usr/bin/env nextflow

nextflow.enable.dsl=2


if (!params.publish_mode) {
	params.publish_mode = "symlink"
}

if (!params.output_dir) {
	params.output_dir = "stag_genome"
}

output_dir = "${params.output_dir}"


process align_marker_genes {
    publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(gene), path(seqs)

	output:
	stdout
	val(gene), emit: gene
	path("${gene}/${gene}.ali"), emit: alignment

	script:

	if (seqs.size() == 2) {
		"""
		echo "PROTEIN"
		mkdir -p ${gene}
		stag align -t $task.cpus -i ${gene}.fna -p ${gene}.faa -x ${params.taxonomy} -a ${params.hmmlib}/${gene}.hmm -o ${gene}/${gene}.ali
		""" 
	} else {
		"""
		echo "NO PROTEIN"
		mkdir -p ${gene}
		stag align -t $task.cpus -i ${gene}.fna -x ${params.taxonomy} -a ${params.hmmlib}/${gene}.hmm -o ${gene}/${gene}.ali
		"""       	
	}

}


process create_marker_dbs {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	tuple val(gene), path(seqs)
	path(alignment)

	output:
	stdout
	val(gene), emit: gene
	path("${gene}/${gene}.stagDB"), emit: stag_db
	path("${gene}/${gene}.cross_val"), emit: cross_val

	script:

	if (seqs.size() == 2) {
		"""
		mkdir -p ${gene}
		stag create_db -t $task.cpus -s ${alignment} -x ${params.taxonomy} -a ${params.hmmlib}/${gene}.hmm -o ${gene}/${gene}.stagDB -C ${gene}/${gene}.cross_val -p ${gene}.faa
		"""
	} else {
		"""
		mkdir -p ${gene}
		stag create_db -t $task.cpus -s ${alignment} -x ${params.taxonomy} -a ${params.hmmlib}/${gene}.hmm -o ${gene}/${gene}.stagDB -C ${gene}/${gene}.cross_val
		"""
	}
	


}


process concat_alignments {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(alignments)

	output:
	stdout
	path("concat_alignment.txt"), emit: alignment

	script:
	"""
	concat_alignment ${alignments} > concat_alignment.txt
	"""

}


process create_genome_db {
	publishDir "$output_dir", mode: params.publish_mode
	
	input:
	path(concat_alignment)

	output:
	stdout
	path("concat_alignment_db"), emit: stag_db
	path("concat_alignment_db.classifiers.dat"), emit: classifiers
	path("concat_alignment_db.taxfunc.dat"), emit: taxfunc
	path("concat_alignment_db.cross.txt"), emit: cv

	script:
	"""
	touch dummy.txt
	stag create_db -t $task.cpus -s ${concat_alignment} -o concat_alignment_db -C concat_alignment_db.cross.txt -x ${params.taxonomy} -a dummy.txt 
	"""
		
}


workflow {
	
	seq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{faa,fna}")
		.map { file ->
			def sample = file.name.replaceAll(/.(faa|fna)$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)

	seq_ch.view()


	align_marker_genes(seq_ch)
	//align_marker_genes.out.view()
	create_marker_dbs(seq_ch, align_marker_genes.out.alignment)

	concat_alignments(align_marker_genes.out.alignment.collect())
	create_genome_db(concat_alignments.out.alignment)
}
/*

	fastq_ch = Channel
    	.fromPath(params.input_dir + "/" + "**.{fastq,fq,fastq.gz,fq.gz}")
        .map { file ->
				def sample = file.name.replaceAll(/.(fastq|fq)(.gz)?$/, "")
				sample = sample.replaceAll(/_R?[12]$/, "")
				return tuple(sample, file)
		}
		.groupTuple()
	//fastq_ch.view()

	bam_ch = Channel
		.fromPath(params.input_dir + "/" + "**.bam")
		.map { file ->
			def sample = file.name.replaceAll(/.bam$/, "")
			return tuple(sample, file)
		}
		.groupTuple(sort: true)
	//bam_ch.view()

	make_dummy_fastqs(fastq_ch)
	make_dummy_bam(bam_ch)
	
	bam2fq(bam_ch)
	fq2bam(fastq_ch)

	if (params.bam_input) {
		count_reads(make_dummy_bam.out.sample, make_dummy_bam.out.bam)
		pathseq(make_dummy_bam.out.sample, make_dummy_bam.out.bam)
		if (bam2fq.out.r2 != null) {
			kraken2_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			motus_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			mtag_extraction_paired(bam2fq.out.sample, bam2fq.out.r1, bam2fq.out.r2)
			mapseq_paired(mtag_extraction_paired.out.sample, mtag_extraction_paired.out.bac_lsu_r1, mtag_extraction_paired.out.bac_ssu_r1, mtag_extraction_paired.out.bac_lsu_r2, mtag_extraction_paired.out.bac_ssu_r2)
			collate_mapseq_paired(mapseq_paired.out.bac_lsu_r1.collect(), mapseq_paired.out.bac_ssu_r1.collect(), mapseq_paired.out.bac_lsu_r2.collect(), mapseq_paired.out.bac_ssu_r2.collect())
		} else {
			kraken2_single(bam2fq.out.sample, bam2fq.out.r1)
			motus_single(bam2fq.out.sample, bam2fq.out.r1)
			mtag_extraction_single(bam2fq.out.sample, bam2fq.out.r1)
			mapseq_single(mtag_extraction_single.out.sample, mtag_extraction_single.out.bac_lsu_r1, mtag_extraction_single.out.bac_ssu_r1)
			collate_mapseq_single(mapseq_single.out.bac_lsu_r1.collect(), mapseq_single.out.bac_ssu_r1.collect())
		}
	} else {
		count_reads(fq2bam.out.sample, fq2bam.out.bam)
		pathseq(fq2bam.out.sample, fq2bam.out.bam)
		if (make_dummy_fastqs.out.r2 != null) {
			kraken2_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
			motus_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
            mtag_extraction_paired(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1, make_dummy_fastqs.out.r2)
			mapseq_paired(mtag_extraction_paired.out.sample, mtag_extraction_paired.out.bac_lsu_r1, mtag_extraction_paired.out.bac_ssu_r1, mtag_extraction_paired.out.bac_lsu_r2, mtag_extraction_paired.out.bac_ssu_r2)
			collate_mapseq_paired(mapseq_paired.out.bac_lsu_r1.collect(), mapseq_paired.out.bac_ssu_r1.collect(), mapseq_paired.out.bac_lsu_r2.collect(), mapseq_paired.out.bac_ssu_r2.collect())
		} else {
			kraken2_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			motus_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			mtag_extraction_single(make_dummy_fastqs.out.sample, make_dummy_fastqs.out.r1)
			mapseq_single(mtag_extraction_single.out.sample, mtag_extraction_single.out.bac_lsu_r1, mtag_extraction_single.out.bac_ssu_r1)
			collate_mapseq_single(mapseq_single.out.bac_lsu_r1.collect(), mapseq_single.out.bac_ssu_r1.collect())
		}
	}


}
*/
