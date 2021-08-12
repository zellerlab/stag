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
	tuple val(gene), path("${gene}/${gene}.ali"), emit: alignment

	script:

	if (seqs.size() == 2) {
		"""
		echo "PROTEIN"
		mkdir -p ${gene}
		touch ${gene}/${gene}.ali
		""" 
	} else {
		"""
		echo "NO PROTEIN"
		mkdir -p ${gene}
		touch ${gene}/${gene}.ali
		"""       	
	}

}


process concat_alignments {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(alignments)

	output:
	stdout
	tuple val("genome"), path("concat_alignment.txt"), emit: alignment

	script:
	"""
	echo $alignments
	touch concat_alignment.txt
	"""

}



process learn_function {

	input:
	tuple val(gene), path(alignment), val(level)

	output:
	stdout
	tuple val(gene), path("${gene}.${level}.lfunc.dat"), emit: lfunc

	script:
	"""
	touch "${gene}.${level}.lfunc.dat"
	"""
}


process train_classifiers {

	input:
	tuple val(gene), path(alignment)

	output:
	stdout
	tuple val(gene), path("${gene}.classifiers.dat"), emit: classifiers

	script:
	"""
	touch "${gene}.classifiers.dat"
	"""
}


process save_db {

	input:
	tuple val(gene), path(lfunc), path(classifiers)
	
	output:
	stdout
	tuple val(gene), path("${gene}.stagDB"), emit: db

	script:
	"""
	touch ${gene}.stagDB
	"""
}


workflow {
	
	/* 
		1. We expect marker gene (+ opt. protein) sequences arranged in --input_dir via the pattern
	       `<marker_gene_id>.(fna|faa)`.

		   The input channel will crawl the whole --input_dir subtree.
	*/
	
	seq_ch = Channel
		.fromPath(params.input_dir + "/" + "**.{faa,fna}")
		.map { file ->
			def gene = file.name.replaceAll(/.(faa|fna)$/, "")
			return tuple(gene, file)
		}
		.groupTuple(sort: true)

	seq_ch.view()

	/*
		2. Align the marker gene sets need to be against gene-specific hmms (provided via config file).
	*/

	align_marker_genes(seq_ch)
	aln_ch = align_marker_genes.out.alignment

	/*
		3. Concatenate the individual alignments into a genome alignment.
	*/

	all_aln_ch = aln_ch 
		.map { rec -> return rec[1] }
		.collect()

	concat_alignments(all_aln_ch)

	/*
		4. Train the classifiers on the individual marker gene alignments as well as on the genome alignment.
	*/	

	clf_input_ch = aln_ch.concat(concat_alignments.out.alignment)
	train_classifiers(clf_input_ch)

	/*
		5. Estimate the learning function on the individual marker gene alignments as well as on the genome alignment.
	*/

	levels = Channel.of(0..6)  // levels could come via params/config

	genome_db_learn_input_ch = concat_alignments.out.alignment.combine(levels)
	learn_input_ch = aln_ch.combine(levels).concat(genome_db_learn_input_ch)
	learn_function(learn_input_ch)

	lf_combine_ch = learn_function.out.lfunc
		.groupTuple(sort: true)

	/*
		6. Combine the classifier trainings and learning function estimates and save the database.
	*/

	lf_clf_combine_ch = lf_combine_ch.join(train_classifiers.out.classifiers)

	save_db(lf_clf_combine_ch)
	save_db.out.db.view()

}
