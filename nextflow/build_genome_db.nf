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
	tuple val(gene), path("${gene}/${gene}.ali"), emit: alignment

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


process concat_alignments {
	publishDir "$output_dir", mode: params.publish_mode

	input:
	path(alignments)

	output:
	tuple val("genome"), path("concat_alignment.txt"), emit: alignment

	script:
	"""
	echo \$(ls $alignments)
    concat_alignment \$(ls ${alignments}) > concat_alignment.txt
	"""

}


process learn_function {

	input:
	tuple val(gene), path(alignment), val(level)

	output:
	tuple val(gene), path("${gene}.${level}.lfunc.dat"), emit: lfunc

	script:
	"""
	learn_function ${params.taxonomy} ${alignment} ${level} -o ${gene} -t $task.cpus
	"""
}


process train_classifiers {

	input:
	tuple val(gene), path(alignment)

	output:
	tuple val(gene), path("${gene}.classifiers.dat"), emit: classifiers

	script:
	"""
	train_classifiers ${params.taxonomy} ${alignment} -o ${gene} -t $task.cpus
	"""
}


process save_db {
	publishDir "${output_dir}/gene_databases", mode: params.publish_mode
	label "save_db"

	input:
	tuple val(gene), path(lfunc), path(classifiers), path(alignment)
	
	output:
	tuple val(gene), path("${gene}.stagDB"), emit: db
	path("${gene}.cross_val")

	script:
	def hmm_file = ("${gene}" == "genome") ? "dummy_hmm.txt" : "${params.hmmlib}/${gene}.hmm";
	"""
	touch dummy_hmm.txt
	touch protein_stuff.txt
	save_db ${params.taxonomy} ${alignment} ${classifiers} ${hmm_file} protein_stuff.txt \$(ls *.lfunc.dat) -o ${gene}
	"""
}


process save_genome_db {
	publishDir "${output_dir}", mode: params.publish_mode
	label "save_db"

	input:
	val(dbname)
	path(databases)
	path(thresholds)
	path(lengths)

	output:
	tuple val(dbname), path("db/${dbname}.stagDB"), emit: db

	script:
	"""
	mkdir -p db/
	join -1 1 -2 1 <(sort -k1,1 ${thresholds}) <(sort -k1,1 ${lengths}) | tr " " "\\t" > threshold_lengths.txt
	stag train_genome -i \$(ls ${databases} | grep -v '^genome' | tr "\\n" ',' | sed "s/,\$//") -T threshold_lengths.txt -C genome.stagDB -o db/${dbname}.stagDB
	"""
}


process extract_alignment_lengths {
	input:
	path(alignments)

	output:
	path("alignment_lengths.txt"), emit: align_lengths

	script:
	"""
	awk -F '\t' -v OFS='\t' '{ print FILENAME,\$1; nextfile }' \$(ls ${alignments}) | sed 's/\\.ali/.stagDB/' > alignment_lengths.txt
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
		2. Align the marker gene sets against gene-specific hmms (provided via config file).
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

	extract_alignment_lengths(all_aln_ch)

	/*
		4. Train the classifiers on the individual marker gene alignments as well as on the genome alignment.
	*/	

	clf_input_ch = aln_ch.concat(concat_alignments.out.alignment)
	clf_input_ch.view()
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

	lf_clf_combine_ch = lf_combine_ch.join(train_classifiers.out.classifiers.join(clf_input_ch))
	lf_clf_combine_ch.view()

	save_db(lf_clf_combine_ch)
	save_db.out.db.view()


	db_ch = save_db.out.db
		.map { gene, db -> return db }
		.collect()
	db_ch.view()

	save_genome_db(
		params.dbname,
		db_ch,
		params.thresholds,
		extract_alignment_lengths.out.align_lengths
	)
}
