#!/usr/bin/env python
import os
import pathlib
import shutil
import subprocess
import sys
import tempfile

from . import __version__ as tool_version

import stag.align as align
import stag.create_db as create_db
import stag.classify as classify
import stag.check_create_db_input_files as check_create_db_input_files
import stag.correct_seq as correct_seq
import stag.unzip_db as unzip_db
import stag.classify_genome as classify_genome
import stag.train_genome as train_genome
import stag.convert_ali as convert_ali

from .helpers import check_file_exists, check_file_doesnt_exists

from stag.classify_genome import validate_genome_files

from stag.handle_args import (
    get_args, handle_error,
    print_menu_align, print_menu_create_db, print_menu_train, print_menu_classify, print_menu_check_input,
    print_menu_correct_seq, print_menu_convert_ali, print_menu_unzip_db, print_menu_train_genome,
    print_menu_classify_genome
    )


def run_test():
    popenCMD = "stag_test"
    child = subprocess.Popen(popenCMD)
    child.communicate()
    rc = child.wait()
    return(rc)


def run_align(args):
    # check that '-i' and '-a' have been provided
    if not args.fasta_input:
        handle_error("missing <seqfile> (-i)", print_menu_align)
    if not args.template_al:
        handle_error("missing <hmmfile>/<cmfile> (-a)", print_menu_align)

    check_file_exists(args.fasta_input, isfasta=True)
    check_file_exists(args.template_al, isfasta=False)

    if args.protein_fasta_input:
        check_file_exists(args.protein_fasta_input, isfasta=True)

    if args.output is None:
        for ali in align.align_generator(
            args.fasta_input, args.protein_fasta_input, args.template_al,
            args.use_cm_align, args.threads, args.verbose, False, args.min_perc_state
        ):
            print(ali)
    else:
        align.align_file(
            args.fasta_input, args.protein_fasta_input, args.template_al, args.use_cm_align, args.threads,
            args.verbose, args.output, args.min_perc_state
        )


def run_create_db(args):
    if not args.aligned_sequences:
        # check that '-s' has been provided (alignment produced by stag align)
        handle_error("missing <aligned_file> (-s)", print_menu_create_db)
    if not args.taxonomy:
        # check that '-x' has been provided (taxonomy file)
        handle_error("missing <taxonomy_file> (-x)", print_menu_create_db)
    if not args.template_al:
        # check that the hmm file is provided
        handle_error("missing <hmmfile>/<cmfile> (-a)", print_menu_create_db)
    if not args.output:
        # check that output is set
        handle_error("missing <output_DB> (-o)", print_menu_create_db)

    check_file_exists(args.aligned_sequences, isfasta=False)
    check_file_exists(args.taxonomy, isfasta=False)
    check_file_exists(args.template_al, isfasta=False)

    if not args.force_rewrite:
        check_file_doesnt_exists(args.output)

    create_db.create_db(
        args.aligned_sequences, args.taxonomy, args.verbose, args.output, args.use_cm_align,
        args.template_al, args.intermediate_cross_val, args.protein_fasta_input,
        args.penalty_logistic, args.solver_logistic, max_iter=args.solver_iterations,
        procs=args.threads
    )


def run_train(args):
    # check that '-i' and '-a' have been provided
    if not args.fasta_input:
        handle_error("missing <seqfile> (-i)", print_menu_train)
    elif not args.template_al:
        handle_error("missing <hmmfile>/<cmfile> (-a)", print_menu_train)
    elif not args.taxonomy:
        # check that '-x' has been provided (taxonomy file)
        handle_error("missing <taxonomy_file> (-x)", print_menu_train)
    elif not args.output:
        # check that output is set
        handle_error("missing <output_DB> (-o)", print_menu_train)

    check_file_exists(args.fasta_input, isfasta=True)
    check_file_exists(args.template_al, isfasta=False)

    if args.protein_fasta_input:
        check_file_exists(args.protein_fasta_input, isfasta=True)

    check_file_exists(args.taxonomy, isfasta=False)

    if not args.force_rewrite:
        check_file_doesnt_exists(args.output)

    # we create a temporary file that will contain the alignments
    al_file = tempfile.NamedTemporaryFile(delete=False, mode="w")
    os.chmod(al_file.name, 0o644)

    with al_file:
        align.align_file(
            args.fasta_input, args.protein_fasta_input, args.template_al, args.use_cm_align,
            args.threads, args.verbose, al_file.name, args.min_perc_state
        )

        create_db.create_db(
            al_file.name, args.taxonomy, args.verbose, args.output, args.use_cm_align,
            args.template_al, args.intermediate_cross_val, args.protein_fasta_input,
            args.penalty_logistic, args.solver_logistic, max_iter=args.solver_iterations,
            procs=args.threads
        )

    # what to do with intermediate alignment -------------------------------
    if not args.intermediate_al:
        # remove it
        os.remove(al_file.name)
    else:
        # save it
        shutil.move(al_file.name, args.intermediate_al)


def run_classify(args):
    # check that '-i' has been provided (alignment produced by stag align)
    if not args.fasta_input and not args.aligned_sequences:
        handle_error("missing <fasta_seqs> (-i) or <aligned_seq> (-s)", print_menu_classify)
    if not args.database:
        # check that '-d' has been provided (taxonomy file)
        handle_error("missing <database> (-d)", print_menu_classify)

    check_file_exists(args.database, isfasta=False)

    if args.fasta_input:
        check_file_exists(args.fasta_input, isfasta=True)

    if args.protein_fasta_input:
        check_file_exists(args.protein_fasta_input, isfasta=True)

    # if -S is provided, we remove the file if it exists, since in the
    # function it appends only
    if args.intermediate_al and os.path.isfile(args.intermediate_al):
        os.remove(args.intermediate_al)

    classify.classify(
        args.database, fasta_input=args.fasta_input, protein_fasta_input=args.protein_fasta_input,
        verbose=args.verbose, threads=args.threads, output=args.output, long_out=args.long_out,
        current_tool_version=tool_version, aligned_sequences=args.aligned_sequences,
        save_ali_to_file=args.intermediate_al, min_perc_state=args.min_perc_state
    )


def run_check_input(args):
    if not args.fasta_input:
        handle_error("missing <fasta_seqs> (-i)", print_menu_check_input)
    if not args.taxonomy:
        handle_error("missing <taxonomy_file> (-x)", print_menu_check_input)
    if not args.template_al:
        handle_error("missing <hmmfile>/<cmfile> (-a)", print_menu_check_input)

    check_create_db_input_files.check_input_files(
        args.fasta_input, args.protein_fasta_input, args.taxonomy,
        args.template_al, args.use_cm_align, args.warning_file_check_input
    )


def run_correct_seq(args):
    # check if the sequences are in correct orientation, if they are not, then
    # take reverse complement. Save to -o all the seqeunces is correct order

    # check that '-i' and '-a' have been provided
    if not args.fasta_input:
        handle_error("missing <seqfile> (-i)", print_menu_correct_seq)
    if not args.template_al:
        handle_error("missing <hmmfile>/<cmfile> (-a)", print_menu_correct_seq)

    check_file_exists(args.fasta_input, isfasta=True)
    check_file_exists(args.template_al, isfasta=False)

    correct_seq.correct_seq(
        args.fasta_input, args.template_al, args.use_cm_align, args.threads, args.verbose,
        args.min_perc_state, args.output
    )


def run_convert_ali(args):
    # check that '-i' and '-o' have been provided
    if not args.fasta_input:
        handle_error("missing <file_in> (-i)", print_menu_convert_ali)
    if not args.output:
        handle_error("missing <file_out> (-o)", print_menu_convert_ali)

    check_file_exists(args.fasta_input, isfasta=False)

    convert_ali.convert_ali(args.fasta_input, args.output, args.verbose)


def run_unzip_db(args):
    # check that '-d' and '-o' have been provided
    if not args.database:
        error = "missing <database> (-d)"
    if not args.output:
        error = "missing <dir_out> (-o)"

    if error:
        handle_error(error, print_menu_unzip_db)

    # check that '-d' is a file
    check_file_exists(args.database, isfasta=False)

    # call function
    unzip_db.unzip_db(args.database, args.verbose, args.output)


def run_train_genome(args):
    # We want to have a database for classifying genomes. It will contains the
    # calssifiers for the seingle genes
    # Input: single gene databases
    # Output: a database file (hdf5) that can be used by stag classify_genome
    if not args.output:
        handle_error("missing <output_DB> (-o)", print_menu_train_genome)
    if not args.fasta_input:
        handle_error("missing <list_gene_DBs> (-i)", print_menu_train_genome)
    if not args.file_thresholds:
        handle_error("missing <gene_thresholds> (-T)", print_menu_train_genome)
    if not args.intermediate_cross_val:
        handle_error("missing <concat_genes_DB> (-C)", print_menu_train_genome)

    train_genome.train_genome(
        args.output, args.fasta_input, args.file_thresholds,
        args.threads, args.verbose, args.intermediate_cross_val
    )


def run_classify_genome(args):
    if not args.database:
        handle_error("missing <database> (-d)", print_menu_classify_genome)
    if not any((args.fasta_input, args.dir_input, args.marker_genes)):
        handle_error("you need to provide at least -i, -D, or -G.", print_menu_classify_genome)
    if sum(map(bool, (args.fasta_input, args.dir_input, args.marker_genes))) != 1:
        handle_error("options -i, -D, and -G are mutually exclusive", print_menu_classify_genome)
    if args.dir_input and not os.path.isdir(args.dir_input):
        handle_error("-D is not a directory.", print_menu_classify_genome)
    if args.marker_genes and not os.path.isfile(args.marker_genes):
        handle_error("-G is not a valid file.", print_menu_classify_genome)
    if not args.output:
        handle_error("missing output directory (-o)", print_menu_classify_genome)

    # find files to classify
    marker_genes, list_files = [], []
    if args.fasta_input:
        check_file_exists(args.fasta_input, isfasta=True)
        list_files.append(args.fasta_input)
    elif args.marker_genes:
        marker_genes.append(args.marker_genes)
    else:
        for f in os.listdir(args.dir_input):
            f = os.path.join(args.dir_input, f)
            if os.path.isfile(f):
                try:
                    with open(f) as _in:
                        if _in.read(1)[0] == ">":
                            list_files.append(f)
                except Exception:
                    if args.verbose > 1:
                        sys.stderr.write("[W::main] Warning: ")
                        sys.stderr.write("Cannot open file: {}\n".format(f))

        if not list_files:
            handle_error("no fasta files found in the directory.", None)
        sys.stderr.write(" Found "+str(len(list_files))+" fasta files\n")

    if os.path.isdir(args.output):
        if args.force_rewrite:
            shutil.rmtree(args.output)
        else:
            handle_error("output directory (-o {}) exists already.".format(args.output), None)

    # create output dir
    try:
        pathlib.Path(args.output).mkdir(exist_ok=True, parents=True)
    except Exception:
        handle_error("creating the output directory (-o {}).".format(args.output), None)

    if list_files:
        validate_genome_files(list_files)

    classify_genome.classify_genome(
        args.database, genome_files=list_files, marker_genes=marker_genes,
        verbose=args.verbose, threads=args.threads,
        output=args.output, long_out=args.long_out, keep_all_genes=args.keep_all_genes
    )


def main(argv=None):

    args = get_args()

    if args.command == 'test':
        run_test()
    else:

        routines = {
            "align": run_align,
            "create_db": run_create_db,
            "train": run_train,
            "classify": run_classify,
            "check_input": run_check_input,
            "correct_seq": run_correct_seq,
            "convert_ali": run_convert_ali,
            "unzip_db": run_unzip_db,
            "train_genome": run_train_genome,
            "classify_genome": run_classify_genome,
        }

        # set defaults for the parameters
        if args.min_perc_state is None:
            args.min_perc_state = 5 if args.command == "correct_seq" else 0

        if args.threads < 1:
            handle_error("number of threads (-t) is less than 1", None)
        if args.min_perc_state < 0 or args.min_perc_state > 100:
            handle_error("-m should be between 0 and 100. It represents the percentage\n"
                         "of internal states covered by the sequence (i.e. the number of features).", None)

        run_routine = routines.get(args.command)
        if run_routine is None:
            raise ValueError(f"Subroutine {args.command} is unknown.")

        run_routine(args)


if __name__ == '__main__':
    main()
