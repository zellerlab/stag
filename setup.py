# coding: utf-8
from setuptools import setup, find_packages
from setuptools.extension import Extension
from distutils.extension import Extension
from codecs import open
from os import path
import glob
import re
import sys


here = path.abspath(path.dirname("__file__"))

with open(path.join(here, "DESCRIPTION.md"), encoding="utf-8") as description:
	description = long_description = description.read()

	name="stag"
	version = [line.strip().split(" ")[-1] for line in open("stag/__init__.py") if line.startswith("__version__")][0]

	if sys.version_info.major != 3:
		raise EnvironmentError("""{toolname} is a python module that requires python3, and is not compatible with python2.""".format(toolname=name))

	setup(
		name=name,
		version=version,
		description="stag - Supervised Taxonomic Assignment of marker Genes", 
		long_description=long_description,
		url="https://github.com/zellerlab/stag",
		author="Alessio Milanese",
		author_email="milanese.alessio@gmail.com",
		license="GPLv3+",
		classifiers=[
			"Development Status :: 4 - Beta",
			"Topic :: Scientific Engineering :: Bio/Informatics",
			"License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
			"Operating System :: POSIX :: Linux",
			"Programming Language :: Python :: 3.7"
		],
		zip_safe=False,
		keywords="",
		packages=find_packages(exclude=["test"]),
		install_requires=[
			"h5py",
			"regex",
			"scikit-learn",
			"numpy",
			"pandas"
		],
		entry_points={
			"console_scripts": [
				"stag=stag.__main__:main",
				"stag_test=stag.stag_test:main",
				"concat_alignment=stag.concat_alignment:main",
				"train_classifiers=stag.train_classifiers:main",
				"learn_function=stag.learn_function:main",
				"save_db=stag.save_db:main",
			],
		},
		scripts=["nextflow/build_genome_db.nf"],
		package_data={
			"stag.test": ["gene.hmm", "sequences.fasta", "sequences.taxonomy"]
		},
		include_package_data=True,
		data_files=[],
	)
