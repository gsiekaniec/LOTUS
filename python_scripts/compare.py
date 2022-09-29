#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import logging
from python_scripts.check_files import verif_input_vcf, verif_output, verif_input_config

def read_config_file(config):
	verif_input_config(config)
	with open(config, 'r') as c:
		for line in c:
			if line[0] != '#':
				print(line)	



def main(args):

	input_files = read_config_file(args.config)
        
	#working_directory = Path(vcf_file_pass).parent.absolute()
        
	genome_file = args.genome
	
	out_stats = args.stats
	out_genes = args.genes
	out_profile = args.profile
	out_indel = args.indel
	enrichment = args.enrichment	


	# Logger configuration

	logger = logging.getLogger('g-LOTUS compare')
	logger.setLevel(logging.DEBUG)

	fh = logging.FileHandler(args.log)
	fh.setLevel(logging.DEBUG)

	formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
	fh.setFormatter(formatter)
	
	logger.addHandler(fh)
