import os
import json
import pandas as pd
import requests
from pathlib import Path

def ToppGene_GEOA(list_genes : list, toppgene_name : str, logger = None, pvalue = 0.05, maxres = 100):

	if len(list_genes) == 0:
		if logger:
			logger.warning('No gene in list !')
	else:

		if logger:
			logger.info('Extraction of the biological process from the list of genes')
		print('Extraction of the biological process from the list of genes...')

		lookup='/tmp/lookup.txt'
		ToppGeneEnrichment='/tmp/ToppGeneEnrichment.txt'

		#############
		#ToppGene API

		if logger:
			logger.info('Get gene id using ToppGene API')
		print('Get gene id using ToppGene API')

		# Lookup for genes and get their Entrez id

		lookup_genes = "\""+str("\",\"".join(list_genes))+"\""
		cURL = r"""curl -H 'Content-Type: text/json' -d '{"Symbols":["""+lookup_genes+r"""]}' -X POST https://toppgene.cchmc.org/API/lookup > """+lookup
		print(f' Request : {cURL}\n')
		response = os.system(cURL)


		# Get and process json answer
		json_data = None
		with open(lookup, 'r') as f:
			json_data = json.loads(f.readlines()[0])
		genes = set()

	
		if json_data['Genes']:
			for l in json_data['Genes']:
				for k, v in l.items():
					if k == 'Entrez':
						genes.add(v)
			genes = list(genes)
			print(f'Genes for enrichment {genes}')

			if logger:
				logger.info(f'ToppGene enrichment from API')


			# Get infos from ToppGene API
		
			cURL = r"""curl -H 'Content-Type: text/json' -d '{"Genes":"""+str(genes)+r""", "Type": "GeneOntologyBiologicalProcess", "PValue": """+str(pvalue)+r""", "MinGenes": 1, "MaxGenes": 1500, "MaxResults": """+str(maxres)+r""", "Correction": "FDR"}' https://toppgene.cchmc.org/API/enrich > """+ToppGeneEnrichment
			print(f' Request ToppGene...\n{cURL}\n')
			response = os.system(cURL)
			print('Done.\n')
			print(response)

			# Get and process json answer
			d = {}
			compt = 0
			json_data = None

			with open(ToppGeneEnrichment, 'r') as f:
				try:
					results_data = f.readlines()[0]
					json_data = json.loads(results_data)
				except json.decoder.JSONDecodeError as jsonerr:
					print('Error while reading the json from ToppGene, possibly to much genes in the request')

			if json_data:
				for l in json_data['Annotations']:
					l['URL']='https://www.ebi.ac.uk/QuickGO/term/'+l['ID']
					for k, v in l.items():
						if k == 'Category' and v == 'GeneOntologyBiologicalProcess':
							compt+=1
							d[compt]=l
	
				if d != {}:
					df = pd.DataFrame.from_dict(d).T
					del df['Category']
					del df['Source']

					#Output files verification

					xlsx = toppgene_name
					if not toppgene_name.endswith('.xlsx'):
						xlsx = Path(toppgene_name).with_suffix('.xlsx')
					tsv = toppgene_name
					if not toppgene_name.endswith('.tsv'):
						tsv = Path(toppgene_name).with_suffix('.tsv')

					if logger:
						logger.info(f'Write ToppGene enrichment in {xlsx} and {tsv}')
					print(f'Write ToppGene enrichment in {xlsx} and {tsv}')
					df.to_excel(xlsx)
					df.to_csv(tsv, sep='\t')
				else:
					if logger:
						logger.info('ToppGene don\'t find biological process corresponding to this gene list')
					print('ToppGene don\'t find biological process corresponding to this gene list')


		else:
			print('No valid Gene !')














