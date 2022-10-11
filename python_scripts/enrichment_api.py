import os
import json
import pandas as pd
import requests
from pathlib import Path

def extract_biological_process(list_genes : list, toppgene_name : str, panther_name : str, logger = None, pvalue = 0.05, maxres = 100):

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
		print('Done.\n')

		# Get and process json answer
		with open(lookup, 'r') as f:
			json_data = json.loads(f.readlines()[0])
		genes = set()
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

		# Get and process json answer
		d = {}
		compt = 0
		with open(ToppGeneEnrichment, 'r') as f:
			json_data = json.loads(f.readlines()[0])
	
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

		##############
		#Pantherdb API
	
		print(f'Request Pantherdb...\n')
		if logger:
			logger.info('Panthere enrichment from API')
		print('Panthere enrichment from API')
		headers = {'accept': 'application/js'}
		genes = ','.join([str(gene) for gene in list_genes])

		print(r""" curl -X POST -H 'accept: application/js' "http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="""+genes+r"""&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR" \n""")
		try:
			response = requests.post("http://pantherdb.org/services/oai/pantherdb/enrich/overrep?geneInputList="+genes+"&organism=9606&annotDataSet=GO%3A0008150&enrichmentTestType=FISHER&correction=FDR", headers=headers, timeout=45)
	
			json_data = json.loads(response.text) # Get json data
			print('Done.\n')
	
			if 'search' in json_data.keys() and 'error' in json_data['search'].keys():			# Search for errors in the answer
				raise ValueError('Something went wrong with the Panther request (POST) !')

			#Output files verification
	
			xlsx = panther_name	
			if not panther_name.endswith('.xlsx'):
				xlsx = Path(panther_name).with_suffix('.xlsx')
			tsv = panther_name
			if not panther_name.endswith('.tsv'):
				tsv = Path(panther_name).with_suffix('.tsv')

			if logger:
				logger.info(f'Write Panthere enrichment in {xlsx} and {tsv}')
			print(f'Write Panthere enrichment in {xlsx} and {tsv}')

			# Data processing

			d = {}
			compt = 0
			nb_unclassified_GO=0
			col = ["Id", "name", "PValue", "FDR", "Genes number"]
			for l in json_data['results']['result']:
				if float(l['pValue']) < pvalue and int(l['number_in_list']) != 0:
					compt += 1
					if l["term"]["label"] != 'UNCLASSIFIED':
						new = [l["term"]["id"], l["term"]["label"], l["pValue"], l["fdr"], l["number_in_list"]]
						d[compt]= new
					elif logger:
						nb_unclassified_GO+=1		
			if logger:
				logger.warning(f'{nb_unclassified_GO} unclassified gene onthology !')
			df = pd.DataFrame.from_dict(d, orient='index',  columns=col)
			df['Genes number'] = df['Genes number'].astype(int)
			df = df.sort_values(by='Genes number', ascending=False)
			df = df.sort_values(by=['PValue','FDR'])
			df.index = [i for i in range(len(df.index))]
		
			df.to_excel(xlsx)
			df.to_csv(tsv, sep='\t')
	
		except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as err:
			print ('Server taking too long. Try again later')
		except ValueError as valerr:
			print('Error in the post request or the API has changed !')
	


#extract_biological_process(['BRCA1','BRCA2'], 'yt', 'yp')





