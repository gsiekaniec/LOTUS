import os
import json
import pandas as pd
import requests
from pathlib import Path


def Panther_GEOA(list_genes : list, panther_name : str, logger = None, pvalue = 0.05, maxres = 100):

	###############
	#Pantherdb API

	if len(list_genes) == 0:
		if logger:
			logger.warning('No gene in list !')
	else:

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

			###HERE

			if 'search' in json_data.keys() and 'error' in json_data['search'].keys():                      # Search for errors in the answer
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


























