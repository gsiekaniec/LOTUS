3#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from matplotlib import pyplot as plt
import pandas as pd
import re
import numpy as np
import pickle as pk
import gzip
from pathlib import Path
from collections import OrderedDict
from matplotlib.ticker import FixedLocator, FixedFormatter, FuncFormatter
from itertools import chain

COLOR_B1 = "darkorange"
COLOR_B2 = "black"
COLOR_SAMPLE = "cornflowerblue"
WIDTH = 20
LENGTH = 40
SIZE_B1 = 1.3
SIZE_B2 = 1.3
SIZE_SAMPLE = 1.3


COLOR_G_STAINING = {
	'gpos1': "#FDFDFD",
	'gpos2': "#FBFBFB",
	'gpos3': "#F8F8F8",
	'gpos4': "#F6F6F6",
	'gpos5': "#F3F3F3",
	'gpos6': "#F1F1F1",
	'gpos7': "#EEEEEE",
	'gpos8': "#ECECEC",
	'gpos9': "#E9E9E9",
	'gpos10': "#E6E6E6",
	'gpos11': "#E4E4E4",
	'gpos12': "#E1E1E1",
	'gpos13': "#DFDFDF",
	'gpos14': "#DCDCDC",
	'gpos15': "#DADADA",
	'gpos16': "#D7D7D7",
	'gpos17': "#D4D4D4",
	'gpos18': "#D2D2D2",
	'gpos19': "#CFCFCF",
	'gpos20': "#CDCDCD",
	'gpos21': "#CACACA",
	'gpos22': "#C8C8C8",
	'gpos23': "#C5C5C5",
	'gpos24': "#C3C3C3",
	'gpos25': "#C0C0C0",
	'gpos26': "#BDBDBD",
	'gpos27': "#BBBBBB",
	'gpos28': "#B8B8B8",
	'gpos29': "#B6B6B6",
	'gpos30': "#B3B3B3",
	'gpos31': "#B1B1B1",
	'gpos32': "#AEAEAE",
	'gpos33': "#ACACAC",
	'gpos34': "#A9A9A9",
	'gpos35': "#A6A6A6",
	'gpos36': "#A4A4A4",
	'gpos37': "#A1A1A1",
	'gpos38': "#9F9F9F",
	'gpos39': "#9C9C9C",
	'gpos40': "#9A9A9A",
	'gpos41': "#979797",
	'gpos42': "#949494",
	'gpos43': "#929292",
	'gpos44': "#8F8F8F",
	'gpos45': "#8D8D8D",
	'gpos46': "#8A8A8A",
	'gpos47': "#888888",
	'gpos48': "#858585",
	'gpos49': "#838383",
	'gpos50': "#808080",
	'gpos51': "#7D7D7D",
	'gpos52': "#7B7B7B",
	'gpos53': "#787878",
	'gpos54': "#767676",
	'gpos55': "#737373",
	'gpos56': "#717171",
	'gpos57': "#6E6E6E",
	'gpos58': "#6C6C6C",
	'gpos59': "#696969",
	'gpos60': "#666666",
	'gpos61': "#646464",
	'gpos62': "#616161",
	'gpos63': "#5F5F5F",
	'gpos64': "#5C5C5C",
	'gpos65': "#5A5A5A",
	'gpos66': "#575757",
	'gpos67': "#545454",
	'gpos68': "#525252",
	'gpos69': "#4F4F4F",
	'gpos70': "#4D4D4D",
	'gpos71': "#4A4A4A",
	'gpos72': "#484848",
	'gpos73': "#454545",
	'gpos74': "#434343",
	'gpos75': "#404040",
	'gpos76': "#3D3D3D",
	'gpos77': "#3B3B3B",
	'gpos78': "#383838",
	'gpos79': "#363636",
	'gpos80': "#333333",
	'gpos81': "#313131",
	'gpos82': "#2E2E2E",
	'gpos83': "#2C2C2C",
	'gpos84': "#292929",
	'gpos85': "#262626",
	'gpos86': "#242424",
	'gpos87': "#212121",
	'gpos88': "#1F1F1F",
	'gpos89': "#1C1C1C",
	'gpos90': "#1A1A1A",
	'gpos91': "#171717",
	'gpos92': "#141414",
	'gpos93': "#121212",
	'gpos94': "#0F0F0F",
	'gpos95': "#0D0D0D",
	'gpos96': "#0A0A0A",
	'gpos97': "#080808",
	'gpos98': "#050505",
	'gpos99': "#030303",
	'gpos100': "#000000",
	'gneg': "#FFFFFF",
	'acen': "#660033",
	'gvar': "#660099",
	'stalk': "#6600CC"
}


def extract_num(s, p, ret=0):
	search = p.search(s)
	if search:
		return int(search.groups()[0])
	else:
		return ret


def count_impacted_genes(pos, span, positions):
	counted_genes = 0
	pos1 = pos - (0.5 * span)
	pos2 = pos + (0.5 * span)
	for (start, end) in positions:
		if (pos1 <= start and pos2 >= start) or (pos1 <= end and pos2 >= end):
			counted_genes += 1
	return counted_genes


def get_value_for_chr(position, chr, length, span, step):
	positions = []
	counter = []
	for pos in range(int(span / 2), length - (int(span / 2)), step):
		positions.append(pos)
		counter.append((count_impacted_genes(pos, span, position[chr])))
	# add beginning of the plot
	positions = [i for i in range(0, int(span / 2), step)] + positions
	counter = [counter[0] for i in range(0, int(span / 2), step)] + counter
	# add ending of the plot
	positions = positions + [i for i in range(length - (int(span / 2)), length+step, step)]
	counter = counter + [counter[-1] for i in range(length - (int(span / 2)), length+step, step)]

	return positions, counter


def get_plot(position, chr, chromosomes_size, span, step):
	if chr in position.keys():
		positions, counter = get_value_for_chr(position, chr, chromosomes_size[chr], span, step)
		return positions, counter
	else:
		return [i for i in range(0, chromosomes_size[chr], step)],[0 for i in range(0, chromosomes_size[chr], step)]


def scientific(x, pos):
	'''Return a scientific notation of x'''
	return '%.2e' % x


def cross_product(maxi, maxi2, y):
	if maxi != 0:
		y2 = [(maxi2 * value)/maxi for value in y]
	else:
		y2 = [0 for value in y]
	return y2


def sorted_nicely(l):
    """ Sort the given iterable in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    return sorted(l, key = alphanum_key)


def create_plot(cytobands_infos, curves_infos_b1, curves_infos_b2, curves_infos_sample, chromosomes_size, counts, plot_size, plot_name, title = '', logger = None):

	if not str(plot_name).endswith('.svg'):
		plot_name = Path(plot_name).with_suffix('.svg')

	if logger:
                logger.info('- Creation of the chromosomes plot -')
	print(f'Creation of the chromosomes plot...')
	fig = plt.figure(figsize=(WIDTH, LENGTH))
	ax = fig.add_subplot()
	ax.tick_params(axis='both', which='major', labelsize=10.5, right=True, labelright=True, left=False, labelleft=False)
	ax.tick_params(axis='both', which='minor', labelsize=8.5)
	ax.tick_params(axis='x', labelsize=9)

	max_chr_size = max(chromosomes_size.values())
	chromosome_windows_size = 17
	chromosome_windows_low = list(range(chromosome_windows_size*24,-chromosome_windows_size,-chromosome_windows_size))
	bar_size = 1.1
	x_ticks = [i for i in range(0, max_chr_size, int(max_chr_size / 12))]
	pos_ticks2_basis = [0,plot_size/2,plot_size]
	y_ticks2 = []
	pos_y_ticks = []
	pos_y_ticks2 = []

	# Graph creation for each chromosome
	for num, chromosome in enumerate(sorted_nicely(chromosomes_size.keys())):

		bar_pos = chromosome_windows_low[num]+(chromosome_windows_size/2)-(bar_size/2)
		cyto = cytobands_infos[chromosome][0]
		color_cyto = cytobands_infos[chromosome][1]

		##################################################
		# Creation of the chromosomes with their cytobands

		ax.hlines(y=chromosome_windows_low[num], xmin=0 - (0.01 * max_chr_size), xmax=max_chr_size + (0.01 * max_chr_size), color='k', linestyle='--', lw=0.3)
		ax.hlines(y=bar_pos+(bar_size/2), xmin=0 - (0.01 * max_chr_size), xmax=0, color='k', linestyle='-', lw=0.5)
		ax.hlines(y=bar_pos+(bar_size/2), xmin=chromosomes_size[chromosome], xmax=max_chr_size + (0.01 * max_chr_size), color='k', linestyle='-', lw=0.5)
		ax.broken_barh([(0, int(cyto[-1][0] + cyto[-1][1]))], (bar_pos, bar_size), linewidth=2, edgecolor=['black'])
		ax.broken_barh(cyto, (bar_pos, bar_size), facecolors=color_cyto)
		pos_y_ticks.append(bar_pos+(0.5*bar_size))

		#############################
		# Curves creation

		max_y = counts['maximum'][chromosome]
		min_y = counts['minimum'][chromosome]

		# B1
		y = [y+bar_pos+bar_size+0.5 for y in curves_infos_b1[chromosome][1]]
		x = curves_infos_b1[chromosome][0]
		plt.plot(x, y, color=COLOR_B1, linewidth=SIZE_B1, label="Variants number (First time point)")
		if max_y != 0:
			for i in pos_ticks2_basis[::-1]:
				pos_y_ticks2.append(bar_pos+bar_size+0.5+i)
			for i in cross_product(plot_size, max_y, pos_ticks2_basis[::-1]):
				y_ticks2.append(round(i,1))
		else:
			pos_y_ticks2.append(bar_pos + bar_size + 0.5 + 0)
			y_ticks2.append(0)

		# B2
		y = [y + bar_pos+ bar_size+0.5 for y in curves_infos_b2[chromosome][1]]
		x = curves_infos_b2[chromosome][0]
		plt.plot(x, y, color=COLOR_B2, linewidth=SIZE_B2, label="Variants number (Second time point)")

		# SAMPLE
		y = [bar_pos-y-0.5 for y in curves_infos_sample[chromosome][1]]
		x = curves_infos_sample[chromosome][0]
		plt.plot(x, y, color=COLOR_SAMPLE, linewidth=SIZE_SAMPLE, label="Samples number")
		if min_y != 0:
			for i in pos_ticks2_basis[::-1]:
				pos_y_ticks2.append(bar_pos-0.5-i)
			for i in cross_product(plot_size, min_y, pos_ticks2_basis[::-1]):
				if i != 0:
					y_ticks2.append(round(-i,1))
				else:
					y_ticks2.append(0)
		else:
			pos_y_ticks2.append(bar_pos - 0.5 + 0)
			y_ticks2.append(0)

	# Create the legend
	
	handles, labels = plt.gca().get_legend_handles_labels()
	by_label = dict(zip(labels, handles))
	plt.legend(by_label.values(), by_label.keys(), title='Union number for all sample', loc='upper center', bbox_to_anchor=(0.5, -0.02))

	# Put label and minor and major ticks

	ax.yaxis.set_major_locator(FixedLocator(pos_y_ticks))
	ax.yaxis.set_major_formatter(FixedFormatter(list(sorted_nicely(chromosomes_size.keys()))))

	ax.set_xlabel('Chromosome position')
	ax.set_xticks(x_ticks)
	ax.set_xticklabels(x_ticks)
	ax.yaxis.set_minor_locator(FixedLocator(pos_y_ticks2))
	ax.yaxis.set_minor_formatter(FixedFormatter(y_ticks2))

	# Set the min/max limits of x and f axis
	ax.set_ylim(0+chromosome_windows_size, chromosome_windows_size*24+chromosome_windows_size)
	ax.set_xlim(0 - (0.01 * max_chr_size), max_chr_size + (0.01 * max_chr_size))

	# Set the scientific notation for x axis
	scientific_formatter = FuncFormatter(scientific)
	ax.xaxis.set_major_formatter(scientific_formatter)

	plt.title(title)

	plt.savefig(plot_name, format='svg', bbox_inches='tight')
	
	if not str(plot_name).endswith('.png'):
		plot_name = Path(plot_name).with_suffix('.png')
	plt.savefig(plot_name, format='png', bbox_inches='tight')

	plt.close()

	if logger:
                logger.info('- Done -')
	print('Done.')


def create_chromosomes_plot(gene_positions_b1 : dict, gene_positions_b2 : dict, gene_positions_sample : dict, cytoband_file : str, plot_name : str, step_on : int = 500000, logger = None) :

	span_on = step_on*2

	print(f'Load cytoband infos from hg38_cytoband.tsv')
	cytoband = pd.read_csv(cytoband_file, sep='\t')
	cytoband = cytoband.dropna()
	pattern = re.compile(r'chr(\d+)')
	chromosomes = sorted(list(set(cytoband['#chrom'])), key=lambda x: extract_num(x, pattern, float('inf')))

	# Initialisation of needed variables

	chromosomes_size = {}
	counts = {'maximum':{}, 'minimum':{}}		# dictionnary containing maximum/minimum number of gene counts for chromosomes positions (used to get the size of the plot)
	cytobands_infos = {}

	curves_infos_b1 = {}
	curves_infos_b2 = {}
	curves_infos_sample = {}

	print(f'Get informations from cytoband and genes positions dictionnary for each chromosome to create the plot...')
	if logger:
		logger.info(f'- Get informations from cytoband and genes positions dictionnary for each chromosome to create the plot -')		
	for num, chromosome in enumerate(chromosomes):

		size_cyto = list(cytoband.loc[cytoband['#chrom'] == chromosome]['chromEnd']-cytoband.loc[cytoband['#chrom'] == chromosome]['chromStart'])
		start_cyto = list(cytoband.loc[cytoband['#chrom'] == chromosome]['chromStart'])
		cyto = list(zip(start_cyto, size_cyto))
		color_cyto = [COLOR_G_STAINING[giestain] for giestain in list(cytoband.loc[cytoband['#chrom'] == chromosome]['gieStain'])]

		cytobands_infos[chromosome] = (cyto,color_cyto)
		chr_size = int(cyto[-1][0] + cyto[-1][1])
		chromosomes_size[chromosome] = chr_size

		x1, y1 = get_plot(gene_positions_b1, chromosome, chromosomes_size, span_on, step_on)
		x2, y2 = get_plot(gene_positions_b2, chromosome, chromosomes_size, span_on, step_on)
		x3, y3 = get_plot(gene_positions_sample, chromosome, chromosomes_size, span_on, step_on)

		counts['maximum'][chromosome] = max(y1+y2)
		counts['minimum'][chromosome] = -max(y3)

		plot_size = 5
		y1 = cross_product(max(y1+y2), plot_size, y1)
		y2 = cross_product(max(y1 + y2), plot_size, y2)
		y3 = cross_product(max(y3), plot_size, y3)

		curves_infos_b1[chromosome] = (x1, y1)
		curves_infos_b2[chromosome] = (x2, y2)
		curves_infos_sample[chromosome] = (x3, y3)

	create_plot(cytobands_infos, curves_infos_b1, curves_infos_b2, curves_infos_sample, chromosomes_size, counts, plot_size, plot_name, '', logger)



###########
# TO TEST #
###########
'''
gene_positions = {'chr1' : [(125600,136700),(145600,150500),(20000000,20003000),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr17' :  [(43044295, 43170245), (43044295, 43170245)]}
gene_positions2 = {'chr1' : [(125600,136700),(20000000,20003000),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr19' :  [(43044295, 43170245), (43044295, 43170245)]}
gene_positions3 = {'chr1' : [(125600,136700),(20000000,20003000),(20000000,20003000),(20000000,20003000),(20000000,20003000),(20000000,20003000),(20000050,20003050),(20003020,20003025),(20003100,20003150)], 'chr2' : [(125600000,125600050)], 'chr19' :  [(43044295, 43170245), (43044295, 43170245)]}
cytoband_file = '../LOTUS_external_files/hg38_cytoband.tsv'
create_chromosomes_plot(gene_positions, gene_positions2, gene_positions3, cytoband_file, 'Imnotaplot.svg')
'''








