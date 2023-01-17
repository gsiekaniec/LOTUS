#!/usr/bin/env python3
# coding: utf-8

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOTUS : LOngiTUdinal comparative genomic Study
#   Authors: G. Siekaniec, W. Gouraud
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os 
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/matplot_configs/"
import argparse 
import matplotlib
import logging
matplotlib.use('agg')
import python_scripts.filter
import python_scripts.summarise
import python_scripts.compare
import python_scripts.merge


__version__ = '0.0.1'

lotus_ascii = r'''
-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------

                                        .
                                      .. ..
                                     .     .
                                    .       .
          ...           . ..       .         .       ....           ...
          .    ..      ..   ..     .          .    .    .      ...   ..
          .       ..   .      ..  .           .. .      .   ..       .
           .         ...        . .           ...       ....         .
           ..          .         ..           .         ..          .
 .    ..... ..         .           .         .          ..         .  ....   ..
  .          ...       .            .       .           .        ..          .
   ..             ..   ..            .     .            .    ..             .
     .               .. .            ..   .            ....               ..
      .                 ..            .  .            ..                 ..
       ..                  ..         .. .          ..                  .
         .                   ..        ...       ..                   ..
           .                    .      ..      ..                   ..
             ..                   .    ..    ..                  ..
                ...                 .  ..  ..                ...
                      ....            . ...           ....
                            ......   .......   .....
                            _       _
                           | |     | |
                           | | ___ | |_ _   _ ___
                           | |/ _ \| __| | | / __|
                           | | (_) | |_| |_| \__ \
                           |_|\___/ \__|\__,_|___/
                         
                         
'''

ascii_help= r'''

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_summarise = r'''                                                      _
                                                     (_)         
             ___ _   _ _ __ ___  _ __ ___   __ _ _ __ _ ___  ___ 
            / __| | | | '_ ` _ \| '_ ` _ \ / _` | '__| / __|/ _ \
            \__ \ |_| | | | | | | | | | | | (_| | |  | \__ \  __/
            |___/\__,_|_| |_| |_|_| |_| |_|\__,_|_|  |_|___/\___|                                                         

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_filter= r'''                              __ _ _ _            
                             / _(_) | |           
                            | |_ _| | |_ ___ _ __ 
                            |  _| | | __/ _ \ '__|
                            | | | | | ||  __/ |   
                            |_| |_|_|\__\___|_|

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_compare= r'''                  ___ ___  _ __ ___  _ __   __ _ _ __ ___ 
                 / __/ _ \| '_ ` _ \| '_ \ / _` | '__/ _ \
                | (_| (_) | | | | | | |_) | (_| | | |  __/
                 \___\___/|_| |_| |_| .__/ \__,_|_|  \___|
                                    | |                   
                                    |_|                  

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_merge= r'''                        _ __ ___   ___ _ __ __ _  ___ 
                       | '_ ` _ \ / _ \ '__/ _` |/ _ \
                       | | | | | |  __/ | | (_| |  __/
                       |_| |_| |_|\___|_|  \__, |\___|
                                            __/ |     
                                           |___/

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''


if __name__ == '__main__':

    print(lotus_ascii)
    
    #Global parser####
    
    parser = argparse.ArgumentParser()
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'
    parser.add_argument('--version', action='version',
        version=f'LOTUS v{__version__}',
        help='Display LOTUS version')
    parser.add_argument('--log', '-l', dest='log', metavar='LOG_FILE',
    default='LOTUS.log',
    help='Log file name.')
    
    ###########
    #Subparser#
    
    subparsers = parser.add_subparsers(help='Functions')
    
    ###############
    #Filter parser#
    
    parser_filter = subparsers.add_parser('filter',help='Simple filter on the vcf file from Funcotator using multiple informations to keep only trustworthy somatic variants.')
    
    required_filter = parser_filter.add_argument_group('required arguments')
    optional_filter = parser_filter.add_argument_group('optional arguments')

    required_filter.add_argument('--vcf', '-v', dest='vcf', metavar='IN_FILE', type=str,
    required=True, 
    help='Result vcf file from Funcotator output.'
    )

    optional_filter.add_argument('--output', '-o', dest='out', metavar='OUT_FILE', 
    default='output.filtered.vcf',
    help='Filtered vcf file. Default = "output.filtered.vcf".')

    optional_filter.add_argument('--working-method', '-wm', dest='working_method', metavar='WORKING_METHOD' ,
    default='InMemory', choices=['InMemory', 'Direct'],
    help='"InMemory" (default) loads the vcf file in memory into a list (more speed but higher memory consumption) or "Direct" reads and modifies the vcf file on the fly (slow speed but low memory consumption).')

    optional_filter.add_argument('--MBQ', dest='mbq', metavar='MBQ',
    default=20, type=int,
    help='Min median base variant quality for variant. Default = 20.')
   
    optional_filter.add_argument('--DP', dest='dp', metavar='DP',
    default=10, type=int,
    help='Min variant coverage. Default = 10.')    

    optional_filter.add_argument('--AF', dest='af', metavar='AF',
    default=0.1, type=float,
    help='Min fractions of variant in the tumor. Default = 0,1.')

    optional_filter.add_argument('--AD', dest='ad', metavar='AD',
    default=5, type=int,
    help='Min variant depths. Default = 5.')

    optional_filter.add_argument('--POPAF', dest='popaf', metavar='POPAF',
    default=0.00001, type=float,
    help='Max population (often GnomAD) variant frequencies. Default = 0,00001.')

    optional_filter.add_argument('--unpaired', dest='paired',
    action='store_false',
    help='Argument to use if the reads used are unpaired (single end). Default = True.'
    )

    parser_filter.set_defaults(parser_filter=True, parser_summarise=False, parser_compare=False, parser_merge=False)

    ##################
    #Summarise parser#

    parser_summarise = subparsers.add_parser('summarise', help='Allows to extract a lot of statistics from a vcf file.')
    
    required_summarise = parser_summarise.add_argument_group('Required argument')
    optional_summarise = parser_summarise.add_argument_group('Optional argument')
    
    optional_summarise.add_argument('--vcf', '-v', dest='vcf', metavar='IN_FILE', 
    default=None, 
    help='Vcf file from filter output (.filtered.vcf).'
    )

    required_summarise.add_argument('--vcf_pass', '-vp', dest='vcf_pass', metavar='IN_FILE2',
    required=True,
    help='Vcf file containing variants that pass filter (.filtered.pass.vcf).'
    )

    required_summarise.add_argument('--genome', '-g', dest='genome', metavar='GENOME_FILE', type=str,
    required=True,
    help='Genome fasta file (allowed extensions : .fasta, .fa, .fan) or pickle (.pk, .pickle) file created after the first run.'
    )

    optional_summarise.add_argument('--statistics', '-s', dest='stats', metavar='OUTPUT_STATS', 
    default='stats.txt', 
    help='Output statistics file. Default = "stats.txt".'
    )

    optional_summarise.add_argument('--genes', '-genes', dest='genes', metavar='OUTPUT_GENES',
    default='genes.txt',
    help='Output file containing genes impacted by variants. Default = "genes.txt".'
    )

    optional_summarise.add_argument('--profile', '-p', dest='profile',  metavar='OUTPUT_PROFILE', 
    default='profile.svg',
    help='SVG file that shows the mutations profile of the vcf file. Default = "profil.svg".'
    )

    optional_summarise.add_argument('--indel', '-i', dest='indel',  metavar='OUTPUT_INDEL',
    default='indel.svg',
    help='SVG file that shows the indel mutations size of the vcf file. Default = "indel.svg".'
    )

    optional_summarise.add_argument('--enrichment', dest='enrichment',
    action='store_true',
    help='Did the GO enrichment analysis on the genes list using ToppGene and Panther and returns the biological processes (works if the APIs are not down). Default = False.'
    )

    parser_summarise.set_defaults(parser_filter=False, parser_summarise=True, parser_compare=False, parser_merge=False)

    ################
    #Compare parser#

    parser_compare = subparsers.add_parser('compare', help='Compare multiple vcf files longitudinally')
    
    required_compare = parser_compare.add_argument_group('Required argument')
    optional_compare = parser_compare.add_argument_group('Optional argument')
    
    required_compare.add_argument('--config', '-c', dest='config', metavar='CONFIG_FILE',
    required=True,
    help='Configuration file containing path to vcf file (filtered.vcf and pass.vcf file from LOTUS filter) and tsv files for indel and snp from LOTUS summarise. Example available on github: https://github.com/gsiekaniec/LOTUS/blob/main/example_config.txt'
    )

    required_compare.add_argument('--gff3', '-gff3', dest='gff3', metavar='GFF3_fILE',
    required=True,
    help='Gff3 file. This file can be found here : https://ftp.ensembl.org/pub/release-108/gff3/homo_sapiens/ or in the LOTUS github.'
    )

    optional_compare.add_argument('--output', '-o', dest='out', metavar='OUTPUT', 
    default='genes.xlsx', 
    help='Excel file containing the genes specific to the first or second biopsy. Default = "genes.xlsx" wich give "{vcf1}_{vcf2}_genes.tsv/.xlsx".'
    )
    
    optional_compare.add_argument('--profile', '-p', dest='profile', metavar='PROFILE',
    default='profile.svg',
    help='SVG file that shows the comparison between mutations profiles of the two vcf file. Default = "profile.svg" wich give "{vcf1}_{vcf2}_profile.svg".'
    )

    optional_compare.add_argument('--indel', '-i', dest='indel',  metavar='OUTPUT_INDEL',
    default='indel.svg',
    help='SVG file that shows the indel mutations size of the vcf file. Default = "indel.svg" wich give "{vcf1}_{vcf2}_indel.svg".'
    )

    optional_compare.add_argument('--enrichment', dest='enrichment',
    action='store_true',
    help='Did the GO enrichment analysis on the genes list using ToppGene and Panther and returns the biological processes (works if the APIs are not down). Default = False.'
    )

    optional_compare.add_argument('--pickle_gff3', dest='pk_gff3',
    action='store_true',
    help='Did the gff3 file given is a pickle file from previous lauch. Default = False.'
    )

    optional_compare.add_argument('--additional_gene_information', dest='agi',
    action='store_true',
    help='Add gene informations using the LOTUS file containing information from tumorspecific database (CancerHotSpot, CIViC, COSMIC, DoCM, IntOGen and TSGene 2.0). Default = False.'
    )

    parser_compare.set_defaults(parser_filter=False, parser_summarise=False, parser_compare=True, parser_merge=False)


    #Merge parser

    parser_merge = subparsers.add_parser('merge', help='Merging results to find the genes impacted in all patient.')

    required_merge = parser_merge.add_argument_group('Required argument')
    optional_merge = parser_merge.add_argument_group('Optional argument')

    required_merge.add_argument('--config', '-c', dest='conf',  metavar='CONFIGURATION_FILE',
    required=True,
    help='Configuration file containing genes list from all patients. Merged patients results.')

    optional_merge.add_argument('--output', '-o', dest='output', metavar='OUTPUT',
    default='union.xlsx',
    help='Output file name. Default = union.xlsx.'
    )

    optional_merge.add_argument('--cytoband', '-cyto', dest='cyto', metavar='CYTOBAND_FILE',
    default=None,
    help='Human cytoband file for the corresponding genome version. This file can be download here : https://genome.ucsc.edu/cgi-bin/hgTables or in the LOTUS github (for hg38). Default = 50000.'
    )

    optional_merge.add_argument('--chromosome-step', '-step', dest='step', metavar='STEP',
    default=500000, type=int,
    help='Frame used for counting the number of genes along the chromosomes. Default = 500000.'
    )

    optional_merge.add_argument('--chromosomes_output', '-co', dest='chromosomes_output', metavar='CHROMOSOMES_OUTPUT',
    default='chromosomes.svg',
    help='Output file name for the chromosomes plot. Default = chromosomes.svg'
    )

    optional_merge.add_argument('--upset', '-u', dest='upset', metavar='UPSET_OUTPUT',
    default=None,
    help='Output name for upset plot. Default = None.'
    )

    optional_merge.add_argument('--weakness_threshold', '-w', dest='threshold', metavar='WEAKNESS_THRESHOLD',
    default=100, type=int,
    help='Mean weakness threshold to save take a gene into account. Default = 100.'
    )
	
    optional_merge.add_argument('--min_subset_size', '-minsb', dest='min_subset_size', metavar='MIN_SUBSET_SIZE',
    default=1, type=int,
    help='Minimum size of a subset (nb of genes by subset) to be shown in the UpSetPlot. All subsets with a size smaller than this threshold will be omitted from plotting. Default = 1.'
    )

    optional_merge.add_argument('--max_subset_size', '-maxsb', dest='max_subset_size', metavar='MAX_SUBSET_SIZE',
    default=0, type=int,
    help='Maximum size of a subset (nb of genes by subset) to be shown in the UpSetPlot. All subsets with a size greater than this threshold will be omitted from plotting. Default = 0 (take the maximum possible value).'
    )

    optional_merge.add_argument('--min_degree', '-mind', dest='min_degree', metavar='MIN_DEGREE',
    default=1, type=int,
    help='Minimum degree of a subset (nb of patients) to be shown in the UpSetPlot. Default = 1.'
    )

    optional_merge.add_argument('--max_degree', '-maxd', dest='max_degree', metavar='MAX_DEGREE',
    default=0, type=int,
    help='Maximum degree of a subset (nb of patients) to be shown in the UpSetPlot. Default = 0 (take the maximum possible value).'
    )

    optional_merge.add_argument('--additional_gene_information', dest='agi',
    action='store_true',
    help='Add gene informations using the LOTUS file containing information from tumorspecific database (CancerHotSpot, CIViC, COSMIC, DoCM, IntOGen and TSGene 2.0). Default = False.'
    )

    optional_merge.add_argument('--enrichment', dest='enrichment',
    action='store_true',
    help='Did the GO enrichment analysis on the genes list using ToppGene and Panther and returns the biological processes (works if the APIs are not down). Default = False.'
    )

    parser_merge.set_defaults(parser_filter=False, parser_summarise=False, parser_compare=False, parser_merge=True)


    #End parser#######

    args = parser.parse_args()

    # Logger configuration    

    logger = logging.getLogger('LOTUS main')
    logger.setLevel(logging.DEBUG)
    
    fh = logging.FileHandler(args.log)
    fh.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)    

    logger.addHandler(fh)

    logger.info(f'---------------- LOTUS v{__version__} ----------------')

    if args != argparse.Namespace(log=args.log):
        if args.parser_filter:
            print(ascii_filter)
            python_scripts.filter.main(args)
        elif args.parser_summarise:
            print(ascii_summarise)
            python_scripts.summarise.main(args)
        elif args.parser_compare:
            print(ascii_compare)
            python_scripts.compare.main(args)
        elif args.parser_merge:
            print(ascii_merge)
            python_scripts.merge.main(args)
        logger.info(f'---------------- LOTUS closes ----------------')
    else :
        parser.print_help()
        print(ascii_help)
        logger.info(f'---------------- LOTUS help ----------------')
