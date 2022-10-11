#!/usr/bin/env python3
# coding: utf-8

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOTUS : LOngitudinal TUmors Study
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
                         __ _      | | ___ | |_ _   _ ___
                        / _` |  _  | |/ _ \| __| | | / __|
                       | (_| | |_| | | (_) | |_| |_| \__ \
                        \__, |     |_|\___/ \__|\__,_|___/
                         __/ |
                        |___/ 
'''

ascii_help= r'''

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_summarise = r'''                                                         _
                                                        (_)         
                ___ _   _ _ __ ___  _ __ ___   __ _ _ __ _ ___  ___ 
               / __| | | | '_ ` _ \| '_ ` _ \ / _` | '__| / __|/ _ \
               \__ \ |_| | | | | | | | | | | | (_| | |  | \__ \  __/
               |___/\__,_|_| |_| |_|_| |_| |_|\__,_|_|  |_|___/\___|                                                         

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_filter= r'''                                 __ _ _ _            
                                / _(_) | |           
                               | |_ _| | |_ ___ _ __ 
                               |  _| | | __/ _ \ '__|
                               | | | | | ||  __/ |   
                               |_| |_|_|\__\___|_|

-----------------------------------------------------------------------------------

-----------------------------------------------------------------------------------
'''

ascii_compare= r'''                     ___ ___  _ __ ___  _ __   __ _ _ __ ___ 
                    / __/ _ \| '_ ` _ \| '_ \ / _` | '__/ _ \
                   | (_| (_) | | | | | | |_) | (_| | | |  __/
                    \___\___/|_| |_| |_| .__/ \__,_|_|  \___|
                                       | |                   
                                       |_|                  

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
    
    parser_filter = subparsers.add_parser('filter',help='Simple filter on the vcf file from Mutect2, Funcotator... using multiple informations to keep only trustworthy somatic variants.')
    
    required_filter = parser_filter.add_argument_group('required arguments')
    optional_filter = parser_filter.add_argument_group('optional arguments')

    required_filter.add_argument('--vcf', '-v', dest='vcf', metavar='IN_FILE', type=str,
    required=True, 
    help='Result vcf file from Mutect2, FilterMutectCalls or Funcotator output.'
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

    parser_filter.set_defaults(parser_filter=True, parser_summarise=False, parser_compare=False)

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

    parser_summarise.set_defaults(parser_filter=False, parser_summarise=True, parser_compare=False)

    ################
    #Compare parser#

    parser_compare = subparsers.add_parser('compare', help='Compare two vcf files longitudinally (vcf1 = first in time)')
    
    required_compare = parser_compare.add_argument_group('Required argument')
    optional_compare = parser_compare.add_argument_group('Optional argument')
    
    required_compare.add_argument('--config', '-c', dest='config', metavar='CONFIG_FILE',
    required=True,
    help='Configuration file containing path to vcf file (filtered.vcf and pass.vcf file from g-LOTUS filter) and tsv files for indel and snp from g-LOTUS summarise. Example available on github: https://github.com/gsiekaniec/g-LOTUS/blob/main/example_config.txt'
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

    optional_compare.add_argument('--statistics', '-s', dest='stats', metavar='OUTPUT_STATS',
    default='stats.txt',
    help='Output statistics file. Default = "stats.txt" wich give "{vcf1}_{vcf2}_stats.txt".'
    )

    parser_compare.set_defaults(parser_filter=False, parser_summarise=False, parser_compare=True)


    #End parser#######

    args = parser.parse_args()

    # Logger configuration    

    logger = logging.getLogger('g-LOTUS main')
    logger.setLevel(logging.DEBUG)
    
    fh = logging.FileHandler(args.log)
    fh.setLevel(logging.DEBUG)
    
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    fh.setFormatter(formatter)    

    logger.addHandler(fh)

    logger.info(f'---------------- g-LOTUS v{__version__} ----------------')

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
        logger.info(f'---------------- g-LOTUS closes ----------------')
    else :
        parser.print_help()
        print(ascii_help)
        logger.info(f'---------------- g-LOTUS help ----------------')
