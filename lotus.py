#!/usr/bin/env python3
# coding: utf-8

# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
#   LOTUS : LOngitudinal TUmors Study
#   Authors: G. Siekaniec
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #

import os 
os.environ['MPLCONFIGDIR'] = os.getcwd() + "/matplot_configs/"
import argparse 
import matplotlib
matplotlib.use('agg')
import python_scripts.filter
import python_scripts.summary
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
 .    ..... ..         .           .         .          ..         .  ....    ..
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

ascii_summary= r'''                  ___ _   _ _ __ ___  _ __ ___   __ _ _ __ _   _ 
                 / __| | | | '_ ` _ \| '_ ` _ \ / _` | '__| | | |
                 \__ \ |_| | | | | | | | | | | | (_| | |  | |_| |
                 |___/\__,_|_| |_| |_|_| |_| |_|\__,_|_|   \__, |
                                                            __/ |
                                                           |___/ 

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


if __name__ == '__main__':

    print(lotus_ascii)
    
    #Global parser####
    
    parser = argparse.ArgumentParser()
    parser._positionals.title = 'Subcommands'
    parser._optionals.title = 'Global arguments'
    parser.add_argument('--version', action='version',
        version=f'LOTUS v{__version__}',
        help='Display LOTUS version'
    )
    
    #Subparser
    
    subparsers = parser.add_subparsers(help='Functions')
    
    #Filter parser
    
    parser_filter = subparsers.add_parser('filter',help='Simple filter on the vcf file from Mutect2, Funcotator... using multiple informations to keep only trustworthy somatic variants.')
    
    required_filter = parser_filter.add_argument_group('required arguments')
    optional_filter = parser_filter.add_argument_group('optional arguments')

    required_filter.add_argument('--vcf', '-v', dest='vcf', metavar='IN_FILE', type=str,
    required=True, 
    help='Result vcf file from Mutect2, FilterMutectCalls or Funcotator output.'
    )

    optional_filter.add_argument('--output', '-o', dest='out', metavar='OUT_FILE' , 
    default='filtered.vcf',
    help='Filtered vcf file.')

    parser_filter.set_defaults(parser_filter=True, parser_summary=False, parser_compare=False, parser_merge=False)

    #Summary parser

    parser_summary = subparsers.add_parser('summary', help='Allows to extract a lot of statistics from a vcf file.')
    
    required_summary = parser_summary.add_argument_group('Required argument')
    optional_summary = parser_summary.add_argument_group('Optional argument')
    
    required_summary.add_argument('--vcf', '-v', dest='vcf', metavar='IN_FILE', 
    required=True, 
    help='Vcf file from Mutect2, FilterMutectCalls or Funcotator output.'
    )

    required_summary.add_argument('--genome', '-g', dest='genome', metavar='GENOME_FILE', type=str,
    required=True,
    help='Genome fasta file (allowed extensions : .fasta, .fa, .fan) or pickle (.pk, .pickle) file created after the first run.'
    )

    optional_summary.add_argument('--statistics', '-s', dest='stats', metavar='OUTPUT_STATS', 
    default='stats.txt', 
    help='Output statistics file.'
    )

    optional_summary.add_argument('--genes', '-genes', dest='genes', metavar='OUTPUT_GENES',
    default='genes.txt',
    help='Output file containing genes impacted by variants.'
    )

    optional_summary.add_argument('--profil', '-p', dest='profil',  metavar='OUTPUT_PROFIL', 
    default='profil.svg', 
    help='SVG file that shows the mutations profile of the vcf file.'
    )
 
    optional_summary.add_argument('--pass-only', '-pass', dest='pass_only',
    action='store_true', 
    help='Vcf contain only variant that pass filters. If False, also count the number of variants that do not pass the germline, PON and germline+PON filters and the total number of variants.  If true, consider all variants of the vcf file as passing the filters.'
    )

    parser_summary.set_defaults(parser_filter=False, parser_summary=True, parser_compare=False, parser_merge=False)

    #Compare parser

    parser_compare = subparsers.add_parser('compare', help='Compare two vcf files longitudinally (vcf1 = first in time)')
    
    required_compare = parser_compare.add_argument_group('Required argument')
    optional_compare = parser_compare.add_argument_group('Optional argument')
    
    required_compare.add_argument('--vcf1', '-v1', dest='vcf1', 
    required=True, 
    help='Vcf file from Mutect2, FilterMutectCalls or Funcotator output corresponding to the first biopsy.'
    )

    required_compare.add_argument('--vcf2', '-v2', dest='vcf2',
    required=True,
    help='Vcf file from Mutect2, FilterMutectCalls or Funcotator output corresponding to the second biopsy.'
    )

    optional_compare.add_argument('--functional', '-f', dest='functional', metavar='FUNCTIONAL',
    default=True, 
    help='Consider only mutations that may have a functional impact.'
    )

    optional_compare.add_argument('--output', '-o', dest='out', metavar='OUTPUT', 
    default='stats_comparison.xlsx', 
    help='Excel file containing the genes specific to the first or second biopsy.'
    )
    
    optional_compare.add_argument('--profil', '-p', dest='profil', metavar='PROFIL',
    default='profil_comparison.svg',
    help='SVG file that shows the comparison between mutations profiles of the two vcf file.'
    )    

    parser_compare.set_defaults(parser_filter=False, parser_summary=False, parser_compare=True, parser_merge=False)

    #Merge parser

    parser_merge = subparsers.add_parser('merge', help='Merge ')
    
    required_merge = parser_merge.add_argument_group('Required argument')
    optional_merge = parser_merge.add_argument_group('Optional argument')
    
    required_merge.add_argument('--repertory', '-r', dest='repertory',  metavar='COMPARISON_REPERTORY', 
    required=True,          
    help='Repertory containing multiple stats from comparison between biopsy. Merged patients results to search for biological process affected by the mutated genes.') 
    
    optional_merge.add_argument('--output_repertory', '-out', dest='out', 
    default='.', 
    help='Output repertory.'
    )
    
    parser_compare.set_defaults(parser_filter=False, parser_summary=False, parser_compare=False, parser_merge=True)



    #End parser#######

    args = parser.parse_args()
    
    if args != argparse.Namespace(): 
        if args.parser_filter:
            print(ascii_filter)
            python_scripts.filter.main(args)
        elif args.parser_summary:
            print(ascii_summary)
            python_scripts.summary.main(args)
        elif args.parser_compare:
            python_scripts.compare.main(args)
        elif args.parser_merge:
            python_scripts.merge.main(args)
    else :
        parser.print_help()
