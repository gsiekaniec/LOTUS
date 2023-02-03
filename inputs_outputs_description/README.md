# Description of the inputs/outputs needed to understand LOTUS

## Filter

#### Intputs

- **Raw vcf file**

Raw vcf file from Funcotator

- **Filters**

median base variant quality (MBQ) ≥ 20; variant coverage (DP) ≥ 10; fractions of variant in the tumour (AF) ≥ 0.1; variant depths (AD) ≥ 5; population variant frequencies, often from gnomAD (10-POPAF) < 1e-5; at least one paired read to support a variant (MFRL) > 0; a variant must be in a coding region (RNA or protein) and not silent (Funcotation information).

The basic filters of LOTUS are the following:
  - Median base variant quality (MBQ) ≥ 20.
  - Variant coverage (DP) ≥ 10.
  - Fractions of variant in the tumour (AF) ≥ 0.1.
  - Variant depths (AD) ≥ 5.
  - Population variant frequencies (10<sup>-POPAF</sup>) < 1e-5.
  - At least one paired read to support a variant (MFRL) > 0. It is possible to skip this filter by using the ```--unpaired``` option which indicates that the original reads are not paired.
  - A variant must be in a coding region (RNA or protein) and non-silent (Funcotator information). The list of Funcotator annotations allowing to keep a variant is the following: X, Y, Z, AA.

These filters can be modified in order to filter the varaint in a more or less stringent way.

### Outputs

- **Filtered vcf file** 

The ```filtered.vcf``` file is a vcf file that contains the variants of the original vcf file with an annotation indicating the filters not passed in the *Info field*.

- **Passed vcf file** 

The ```passed.vcf``` file is a vcf file that contains only the variants passing the filters. 


## Summarise

#### Intputs


### Outputs


## Compare

#### Intputs


### Outputs


## Merge

#### Intputs

- **Configuration file**

The only file needed for the *merge* module is the configuration file (```.txt```).

This configuration file contains the list of genes (```MutatedGenes.tsv``` or ```MutatedGenes.xlsx``` from the *compare* step) for all samples, one file per line (either xlsx or tsv). For example:

``` 
sample1_TP1_TP2.MutatedGenes.xlsx
sample2_TP1_TP2.MutatedGenes.tsv
sample3_TP1_TP2.MutatedGenes.tsv
```

- **Cytoband file**

To create the graph showing the impacted genes on the cytoband maps of the chromosomes (```.tsv```), the file containing the position of these cytobands must be provided. 

For the hg38 version of the human genome, this file can be found [here](https://genome.ucsc.edu/cgi-bin/hgTables) or provided with [LOTUS](https://github.com/gsiekaniec/LOTUS/blob/main/LOTUS_external_files/hg38_cytoband.tsv).


### Outputs

- **union.MutatedGenes.tsv|.xlsx file** 

This file contains the list of common impacted genes for all samples given in the configuration file. For each impacted gene a large amount of information can be found such as :

XXX


## Gene Ontology Enrichment Analysis

In addition to the outputs presented above, the three modules *summarise*, *compare* and *merge* offer the possibility to request a GOEA based on biological process terms from the output gene lists using the ```--enrichment``` option. 

These GOEAs are performed using the ToppGene[^1] and PANTHER[^2] APIs and output in files suffixed with ```Panther_enrichment.tsv|.xlsx``` and ```ToppGene_enrichment.tsv|.xlsx```.

The information contained in the columns of these files is as follows:

XXX


[^1]: [Chen,J. et al. (2009) ToppGene Suite for gene list enrichment analysis and candidate gene prioritization. Nucleic Acids Res., 37, W305-11.)](https://academic.oup.com/nar/article/37/suppl_2/W305/1149611?login=true)
[^2]: [Mi,H. et al. (2019) Protocol Update for large-scale genome and gene function analysis with the PANTHER classification system (v.14.0). Nat. Protoc., 14, 703-721.)](https://www.nature.com/articles/s41596-019-0128-8)


