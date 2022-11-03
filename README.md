
<p align="center">
  <img width="500" height="500" src="img/g-LOTUS.gif">
</p>

----

# Python packages

:file_folder: Packages from the [Python Standard Library](https://docs.python.org/3/library/) used:

  - [collections](https://docs.python.org/3/library/collections.html)
  - [copy](https://docs.python.org/3/library/copy.html)
  - [itertools](https://docs.python.org/3/library/itertools.html)
  - [logging](https://docs.python.org/3/library/logging.html)
  - [os](https://docs.python.org/3/library/os.html)
  - [pathlib](https://docs.python.org/3/library/pathlib.html)
  - [pickle](https://docs.python.org/3/library/pickle.html)
  - [sys](https://docs.python.org/3/library/sys.html)
  - [uuid](https://docs.python.org/3/library/uuid.html)
  - [warnings](https://docs.python.org/3/library/warnings.html)
  
:file_folder: Required python packages to run g-LOTUS:
  
  - [matplotlib](https://matplotlib.org/)
  - [more_itertool](https://more-itertools.readthedocs.io/en/stable/)
  - [numpy](https://numpy.org/)
  - [pandas](https://pandas.pydata.org/)
  - [pyfastx](https://pyfastx.readthedocs.io/en/latest/)
  - [pytest](https://docs.pytest.org/en/7.2.x/)
  - [tqdm](https://tqdm.github.io/)
  
  
----

# g-LOTUS informations

g-LOTUS is composed of the following four modules to process vcf files in GATK output (annotated with Funcotator):

<p align="center">
  <a href="https://github.com/gsiekaniec/g-LOTUS/blob/main/README.md#filter">
    <img width="200" height="50" src="img/filter.png">
  </a>
  <a href="https://github.com/gsiekaniec/g-LOTUS/blob/main/README.md#summarise">
    <img width="200" height="50" src="img/summarise.png">
  </a>
  <a href="https://github.com/gsiekaniec/g-LOTUS/blob/main/README.md#compare">
    <img width="200" height="50" src="img/compare.png">
  </a>
  <a href="https://github.com/gsiekaniec/g-LOTUS/blob/main/README.md#merge">
    <img width="200" height="50" src="img/merge.png">
  </a>
</p>

----

## Preliminary steps





----

## Filter

Simple filter on the vcf file from Funcotator using multiple informations to keep only trustworthy somatic variants.

<details><summary>Parameters</summary>

| Parameters | Description | Default |
|----------|:-------------:|:-------------:|
| --vcf, -v | Result vcf file from Funcotator output. |  |
| --output, -o | Filtered vcf file. | output.filtered.vcf |
| --working-method, -w | "InMemory" (default) loads the vcf file in memory into a list (more speed but higher memory consumption) or "Direct" reads and modifies the vcf file on the fly (slow speed but low memory consumption). | InMemory |
| --MBQ | Minimum median base variant quality for variant. | 20 |
| --DP | Minimum variant coverage. | 10 |
| --AF | Minimum fractions of variant in the tumor. | 0.1 |
| --AD | Minimum variant depths. | 5 |
| --POPAF | Maximum population (often GnomAD) variant frequencies. | 0.00001 |

</details>

----

## Summarise

Allows to extract a lot of statistics from a vcf file.

<details><summary>Parameters</summary>

| Parameters | Description | Default |
|----------|:-------------:|:-------------:|
| --vcf, -v | Vcf file containing variants that pass filter (*.filtered.pass.vcf). | None |
| --vcf_pass, -vp | Vcf file containing variants that pass filter (*.filtered.pass.vcf). |  |
| --genome, -g | Genome fasta file (allowed extensions : .fasta, .fa, .fan) or pickle (.pk, .pickle) file created after a first run. |  |
| --statistics, -s | Output statistics file. | stats.txt |
| --genes, -genes | Output file containing genes impacted by variants. | genes.txt |
| --profile, p | SVG file that shows the mutations profile of the vcf file. | profil.svg |
| --indel, -i | SVG file that shows the indel mutations size of the vcf file. | indel.svg |
| --enrichment | Did the GO enrichment analysis on the genes list using ToppGene and Panther and returns the biological processes (works if the APIs are not down). | False |

</details>

----

## Compare

Compare multiple vcf files longitudinally.

<details><summary>Parameters</summary>

| Parameters | Description | Default |
|----------|:-------------:|:-------------:|
| --config, -c | Configuration file containing path to vcf file (filtered.vcf and pass.vcf file from g-LOTUS filter) and tsv files for indel and snp from g-LOTUS summarise. Example available [here](https://github.com/gsiekaniec/g-LOTUS/blob/main/example_config.txt). |  |
| --output, -o | Excel file containing the genes specific to the first or second biopsy. | "genes.xlsx" wich give "{vcf1}_{vcf2}_genes.tsv/.xlsx". |
| --profile, -p | SVG file that shows the comparison between mutations profiles of the two vcf file. | "profile.svg" wich give "{vcf1}_{vcf2}_profile.svg". |
| --indel, -i | SVG file that shows the indel mutations size of the vcf file. | "indel.svg" wich give "{vcf1}_{vcf2}_indel.svg". |
| --enrichment | Did the GO enrichment analysis on the genes list using ToppGene and Panther and returns the biological processes (works if the APIs are not down). | False |

</details>

----

## Merge

Merging results to find the genes impacted in all patient.

<details><summary>Parameters</summary>

| Parameters | Description | Default |
|----------|:-------------:|:-------------:|
| --config, -c | Configuration file containing genes list from all patients. Merged patients results. |  |
| --output, -o | Output file name. | union.xlsx |
| --upset, -u | Output name for upset plot. | upset_plot.svg |
| --weakness_threshold, -wt | Mean weakness threshold to save take a gene into account. | 100 |
| --min_subset_size, -minsb | Minimum size of a subset (nb of genes by subset) to be shown in the UpSetPlot. All subsets with a size smaller than this threshold will be omitted from plotting. | 1 |
| --max_subset_size, -maxsb | Maximum size of a subset (nb of genes by subset) to be shown in the UpSetPlot. All subsets with a size greater than this threshold will be omitted from plotting. | 0 |
| --min_degree, -mind | Minimum degree of a subset (nb of patients) to be shown in the UpSetPlot. | 1 |
| --max_degree, -maxd | Maximum degree of a subset (nb of patients) to be shown in the UpSetPlot. | 0 |
  
</details>


----


# Tests

To run tests:

``` 
python -m py.test tests
```
