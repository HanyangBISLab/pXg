# pXg: proteomics X genomics
<p align="center">
	<img src="https://github.com/HanyangBISLab/pXg/blob/main/img/logo.svg"/>
</p>
Thank you, <a href="https://www.bcm.edu/people-search/jong-min-choi-19517" target="_blank">Dr. Jongmin Choi</a>, for designing the logo!

---
- [About pXg](#about-pxg)
- [Usage](#usage)
  - [Input](#input)
  - [Output](#output)
      - [pXg result](#pxg-result)
      - [Unknown sequences](#unknown-sequences)
  - [Amino acid variant table format](#amino-acid-variant-table-format)
  - [Command-line interface](#command-line-interface)
    - [List of parameters](#list-of-parameters)
    - [Basic command](#basic-command)
- [Tutorial](#tutorial)
  - [RNA-Seq alignment](#rna-seq-alignment)
  - [SAM/BAM preparation](#sorted-sam/bam-preparation)
  - [Toy example](#toy-example)
  - [Run pXg](#run-pxg)
  - [Run Percolator using the pXg results](#run-percolator-using-the-pxg-results)
  - [IGV viewer](#igv-viewer)
- [TODO](#todo)
  - [GTF Export](#gtf-export)
- [Citation](#citation)
---

## About pXg

pXg (proteomics X genomics), a software tool that enables the reliable identification of both canonical and noncanonical MHC-I-associated peptides (MAPs) from de novo peptide sequencing by utilizing RNA-Seq data.
<br>

## Usage
pXg can be integrated with any search engines such as PEAKS and pNovo3.
It was developed for the reliable identification of noncanonical MAPs from de novo peptide sequencing; however, it can also be used to capture the number of reads mapped to each peptide sequence.
### Input
|Input    | Description    | Format    | Mandatory   |
| :---:   | :---:       | :---:     | :---:       |
| Searh result       | A list of PSMs identified from a search engine (e.g. <a href="https://www.bioinfor.com/peaks-studio/" target="_blank">PEAKS</a>, <a href="http://pfind.org/software/pNovo/index.html" target="_blank">pNovo3</a>, <a href="https://github.com/Noble-Lab/casanovo" target="_blank">Casanovo</a>)     | TSV or CSV | Yes   |
| Gene annotation    | It must be the same file used in the read alignment (e.g. <a href="https://www.gencodegenes.org/" target="_blank">Gencode</a>, <a href="https://ensemblgenomes.org/" target="_blank">Ensembl</a>)       | GTF        | Yes   |
| RNA-Seq reads      | Mapped and unmapped RNA-Seq reads. The file must be sorted by coordinates. Multiple SAM/BAM files should be separated by comma (,) | SAM/BAM        | Yes   |
| Protein sequences  | Canonical and contaminant protein sequences (e.g. UniProt)                        | Fasta      | No    |

*pXg is not applicable to the flat formatted output in pNovo3. A user must convert the flat format to CSV or TSV.<br>
*Since version 2.3.0, pXg can support multiple SAM/BAM files. "Reads" column indicates sum of reads from multiple SAM/BAM files. Reads in each SAM/BAM file is appended to the last columns. <br>

### Output
|Output    | Description    | Format   | Mandatory   |
| :---:   | :---:       | :---:     | :---:       |
| pXg result                | This is a main output file and contains a list of identification as TSV format         | TSV         | Yes   |
| pXg result for Percolator | This is a main output file and contains a list of identification as PIN format         | PIN         | Yes   |
| Matched reads*             | Matched reads to peptides passing all filters                            | SAM         | No    |

*Although the pXg result contains PSM information with corresponding RNA-Seq counts, it is not suitable for visualization. <br>
 Two output files (matched reads and peptides) are available for direct use in <a href="https://software.broadinstitute.org/software/igv/" target="_blank">IGV</a>, making visualization easier. <br>

#### pXg Result
|Field    | Description    | Value   |
| :---:   | :---:       | :---:     |
| SpecID | Identifier of a spectrum | String         |
| GenomicID | Identifier of genomic sequence | Integer         |
| Label | Target (1) and decoy (-1) labels | 1\|-1         |
| nDeltaScore | Normalized difference between main scores of current rank and top-rank peptides | Float       |
| Rank | Rank of candidate peptides | Integer       |
| GenomicLociCount | The number of genomic locations | Integer       |
| AminoAcidVariant | an amino acid substitution in a format [position]:[original amino acid]>[altered amino acid] | String      |
| InferredPeptide | Translated nucleotide sequence with a PTM annotation | String       |
| InferredSequence | Translated nucleotide sequence without a PTM annotation | String       |
| GenomicLoci | Genomic location of the peptide | String       |
| Strand | Strand of matched sequence | +\|-       |
| ObservedLeftFlankNucleotide | Nucleotide sequence of the left flank of the peptide | String       |
| ObservedNucleotide | Nucleotide sequence of the peptide | String       |
| ObservedRightFlankNucleotide | Nucleotide sequence of the right flank of the peptide | String       |
| ReferenceLeftFlankNucleotide | Reference nucleotide sequence of the left flank of the peptide | String       |
| ReferenceNucleotide | Reference nucleotide sequence of the peptide | String       |
| ReferenceRightFlankNucleotide | Reference nucleotide sequence of the right flank of the peptide | String       |
| Mutations | Genomic information of mutations in the peptide | String       |
| MutationStatus | Indication of alteration caused by the mutations | Altered\|Same    |
| TranscriptIDs | Matched transcript IDs | String    |
| GeneIDs | Matched gene IDs | String    |
| GeneIDCount | The number of matched gene IDs | Integer    |
| GeneNames | Matched gene names | String    |
| GeneNameCount | The number of matched gene names | Integer    |
| PercentFullDistance | Proportion of start genomic loci in the longest transcripts (exons + introns) | Float    |
| PercentExonDistance | Proportion of start genomic loci in the longest transcripts (exons) |  Float    |
| PercentCDSDistance | Proportion of start genomic loci in the longest transcripts (CDSs) |  Float    |
| FromCDSStartSite | Distance from the start site |  String    |
| FromCDSStopSite | Distance from the stop site |  String    |
| Events | Type of identified feature |  String    |
| EventCount | The number of events |  Integer    |
| FastaIDs | Matched identifiers in a given fasta sequences |  String   |
| FastaIDCount | The number of FastaIDs |  Integer   |
| Reads | Geometirc mean of matched reads (or RPHM) from all SAM/BAM files |  Float   |
| MeanQScore | Mean of Phred scores |  Float   |
| IsReference | Reference (true) or non-reference (false) status |  true\|false  |
| SAM/BAM file name | The number of matched reads (or RPHM) in each SAM/BAM file |  Float  |


#### Unknown sequences
 Unknown sequences include sequence information from "unknown" events. The header line begins with ">[PEPTIDE]". Following the header line is the matched read information, which includes the sequence identifier, genomic location (if available), full sequence, and matched sequence.

### Amino acid variant table format
|aaRNA    | aaPeptide|
| :---:   | :---:    |
| W | F |
| W | M |
Single amino acid variants (SAAVs) can arise after translation and therefore cannot always be detected at the RNA level.<br>
pXg allows for the consideration of SAAVs based on de novo peptide sequencing results.<br>
As an example, we illustrate two SAAVs: W→F and W→M.<br>


### Command-line interface
#### List of Parameters
|Option    | Description    | Value   | Mandatory   |
| :---:   | :---:       | :---:     | :---:     |
| gtf       | GTF file path. We recommand to use the same gtf corresponding to alignment | String |Yes   |
| bam       | SAM/BAM file path. The file must be sorted by coordinate. Multiple SAM/BAM files should be separated by comma (,) | String |Yes   |
| psm       | PSM file path. It is expected that the psm file is derived from proteomics search by de novo or database search engine |String |Yes   |
| identifier_index       | PSM identifier indicies (one-based). One or more indicies can be specified by comma separated. ex> 3,5,7 |Integer|Yes   |
| peptide_index       | Peptide index in the psm file |Integer|Yes   |
| charge_index     | Charge state index in the psm file |Integer|Yes   |
| output         | Base output name of pXg |String|Yes   |
| comment     | Specify the starting characters of comment lines to be ignored during processing. Lines beginning with these characters will be skipped. The default value is #\|@\|%\|MTD |String|No   |
| min_score            | Specify the minimum score threshold for peptide-spectrum matches (PSMs) to be included. The default value is 0 |float|No   |
| count            | Specify which reads to consider for counting. The default is primary |primary\|all |No   |
| sep            | Specify the column separator. Possible values are csv or tsv. Default is tsv |tsv\|csv|No   |
| mode           | Specify strandedness (default is auto). auto: auto-detection. only available in paired-ends. fr: first-forward second-reverse. rf: first-reverse second-forward. r: reverse single end. f: forward single end. none: non-strandedness |f\|r\|fr\|rf\|auto\|none|No   |
| add_index  | Specify the indices for additional features to generate PIN file. Several features can be added by comma separator. ex> 5,6,7|Integer|No  |
| il_equivalent           | Controls whether pXg treats isoleucine (I) and leucine (L) as the same/equivalent with respect to a peptide identification. Default is true |true\|false|No   |
| lengths        | Range of peptide length to consider. Default is 8-15. You can write in this way (min-max, both inclusive) : 8-13 |Integer|No   |
| fasta     | Canonical sequence database to report conservative assignment of noncanonical PSMs |String|No   |
| rank           | How many candidates will be considered per a scan. Default is 100 (in other words, use all ranked candidates) |Integer|No   |
| aa_variant           | File path of amino acid variant table |String|No   |
| output_sam        | Report matched reads as SAM format (true or false). Default is false |true\|false|No   |
| output_canonical  | Report caonical peptides in the out_sam file (true or false). Default is true |true\|false|No   |
| output_noncanonical| Report noncaonical peptides in the out_sam file (true or false). Default is true |true\|false|No   |
| penalty_mutation   | Penalty per a mutation. Default is 1 |Float|No   |
| penalty_alternative_splicing         | Penalty for alternative splicing. Default is 10 |Float|No   |
| penalty_5utr       | Penalty for 5`-UTR. Default is 20 |Float|No   |
| penalty_3utr       | Penalty for 3`-UTR. Default is 20 |Float|No   |
| penalty_ncrna      | Penalty for noncoding RNA. Default is 20 |Float|No   |
| penalty_frameshift         | Penalty for frame shift. Default is 20 |Float|No   |
| penalty_intron_retention         | Penalty for intron region. Default is 30 |Float|No   |
| penalty_intergenic_region        | Penalty for intergenic region. Default is 30 |Float|No   |
| penalty_asrna      | Penalty for antisense RNA. Default is 30 |Float|No   |
| penalty_softclip      | Penalty for softclip reads. Default is 50 |Float|No   |
| penalty_unknown    | Penalty for unmapped reads. Default is 100 |Float|No   |
| gtf_partition_size*       | The size of treating genomic region at once. Default is 5000000 |Integer|No   |
| sam_partition_size*       | The size of treating number of reads at once. Default is 1000000 |Integer|No   |
| threads*                  | The number of threads. Default is 4|Integer|No   |

*size parameters can effect memory usage and time. If your machine does not have enough memory, then decrease those values.

#### Basic command
```bash
java -Xmx30G -jar pXg.jar \
--gtf [gene annotation file path] \
--bam [sorted SAM/BAM file path] \
--psm [de novo result file path] \
--fasta [protein sequence fasta file paht] \
--identifier_index [index of file name column] \
--charge_index [index of chage state column] \
--peptide_index [index of peptide column] \
--score_index [index of search score column] \
--output [base output file name]
```

## Tutorial
This tutorial aims to understand how to run pXg and estimate FDR from the result. It contains 1) running <a href="https://github.com/alexdobin/STAR" target="_blank">STAR2</a> aligner with 2-pass parameter, 2) preparing SAM file from the alignment, 3) running pXg and 4) several post-processing including <a href="https://github.com/percolator/percolator" target="_blank">Percolator</a>, merging pXg result with the result of Percolator and estimating separated FDR.
Note that it neither contains how to run de novo peptide sequencing engines such as <a href="https://www.bioinfor.com/peaks-studio/" target="_blank">PEAKS</a>, <a href="http://pfind.org/software/pNovo/index.html" target="_blank">pNovo3</a> and <a href="https://github.com/Noble-Lab/casanovo" target="_blank">Casanovo</a> AND how to create deep learning based features.

### RNA-Seq alignment
We recommand to align fastq files using STAR2 with The Cancer Genome Atlas (TCGA) <a href="https://docs.gdc.cancer.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline" target="_blank">two-pass alignment option</a>. 

### Sorted SAM/BAM preparation
Once you get the aligned BAM or SAM file, you MUST sort the file by chromosomal coordinates.

We provide a code for preprocessing SAM file using <a href="http://www.htslib.org/" target="_blank">SAMtools</a> below:
```bash
samtools sort -o in.sorted.bam in.bam -@ 8
samtools index in.sorted.bam -@ 8
```

The "in.sorted.bam" is used for pXg input.

### Toy example
In this tutorial, toy datasets including 1) de novo results, 2) in.sorted.sam, 3) gene annotation (GTF) and 4) protein sequence fasta file are provided in the <a href="https://github.com/progistar/pXg/tree/main/tutorial" target="_blank">tutorial</a> folder so that a user can try to run the pXg pipeline. 

### Run pXg
Using the toy datasets, you can run the pXg pipline using following command: <br>
```bash
java -Xmx2G -jar pXg.v2.0.1.jar \
--gtf toy.gtf \
--sam toy.sorted.sam \
--psm toy.psm.csv \
--fasta toy.fasta \
--output toy \
--identifier_index 2,5 \
--peptide_index 4 \
--score_index 8 \
--charge_index 11 \
--add_index 15 \
--sep csv \
--mode none \
--threads 2
```
This may take about 2 mins.<br>

Note that the memory option "-Xmx50G" depends on the size of BAM file. In our experience, "-Xmx30G" is enough to deal with ~20G file. 

### Run Percolator using the pXg results
Once you get the pXg result, you can add more features such as spectral similarity and delta retention time described in our manuscript. Without the additional features, still it is possible to run <a href="https://github.com/percolator/percolator" target="_blank">Percolator</a> and estimate FDR from the pXg results.<br>
We recommand to use Percolator version >= v3.06.1 because former versions have an issue to print proteinIds.<br>
Post processing codes are also provided in the tutorial folder (post_process.ipynb).

### IGV viewer
<img src="https://github.com/progistar/pXg/blob/main/img/toy.viewer.png"/>
When pXg finishes identifying peptides, the resulting GTF and SAM files are immediately available in the <a href="https://software.broadinstitute.org/software/igv/" target="_blank">IGV viewer</a>.

## TODO
### GTF Export
* Export GTF format from pXg result.

## Citation
<a href="https://www.mcponline.org/article/S1535-9476(24)00033-1/" target="_blank">pXg: Comprehensive Identification of Noncanonical MHC-I–Associated Peptides From De Novo Peptide Sequencing Using RNA-Seq Reads. Seunghyuk Choi and Eunok Paek, Molecular & Cellular Proteomics 2024.</a><br>


