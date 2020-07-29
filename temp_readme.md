### Prioritization and gene mode arguments
| Argument           | Choices      | Description |
|:---------          |:-----------: |:------------|
|-p, --prioritize | | Turn on prioritization for non-coding variants.<br>VC_VARCLASS will be added to INFO |
|-g, --gene_mode | annotate, filter | Activate gene based annotation/filter |
|-t, --gene_type | controlled,closest,both | Which genes to consider for NC regions when running in gene mode | 
|-l, --genes_list | | List of genes of interest<br>Can be a comma-separated list or a file with one gene per line |
|--impact | HIGH, MODERATE,<br>LOW,MODIFIER | Only report NC vars if one of the controlled genes is already affected by a coding variants with at least this impact (require SnpEFF ANN field in the input VCF) |

### Additional arguments
| Argument           | Description |
|:---------          |:------------|
|--vcfanno | Full path to vcfanno executable |
|--bed_dir | Directory containing RegDB bed files |
|--AF_file | gnomAD VCF file |
|--scores_dir | Directory containing prediction scores tables 
|--separate_fields | Make multiple INFO fields instead of a single NC_ANNO annotation |
| --logfile | Filename for the log file |
| --threads | Number of threads for vcfanno annotation |
| -w, --overwrite | Overwrite output file if already exists |