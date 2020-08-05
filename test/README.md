# GREEN-VARAN tools test

If everything is configured properly, you should be able to run the following tests from the tool main folder and see the corresponding results in test/out

## GREEN-VARAN
### Basic annotation
Annotate with regulatory region information, and 3 prediction scores

```
python GREEN-VARAN.py -i test/VCF/test.snpEff.vcf.gz -o test/out/test_standard.vcf -b GRCh38 -m annotate -s ReMM -s NCBoost -s LinSight --threads 4
```

### Annotate variants of interest based on a gene list
Annotate with regulatory region information, and 3 prediction scores. Now marks as NC_VOI=1 regulatory vars potentially affecting interesting genes from a list. Only genes experimentally linked to regions are considered (-t controlled)

```
python GREEN-VARAN.py -i test/VCF/test.snpEff.vcf.gz -o test/out/test_VOI.vcf -b GRCh38 -m annotate -s ReMM -s NCBoost -s LinSight -g annotate -t controlled -l test/VCF/genes_list_example.txt --threads 4
```

### Prioritize variants
Run in prioritize mode. This mode uses regulatory regions information and 3 prediction scores to priotize variants with likely impact on regulation

```
python GREEN-VARAN.py -i test/VCF/test.snpEff.vcf.gz -o test/out/test_prioritize.vcf -b GRCh38 -m annotate -p --threads 4
```

## SV_annotation
Using default configuration and supporting files the input VCF is annotated with:
- population AF from gnomAD and 1000G
- overlapping genes
- overlapping CDS
- overlapping GREEN-DB regions with controlled genes and constraint values

```
python SV_annotation.py -i test/VCF/GRCh38.test.SV.vcf.gz -o test/out/SV_test.vcf -c SV_annotation.json -b GRCh38
```

## GREEN-DB query
### Region IDs list input
We use a file containing a list of region IDs, one per line

```
GREEN-DB_query.py -r test/query/regions_IDs.list -o test/out/query_regionIDs -b GRCh38
```

### Variants table input 
We use a table containing some example variants regulating SOX10 gene.

```
GREEN-DB_query.py -t test/query/SOX10_example.tsv -o test/out/query_SOX10 -b GRCh37
```