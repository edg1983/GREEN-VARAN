# GREEN-VARAN tools test

If everything is configured properly and you have annotations files downloaded in the resources folder, you should be able to run the following tests from the tool main folder and see the corresponding results in test/out

## GREEN-VARAN
### Basic annotation
Annotate regulatory regions information and update snpEff ANN field. Also add greenvaran_VOI for some genes of intereset. You will say a warning about prioritization, that's fine

```
greenvaran smallvars \
    -i test/VCF/GRCh38.test.smallvars.vcf.gz \
    -o test/out/test_smallvars.vcf \
    --db resources/GRCh38/GRCh38_GREEN-DB.bed.gz \
    --dbschema config/greendb_schema_v2.5.json \
    --config config/prioritize_smallvars.json \
    --genes test/VCF/genes_list_example.txt 
```

### Use the workflow to annotate and prioritize variants
Using the Nextflow workflow it is possible to automatically add all needed annotations and then prioritize regulatory variatns using greenvaran. Note that this will run on local machine using 10 threads for annotations.

```
nextflow workflow/main.nf -prfile local \
    --input test/VCF/GRCh38.test.smallvars.vcf.gz \
    --out test/out/workflow \
    --build GRCh38 \
    --scores best \
    --regions best \
    --AF \
    --greenvaran_config config/prioritize_smallvars.json \
    --greenvaran_dbschema config/greendb_schema_v2.5.json
```

## GREEN-DB query
### Region IDs list input
We use a file containing a list of region IDs, one per line

```
greendb_query -r test/query/regions_IDs.list -o test/out/query_regionIDs -b GRCh38
```

### Variants table input 
We use a table containing some example variants regulating SOX10 gene.

```
greendb_query -t test/query/SOX10_example.tsv -o test/out/query_SOX10 -b GRCh37
```