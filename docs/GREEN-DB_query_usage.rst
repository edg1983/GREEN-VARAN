greendb_query tool usage
========================

greendb_query assists in quering the GREEN-DB database.
Given a list of region IDs, variant IDs or a table or variant and relevant regions, the tool generates a set of tables
containing detailed information on the regions of interest, overlap with additional supporting regions
(TFBS, DNase HS peaks, UCNE, dbSuper), gene-region connections, tissue of activity and associated phenotypes.

.. code-block:: bash

    greendb_query [-h] (-v VARIDS | -r REGIDS | -t TABLE) -o OUTPREFIX -g
                         {GRCh37,GRCh38} --db GREENDB [--logfile LOGFILE]

Possible inputs
~~~~~~~~~~~~~~~
The tools allows to query GREEN-DB using 3 different type of inputs.
Only one type of input can be specified.

1. List of regions (-r)
#######################
If you are simply interested in detailed information on a list of regions, you can use the -r input.
This argument accepts a comma-separated list of regions (like ``ID1,ID2``) or a text file with one region ID per line.

2. VCF file (-v)
################
If you have a small list of variants for which you want to extract overalpping regulatory regions, you can 
input a them as a comma-separated list of variant IDs (like ``var1,var2``) or a text file with one variant ID per line
A variant ID has the format ``chrom_pos_ref_alt``

3. Variant-regions table (-t)
#############################
If you have a list of variants of interest for which you know the relevant GREEN-DB region IDs
you can query the DB directly providing a tab separated text file with **no header** and 2 columns:

- column 1: variant ID in the format chrom_pos_ref_alt
- column 2: comma-separated list of region IDs overlapping the variant

This table can be generated automatically from a VCF annotated with greenvaran by using ``greenvaran querytab``


Output tables
~~~~~~~~~~~~~
The tool will generate 6 tables with the provided prefix. Some table may be empty if the corresponding information is missing. 
Output tables structure is described below

regions
#######
Details on the regions of interest

| **1. regionID**: GREEN-DB region ID
| **2-4. chrom, start, stop**: genomic location of the region
| **5. type**: region type as extracted from the source dataset
| **6. std_type**: one of the 5 main region types (enhancer, promoter, silencer, bivalent, insulator)
| **7. DB_source**: comma-separated list of sources supporting the region
| **8. PhyloP100_median**: median PhyloP100 conservation value across the region
| **9. constraint_pct**: constraint metric. range 0-1 with higher values equals more intolerant to variants
| **10. controlled_gene**: comma-separated list of gene symbols for controlled genes with experimental support
| **11-13. closestGene_symbol, _ensg, _dist**: symbol, ensembl IDs and distance for the closeset gene
| **14. cell_or_tissues**: comma-seprated list of cell types and tissues where the region is active
| **15. detection_method**: comma-separated list of methods supporting this regions
| **16. phenotype**: comma-separated list of phenotypes eventually associated to this region

gene_details
############
Details on the controlled genes, reporting the tissue where the gene-region interaction is detected

| **1. regionID**: GREEN-DB region ID
| **2-4. chrom, start, stop**: genomic location of the region
| **5. std_type**: one of the 5 main region types (enhancer, promoter, silencer, bivalent, insulator)
| **6. controlled_gene**: gene symbol for controlled gene
| **7. detection method**: method supporting this interaction
| **8. tissue_of_interaction**: comma-separated list of cell types and tissues where this region-gene interaction is detected
| **9. same_TAD**: 0/1 value indicating if the reported interaction occurs in the same TAD according to TADKB

pheno_details
#############
Details on the phenotypes potentially associated with the regions of interest

| **1. regionID**: GREEN-DB region ID
| **2-4. chrom, start, stop**: genomic location of the region
| **5. std_type**: one of the 5 main region types (enhancer, promoter, silencer, bivalent, insulator)
| **6. phenotype**: phenotype eventually associated to this region
| **7. detection method**: method supporting this association. Note that when the method is GENE2HPO this means that the phenotype is inferred from HPOs associated to the controlled gene(s)
| **8. DB source**: source supporting this association


DNase, dbSuper, TFBS, UCNE
##########################
For each of the 4 functional elements a table is generated with details on each element overlapping the region(s) / variant(s) of interest.

| **1. regionID**: GREEN-DB region ID
| **2-4. dataset_chrom, _start, _stop**: genomic location of the functional element
| **5. dataset_ID**: database ID of the functional element
| **6. dataset_cell_or_tissue**: comma-separated list of cell types and tissues where the element is detected
|   cell and tissue information is not available for UCNE

Variant(s) of interest
######################
When the input contains variants of interest (-t, -v), an additional column is added to all tables.
A region or element is reported in the output only if it overlaps with one of the variants.

| **var_id**: comma-separated list of variant ID(s) (chrom_pos_ref_alt) of the variant(s) overlapping this feature 

Arguments list
~~~~~~~~~~~~~~
-v VARID, --vcf VARID
    | Comma separated list of variant IDs or file with a list of variant IDs
-r REGIDS, --regIDs REGIDS
    | Comma separated list of region IDs or file with a list of region IDs
-t TABLE, --table TABLE
    | Tab-separated file with
    | col1 (chr_pos_ref_alt)
    | col2 comma-separated list of region IDs
-o OUTPREFIX, --outprefix OUTPREFIX
    Prefix for output files
-g BUILD, --genome BUILD
    | Possible values: ``{GRCh37,GRCh38}``
    | Genome build for the query
--db GREENDB
    | Location of the GREEN-DB SQLite database file (.db)
--logfile LOGFILE
    | Custom location for the log file