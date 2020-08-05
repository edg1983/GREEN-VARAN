Download resources
==================

To perform annotations GREEN-VARAN tools require a series of supporting files.
By default, various resources are expected in the ``resources`` folder within the tool folder.
The expected folder structure in resources is as follow

.. code-block:: bash

    .
    |-- SQlite
    |   `-- RegulatoryRegions.db
    |-- SV_annotations
    |   `-- BED files used for SV annotations
    |-- bed_files
    |   `-- GREEN-DB BED files
    |-- scores
    |   `-- scores annotations
    `-- AF
        `-- population AF annotations

When you clone the GREEN-VARAN repository a prepare_resources.sh script is provided in resources folder
to automatically download supporting files and setup the resources folder.

Otherwise, single files are available for download from Zenodo repository

.. csv-table::
    :header: "Annotation","Sub-folder","Link"
    :widths: 20,20,60

    GRCh37_CADD,scores,https://zenodo.org/record/3956385/GRCh37_CADD.tsv.gz
    GRCh37_CADD,scores,https://zenodo.org/record/3956385/GRCh37_CADD.tsv.gz.csi
    GRCh37_DANN,scores,https://zenodo.org/record/3957486/GRCh37_DANN.tsv.gz
    GRCh37_DANN,scores,https://zenodo.org/record/3957486/GRCh37_DANN.tsv.gz.csi
    GRCh37_ExPECTO,scores,https://zenodo.org/record/3956168/GRCh37_ExPECTO.tsv.gz
    GRCh37_ExPECTO,scores,https://zenodo.org/record/3956168/GRCh37_ExPECTO.tsv.gz.csi
    GRCh37_FIRE,scores,https://zenodo.org/record/3957356/GRCh37_FIRE.tsv.gz
    GRCh37_FIRE,scores,https://zenodo.org/record/3957356/GRCh37_FIRE.tsv.gz.csi
    GRCh37_LinSight,scores,https://zenodo.org/record/3956168/GRCh37_LinSight.bed.gz
    GRCh37_LinSight,scores,https://zenodo.org/record/3956168/GRCh37_LinSight.bed.gz.csi
    GRCh37_NCBoost,scores,https://zenodo.org/record/3956168/GRCh37_NCBoost.tsv.gz
    GRCh37_NCBoost,scores,https://zenodo.org/record/3956168/GRCh37_NCBoost.tsv.gz.csi
    GRCh37_ReMM,scores,https://zenodo.org/record/3956168/GRCh37_ReMM.tsv.gz
    GRCh37_ReMM,scores,https://zenodo.org/record/3956168/GRCh37_ReMM.tsv.gz.csi
    GRCh37_gnomAD,AF,https://zenodo.org/record/3957637/GRCh37_gnomad.genomes.vcf.gz
    GRCh37_gnomAD,AF,https://zenodo.org/record/3957637/GRCh37_gnomad.genomes.vcf.gz.csi
    GRCh38_CADD,scores,https://zenodo.org/record/3956227/GRCh38_CADD.tsv.gz
    GRCh38_CADD,scores,https://zenodo.org/record/3956227/GRCh38_CADD.tsv.gz.csi
    GRCh38_DANN,scores,https://zenodo.org/record/3957428/GRCh38_DANN.tsv.gz
    GRCh38_DANN,scores,https://zenodo.org/record/3957428/GRCh38_DANN.tsv.gz.csi
    GRCh38_ExPECTO,scores,https://zenodo.org/record/3955933/GRCh38_ExPECTO.tsv.gz
    GRCh38_ExPECTO,scores,https://zenodo.org/record/3955933/GRCh38_ExPECTO.tsv.gz.csi
    GRCh38_FIRE,scores,https://zenodo.org/record/3957216/GRCh38_FIRE.tsv.gz
    GRCh38_FIRE,scores,https://zenodo.org/record/3957216/GRCh38_FIRE.tsv.gz.csi
    GRCh38_LinSight,scores,https://zenodo.org/record/3955933/GRCh38_LinSight.bed.gz
    GRCh38_LinSight,scores,https://zenodo.org/record/3955933/GRCh38_LinSight.bed.gz.csi
    GRCh38_NCBoost,scores,https://zenodo.org/record/3955933/GRCh38_NCBoost.tsv.gz
    GRCh38_NCBoost,scores,https://zenodo.org/record/3955933/GRCh38_NCBoost.tsv.gz.csi
    GRCh38_ReMM,scores,https://zenodo.org/record/3955933/GRCh38_ReMM.tsv.gz
    GRCh38_ReMM,scores,https://zenodo.org/record/3955933/GRCh38_ReMM.tsv.gz.csi
    GRCh38_gnomAD,AF,https://zenodo.org/record/3957637/GRCh38_gnomad.genomes.vcf.gz
    GRCh38_gnomAD,AF,https://zenodo.org/record/3957637/GRCh38_gnomad.genomes.vcf.gz.csi
    SV_annotations,SV_annotations,https://zenodo.org/record/3970785/files/SV_annotations.tar.gz
    bed_files,bed_files,none