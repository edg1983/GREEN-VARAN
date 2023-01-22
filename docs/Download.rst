Download resources
==================

greenvaran tool
~~~~~~~~~~~~~~~
The greenvaran annotation tool only need the GREEN-DB BED file and index for your genome build available from
https://zenodo.org/record/5636209


greendb_query tool
~~~~~~~~~~~~~~~~~~
The greendb query tool only need the GREEN-DB SQlite file (.db.gz) available from https://zenodo.org/record/5636209
Remember to decompress this before use


GREEN-VARAN workflow
~~~~~~~~~~~~~~~~~~~~

To perform annotations GREEN-VARAN Nextflow workflow requires a series of supporting files.
By default, various resources are expected in the ``resources`` folder within the main tool folder.
If you pass an alternative resource folder using ``--resource_folder`` option, the same structure is expected in this folder
The expected folder structure is as follows and the expected file names are those listed in the Zenodo repository table below

.. code-block:: bash

    .
    |-- SQlite
    |   `-- GREEN-DB_v2.5.db
    |-- GRCh37
    |   `-- BED / TSV files used for GRCh37 genome build
    `-- GRCh38
        `-- BED / TSV files used for GRCh38 genome build

When you clone the GREEN-VARAN repository you can use the Nextflow workflow ``workflow/download.nf`` to download files and prepare the resource folder.
Use the ``--list_data`` option to see the full list of available resource and the expected path for each one.

Otherwise, single files are available for download from Zenodo repository

.. csv-table::
    :header: "Annotation","Category","File"
    :widths: 20,20,60

    GRCh37_CADD,scores,https://zenodo.org/record/3956385/files/GRCh37_CADD.tsv.gz
    GRCh37_CADD,scores,https://zenodo.org/record/3956385/files/GRCh37_CADD.tsv.gz.csi
    GRCh37_DANN,scores,https://zenodo.org/record/3957486/files/GRCh37_DANN.tsv.gz
    GRCh37_DANN,scores,https://zenodo.org/record/3957486/files/GRCh37_DANN.tsv.gz.csi
    GRCh37_ExPECTO,scores,https://zenodo.org/record/3956168/files/GRCh37_ExPECTO.tsv.gz
    GRCh37_ExPECTO,scores,https://zenodo.org/record/3956168/files/GRCh37_ExPECTO.tsv.gz.csi
    GRCh37_FIRE,scores,https://zenodo.org/record/3957356/files/GRCh37_FIRE.tsv.gz
    GRCh37_FIRE,scores,https://zenodo.org/record/3957356/files/GRCh37_FIRE.tsv.gz.csi
    GRCh37_LinSight,scores,https://zenodo.org/record/3956168/files/GRCh37_LinSight.bed.gz
    GRCh37_LinSight,scores,https://zenodo.org/record/3956168/files/GRCh37_LinSight.bed.gz.csi
    GRCh37_NCBoost,scores,https://zenodo.org/record/3956168/files/GRCh37_NCBoost.tsv.gz
    GRCh37_NCBoost,scores,https://zenodo.org/record/3956168/files/GRCh37_NCBoost.tsv.gz.csi
    GRCh37_ReMM,scores,https://zenodo.org/record/3956168/files/GRCh37_ReMM.tsv.gz
    GRCh37_ReMM,scores,https://zenodo.org/record/3956168/files/GRCh37_ReMM.tsv.gz.csi
    GRCh37_PhyloP100,scores,https://zenodo.org/record/3973181/files/GRCh37_PhyloP100.bed.gz
    GRCh37_PhyloP100,scores,https://zenodo.org/record/3973181/files/GRCh37_PhyloP100.bed.gz.csi
    GRCh37_Eigen,scores,https://zenodo.org/record/3982095/files/GRCh37_Eigen.tsv.gz
    GRCh37_Eigen,scores,https://zenodo.org/record/3982095/files/GRCh37_Eigen.tsv.gz.csi
    GRCh37_FATHMM_XF,scores,https://zenodo.org/record/3982392/files/GRCh37_FATHMM-XF_NC.tsv.gz
    GRCh37_FATHMM_XF,scores,https://zenodo.org/record/3982392/files/GRCh37_FATHMM-XF_NC.tsv.gz.csi
    GRCh37_FATHMM_MKL,scores,https://zenodo.org/record/3981113/files/GRCh37_FATHMM-MKL_NC.tsv.gz
    GRCh37_FATHMM_MKL,scores,https://zenodo.org/record/3981113/files/GRCh37_FATHMM-MKL_NC.tsv.gz.csi
    GRCh37_GWAVA,scores,https://zenodo.org/record/3956168/files/GRCh37_gwava.bed.gz
    GRCh37_GWAVA,scores,https://zenodo.org/record/3956168/files/GRCh37_gwava.bed.gz.csi
    GRCh37_gnomAD,AF,https://zenodo.org/record/3957637/files/GRCh37_gnomad.genomes.vcf.gz
    GRCh37_gnomAD,AF,https://zenodo.org/record/3957637/files/GRCh37_gnomad.genomes.vcf.gz.csi
    GRCh37_ncER,scores,https://zenodo.org/record/5636163/files/GRCh37_ncER_perc.bed.gz
    GRCh37_ncER,scores,https://zenodo.org/record/5636163/files/GRCh37_ncER_perc.bed.gz.csi
    GRCh38_CADD,scores,https://zenodo.org/record/3956227/files/GRCh38_CADD.tsv.gz
    GRCh38_CADD,scores,https://zenodo.org/record/3956227/files/GRCh38_CADD.tsv.gz.csi
    GRCh38_DANN,scores,https://zenodo.org/record/3957428/files/GRCh38_DANN.tsv.gz
    GRCh38_DANN,scores,https://zenodo.org/record/3957428/files/GRCh38_DANN.tsv.gz.csi
    GRCh38_ExPECTO,scores,https://zenodo.org/record/3955933/files/GRCh38_ExPECTO.tsv.gz
    GRCh38_ExPECTO,scores,https://zenodo.org/record/3955933/files/GRCh38_ExPECTO.tsv.gz.csi
    GRCh38_FIRE,scores,https://zenodo.org/record/3957216/files/GRCh38_FIRE.tsv.gz
    GRCh38_FIRE,scores,https://zenodo.org/record/3957216/files/GRCh38_FIRE.tsv.gz.csi
    GRCh38_LinSight,scores,https://zenodo.org/record/3955933/files/GRCh38_LinSight.bed.gz
    GRCh38_LinSight,scores,https://zenodo.org/record/3955933/files/GRCh38_LinSight.bed.gz.csi
    GRCh38_NCBoost,scores,https://zenodo.org/record/3955933/files/GRCh38_NCBoost.tsv.gz
    GRCh38_NCBoost,scores,https://zenodo.org/record/3955933/files/GRCh38_NCBoost.tsv.gz.csi
    GRCh38_ReMM,scores,https://zenodo.org/record/3955933/files/GRCh38_ReMM.tsv.gz
    GRCh38_ReMM,scores,https://zenodo.org/record/3955933/files/GRCh38_ReMM.tsv.gz.csi
    GRCh38_PhyloP100,scores,https://zenodo.org/record/3973181/files/GRCh38_PhyloP100.bed.gz
    GRCh38_PhyloP100,scores,https://zenodo.org/record/3973181/files/GRCh38_PhyloP100.bed.gz.csi
    GRCh38_Eigen,scores,https://zenodo.org/record/3982182/files/GRCh38_Eigen.tsv.gz
    GRCh38_Eigen,scores,https://zenodo.org/record/3982182/files/GRCh38_Eigen.tsv.gz.csi
    GRCh38_FATHMM_XF,scores,https://zenodo.org/record/3982484/files/GRCh38_FATHMM-XF_NC.tsv.gz
    GRCh38_FATHMM_XF,scores,https://zenodo.org/record/3982484/files/GRCh38_FATHMM-XF_NC.tsv.gz.csi
    GRCh38_FATHMM_MKL,scores,https://zenodo.org/record/3981121/files/GRCh38_FATHMM-MKL_NC.tsv.gz
    GRCh38_FATHMM_MKL,scores,https://zenodo.org/record/3981121/files/GRCh38_FATHMM-MKL_NC.tsv.gz.csi
    GRCh38_GWAVA,scores,https://zenodo.org/record/3955933/files/GRCh38_gwava.bed.gz
    GRCh38_GWAVA,scores,https://zenodo.org/record/3955933/files/GRCh38_gwava.bed.gz.csi
    GRCh38_ncER,scores,https://zenodo.org/record/5636163/files/GRCh38_ncER_perc.bed.gz
    GRCh38_ncER,scores,https://zenodo.org/record/5636163/files/GRCh38_ncER_perc.bed.gz.csi
    GRCh37_TFBS,regions,https://zenodo.org/record/5705936/files/GRCh37_TFBS.merged.bed.gz
    GRCh37_TFBS,regions,https://zenodo.org/record/5705936/files/GRCh37_TFBS.merged.bed.gz.csi
    GRCh37_DNase,regions,https://zenodo.org/record/5705936/files/GRCh37_DNase.merged.bed.gz
    GRCh37_DNase,regions,https://zenodo.org/record/5705936/files/GRCh37_DNase.merged.bed.gz.csi
    GRCh37_UCNE,regions,https://zenodo.org/record/5705936/files/GRCh37_UCNE.bed.gz
    GRCh37_UCNE,regions,https://zenodo.org/record/5705936/files/GRCh37_UCNE.bed.gz.csi
    GRCh37_dbSuper,regions,https://zenodo.org/record/5705936/files/GRCh37_dbSuper.bed.gz
    GRCh37_dbSuper,regions,https://zenodo.org/record/5705936/files/GGRCh37_dbSuper.bed.gz.csi
    GRCh37_TAD,regions,https://zenodo.org/record/5705936/files/GRCh37_TAD.bed.gz
    GRCh37_TAD,regions,https://zenodo.org/record/5705936/files/GRCh37_TAD.bed.gz.csi
    GRCh38_TFBS,regions,https://zenodo.org/record/5705936/files/GRCh38_TFBS.merged.bed.gz
    GRCh38_TFBS,regions,https://zenodo.org/record/5705936/files/GRCh38_TFBS.merged.bed.gz.csi
    GRCh38_DNase,regions,https://zenodo.org/record/5705936/files/GRCh38_DNase.merged.bed.gz
    GRCh38_DNase,regions,https://zenodo.org/record/5705936/files/GRCh38_DNase.merged.bed.gz.csi
    GRCh38_UCNE,regions,https://zenodo.org/record/5705936/files/GRCh38_UCNE.bed.gz
    GRCh38_UCNE,regions,https://zenodo.org/record/5705936/files/GRCh38_UCNE.bed.gz.csi
    GRCh38_dbSuper,regions,https://zenodo.org/record/5705936/files/GRCh38_dbSuper.bed.gz
    GRCh38_dbSuper,regions,https://zenodo.org/record/5705936/files/GRCh38_dbSuper.bed.gz.csi
    GRCh38_TAD,regions,https://zenodo.org/record/5705936/files/GRCh38_TAD.bed.gz
    GRCh38_TAD,regions,https://zenodo.org/record/5705936/files/GRCh38_TAD.bed.gz.csi
    GRCh38_gnomAD,AF,https://zenodo.org/record/3957637/files/GRCh38_gnomad.genomes.vcf.gz
    GRCh38_gnomAD,AF,https://zenodo.org/record/3957637/files/GRCh38_gnomad.genomes.vcf.gz.csi
    SV_annotations,SV_annotations,https://zenodo.org/record/3970785/files/SV_annotations.tar.gz
    GRCh37_GREENDB_bed,GREENDB_bed,https://zenodo.org/record/5636209/files/GRCh37_GREEN-DB.bed.gz
    GRCh37_GREENDB_bed,GREENDB_bed,https://zenodo.org/record/5636209/files/GRCh37_GREEN-DB.bed.gz.csi
    GRCh38_GREENDB_bed,GREENDB_bed,https://zenodo.org/record/5636209/files/GRCh38_GREEN-DB.bed.gz
    GRCh38_GREENDB_bed,GREENDB_bed,https://zenodo.org/record/5636209/files/GRCh38_GREEN-DB.bed.gz.csi
    GREENDB_sqlite,GREENDB_sqlite,https://zenodo.org/record/5636209/files/GREEN-DB_v2.5.db.gz