# GREEN-VARAN workflow

To perform small variants prioritization as described in the GREEN-DB manuscript, GREEN-VARAN need some annotations to be already present in your input VCF (see "Prioritization of small variants" in main README)

This Nextflow workflow automate the whole process annotating additional information and then performing greenevaran annotation. The workflow is tested on Nextflow >=v20.10

## Usage

The typical usage scenario start with a VCF file already containing gene consequences annotations from SnpEff or bcftools. Then from the GREEN-VARAN tool main folder you can perform all annotations using the following command.

This will add a minimum set of information to you VCF including:
- population allele AF from gnomAD genomes v3.1.1 (GRCh38) or v2.1.1 (GRCh37)
- functional regions overlaps for TFBS, DNase peaks and UCNE
- prediction score values for ncER, FATHMM, ReMM
- GREEN-DB information on regulatory variants with prioritization levels

```    
nextflow workflow/main.nf \
    -profile local \
    --input input_file.vcf.gz \
    --build GRCh38 \
    --out results \
    --scores best \
    --regions best \
    --AF \
    --greenvaran_config config/prioritize_smallvars.json \
    --greenvaran_dbschema config/greendb_schema_v2.5.json
```

If requested annotation files are missing, they will be automatically downloaded in the default location (`resources` folder within the main GREEN-VARAN folder) 

Note that `--input` can accept multiple vcf.gz files using a pattern like `inputdir/*.vcf.gz`

### Add additional custom annotations

If you have additional custom annotation you want to add to your VCF before greenvaran processing they can be configured in a .toml 
and then you can pass this file to the workflow using `--anno_toml`.

A toml file is a annotation configuration file used by the vcfanno tool and is described in [vcfanno repository](https://github.com/brentp/vcfanno)

A minimal example is reported below

```
[[annotation]]
file="ExAC.vcf" #source file
fields = ["AF", "AF_nfe"] #INFO fields to be extracted from source
ops=["self", "max"] #How to treat source values
names=["exac_af", "exac_af_nfe_max"] #names used in the annotated file

[[annotation]]
file="regions_score.bed.gz"
columns = [4, 5] #When using a BED or TSV files you can refer to values by col index
names=["regions_ids", "score_max"]
ops=["uniq","max"]
```

## Resources

To perform annotations GREEN-VARAN Nextflow workflow requires a series of supporting files.
By default, various resources are expected in the ``resources`` folder within the main tool folder.
You pass an alternative resource folder using ``--resource_folder`` option, bug the same structure is expected in this folder

The expected folder structure is as follows

.. code-block::

    .
    |-- SQlite
    |   `-- GREEN-DB_v2.5.db
    |-- GRCh37
    |   `-- BED / TSV files used for GRCh37 genome build
    `-- GRCh38
        `-- BED / TSV files used for GRCh38 genome build

Use the `--list_data` option to see the full list of available resources and the expected path for each one.

### Automated download

A supporting workflow is provided to automate data download for all resources included in the GREEN-DB collection. 
You can list the available resources and their resulting download location using

`nextflow workflow/download.nf --list_data`

The reccomended set of annotations can be downloaded to the default location using the following command or
you can set an alternative resource folder using `--resource_folder` option

```
nextflow workflow/download.nf \
    -profile local \
    --scores best \
    --regions best \
    --AF \
    --db 
```

Otherwise, single files are available for download from Zenodo repository and all file locations are listed in 
the ``GREENDB_collection.txt`` file under resources folder.

## Workflow configuration

The workflow has pre-configured profiles for most popular schedulers (sge, lsf, slurm) and also a local profile (local).
These profiles determine how many download jobs can be submitted concurrently and the number of threads used for annotation.

You can activate the desired profile using ``-profile`` argument when launching the workflow

**NB.** You need to update the queue name parameter to reflect your local settings, see how to edit the config below

For more details on how to configure the workflow please refer to the [documentation](https://green-varan.readthedocs.io/en/latest)


## Available parameters for main workflow

The complete description of availabel parameters can be obtain with `--help` or in the the [documentation](https://green-varan.readthedocs.io/en/latest)