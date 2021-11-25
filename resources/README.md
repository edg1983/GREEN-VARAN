# Resources folder

This folder store resources needed for annotation using GREEN-VARAN tools. By default GREEN-VARAN Nextflow workflow searches for files in "resources" folder within the main tool folder. The expected folder structure is reported below.
If you prefer to save the resource files in another location, a custom location can eventually be defined using `--resource_folder` option. Note that even when a custom resource folder is specified, the workflow expect the same sub-folders structure.

## Default sub-folder structure
These are the sub-folders expected in the main resources folder

```
.
|-- SQlite
|   `-- GREEN-DB_v2.5.db
|-- GRCh37
|   `-- BED / TSV files used for GRCh37 genome build
`-- GRCh38
    `-- BED / TSV files used for GRCh38 genome build

```

## Minimum set needed for prioritization

The minimum set of annotations needed to reproduce the four level prioritization described in our paper are:
- population AF (from gnomAD)
- flags marking the overlap with TFBS, DNase peaks, UCNE
- prediction scores from ncER, FATHMM-MKL, ReMM

## Prepare the resources folder

If you are fine with the default location and you have Nextflow >= v20.10 installed you can downloaded the essential datasets for GRCh38 or GRCh37 using the following command (adjust `--build` as needed) from the GREEN-VARAN main folder:

```
nextflow workflow/download.nf -profile local \
    --scores best \
    --regions best \
    --db \
    --AF \
    --build GRCh38
```

The script will automatically download the latest version of resource files and required sub-folders will be created in the resources folder.

You can specify a different resources main folder using `--resource_folder`.

To see the complete list of datasets available for download use 

```
nextflow workflow/download.nf --list_data
```

## Manually download the resource files

GREEN-DB BED / SQlite files and pre-processed files for additional regions and scores are available for download from Zenodo (see GREENDB_collection.txt).
