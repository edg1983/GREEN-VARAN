# Resources folder

This folder contains resources needed for annotation using GREEN-VARAN tools.
If you prefer to save the resource files in another location, a custom locations can eventually be defined in the configuration file or through command-line options

## Default folder structure

```
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
```

## Prepare the resources folder

If you are fine with the default location just run

```
bash prepare_resources.sh
```

The script will automatically read resources from GREEN-VARAN_resources.txt and download the latest version of resource files. Required sub-folders will be created in the current folder.

Additional arguments are available to customize resources file name and output dir.

```
bash prepare_resources.sh
    [-r,--resource RESOURCES_FILE]
    [-n,--name RESOURCE_NAME] (default all = download all)
    [-o,--out OUTPUT_DIR] (default to current dir)
    [-l,--list] (list available resources)
    [-h,--help] (get help message)
```


## Download the resource files

RegDB and pre-assembled archives containing these resources are available for download from Zenodo (see GREEN-VARAN_resources.txt)
