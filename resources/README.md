# Resources folder

This folder store resources needed for annotation using GREEN-VARAN tools. By default GREEN-VARAN tool search for files in "resources" folder within the tool folder. The expected folder structure is reported below.
If you prefer to save the resource files in another location, a custom location can eventually be defined in the configuration file (RES_DIR) or through command-line options when running each tool.

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
prepare_resources.sh
    [-r,--resource RESOURCES_FILE]  GREEN-VARAN_resources.txt (default)
    [-n,--name RESOURCE_NAME]       all = download all (default)
                                    scores = download all scores
                                    bed_files = download all bed_files
                                    AF = download AF files
                                    SV_annotation = download all SV annotations
                                    specific_name (see --list)
    [-o,--out OUTPUT_DIR]           current directory (default)
                                    folder where to store downloaded files
                                    expected sub-folders will be created within this folder
    [-l,--list]                     list available resources
    [-h,--help]                     get help message
```


## Download the resource files

RegDB and pre-assembled archives containing these resources are available for download from Zenodo (see GREEN-VARAN_resources.txt)
