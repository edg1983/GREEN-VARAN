SV_annotation tool usage
========================

SV_annotation.py allows annotation of structural variants VCF based on overlap with known regions in external bed files.
Overlap threshold is configurable and the resources provided with GREEN-VARAN allow annotation of population AF
from gnomAD and 1000G, overlapping genes and GREEN-DB information

.. code-block:: bash
    
    SV_annotation.py [-h] -i INPUTVCF 
                  -o OUT [-t TMPDIR] -b {GRCh37,GRCh38}
                  -c CONFIG [-k] [--logfile LOGFILE]

How it works
~~~~~~~~~~~~

SV annotation requires a standard VCF as input and output an uncompressed VCF with annotations.
Annotation files, overlap thresholds and other parameters can be configured modifying the configuration json file.
The input VCF must contain the following information:

- a unique variant ID in ID column
- SVTYPE and END fields in the INFO column. 
    | SVTYPE must follow the standard types definition: DEL, DUP, INS, INV, BND
    | The exact name of the INFO fields containing SVTYPE and END informations can be configured in the config file 

The standard settings will work directly on the output of popular SV caller like Lumpy, Manta and CANVAS

**Deletions, duplications and inversion**
SV_annotation annotate deletions, duplications and inversion by overlap with the set of known intervals provided in configuration.
The tool uses different overlap strategies for the different type of dataset provided in the configuration file (see the configuration section).
Considering a structural variants (A), an annotation region (B) and a given overlap threshold (T):

- AF datasets
    region A is annotated with AF from region B only if at least T fraction of A is overlapped by B
- custom datasets
    region A is annotated with information from region B only if at least T fraction of B is overlapped by A
- genes datasets
    work same as custom datasets

**Insertions**
Since the exact nature of an insertion is difficult to be determined, the tool will only try to annotate the break point where the insertion occurred.
Given that this is a single base position, it may be difficult that it overlaps exactly across different dataset.
As results, most insertions do not get an AF annotation, while you get information on which regions from custom dataset and genes are interrupted by the insertion

**BND**
Generic break-point, usually defined as BND in the VCF are not annotated.

Output annotations fields
#########################
Based on the dataset name and field names provided in the configuration file, the tool will add a single field
like ``DATASET_FIELDNAME=value`` for each configured annotation.

The configuration file
~~~~~~~~~~~~~~~~~~~~~~
A default configuration fils (SV_annotation.json) is provided in the GREEN-VARAN repository.
The configuration file is a standard json file oragnized as follows:

.. code-block:: json

    {
        "BEDTOOLS": "bedtools",
        "SVTYPE": "SVTYPE",
        "END": "END",
        "RES_DIR": "resources/SV_annotations",

        "overlap": {
            "AF_datasets": ["0.75", "FALSE"],
            "custom_datasets": ["0.10", "FALSE"],
            "genes": ["10e-9","FALSE"]
        },

        "AF_datasets": {
            "genome_build": {
                "Dataset1": [
                    ["INS", "INS_file.bed", "5", "AF", "Float"],
                    ["DEL", "DEL_file.bed", "5", "AF", "Float"],
                    ["DUP", "DUP_file.bed", "5", "AF", "Float"],
                    ["INV", "INV_file.bed", "5", "AF", "Float"]
                ]
            }
        },
        "custom_datasets": {
            "genome_build": {
                "Dataset1": [
                    ["INS", "INS_file.bed", "5", "AF", "Float"],
                    ["DEL", "DEL_file.bed", "5", "AF", "Float"],
                    ["DUP", "DUP_file.bed", "5", "AF", "Float"],
                    ["INV", "INV_file.bed", "5", "AF", "Float"]
                ]
            }
        },
        "genes": {
            "genome_build": {
                "gene": [
                    ["INS", "INS_file.bed", "5", "AF", "Float"],
                    ["DEL", "DEL_file.bed", "5", "AF", "Float"],
                    ["DUP", "DUP_file.bed", "5", "AF", "Float"],
                    ["INV", "INV_file.bed", "5", "AF", "Float"]
                ],
                "CDS": [
                    ["INS", "INS_file.bed", "5", "AF", "Float"],
                    ["DEL", "DEL_file.bed", "5", "AF", "Float"],
                    ["DUP", "DUP_file.bed", "5", "AF", "Float"],
                    ["INV", "INV_file.bed", "5", "AF", "Float"]
                ]
            }
        }
    }

Header
######
.. code-block::

    {
        "BEDTOOLS": "bedtools",
        "SVTYPE": "SVTYPE",
        "END": "END",
        "RES_DIR": "resources/SV_annotations",

These tags at the beginning of the file defines the location of bedtools executable and the 
exact INFO field names for SVTYPE and END. RES_DIR defines the folder containing the annotation files provided 
in the subsequent section. This folder is added before the file names and can be left empty when each file 
is provided in a different location.

Overlap
#######
.. code-block::

    "overlap": {
        "AF_datasets": ["0.75", "FALSE"],
        "custom_datasets": ["0.10", "FALSE"],
        "genes": ["10e-9","FALSE"]
    },

The ``overlap`` block defines the thresholds for each of the files types. 
Only the three annotation types defined above are accepted. The first value define the fraction of overlap
and the second value can be TRUE or FALSE and set if overlap must be reciprocal.

Datasets
########
.. code-block::

    "dataset_type": {
        "genome_build": {
            "Dataset1": [
                ["INS", "INS_file.bed", "5", "AF", "Float"],
                ["DEL", "DEL_file.bed", "5", "AF", "Float"],
                ["DUP", "DUP_file.bed", "5", "AF", "Float"],
                ["INV", "INV_file.bed", "5", "AF", "Float"]
            ]
        }
    }

For each accepted dataset type (AF_datasets, custom_datasets, genes) you can define a set of data sources 
for each genome build (like GRCh37, GRCh38). Within the genome_build block you define a dataset name which must
contain 4 files definition, one for each variant type (INS, DEL, DUP, INV).
Each data source contains the following setting:

- variant type
    Must be one of INS, DEL, DUP, INV
- annotation bed file location
    BED like files must be provided as input. First 3 columns are chrom,start,end
- comma-separated list of column numbers from which extract annotations
    For example to get values from column 4 and 5 use ``"4,5"``
- comma-separated list of field names to be used in INFO field
    The final field generated in the INFO output will be equal to ``dataset_fieldname``
- data type of annotation according to VCF standard (String, Integer, Float)
    single value expected, use String if you want to extract mixed values

Arguments list
~~~~~~~~~~~~~~
Mandatory Arguments
###################
-h, --help
    | Shows help message and exit
-i INPUTVCF, --inputvcf INPUTVCF
    | Input vcf[.gz] file
-o OUT, --out OUT
    | VCF output file (at the moment only support plain VCF output)
-b BUILD, --build BUILD 
    | Possible values: ``{GRCh37,GRCh38}``
    | Specify the genome build of input VCF
-c CONFIG_FILE, --config CONFIG_FILE
    | Configuration file (json)

Additional Arguments
####################
-t TMPDIR, --tmpdir TMPDIR
    | Location of temporaty folder to store temp files
    | By default a tmp folder will be created in the working directory
-k, --keeptmp
    | Set to keep temporary files
--logfile LOGFILE
    | Specify alternative location for the log file
