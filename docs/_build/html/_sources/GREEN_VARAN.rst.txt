GREEN-VARAN tool set
====================

**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**

Annotate non-coding regulatory variants in a VCF with a combination of

- our GREEN-DB collection of regulatory regions
- non-coding variant prediction scores (NCBoost, FIRE, LinSight, ReMM, CADD, DANN, ExPecto)
- AF from gnomAD genomes
- conservation with PhyloP100

Installation
~~~~~~~~~~~~

1. Get the tools from the repository
####################################

``git clone https://github.com/edg1983/GREEN-VARAN.git``

2. Download the supporting files
################################

GREEN-DB files and a set of additional files are needed for annotation (see Downloads)
Supporting files are supposed to be in resources folder, but location can be configured in Configuration file 
or changed providing options on the command line. 
Supporting files can be obtained from Zenodo (see Download):

- GREEN-DB bed files and sqlite db file
- processed VCF from gnomAD genomes
- processed tables of prediction scores and phyloP100
- bed files for SV annotations

Requirements
~~~~~~~~~~~~
All tools are tested with Python >=3.4

**Python libraries**

+----------------+-------------------+----------------+
| GREEN-VARAN    | GREEN-DB_query    | SV_annotation  |
+================+===================+================+
| cyvcf2 >= 0.20 | cyvcf2 >= 0.20    | cyvcf2 >= 0.20 |
+----------------+-------------------+----------------+
| pandas >= 0.25 | pandas >= 0.25    |                |
+----------------+-------------------+----------------+
|                | sqlite3 >= 2.6.0  |                |
+----------------+-------------------+----------------+

**External software and libs**

- bedtools >= 2.27
- htslib >= 1.10
- vcfanno >= 0.3.0 (pre-compiled binary is included in GREEN-VARAN release)

Singularity
~~~~~~~~~~~
GREEN-VARAN tool set is also provided as Singularity image (tested on singularity >= 3.2). 
A Singularity recipe is included in the repository or you can download a pre-compiled image from zenodo(LINK).

Usage
#####

The image contains all 3 GREEN-VARAN tools:

- GREEN-VARAN: annotation of small variants VCF
- SV_annotation: annotation of SV VCF
- GREEN-DB_query: query the GREEN-DB for detailed information

To run one of the tool use:

.. code-block:: bash

    standard run
    singularity run \
    --bind resources_folder:/opt/GREEN_VARAN/resources \
    GREEN-VARAN.sif tool_name [tool arguments]


**NB.** The host resources_folder must contain the standard subfolders and files expected by GREEN-VARAN.

Bind specific folders for resources or data
###########################################
If you have stored resources in other locations you have to bind them manually into the container and 
then pass the mounted path to the tool with the corresponding argument (see below)
By default you are expected to read/write from the present working directory. This means that

- all input files, including eventual configuration file, must be present in the working directory
- output files and eventual tmp folder are created in the present working directory

Example if you have input/output/resources in other folders

.. code-block:: bash

    singularity run \
    --bind resources_folder:/opt/GREEN_VARAN/resources \
    --bind input_folder:/input \
    --bind output_fodler:/output \
    --bind bed_files:/bed_files \
    GREEN-VARAN.sif \
    GREEN-VARAN -i /input/input.vcf \
    -o /output/output.vcf \
    --bed_dir /bed_files \ 
    [tool arguments]

Single tools usage
~~~~~~~~~~~~~~~~~~
The GREEN-VARAN tool set includes 3 tools to annotate variants and interact with GREEN-DB.

1. GREEN_VARAN.py
    Perform annotation on small variants VCF. Provides also abiliy to filter for genes of interest,
    select variants on genes already affected by a coding variants and a prioritization function for
    regulatory variants
2. SV_annotation.py
    This tool allow annotation of SV VCF based on overlap with known regions in external bed files.
    Overlap threshold is configurable and the provided resources allow annotation of population AF
    from gnomAD and 1000G, genes overlap and GREEN-DB
3. GREEN-DB_query.py
    This tool assist in quering the GREEN-DB. Given a list of region IDs or a VCF annotated by GREEN-VARAN
    the tool generates a set of tables containing detailed information on the regions of interest

For detailed instruction on the single tools usage please refer to the corresponding page

.. toctree::
   :maxdepth: 2
   
   GREEN_DB
   GREEN_VARAN
   Download

