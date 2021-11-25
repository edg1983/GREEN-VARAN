GREEN-VARAN tool set
====================

**Genomic Regulatory Elements ENcyclopedia VARiant ANnotation**

Annotate variants in a VCF using GREEN-DB to provide information on non-coding regualtory variants and the controlled genes.
Additionally perform prioritization summing up evidences of regulatory impact from GREENDB, population AF, functional regions and prediction scores

Installation
~~~~~~~~~~~~

1. Get the tool binary from the repository
##########################################

The easiest way to run GREEN-VARAN is to download the pre-compiled binaries from the latest release at https://github.com/edg1983/GREEN-VARAN

2. Compile the tool
###################

Alternatively, you can clone the repository 
``git clone https://github.com/edg1983/GREEN-VARAN.git``

And then compile the greenvaran using Nim compiler (https://nim-lang.org/). 
GREEN-VARAN requires
- nim >= 0.10
- hts-nim >= 0.3.4
- argparse 0.10.1 

If you have Singularity installed, you can use the script ``nim_compile.sh`` to create a static binary with no dependencies 
This uses musl-hts-nim as described in hts-nim repository (see https://github.com/brentp/hts-nim#static-binary-with-singularity)

Get GREEN-DB files
~~~~~~~~~~~~~~~~~~

To perform annotations with GREEN-VARAN you will need the GREEN-DB bed files for your genome build.
You can download the GREEN-DB BED file for GRCh37 or GRCh38 from https://zenodo.org/record/5636209

The complete SQLite database is also available from the same repository

GREEN-VARAN Nextflow workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
We also provide a Nextflow workflow that can be used to automate VCF annotation and resource download.
Given a small variants VCF annotated for gene consequences using snpEff or bcftools the workflow can be used to
- automatically add functional regions annotations and non-coding prediction scores
- perform greenvaran prioritization

Missing datasets will be downloaded automatically during the process.
See the dedicated page for more usage information


Singularity
~~~~~~~~~~~
The tool binaries should work on most linux based system. In case you have any issue, we also provdie GREEN-VARAN as Singularity image (tested on singularity >= 3.2). 
A Singularity recipe is included in the repository or you can pull the image from Singularity Library using

``singularity pull library://edg1983/greenvaran/greenvaran:latest``

See GREEN-VARAN usage for more details
Usage
#####

The image contains both greenvaran and greendb_query tools.
The general usage is:

.. code-block:: bash

    singularity run \
    greenvaran.sif \
    tool_name [tool arguments]

Bind specific folders for resources or data
###########################################

The tool need access to input VCF file, required GREEN-DB bed file and config files so remember to bind the corresponding locations in the container 

See the following example where we use the current working directory for input/output, while other files are located
in the default config / resources folder within greenvaran folder. In the example we use GRCh38 genome build

.. code-block:: bash

    singularity run \
    --bind /greenvaran_path/resources/GRCh38:/db_files \
    --bind /greenvaran_path/config:/config_files \
    --bind ${PWD}:/data \
    greenvaran.sif \
    greenvaran -i /data/input.vcf.gz \
    -o /data/output.vcf.gz \
    --db /db_files/GRCh38_GREEN-DB.bed.gz \
    --dbschema /config_files/greendb_schema_v2.5.json
    --config /config_files/prioritize_smallvars.json
    [additional tool arguments]

Single tools usage
~~~~~~~~~~~~~~~~~~
The GREEN-VARAN tool set includes 2 main tools to annotate variants and interact with GREEN-DB.

1. greenvaran
   Perform annotation on small variants or structural variants VCF. 
   Provides prioritization of regulatory variants summing up evidences of impact from GREENDB, population AF, functional regions and prediction scores.
   Variants can also be tagged based on a list of genes of interest.
   Finally, the tool can update standard gene consequence in ANN or BCQS fields to reflect regulated genes.
2. greendb_query
   Assists in quering the GREEN-DB. Given a list of region IDs, a list of variants or a table of variants and relevant GREENDB regions
   the tool generates a set of tables containing detailed information on the regions of interest, region-gene connections, functional regions and tissues.

For detailed instruction on the single tools usage please refer to the corresponding page

.. toctree::
   :maxdepth: 2
   
   GREEN-VARAN_usage
   GREEN-DB_query_usage
   nextflow_workflow

