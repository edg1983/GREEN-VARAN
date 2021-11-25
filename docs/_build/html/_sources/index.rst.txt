.. GREEN-DB and GREEN-VARAN documentation master file, created by
   sphinx-quickstart on Mon Aug  3 15:14:10 2020.
   
GREEN: Genomic Regulatory Elements ENcyclopedia
===============================================

Welcome to the house of GREEN-DB and GREEN-VARAN!
This documentation describes the resources part of the Genomic Regulatory Elements Encyclopedia

The GREEN project is made by 2 main components:

1. The GREEN-DB collection
~~~~~~~~~~~~~~~~~~~~~~~~~~
The collection includes information useful for the annotation of non-coding variants in regulatory regions 

- a database (GREEN-DB) containing ~2.4M regulatory regions in the human genome with information on controlled gene(s) and tissue(s) of activity
- pre-processed indexed BED files representing functional genomic signals (TFBS, DNase peaks, UCNE, TADs)
- pre-processed indexed tables for 12 non-coding variant impact prediction scores and PhyloP100 conservation 

The GREEN-DB files can be downloaded from Zenodo: https://zenodo.org/record/5636209
All the additional pre-preprocessed datasets are also available from Zenodo, see the Download section.

1. The GREEN-VARAN tool set
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This include tools and workflows that can be used to interact with information in the GREEN-DB and annotate VCF files

- annotate small variants or structural variants with regulatory impact information, including possibly controlled genes
- add additional annotations on functional elements and non-coding prediction scores
- prioritize small variants for possible regulatory impact
- given a list of variants or regions, query the GREEN-DB for detailed information

Available from GitHub: https://github.com/edg1983/GREEN-VARAN

3. The prioritization workflow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#TODO

GREEN-DB and GREEN-VARAN are described in our preprint (https://doi.org/10.1101/2020.09.17.301960)

The Download section lists locations to download the GREEN-DB and other resource files for annotation

.. toctree::
   :maxdepth: 2
   :caption: Contents:
   
   GREEN_DB
   GREEN_VARAN
   nextflow_workflow
   Download
   how_to_cite


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
