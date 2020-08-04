### RESOURCES LOCATIONS ###
# Modify values according to file locations
# {base_dir} is translated to the directory where python script is located

VCF_ANNO = "{base_dir}/vcfanno"
BCFTOOLS = "bcftools"

RESOURCE_DIR = "resources"
BED_DIR = '{base_dir}/' + RESOURCE_DIR + '/bed_files'
SCORES_DIR = '{base_dir}/' + RESOURCE_DIR + '/scores'
AF_FILE='{base_dir}/' + RESOURCE_DIR + '/AF/{build}_gnomad.genomes.vcf.gz'

### ANNOTATION CONFIGURATION ###
# You are not supposed to modify this part
NC_ANNO_DATATYPES = {
    'NC_support': 'Float',
    'NC_regionID': 'String',
    'NC_region_type': 'String',
    'NC_constraint': 'Float',
    'NC_methods': 'Integer',
    'NC_genes': 'String',
    'NC_closestGene' : 'String',
    'NC_closestGene_dist' : 'Integer',
    'NC_closestProt' : 'String',
    'NC_closestProt_dist' : 'Integer',
    'NC_tolerant_P' : 'Float',
    'NC_tolerant_label' : 'String',
    'NC_median_PhyloP100' : 'Float'}
NC_ANNO_TAGS = NC_ANNO_DATATYPES.keys()

HEADER_VOI = '##INFO=<ID=NC_VOI,Number=1,Type=Integer,Description="This is a variant of interest from NC annotation"'
HEADER_CLASS = '##INFO=<ID=NC_VARCLASS,Number=1,Type=Integer,Description="Variant class based on GREEN-VARAN prioritization. Higher values means higher probability of regulatory impact"'
HEADER_SINGLE = '##INFO=<ID=NC_ANNO,Number=1,Type=String,Description="' + '|'.join(NC_ANNO_TAGS) + '">'
HEADER_SEPARATED = []
for tag, datatype in NC_ANNO_DATATYPES.items():
    HEADER_SEPARATED.append('##INFO=<ID={tag},Number=1,Type={datatype},Description="Non-coding region annotation: {tag}">'.format(tag=tag,datatype=datatype))

TOML = '''[[annotation]]
file="{bed_dir}/{build}_GREEN-DB.bed.gz"
names=["TMPNC_regionID","TMPNC_region_type","TMPNC_constraint","TMPNC_median_PhyloP100","TMPNC_closestGene","TMPNC_closestGene_dist", "TMPNC_closestProt", "TMPNC_closestProt_dist","TMPNC_genes","TMPNC_methods"]
ops=["self","self","self","max","uniq","uniq","uniq","uniq","uniq","max"]
columns=[4,5,7,8,9,10,11,12,13,14]

[[annotation]]
file="{bed_dir}/{build}_TFBS.merged.bed.gz"
names=["NC_TFname"]
ops=["uniq"]
columns=[4]

[[annotation]]
file="{bed_dir}/{build}_UCNE.sorted.bed.gz"
names=["NC_UCNE"]
columns=[3]
ops = ["flag"]

[[annotation]]
file="{bed_dir}/{build}_dbSuper.sorted.bed.gz"
names=["NC_dbSUPER"]
columns=[3]
ops = ["flag"]

[[annotation]]
file="{bed_dir}/{build}_DNase.merged.bed.gz"
names=["NC_DNase"]
columns=[3]
ops = ["flag"]

[[annotation]]
file="{bed_dir}/{build}_LoF_tolerance.sorted.bed.gz"
names=["TMPNC_tolerant_P","TMPNC_tolerant_label"]
ops=["max","uniq"]
columns=[5,6]
'''

#[[annotation]]
#file="{bed_dir}/{build}_SegWey.sorted.bed.gz"
#names=["NC_SegWey_mean", "NC_SegWey_mean_pct","NC_SegWey_sum","NC_SegWey_sum_pct"]
#ops=["max","max","max","max"]
#columns=[7,8,9,10]

TOML_SCORES = {
'FIRE' : '''[[annotation]]
file="{scores_dir}/{build}_FIRE.tsv.gz"
names=["NC_FIRE"]
ops=["max"]
columns=[5]''',

'LinSight': '''[[annotation]]
file="{scores_dir}/{build}_LinSight.bed.gz"
names=["NC_LinSight"]
ops=["max"]
columns=[4]''',

'ExPECTO': '''[[annotation]]
file="{scores_dir}/{build}_ExPECTO.tsv.gz"
names=["NC_ExPECTO"]
ops=["self"]
columns=[8]''',

'ReMM': '''[[annotation]]
file="{scores_dir}/{build}_ReMM.tsv.gz"
names=["NC_ReMM"]
ops=["max"]
columns=[3]''',

'NCBoost': '''[[annotation]]
file="{scores_dir}/{build}_NCBoost.tsv.gz"
names=["NC_NCBoost"]
ops=["max"]
columns=[5]''',

'CADD': '''[[annotation]]
file="{scores_dir}/{build}_CADD.tsv.gz"
names=["NC_CADD"]
ops=["max"]
columns=[5]''',

'DANN': '''[[annotation]]
file="{scores_dir}/{build}_DANN.tsv.gz"
names=["NC_DANN"]
ops=["max"]
columns=[5]'''
}

TOML_AF = '''[[annotation]]
file="{AF_file}"
names=["gnomAD_AF"{pop_tag}]
ops=["self","self"]
fields=["AF"{pop_field}]
'''