'''
SV annotation
Take a SV VCF and annotate it using fields from other VCF or BED files
Can add gene based annotations as well using a GFF file

Author: Edoardo Giacopuzzi
Email: edoardo.giacopuzzi@well.ox.ac.uk
'''

import json
import sys, os, time
import gzip
import argparse
import subprocess
import tempfile
import re
import logging
import pandas as pd
from shutil import which, rmtree
from datetime import datetime, timedelta
from collections import OrderedDict
from cyvcf2 import VCF

COL_OPERATIONS = {
    'String': 'unique',
    'Float': 'max',
    'Integer': 'max'
}

class INFO():
	def __init__(self, info_field):
		exp = re.compile('(.+)=(.+)')
		self.infos = OrderedDict()
		info_tags = info_field.split(";")
		for value in info_tags:
			m = exp.match(value)
			if m:
				self.infos[m.group(1)] = m.group(2)
			else:
				self.infos[value] = True
	
	def selectINFO(self, info_keys):
		if isinstance(info_keys, list):
			self.infos = { key: self.infos[key] for key in info_keys }
		else:
			self.infos = { info_keys : self.infos[info_keys]}
	
	def dropINFO(self, info_keys):
		if isinstance(info_keys, list):
			for key in info_keys: self.infos.pop(key)
		else:
			self.infos.pop(key)
    
	def addINFO(self, info_key, value):
		self.infos[info_key] = value
    
	def updateINFO(self, info_key, value):
		if info_key in self.infos.keys():
			self.infos[info_key] = value
		else:
			self.addINFO(info_key, value)
	
	def getINFO(self,format):
		if format == "dict":
			return self.infos
		elif format == "string":
			output = [] 
			for key, value in self.infos.items():
				output.append(key + "=" + str(value))
		return ';'.join(output)

def get_stdout(cmd, shell=False):
	proc = subprocess.Popen(cmd, shell=shell,
		stdout = subprocess.PIPE,
		stderr = subprocess.DEVNULL)
    
	for line in iter(proc.stdout.readline, b''):
		yield line.decode('utf-8')

def tokenize(line,sep):
	line = line.rstrip('\n')
	line = line.split(sep)
	return line

def reciprocalOpt(option):
    if option == "TRUE": 
        reciprocal = ['-r'] 
    else:
        reciprocal = []
    return reciprocal

def computeOverlap(dataset, group, a_file, b_file, overlap, reciprocal, cols_idx, fields, tmp_dir):
    filename = tmp_dir + "/SVannot.tmp_" + dataset + "_" + group
    out_file = open(filename, "w+")
    header = ['ID'] + [dataset + "_" + x for x in fields.split(",")]
    out_file.write("\t".join(header) + "\n")
    cols_idx = [3] + [int(x) + 3 for x in cols_idx.split(",")]
    if group == "INS": overlap[1] = '10e-9'
    for line in get_stdout(bedtools + overlap + reciprocal + ['-a',a_file,'-b',b_file]):
        line = tokenize(line, "\t")
        out_file.write("\t".join([line[i] for i in cols_idx]) + "\n")  
    out_file.close()    

    mydf = pd.read_csv(filename, sep="\t", index_col="ID")
    return mydf

def setDataType(dataset, fields, data_types):
    tags = [dataset + "_" + x for x in fields.split(",")]
    data_types = data_types.split(",")
    if len(data_types) == 1:
        tags_dict = {k:data_types[0] for k in tags}
    else:
        tags_dict = {k:v for k,v in zip(tags, data_types)}
    
    return tags_dict

def readFileConfig(file_config, res_dir):
    group = file_config[0]
    b_file = res_dir + "/" + file_config[1]
    cols_idx = file_config[2]
    fields = file_config[3]
    data_types = file_config[4]
    return group,b_file,cols_idx,fields,data_types

def flat(input_list):
    flat_list = [item for sublist in input_list if isinstance(sublist,list) for item in sublist]
    return flat_list

#check file exists and trigger proper actions
def checkfile(files, operation="open", mode="stop", logger=None):
    if isinstance(files, str): files = [files]
    
    outfiles = []
    for myfile in files:
        if operation == "write":
            if os.path.isfile(myfile):
                if mode == "rename":
                    if myfile.endswith(".vcf.gz"):
                        filebase = myfile.replace(r'.vcf.gz$','')
                        extension = ".vcf.gz"
                    elif myfile.endswith(".vcf"):
                        filebase = myfile.replace(r'.vcf$','')
                        extension = ".vcf"
                    newfilename = filebase + '.' + now("") + extension
                    logger.warning("%s already exists", myfile)
                    logger.warning("Saving to: %s", newfilename)
                    outfiles.append(newfilename)
                elif mode == "overwrite":
                    logger.warning("%s already exists and will be overwritten!", myfile)
                    os.remove(myfile)
                    outfiles.append(myfile)	
                elif mode == "stop":
                    logger.critical("%s already exists! Exiting now!", myfile)
                    sys.exit()
            else:
                outfiles.append(myfile)
        elif operation == "open":
            if not os.path.isfile(myfile):
                logger.critical("Unable to open file %s", myfile)
                sys.exit()
    
    if len(outfiles) == 1: 
        return outfiles[0]
    else:
        return outfiles

#create string with current time (short) or date+time(long)
def now(sep=":", string_format="short"):
    now = datetime.now()
    if string_format == "short": 
        current_time = now.strftime("%H{sep}%M{sep}%S".format(sep=sep))
    elif string_format == "long":
        current_time = now.strftime("%Y%m%d_%H{sep}%M{sep}%S".format(sep=sep)) 
    return current_time

#arguments
parser = argparse.ArgumentParser(description='Script to add annotSV annotations to VCF file')
parser.add_argument("-i", "--inputvcf", help="Input VCF file to be annotated", action="store", required=True)
parser.add_argument("-o", "--out", help="Output VCF file", action="store", required=True)
parser.add_argument("-t", "--tmpdir", help="Folder to store temp files", action="store", required=False)
parser.add_argument("-b", "--build", help="Genome build", action="store", choices=["GRCh37", "GRCh38"], required=True)
parser.add_argument("-c", "--config", help="Config file (json)", action="store", required=False, default="SV_annotation.json")
parser.add_argument("-k", "--keeptmp", help="Set to keep tmp files", action="store_true", required=False)
parser.add_argument("--logfile", help="Log file", action="store", required=False)
args = parser.parse_args()

### LOGGER ###
if args.logfile:
    logfile = args.logfile
else:
    logfile = args.out + "_SVANNOT_" + now("","long") + ".log"

#config log for file
logging.basicConfig(filename=logfile, filemode='w', level=logging.INFO, \
    format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
logger=logging.getLogger('NCANNO')

#config log for stdout
console = logging.StreamHandler()
console.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s: %(levelname)-8s %(message)s', datefmt = '%y-%m-%d %H:%M:%S')
console.setFormatter(formatter)
logger.addHandler(console)

### Initial configuration ###
logger.info("### SV annotation started ###")
start_time = time.time()

input_file = args.inputvcf
build = args.build
json_file = args.config
if args.tmpdir:
    tmp_folder = args.tmpdir
else:
    tmp_folder = "tmp" + now("","long")
os.makedirs(tmp_folder, exist_ok=True)
out_filename = checkfile(args.out,"write","rename",logger)

logger.info("Input file: %s", input_file)
logger.info("Output file: %s", out_filename)
logger.info("Config file: %s", json_file)
logger.info("Genome build: %s", build)

checkfile([input_file, json_file], logger=logger)

with open(json_file) as json_config:
    config = json.load(json_config)
res_dir = config['RES_DIR']
bedtools = config['BEDTOOLS']

if not (os.access(bedtools, os.X_OK) or which(bedtools)):
    logger.critical("Unable to execute bedtools from %s", bedtools)
    sys.exit()
bedtools = [bedtools,'intersect','-wa','-wb']

tmp_files = {}
tmp_files['INS'] = tmp_folder + "/SVannot.tmp_INS.bed"
tmp_files['DEL'] = tmp_folder + "/SVannot.tmp_DEL.bed"
tmp_files['DUP'] = tmp_folder + "/SVannot.tmp_DUP.bed"
tmp_files['INV'] = tmp_folder + "/SVannot.tmp_INV.bed"

#Read input VCF and create temp bed files
logger.info("### READING INPUT VCF ###")
vcf = VCF(input_file)
skipped = 0
out = {}
counts = {'INS': 0, 'DEL': 0, 'DUP': 0, 'INV': 0}

out['INS'] = open(tmp_files['INS'], "w+")
out['DEL'] = open(tmp_files['DEL'], "w+")
out['DUP'] = open(tmp_files['DUP'], "w+")
out['INV'] = open(tmp_files['INV'], "w+")
tot_vars = 0
for v in vcf:
    tot_vars = tot_vars + 1
    try:
        out[v.INFO.get(config['SVTYPE'])].write("{chrom}\t{start}\t{stop}\t{ID}\n".format(chrom=v.CHROM,start=v.POS,stop=v.INFO.get(config['END']),ID=v.ID))
        counts[v.INFO.get(config['SVTYPE'])] += 1
    except:
        skipped += 1
out['INS'].close()
out['DEL'].close()
out['DUP'].close()
out['INV'].close()
logger.info("Variants loaded from input file")
for key, value in counts.items(): 
    logger.info("%s: %d", key, value)
logger.info("Skipped vars: %d", skipped)

#Run bedtools intersect on the single datasets
logger.info("### ANNOTATE VARIANTS ###")
annotations = []
tags_dataTypes = {}

logger.info("1. Processing AF datasets")
logger.info("Overlap: %s", config['overlap']['AF_datasets'][0])
logger.info("Reciprocal: %s", config['overlap']['AF_datasets'][1])
overlap = ['-f', config['overlap']['AF_datasets'][0]]
reciprocal = reciprocalOpt(config['overlap']['AF_datasets'][1]) 

for dataset, items in config['AF_datasets'][build].items():
    logger.info("Annotation for dataset %s", dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config,res_dir=res_dir)
        a_file = tmp_files[group]
        logger.info("%s: %s", group, b_file)
        checkfile([b_file], logger=logger)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

logger.info("2. Processing custom datasets")
logger.info("Overlap: %s", config['overlap']['custom_datasets'][0])
logger.info("Reciprocal: %s", config['overlap']['custom_datasets'][1])
overlap = ['-F', config['overlap']['custom_datasets'][0]]
reciprocal = reciprocalOpt(config['overlap']['custom_datasets'][1]) 

for dataset, items in config['custom_datasets'][build].items():
    logger.info("Annotation for dataset %s", dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config, res_dir=res_dir)
        a_file = tmp_files[group]
        logger.info("%s: %s", group, b_file)
        checkfile([b_file], logger=logger)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

#full_annots = pd.concat(annotations, axis=0, join='outer', sort=False)
#full_annots.fillna(0, inplace=True)

#col_operations = {}
#for col in full_annots.columns:
#    col_operations[col] = COL_OPERATIONS[tags_dataTypes[col]]

#full_annots = full_annots.groupby('ID').agg(col_operations)
#full_annots.index = full_annots.index.map(str)

logger.info("3. Processing genes datasets")
logger.info("Overlap: %s", config['overlap']['genes'][0])
logger.info("Reciprocal: %s", config['overlap']['genes'][1])
overlap = ['-F', config['overlap']['genes'][0]]
reciprocal = reciprocalOpt(config['overlap']['genes'][1]) 

for dataset, items in config['genes'][build].items():
    logger.info("Annotation for dataset %s", dataset)
    for file_config in items:
        group,b_file,cols_idx,fields,data_types = readFileConfig(file_config, res_dir=res_dir)
        a_file = tmp_files[group]
        logger.info("%s: %s", group, b_file)
        checkfile([b_file], logger=logger)

        annotations.append(computeOverlap(dataset,group,a_file,b_file,overlap,reciprocal,cols_idx,fields,tmp_folder))
        tags_dataTypes.update(setDataType(dataset,fields,data_types))

logger.info("### MERGING ANNOTATIONS ###")
full_annots = pd.concat(annotations, axis=0, join='outer', sort=False)
full_annots.fillna(0, inplace=True)

col_operations = {}
for col in full_annots.columns:
    col_operations[col] = COL_OPERATIONS[tags_dataTypes[col]]

full_annots = full_annots.groupby('ID').agg(col_operations)
full_annots.index = full_annots.index.map(str)

#Add the annotations and write new VCF
step_time = time.time()
start_annot_time = time.time()
logger.info("### ADD ANNOTATIONS TO OUTPUT VCF ###")
logger.info("Output file: %s", out_filename)
outfile = open(out_filename, "w+")
if args.inputvcf.endswith("vcf.gz"):
    vcf = gzip.open(args.inputvcf,"rt")
elif args.inputvcf.endswith("vcf"):
    vcf = open(args.inputvcf,"r")

line = vcf.readline()
while line.startswith('##'):
    outfile.write(line)
    line = vcf.readline()

#Add INFO tags descriptions to header
for tag, data_type in tags_dataTypes.items():
    header_line = '##INFO=<ID={tag},Number=1,Type={datatype},Description="SV annotation">'.format(tag=tag,datatype=data_type)
    outfile.write(header_line + "\n")

while line.startswith('#'):
    outfile.write(line)        
    line = vcf.readline()

#Annotate vars
n = 0
annotated = 0
while line:
    n += 1
    is_annotated = 0
    line = tokenize(line,"\t")
    ID = line[2]
    infos = INFO(line[7])
    for col in full_annots.columns:
        try:
            newannot = full_annots.loc[ID,col]
            if isinstance(newannot, (float, int)):
                newannot = str(round(newannot,4))
                infos.addINFO(col, newannot)
                is_annotated = 1
            else:
                newannot = [str(x) for x in newannot if x != 0]
                if len(newannot) > 0:
                    newannot = [x.split(",") for x in newannot]
                    newannot = set(flat(newannot))
                    newannot = ",".join(newannot)
                    infos.addINFO(col, newannot)
                    is_annotated = 1
            newline = line[0:7]
            newline.append(infos.getINFO("string"))
            newline.extend(line[8:])   
        except:
            newline = line
    if is_annotated == 1: 
        annotated += 1

    outfile.write("\t".join(newline) + "\n")
    
    #Counter and logging
    if n % 5000 == 0:
        elapsed_time = time.time() - step_time
        step_time = time.time()
        batch_time = time.time() - start_annot_time
        total_time = round(int(tot_vars) / (n/int(batch_time)))
        eta = timedelta(seconds=(total_time - elapsed_time))
        logger.info("Last batch %s sec: %d variants processed - %s%% - ETA %s", round(elapsed_time,3), n, round((n/tot_vars) * 100,2), eta)

    line = vcf.readline()

outfile.close()
logger.info("### ANNOTATION FINISHED ###")

if not args.keeptmp:
    logger.info("Removing temp files")
    rmtree(tmp_folder)
    #tmp_files = os.listdir(tmp_folder)
    #for myfile in tmp_files:
    #    if myfile.startswith('SVannot.tmp_'):
    #        os.remove(tmp_folder + "/" + myfile)

end_time = time.time() - start_time
end_time = str(timedelta(seconds=round(end_time)))
logger.info("Annotation completed in %s", end_time)
logger.info("Variants annotated: %s",annotated)