'''
Regulatory regions annotator
Annotate non-coding regions from the RegulatoryRegions.db
Regulatory variant impact predictions from LINSight, FIRE, ExPECTO, NCBoost, ReMM, CADD, DANN
Can add var AF from gnomAD and conservation from PhyloP100
The tool runs vcfanno under the cover for fast annotation
Author: Edoardo Giacopuzzi
Nov 2019
'''

import argparse, os, re, sys, subprocess, time
from shutil import which
from collections import OrderedDict
from datetime import datetime,timedelta
import sqlite3
import logging
from Configuration import *

VERSION=0.1
BASE_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))
INTERVALS = {
    '20000' : 10000,
    '50000' : 25000,
    '150000' : 50000,
    '500000' : 100000,
    '1000000' : 500000
}


#check file exists and trigger proper actions
def checkfile(files, operation="open", mode="stop", logger=None):
    if isinstance(files, str): files = [files]
    
    outfiles = []
    for myfile in files:
        if operation == "write":
            if os.path.isfile(myfile):
                if mode == "rename":
                    filebase, extension = os.path.splitext(myfile)
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

def checkfolder(folders, logger=None):
    if isinstance(folders, str): folders = [folders]
    
    for folder in folders:
        if not os.path.isdir(folder):
            logger.critical("{folder} folder not found!".format(folder=folder))
            sys.exit()
    return True

#create string with current time (short) or date+time(long)
def now(sep=":", string_format="short"):
    now = datetime.now()
    if string_format == "short": 
        current_time = now.strftime("%H{sep}%M{sep}%S".format(sep=sep))
    elif string_format == "long":
        current_time = now.strftime("%Y%m%d_%H{sep}%M{sep}%S".format(sep=sep)) 
    return current_time

#run external command and return exit code and stderr
def run(cmd, shell=False):
    proc = subprocess.Popen(cmd, shell=shell,    
        stdout = subprocess.PIPE,
        stderr = subprocess.PIPE   
    )
    stdout, stderr = proc.communicate()
 
    return proc.returncode, stdout, stderr

#run external command and yield lines from stdout
def get_stdout(command, shell=False):
    program = subprocess.Popen(command, shell=shell, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    for line in iter(program.stdout.readline, b''):
        yield line.decode('utf-8')

#create a database connection to the SQLite database
def create_connection(db_file):
    conn = None
    try:
        conn = sqlite3.Connection(db_file)
        cur = conn.cursor()
        cur.execute('PRAGMA cache_size = 1000000')
        cur.execute('PRAGMA locking_mode = EXCLUSIVE')
        cur.execute('PRAGMA temp_store = 2')
        cur.execute('PRAGMA synchronous = 0')
        cur.execute('PRAGMA journal_mode = OFF')
        cur.execute('PRAGMA threads = 10')
    except:
        conn = None
 
    return conn
 
#prepare query for methods or genes tables 
def prepare_query(table, ids):
    if table == "genes":
        query = "SELECT GROUP_CONCAT(DISTINCT gene_symbol) FROM genes WHERE regionID IN ({regions})".format(regions=','.join(['?']*len(ids)))
    elif table =="methods":
        query = "SELECT COUNT(DISTINCT method) FROM methods WHERE regionID IN ({regions})".format(regions=','.join(['?']*len(ids)))
    return query

#get data from the database
def getData(conn, query, ids):
    cur = conn.cursor()
    cur.execute(query,ids)
    row = cur.fetchone()
    if row:
        return row
    else:
        return 0

#store, update, read from INFO fields
class INFO():
    def __init__(self,info_field):
        exp = re.compile('(.+)=(.+)')
        self.infos = OrderedDict()
        info_tags = info_field.split(";")
        for value in info_tags:
            m = exp.match(value)
            if m:
                self.infos[m.group(1)] = m.group(2)
            else:
                self.infos[value] = True
    
    #Remove a list of tags from INFO
    def dropINFO(self, info_keys):
        if isinstance(info_keys, list):
            for key in info_keys: 
                try:
                    self.infos.pop(key)
                except:
                    pass
        else:
            self.infos.pop(key)

    #add a tag value pair to INFO
    def addINFO(self, info_dict):
        for key, value in info_dict.items():
            self.infos[key] = value
    
    #change the value of a INFO tag
    def updateINFO(self, info_key, value):
        if info_key in self.infos.keys():
            self.infos[info_key] = value
        else:
            self.addINFO({info_key: value})
    
    #Return content of INFO as dict of string
    def buildINFO(self,format):
        if format == "dict":
            return self.infos
        elif format == "string":
            output = [] 
            for key, value in self.infos.items():
                output.append(key + "=" + str(value))
        return ';'.join(output)

#Manage NC_ANNO field configuration. 
#Dictionary of values extracted from VCF line is passed when initialize
class NC_ANNO():
    def __init__(self, values_dict):
        self.data_dict = values_dict
        self.values = {}
    
    #Configure each annotation value as needed (store/binary)
    def configureElement(self,key,tag,mode="store"):
        if mode == "store":
            if key not in self.data_dict.keys():
                self.values[tag] = ""
            else:
                self.values[tag] = self.data_dict[key]
        elif mode == "binary":
            if key not in self.data_dict.keys():
                self.values[tag] = 0
            else:
                self.values[tag] = 1 
    
    #Manually add a specific value to the NC_ANNO field
    def updateValue(self,tag,value):
        self.values[tag] = value

    #buildAnno actually create the string for NC_ANNO field
    def buildAnno(self,order,individual_fields=False,sep='|'):
        if individual_fields:
            out_field = OrderedDict()
            for key in order:
                out_key = "NC_" + key
                out_field[out_key] = self.values[key]
            return out_field
        else:
            out_field = []
            for key in order:
                out_field.append(self.values[key])
            return {'NC_ANNO' : sep.join(str(f) for f in out_field)}
        
###############
## Arguments ##
###############
parser = argparse.ArgumentParser(description='Script to annotate regulatory informations')
parser.add_argument("-i", "--vcf", help="Input vcf[.gz] file", action="store", required=True)
parser.add_argument("-o", "--output", help="VCF output file", action="store", required=True)

parser.add_argument("-b", "--build", help="Genome build", action="store", choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("--regDB", help="RegDB database file", action="store", required=False, default=DB_FILE)
parser.add_argument("--bed_dir", help="Directory containing RegDB bed files", action="store", required=False, default=BED_DIR)
parser.add_argument("--AF_file", help="gnomAD VCF", action="store", required=False, default=AF_FILE)
parser.add_argument("--scores_dir", help="Directory containing prediction scores tables", action="store", required=False, default=SCORES_DIR)

parser.add_argument("-g", "--gene", help="Gene of interest, can be specified multiple times", action="append", required=False)
parser.add_argument("-f", "--allelefreq", help="Add gnomAD AF annotations", action="store_true", required=False)
parser.add_argument("-c", "--scores", help="Add selected prediction score for non-coding vars. Repeat to add multiple scores", choices=['ReMM','FIRE','LinSight','ExPECTO','NCBoost','DANN','CADD'], action="append", required=False)
parser.add_argument("-s", "--allscores", help="Add all prediction score for non-coding vars (ReMM,FIRE,LinSight,ExPECTO,NCBoost)", action="store_true", required=False)
parser.add_argument("-e", "--separate_fields", help="Make multiple fields instead of a single NC_ANNO annotation", action="store_true", required=False)
parser.add_argument("-l", "--logfile", help="Log file", action="store", required=False)
parser.add_argument("-t", "--threads", help="Number of threads for annotation", action="store", required=False, default=4)
parser.add_argument("-w", "--overwrite", help="Set if you want to overwrite output file if already exists", action="store_true", required=False)
args = parser.parse_args()

###########################
## Logging configuration ##
###########################
if args.logfile:
    logfile = args.logfile
else:
    logfile = args.output + "_NCANNO_" + now("","long") + ".log"

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

logger.info("###  NC ANNO started ###")
logger.info("# DETECTED CONFIGURATION #")
logger.info("Python version detected: %s.%s", str(sys.version_info[0]), str(sys.version_info[1]))
if sys.version_info[0] < 3:
    logger.critical("Python version 3 is required")

######################################
## Initial configuration and checks ##
######################################

#Set file paths and scores
vcf_anno = VCF_ANNO.format(base_dir=BASE_DIR)
threads = args.threads
db_file = args.regDB.format(base_dir=BASE_DIR)
scores_dir = args.scores_dir.format(base_dir=BASE_DIR)
bed_dir = args.bed_dir.format(base_dir=BASE_DIR)
AF_file = args.AF_file.format(base_dir=BASE_DIR, build=args.build)
if args.allscores:
    scores = ['ReMM','FIRE','LinSight','ExPECTO','NCBoost','DANN','CADD']
elif args.scores:
    scores = args.scores

#Detected configuration
start_time = time.time()
step_time = time.time()
logger.info("threads: %s", threads)
logger.info("vcfanno path: %s", vcf_anno)
logger.info("Genome build: %s", args.build)
logger.info("Regulatory regions db: %s", db_file)
logger.info("Scores folder: %s", scores_dir)
logger.info("BED files folder: %s", bed_dir)
logger.info("AF file: %s", AF_file)

#Check files and folders 
checkfolder([scores_dir, bed_dir], logger=logger)
checkfile([args.vcf, db_file, AF_file, vcf_anno], logger=logger)
if args.overwrite == True: 
    mode = "overwrite"
else:
    mode = "rename"
out_VCFfile = checkfile(args.output,"write",mode,logger) 

#Open connection to database
anno_db = create_connection(db_file)
if anno_db is None:
    logger.critical("Failed to connect to Regulatory regions db")
    sys.exit()

#Check vcfanno/bcftools are executable
if not (os.access(vcf_anno, os.X_OK) or which(vcf_anno)):
    logger.critical("Unable to execute vcfanno from %s", vcf_anno)
    sys.exit()
if not (os.access(BCFTOOLS, os.X_OK) or which(BCFTOOLS)):
    logger.critical("Unable to execute bcftools from %s", BCFTOOLS)
    sys.exit()

#Generate toml file for vcfanno based on build and script folder
toml_file = args.output + '_' + now('') + '_Annotations.toml'
out = open(toml_file, 'w+')
anno_configuration = TOML.format(bed_dir=bed_dir, build=args.build)

if args.allscores or args.scores: 
    for score in scores:
        anno_configuration = anno_configuration + '\n' + TOML_SCORES[score].format(build=args.build,scores_dir=scores_dir)
    logger.info("Score annotation active") 
    logger.info("Following scores will be added to VCF: %s", ','.join(scores))

if args.allelefreq: 
    anno_configuration = anno_configuration + '\n' + TOML_AF.format(AF_file=AF_file)
    logger.info("gnomAD allele frequency annotation active")
    logger.info("Following fields will be added to VCF: gnomAD_AF, gnomAD_AF_max")

out.write(anno_configuration + '\n')
out.close()

#Get total number of vars in VCF
notused_returnCode, out, notused_err = run([BCFTOOLS,'index','--nrecords',args.vcf])
tot_vars = int(out.decode("utf-8"))
logger.info("\n###############################################################")
logger.info("Annotating %s", args.vcf)
logger.info("%d variants in VCF", tot_vars)

logger.info("Writing annotated VCF: %s", out_VCFfile)
if args.gene:
    logger.info("GENE option is active")
    logger.info("Only variants in regions associated to the following genes will be written")
    logger.info("%s", '; '.join(args.gene))

outVCF = open(out_VCFfile, "w+")

logger.info("Start annotation")
n = 0
interval = 5000
written_vars_count = 0

#########################
## Perform annotations ##
#########################

#Grab annotated lines streaming vcfanno output
os.environ['GOGC'] = "5000"
os.environ['IRELATE_MAX_CHUNK'] = "12000"
os.environ['IRELATE_MAX_GAP'] = "1000"
for line in get_stdout([vcf_anno,"-p",threads,toml_file,args.vcf]):
    
    write_this_var = True

    #Update the header (remove temporary tags and add the NC_ANNO tag(s))
    if line.startswith('##'):
        if not line.startswith('##INFO=<ID=TMPNC_'):
            outVCF.write(line)
    elif line.startswith('#CHROM'):
        if args.separate_fields:
            header_line = '\n'.join(HEADER_SEPARATED)
        else:
            header_line = HEADER_SINGLE
        
        outVCF.write(header_line + '\n')
        outVCF.write(line)
   
    #Processing annotations
    else:
        n += 1
        line = line.rstrip('\n')
        tokens = line.split('\t')
        info = INFO(tokens[7])
        nc_anno = NC_ANNO(dict((k, info.infos[k]) for k in NC_ANNO_OPS.keys() if k in info.infos)) 

        #A region is annotated only if it is in the region table
        if "TMPNC_regionID" in info.infos.keys():

            #Regions type, associated genes and methods: reg_type, reg_genes, reg_methods
            regionIDs = info.infos['TMPNC_regionID'].split(",")

            #Read and convert to unique list region type and subtype
            region_type = set(info.infos['TMPNC_type'].split(','))
            region_subtype = set(info.infos['TMPNC_subtype'].split(','))
            nc_anno.updateValue('region_type', ",".join(region_type))
            nc_anno.updateValue('region_subtype', ",".join(region_subtype))

            #Query the database to get associated genes and supporting methods
            query_genes = prepare_query('genes',regionIDs)
            query_methods = prepare_query('methods',regionIDs)
            reg_genes = getData(anno_db,query_genes,regionIDs)[0]
            if reg_genes is None: 
                nc_anno.updateValue('genes', '')
            else:
                nc_anno.updateValue('genes', reg_genes)

            reg_methods = getData(anno_db,query_methods,regionIDs)[0]
            if reg_methods is None: 
                nc_anno.updateValue('methods', 0)
            else:
                nc_anno.updateValue('methods', reg_methods)
            
            #TFBS informations: is_TFBS, ref_TFname
            try:
                nc_anno.updateValue('TFname', info.infos['TMPNC_TFname'])
                nc_anno.updateValue('is_TFBS', 1)
            except:
                nc_anno.updateValue('TFname', '')
                nc_anno.updateValue('is_TFBS', 0)
            else:
                #Set TFName to empty string if it has missing -99 value in the DB
                if info.infos['TMPNC_TFname'] == "-99": nc_anno.updateValue('TFname', '')

            #Process additional non-coding annotations
            for name, anno_type in NC_ANNO_OPS.items():
                tag = re.sub('TMPNC_','',name)
                nc_anno.configureElement(name,tag,anno_type)

            #Calculate global score
            global_score = [
                nc_anno.values['methods'], 
                nc_anno.values['is_TFBS'], 
                nc_anno.values['DNase'],
                nc_anno.values['UCNE'],
                nc_anno.values['dbSUPER'],
                #int(nc_anno.values['SegWey_top'] != ""),
                #int(nc_anno.values['ENCODE_HMM'] != ""),
                int(nc_anno.values['genes'] != "")
            ]
            nc_anno.updateValue('global_score', sum(int(i) for i in global_score))
            
            #Drop temporary tags and add NC_ANNO formatted
            info.dropINFO(TMP_ANNOTATIONS)
            nc_info_field = nc_anno.buildAnno(TAGS_ORDER,args.separate_fields)
            info.addINFO(nc_info_field)
        else:
            #No field is added if no annotations
            info.dropINFO(TMP_ANNOTATIONS)

        tokens[7] = info.buildINFO("string")
        
        if args.gene:
            if 'genes' in nc_anno.values.keys():
                if set(args.gene).isdisjoint(nc_anno.values['genes'].split(',')): write_this_var = False 
            else:
                write_this_var = False

        if write_this_var == True:
            written_vars_count += 1
            outVCF.write('\t'.join(tokens) + "\n")

        #Counter and logging
        if n % interval == 0:
            interval = INTERVALS.get(str(interval), 5000)
            elapsed_time = time.time() - step_time
            step_time = time.time()
            logger.info("Last batch %s sec: %d variants processed - %s%%", round(elapsed_time,3), n, round((n/tot_vars) * 100,2))

outVCF.close()
anno_db.close()
#os.remove(toml_file)

end_time = time.time() - start_time
end_time = str(timedelta(seconds=round(end_time)))
logger.info("Annotation of %d variants completed in %s", written_vars_count, end_time)

#TODO implement gene filter to add a tag in var INFO 
#TODO Allow tag/filtering of vars in reg regions linked to a gene with a severe var hit

#TODO Develop detail mode to query single variants in the DB for detailed annotations
#This may need to integrate more tables in the release, like SegWey detailed scores
#TODO Develop plotting for detail mode
#TODO Add genome browser integration for detail mode

#TODO Add automated genes selection: 
# custom list of genes of interest
# HPO associated genes
# GADO top genes based on HPOs
# OMIM associated genes (through HPO if HPOs are not specified)
# disease associated genes (through PanelApp)

