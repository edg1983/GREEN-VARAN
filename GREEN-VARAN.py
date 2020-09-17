'''
GREEN-VARAN
Genome Regulatory Elements ENcyclopedia VARiant ANnotation
Annotate non-coding variants with information from GREEN-DB
Regulatory variant impact predictions from LINSight, FIRE, ExPECTO, NCBoost, ReMM, CADD, DANN
Can add var AF from gnomAD and conservation from PhyloP100
The tool runs vcfanno under the cover for fast annotation

Author: Edoardo Giacopuzzi
'''

import argparse, os, re, sys, subprocess, time
from shutil import which
from collections import OrderedDict
from datetime import datetime,timedelta
import logging
from Configuration import *
from cyvcf2 import VCF, Writer

VERSION=1.0
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

def checkfolder(folders, logger=None):
    if isinstance(folders, str): folders = [folders]
    
    for folder in folders:
        if not os.path.isdir(folder):
            logger.critical("{folder} folder not found!".format(folder=folder))
            sys.exit()
    return True

#write to standard file or bgzipped
def writeFile(line,out_file,gzout):
    if gzout:
        print("TODO - Implement gz")
    else:
        out_file.write(line)

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

#Read genes list, either comma-separated or file with list
def readGenes(genes):
    if os.path.isfile(genes):
        gene_list = open(genes).read().splitlines()
    else:
        gene_list = genes.split(",")
    return(gene_list)

#store, update, read from INFO fields
class INFO():
    def __init__(self,info_field):
        exp = re.compile('(.+)=(.*)')
        self.infos = OrderedDict()
        info_tags = info_field.split(";")
        for value in info_tags:
            m = exp.match(value)
            if m:
                self.infos[m.group(1)] = m.group(2)
            else:
                self.infos[value] = 1
    
    #Remove a list of tags from INFO
    def dropINFO(self, info_keys):
        for key in info_keys: 
            try:
                self.infos.pop(key)
            except:
                pass

    #add a tag value pair to INFO
    def addINFO(self, info_dict):
        for key, value in info_dict.items():
            self.infos[key] = value
    
    #Return content of INFO as dict of string
    def buildINFO(self):
        output = [key + "=" + str(value) for key, value in self.infos.items()] 
        return ';'.join(output)

#Manage NC_ANNO field configuration. 
#Dictionary of values extracted from VCF line is passed when initialize
class NC_ANNO():
    def __init__(self, info_dict, keys):
        self.values = {}
        for key in keys:
            tmpkey = "TMP" + key
            try:
                self.values[key] = info_dict[tmpkey]
            except:
                pass

    #Manually add a specific value to the NC_ANNO field
    def updateValue(self,tag,value):
        self.values[tag] = value

    #buildAnno actually create the string for NC_ANNO field
    def buildAnno(self,order,individual_fields=False,sep='|'):
        if individual_fields:
            out_field = OrderedDict((key, self.values[key]) for key in order if key in self.values.keys())
            return out_field
        else:
            out_field = [self.values.get(key,'') for key in order]
            return {'NC_ANNO' : sep.join(str(f) for f in out_field)}

#Set variant class level for prioritization
def prioritize(parameters):
    if (parameters['remm'] >= 0.963 or \
        parameters['fathmm_mkl'] >= 0.908 or \
        parameters['ncboost'] >= 0.238) and \
        parameters['is_RegDB'] == 1 and \
        (parameters['is_TFBS'] == 1 or \
        parameters['is_DNase'] == 1 or \
        parameters['is_dbSuper'] == 1 or \
        parameters['is_UCNE'] == 1) and \
        parameters['constraint_max'] >= 0.9:
        level = 13
    elif (parameters['remm'] >= 0.963 or \
        parameters['fathmm_mkl'] >= 0.908 or \
        parameters['ncboost'] >= 0.238) and \
        parameters['is_RegDB'] == 1 and \
        parameters['constraint_max'] >= 0.9:
        level = 12
    elif (parameters['remm'] >= 0.963 or \
        parameters['fathmm_mkl'] >= 0.908 or \
        parameters['ncboost'] >= 0.238) and \
        parameters['is_RegDB'] == 1:
        level = 11
    elif (parameters['remm'] >= 0.203 or \
        parameters['fathmm_mkl'] >= 0.058 or \
        parameters['ncboost'] >= 0.062) and \
        parameters['is_RegDB'] == 1 and \
        (parameters['is_TFBS'] == 1 or \
        parameters['is_DNase'] == 1 or \
        parameters['is_dbSuper'] == 1 or \
        parameters['is_UCNE'] == 1) and \
        parameters['constraint_max'] >= 0.9:
        level = 10
    elif (parameters['remm'] >= 0.203 or \
        parameters['fathmm_mkl'] >= 0.058 or \
        parameters['ncboost'] >= 0.062) and \
        parameters['is_RegDB'] == 1 and \
        parameters['constraint_max'] >= 0.9:
        level = 9
    elif (parameters['remm'] >= 0.203 or \
        parameters['fathmm_mkl'] >= 0.058 or \
        parameters['ncboost'] >= 0.062) and \
        parameters['is_RegDB'] == 1:
        level = 8
    elif (parameters['remm'] >= 0.963 or \
        parameters['fathmm_mkl'] >= 0.908 or \
        parameters['ncboost'] >= 0.238) and \
        (parameters['is_TFBS'] == 1 or \
        parameters['is_DNase'] == 1 or \
        parameters['is_dbSuper'] == 1 or \
        parameters['is_UCNE'] == 1):
        level = 7
    elif (parameters['remm'] >= 0.203 or \
        parameters['fathmm_mkl'] >= 0.058 or \
        parameters['ncboost'] >= 0.062) and \
        (parameters['is_TFBS'] == 1 or \
        parameters['is_DNase'] == 1 or \
        parameters['is_dbSuper'] == 1 or \
        parameters['is_UCNE'] == 1):
        level = 6
    elif parameters['remm'] < 0.203 and \
        parameters['fathmm_mkl'] < 0.058 and \
        parameters['ncboost'] < 0.062 and \
        parameters['is_RegDB'] == 1 and \
        (parameters['is_TFBS'] == 1 or \
        parameters['is_DNase'] == 1 or \
        parameters['is_dbSuper'] == 1 or \
        parameters['is_UCNE'] == 1) and \
        parameters['constraint_max'] >= 0.9:
        level = 5
    elif parameters['remm'] < 0.203 and \
        parameters['fathmm_mkl'] < 0.058 and \
        parameters['ncboost'] < 0.062 and \
        parameters['is_RegDB'] == 1 and \
        parameters['constraint_max'] >= 0.9:
        level = 4
    elif parameters['remm'] < 0.203 and \
        parameters['fathmm_mkl'] < 0.058 and \
        parameters['ncboost'] < 0.062 and \
        parameters['is_RegDB'] == 1:
        level = 3
    elif (parameters['remm'] >= 0.963 or \
        parameters['fathmm_mkl'] >= 0.908 or \
        parameters['ncboost'] >= 0.238):
        level = 2
    elif (parameters['remm'] >= 0.203 or \
        parameters['fathmm_mkl'] >= 0.058 or \
        parameters['ncboost'] >= 0.062):
        level = 1
    else:
        level = 0
    
    return str(level)

###############
## Arguments ##
###############
parser = argparse.ArgumentParser(description='Script to annotate regulatory informations')
parser.add_argument("-i", "--vcf", help="Input vcf[.gz] file", action="store", required=True)
parser.add_argument("-o", "--output", help="VCF output file", action="store", required=True)
parser.add_argument("-b", "--build", help="Genome build", action="store", choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("-m", "--mode", help="Set running mode", action="store", choices=['annotate','filter_regdb','filter_any'], required=True)
parser.add_argument("-p", "--prioritize", help="Turn on prioritization for non-coding variants", action="store_true", required=False)

parser.add_argument("--vcfanno", help="full path to vcfanno executable", action="store", required=False, default=VCF_ANNO.format(base_dir=BASE_DIR))
parser.add_argument("--bed_dir", help="Directory containing RegDB bed files", action="store", required=False, default=BED_DIR.format(base_dir=BASE_DIR))
parser.add_argument("--AF_file", help="gnomAD VCF", action="store", required=False)
parser.add_argument("--scores_dir", help="Directory containing prediction scores tables", action="store", required=False, default=SCORES_DIR.format(base_dir=BASE_DIR))

parser.add_argument("-g", "--gene_mode", help="Activate gene based annotation/filter", action="store", required=False, choices=['annotate','filter'])
parser.add_argument("-t", "--gene_type", help="Which genes to consider for NC regions", action="store", required=False, default='controlled', choices=['controlled','closest','both'])
parser.add_argument("-l", "--genes_list", help="List of genes of interest, can be comma-separated list or file with one gene per line", action="store", required=False)
parser.add_argument("--impact", help="Only report NC vars if the controlled at least this impact", action="store", required=False, choices=['HIGH', 'MODERATE', 'LOW', 'MODIFIER'])

parser.add_argument("--allelefreq", help="Add gnomAD AF annotations based on global AF or specific pop", action="store", choices=['global','afr','amr','eas','fin','nfe','sas','oth'], required=False)
parser.add_argument("-s", "--scores", help="Add selected prediction score for non-coding vars. Repeat to add multiple scores", choices=SCORES, action="append", required=False)
parser.add_argument("-a", "--allscores", help="Add all prediction scores for non-coding vars (ReMM,FIRE,LinSight,ExPECTO,NCBoost)", action="store_true", required=False)

parser.add_argument("--separate_fields", help="Make multiple fields instead of a single NC_ANNO annotation", action="store_true", required=False)
parser.add_argument("--logfile", help="Log file", action="store", required=False)
parser.add_argument("--threads", help="Number of threads for annotation", action="store", required=False, default="4")
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
                  
print('\033[32m' + '''
 / __)(  _ \(  __)(  __)(  ( \ ___ / )( \ / _\ (  _ \ / _\ (  ( \\
( (_ \ )   / ) _)  ) _) /    /(___)\ \/ //    \ )   //    \/    /
 \___/(__\_)(____)(____)\_)__)      \__/ \_/\_/(__\_)\_/\_/\_)__)                    
                                                             
                                              _.-~`  `~-.
                  _.--~~~---,.__          _.,;; .   -=(@'`\\
               .-`              ``~~~~--~~` ';;;       ____)
            _.'            '.              ';;;;;    '`_.'
         .-~;`               `\           ' ';;;;;__.~`
       .' .'          `'.     |           /  /;''
        \/      .---'' ``)   /'-._____.--'\  \\
       _/|    (`        /  /`              `\ \__
',    `/- \   \      __/  (_                /-\-\-`
  `;'-..___)   |     `/-\-\-`
    `-.       .'
 jgs   `~~~~``
''', '\033[m')

logger.info("###  GREEN-VARAN version %s ###", str(VERSION))
logger.info("# EFFECTIVE CONFIGURATION #")
logger.info("Python version: %s.%s", str(sys.version_info[0]), str(sys.version_info[1]))
if sys.version_info[0] < 3:
    logger.critical("Python version 3 is required")

######################################
## Initial configuration and checks ##
######################################

#Set file paths and scores
vcf_anno = args.vcfanno
threads = args.threads
scores_dir = args.scores_dir
bed_dir = args.bed_dir
if args.AF_file:
    AF_file = args.AF_file
else:
    AF_file = AF_FILE.format(base_dir=BASE_DIR, build=args.build)
if args.allscores:
    scores = set(SCORES)
elif args.scores:
    scores = set(args.scores)
else:
    scores = set()

#Detected configuration
start_time = time.time()
step_time = time.time()
logger.info("threads: %s", threads)
logger.info("Genome build: %s", args.build)
logger.info("vcfanno path: %s", vcf_anno)
logger.info("Scores folder: %s", scores_dir)
logger.info("BED files folder: %s", bed_dir)
logger.info("AF file: %s", AF_file)

#Check files and folders 
checkfolder([scores_dir, bed_dir], logger=logger)
checkfile([args.vcf, AF_file, vcf_anno], logger=logger)
if args.overwrite == True: 
    mode = "overwrite"
else:
    mode = "rename"
if args.output.endswith(".vcf") or args.output.endswith(".vcf.gz"):
    out_VCFfile = checkfile(args.output,"write",mode,logger) 
else:
    logger.critical("Output file extension should be .vcf/.vcf.gz")
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

if args.prioritize:
    logger.info("Variant prioritization active")
    logger.info("Variant class will be addedd to as NC_VARCLASS") 
    if not all(elem in scores for elem in ["ReMM","FATHMM_MKL","NCBoost"]):
        logger.warning("ReMM,FATHMM_MKL,NCBoost required for prioritization. They have been addedd to score annotations")
    scores.update(["ReMM","FATHMM_MKL","NCBoost"])

if len(scores) > 0: 
    for score in scores:
        anno_configuration = anno_configuration + '\n' + TOML_SCORES[score].format(build=args.build,scores_dir=scores_dir)
    logger.info("Score annotation active") 
    logger.info("Following scores will be added to VCF: %s", ','.join(scores))

if args.allelefreq:
    if args.allelefreq == "global":
        pop_field = ""
        pop_tag = ""
    else:
        pop_field = ',"AF_"' + args.allelefreq + '"'
        pop_tag = ',"gnomAD_AF_"' + args.allelefreq + '"'
    anno_configuration = anno_configuration + '\n' + TOML_AF.format(AF_file=AF_file, pop_field=pop_field, pop_tag=pop_tag)
    logger.info("gnomAD allele frequency annotation active")
    logger.info("Following fields will be added to VCF: gnomAD_AF %s", pop_tag)

if args.gene_mode:
    logger.info("GENE option is active in %s mode", args.gene_mode)
    if args.genes_list:
        genes_of_interest = readGenes(args.genes_list)
        logger.info("%d genes loaded from your input list", len(genes_of_interest))
out.write(anno_configuration + '\n')
out.close()

#Get total number of vars in VCF
notused_returnCode, out, notused_err = run([BCFTOOLS,'index','--nrecords',args.vcf])
tot_vars = int(out.decode("utf-8"))
logger.info("\n###############################################################")
logger.info("Input VCF: %s", args.vcf)
logger.info("%d variants in VCF", tot_vars)
logger.info("Output VCF: %s", out_VCFfile)

gene_impacts = { 'HIGH' : set(), 'MODERATE' : set(), 'LOW' : set(), 'MODIFIER' : set()}
if args.impact:
    logger.info("gene impact option is ON: scanning variants impact from ANN field")
    for variant in VCF(args.vcf):
        ANN_field = variant.INFO.get('ANN', False)
        if ANN_field:
            consequences = ANN_field.split(",")
            for c in consequences: 
                snpeff_field = c.split('|')
                gene_impacts[snpeff_field[2]].add(snpeff_field[3])

#########################
## Perform annotations ##
#########################
logger.info("Variants processing started...")
if out_VCFfile.endswith(".vcf"):
    outVCF = open(out_VCFfile, "w+")
    gzout = False
elif out_VCFfile.endswith(".vcf.gz"):
    print("gz not supported")
    sys.exit()
    #print("vcg.gz file detected")
    #outfile = open(out_VCFfile, "wb")
    #outVCF = BGZipWriter(outfile)
    #gzout = True
n = 0
interval = 5000
written_vars_count = 0

#Grab annotated lines streaming vcfanno output
os.environ['GOGC'] = "5000"
os.environ['IRELATE_MAX_CHUNK'] = "12000"
os.environ['IRELATE_MAX_GAP'] = "1000"

for line in get_stdout([vcf_anno,"-p",threads,toml_file,args.vcf]):
    write_this_var = True

    #Update the header (remove temporary tags and add the NC_ANNO tag(s))
    if line.startswith('##'):
        if not line.startswith('##INFO=<ID=TMPNC_'):
            #writeFile(line,outVCF,gzout)
            outVCF.write(line)
    elif line.startswith('#CHROM'):
        if args.separate_fields:
            header_line = '\n'.join(HEADER_SEPARATED)
        else:
            header_line = HEADER_SINGLE
        
        outVCF.write(header_line + '\n')
        if args.gene_mode: outVCF.write(HEADER_VOI + '\n')
        if args.prioritize: outVCF.write(HEADER_CLASS + '\n')
        outVCF.write(line)

    #Processing annotations
    else:
        n += 1
        line = line.rstrip('\n')
        tokens = line.split('\t')
        info = INFO(tokens[7])
        #print("#INFO keys start:", ",".join(info.infos.keys()))
        #A region is annotated only if it is in the region table
        if "TMPNC_regionID" in info.infos.keys():
            nc_anno = NC_ANNO(info.infos, NC_ANNO_TAGS)
            
            #Read and convert to unique list region type, TFBS, reg_genes
            region_type = set(nc_anno.values['NC_region_type'].split(','))
            reg_genes = nc_anno.values.get('NC_genes','')
            reg_genes = set([x for x in reg_genes.split(',') if x != ''])
            TFnames = info.infos.get('NC_TFname','')
            TFnames = set([x for x in TFnames.split(',') if x != ''])

            nc_anno.updateValue('NC_region_type', ",".join(region_type))
            if len(reg_genes) > 0:
                nc_anno.updateValue('NC_genes', ",".join(reg_genes))
            if len(TFnames) > 0:
                info.addINFO({'NC_TFname': ",".join(TFnames)})

            #Set variables for global score and priotitize            
            phyloP = 0
            try:
                if float(nc_anno.values['NC_median_PhyloP100']) > 0: phyloP = nc_anno.values['NC_median_PhyloP100']
            except:
                phyloP = 0

            constraint_values = nc_anno.values['NC_constraint'].split(",")
            constraint_values = [0] + [float(i) for i in constraint_values if i != "NA"]
            constraint_max = max(constraint_values)
            
            parameters = {
                'is_RegDB' : 1,
                'methods' : int(nc_anno.values['NC_methods']),
                'has_gene' : int(nc_anno.values.get('NC_genes',None) is not None),
                'is_TFBS' : int(info.infos.get('NC_TFname', None) is not None),
                'is_DNase' : int(info.infos.get('NC_DNase', None) is not None),
                'is_UCNE' : int(info.infos.get('NC_UCNE', None) is not None),
                'is_dbSuper' : int(info.infos.get('NC_dbSUPER', None) is not None),
                'phyloP' : phyloP,
                'constraint_max' : constraint_max
            }

            #Calculate global regulatory support based on RegDB regions
            global_score = sum([float(parameters[x]) for x in [
                'methods',
                'constraint_max', 
                'is_TFBS', 
                'is_DNase',
                'is_UCNE',
                'is_dbSuper',
                'has_gene',
                'phyloP']] )
            nc_anno.updateValue('NC_support', round(global_score, 2))
            
            #Drop temporary tags and add NC_ANNO formatted
            info.dropINFO(["TMP"+x for x in NC_ANNO_TAGS])
            nc_info_field = nc_anno.buildAnno(NC_ANNO_TAGS,args.separate_fields)

            info.addINFO(nc_info_field)

            #Make a the list of genes of interest if gene option is active
            if args.gene_mode:
                if nc_anno.values['NC_genes']:
                    controlled_genes = nc_anno.values['NC_genes'].split(',')
                else:
                    controlled_genes = []
                closest_genes = nc_anno.values['NC_closestGene'].split(',')
                if args.gene_type == 'controlled':
                    genes = set(controlled_genes)
                elif args.gene_type == 'closest':
                    genes = set(closest_genes)
                else:
                    genes = set(controlled_genes + closest_genes)

                if args.genes_list:
                    genes = genes.intersection(genes_of_interest)

                if args.impact:
                    genes = genes.intersection(gene_impacts[args.impact])

                #Check if genes of interest are among genes associated to the NC region  
                if len(genes) == 0:
                    if args.gene_mode == "filter": write_this_var = False
                    info.addINFO({'NC_VOI' : 0})
                else:
                    if args.gene_mode == "annotate": info.addINFO({'NC_VOI' : 1})

        else:
            #No RegDB fields are added if no annotations
            #but TFBS, DNase, UCNE, dbSuper remain
            info.dropINFO(["TMP"+x for x in NC_ANNO_TAGS])
            
            #Set variable for prioritize
            parameters = {
                'is_RegDB' : 0,
                'has_gene' : 0,
                'is_TFBS' : int(info.infos.get('NC_TFname', None) is not None),
                'is_DNase' : int(info.infos.get('NC_DNase', None) is not None),
                'is_UCNE' : int(info.infos.get('NC_UCNE', None) is not None),
                'is_dbSuper' : int(info.infos.get('NC_dbSUPER', None) is not None),
                'constraint_max' : 0
            }

            #Var is not written if gene_mode filter is active, since we have no gene info
            if args.gene_mode == "filter": write_this_var = False

            #Set var to output based on running mode
            if args.mode == "filter_regdb": write_this_var = False
            if args.mode == "filter_any":
                ann_count = [
                    int(info.infos.get('NC_TFname', 0) != 0), 
                    int(info.infos.get('NC_DNase', 0) != 0),
                    int(info.infos.get('NC_UCNE', 0) != 0),
                    int(info.infos.get('NC_dbSUPER', 0) != 0)
                ]
                if sum(ann_count) == 0: write_this_var = False  

        if args.prioritize:
            parameters['remm'] = float(info.infos.get('NC_ReMM', 0))
            parameters['fathmm_mkl'] = float(info.infos.get('NC_FATHMM_MKL', 0))
            parameters['ncboost'] = float(info.infos.get('NC_NCBoost', 0))
            var_class = prioritize(parameters)
            info.addINFO({'NC_VARCLASS' : var_class})

        #Update INFO with new annotations and write vars
        tokens[7] = info.buildINFO()
        #try:
        #    print("#NCANNO:",nc_anno.values.get("NC_genes","NO GENES"))
        #    print("#INFO:",info.infos.get("TMPNC_genes","NO GENES"))
        #    print("#INFO keys:", ",".join(info.infos.keys()))
        #except:
        #    pass
        #print(tokens[1],tokens[2],tokens[4],tokens[5],tokens[7])
        if write_this_var == True:
            written_vars_count += 1
            outVCF.write('\t'.join(tokens) + "\n")

        #Counter and logging
        if n % interval == 0:
            interval = INTERVALS.get(str(interval), 5000)
            elapsed_time = time.time() - step_time
            step_time = time.time()
            batch_time = time.time() - start_time
            total_time = round(int(tot_vars) / (n/int(batch_time)))
            eta = timedelta(seconds=(total_time - batch_time))
            logger.info("Last batch %s sec: %d variants processed - %s%% - ETA %s", round(elapsed_time,3), n, round((n/tot_vars) * 100,2), eta)

outVCF.close()
#os.remove(toml_file)

end_time = time.time() - start_time
end_time = str(timedelta(seconds=round(end_time)))
logger.info("Annotation completed: %d variants written in %s", written_vars_count, end_time)