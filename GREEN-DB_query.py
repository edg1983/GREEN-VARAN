'''
RegDB query
Retrieve detailed information from RegDB for regions/vars of interes
Input can be
- list of region IDs
- a VCF annotated with GREEN-VARAN
- a tables containing variant ID (chr_pos_ref_alt), comma-separated regionIDs 
Author: Edoardo Giacopuzzi
'''

import argparse, os, re, sys, subprocess, time
from datetime import datetime,timedelta
import sqlite3
import logging
import pandas as pd
from cyvcf2 import VCF

VERSION=1.0
pd.set_option('mode.chained_assignment',None)
BASE_DIR = os.path.dirname(os.path.realpath(sys.argv[0]))
RESOURCE_DIR = "resources"
REG_DB="SQlite/RegulatoryRegions.db"
DB_FILE = "{base_dir}/{resource_dir}/{reg_db}".format(base_dir=BASE_DIR,resource_dir=RESOURCE_DIR,reg_db=REG_DB)
TABLES = ["regions","gene_details","pheno_details", "DNase","dbSuper","TFBS","UCNE"]
TAB_HEADERS = {
    'regions' : ['regionID','chrom','start','stop','type','std_type','DB_source','PhyloP100_median','constraint_pct','controlled_gene','closestGene_symbol', 'closestGene_ensg', 'closestGene_dist','cell_or_tissues','detection_method','phenotype'],
    'gene_details' : ['regionID','chrom','start','stop','std_type','controlled_gene','detection_method','tissue_of_interaction'],
    'pheno_details' : ['regionID','chrom','start','stop','std_type','phenotype','detection_method','DB_source'],
    'DNase' : ['regionID','DNase_chrom','DNase_start','DNase_stop','DNase_ID','DNase_cell_or_tissue'],
    'dbSuper' : ['regionID','dbSuper_chrom','dbSuper_start','dbSuper_stop','dbSuper_ID','dbSuper_cell_or_tissue'],
    'TFBS' : ['regionID','TFBS_chrom','TFBS_start','TFBS_stop','TF_name','TFBS_cell_or_tissue'],
    'UCNE' : ['regionID','UCNE_chrom','UCNE_start','UCNE_stop','UCNE_ID']
}

def exitNow(start_time):
    end_time = time.time() - start_time
    end_time = str(timedelta(seconds=round(end_time)))
    logger.info("Process completed in %s", end_time)
    sys.exit()

#create string with current time (short) or date+time(long)
def now(sep=":", string_format="short"):
    now = datetime.now()
    if string_format == "short": 
        current_time = now.strftime("%H{sep}%M{sep}%S".format(sep=sep))
    elif string_format == "long":
        current_time = now.strftime("%Y%m%d_%H{sep}%M{sep}%S".format(sep=sep)) 
    return current_time

#Read genes list, either comma-separated or file with list
def readList(ID_list):
    if os.path.isfile(ID_list):
        gene_list = open(ID_list).read().splitlines()
    else:
        gene_list = ID_list.split(",")
    return(gene_list)

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
        cur.execute('PRAGMA threads = 8')
    except:
        conn = None
 
    return conn

#prepare query for methods or genes tables 
def prepare_query(table, build, ids):
    if table in ["DNase", "TFBS"]:
        query = "SELECT r.regionID AS regionID, \
        t.chromosome AS {table}_chrom, \
        t.start AS {table}_start, \
        t.stop AS {table}_stop, \
        t.name AS {table}_name, \
        tt.cell_or_tissue AS {table}_tissue \
        FROM {build}_Regions AS r \
        INNER JOIN {build}_regionID_to_{table} AS link ON r.regionID = link.regionID \
        INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
        INNER JOIN {build}_{table}_tissues AS tt ON t.regionID = tt.regionID \
        WHERE r.regionID IN ({regions})".format(
            regions=','.join(['?']*len(ids)),
            build=build,
            table=table)
    elif table == "dbSuper":
        query = "SELECT r.regionID AS regionID, \
        t.chromosome AS {table}_chrom, \
        t.start AS {table}_start, \
        t.stop AS {table}_stop, \
        t.regionID AS {table}_name, \
        tt.cell_or_tissue AS {table}_tissue \
        FROM {build}_Regions AS r \
        INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID \
        INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
        INNER JOIN {table}_tissues AS tt ON t.regionID = tt.regionID \
        WHERE r.regionID IN ({regions})".format(
            regions=','.join(['?']*len(ids)),
            build=build,
            table=table)
    elif table == "UCNE":
        query = "SELECT r.regionID AS regionID, \
        t.chromosome AS {table}_chrom, \
        t.start AS {table}_start, \
        t.stop AS {table}_stop, \
        t.regionID AS {table}_name \
        FROM {build}_Regions AS r \
        INNER JOIN regionID_to_{table} AS link ON r.regionID = link.regionID \
        INNER JOIN {build}_{table} AS t ON link.link_ID = t.regionID \
        WHERE r.regionID IN ({regions})".format(
            regions=','.join(['?']*len(ids)),
            build=build,
            table=table)
    elif table == "regions":
        query = "SELECT r.regionID, r.chromosome, r.start, r.stop, \
        r.type, r.std_type, r.DB_source, r.PhyloP100_median, r.constrain_pct, \
        GROUP_CONCAT(DISTINCT g.gene_symbol), \
        r.closestGene_symbol, r.closestGene_ensg, r.closestGene_dist, \
        GROUP_CONCAT(DISTINCT t.cell_or_tissue), \
        GROUP_CONCAT(DISTINCT m.method), \
        GROUP_CONCAT(DISTINCT p.phenotype) \
        FROM {build}_Regions AS r \
        LEFT JOIN tissues AS t ON r.regionID = t.regionID \
        LEFT JOIN genes AS g ON r.regionID = g.regionID \
        LEFT JOIN methods AS m ON r.regionID = m.regionID \
        LEFT JOIN phenos AS p ON r.regionID = p.regionID \
        WHERE r.regionID IN ({regions}) \
        GROUP BY r.regionID".format(
            regions=','.join(['?']*len(ids)),
            build=build)
    elif table == "gene_details":
        query = "SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
        g.gene_symbol, m.method, t.cell_or_tissue \
        FROM {build}_Regions AS r \
        LEFT JOIN genes AS g ON r.regionID = g.regionID \
        LEFT JOIN tissues AS t ON g.interactionID = t.regionID \
        LEFT JOIN methods AS m ON g.interactionID = m.regionID \
        WHERE r.regionID IN ({regions})".format(
            regions=','.join(['?']*len(ids)),
            build=build)
    elif table == "pheno_details":
        query = "SELECT r.regionID, r.chromosome, r.start, r.stop, r.std_type, \
        p.phenotype, p.method, p.DB_source \
        FROM {build}_Regions AS r \
        LEFT JOIN phenos AS p ON r.regionID = p.regionID \
        WHERE r.regionID IN ({regions})".format(
            regions=','.join(['?']*len(ids)),
            build=build)
    return query

#get data from the database
def getData(conn, query, ids):
    ids = list(ids)
    cur = conn.cursor()
    cur.execute(query,ids)
    #print(query)
    #print(ids)
    results = cur.fetchall()
    if results:
        return results
    else:
        return False

#check file existance
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

#get results tables from db based on prepared query
def get_result_tables(db, build, region_ids, var_id=False, logger=None):
    result_tables = dict()
    for table in TABLES:            
        query = prepare_query(table, build, region_ids)
        records = getData(db,query,region_ids)
        if records:
            dataframe = pd.DataFrame.from_records(getData(db,query,region_ids))
            dataframe.columns = TAB_HEADERS[table]
            dataframe.iloc[:,2] = pd.to_numeric(dataframe.iloc[:,2])
            dataframe.iloc[:,3] = pd.to_numeric(dataframe.iloc[:,3])
            if var_id: dataframe['var_id'] = "VAR_ID"
            result_tables[table] = dataframe
        #else:
        #    result_tables[table] = pd.DataFrame
        if logger: logger.info("Query completed for %s table", table)
    return result_tables

#Process result tables add variant ids
def fillVarID(var_id, query_tables):
    var_tables = dict()
    for table in TABLES:
        if table in query_tables.keys():
            if not query_tables[table].empty: 
                df = query_tables[table]
                variant_fields = var_id.split("_")
                variant_chrom = variant_fields[0]
                variant_pos = int(variant_fields[1])
                df.loc[(df.iloc[:,1] == variant_chrom) & (df.iloc[:,2] <= variant_pos) & (df.iloc[:,3] >= variant_pos), 'var_id'] += "," + var_id
                var_tables[table] = df
    return var_tables

#Save tables to disk
def saveTables(tables,out_tables,logger):
    for table in tables:
        outfile = out_prefix + "." + table + ".tsv"
        if table in out_tables.keys():
            if not out_tables[table].empty:
                out_tables[table] = out_tables[table].drop_duplicates()
                out_tables[table].to_csv(outfile,sep="\t",index=False)
                logger.info("%s table saved to: %s", table, outfile)
            else:
                logger.info("No info found for %s table", table)
        else:
            logger.info("No info found for %s table", table)

###############
## Arguments ##
###############
parser = argparse.ArgumentParser(description='Script to annotate regulatory informations')
inputs = parser.add_mutually_exclusive_group(required=True)
inputs.add_argument("-v", "--vcf", help="Input vcf[.gz] file. Must be annotated with GREEN-VARAN", action="store")
inputs.add_argument("-r", "--regIDs", help="Comma separated list of region IDs or file with a list of region IDs", action="store")
inputs.add_argument("-t", "--table", help="Tab-separated file with col1: variant (chr_pos_ref_alt); col2: comma-separated list of region IDs", action="store")

parser.add_argument("-o", "--outprefix", help="Prefix for output files", action="store", required=True)

parser.add_argument("-b", "--build", help="Genome build", action="store", choices=['GRCh37', 'GRCh38'], required=True)
parser.add_argument("--regDB", help="RegDB database file", action="store", required=False, default=DB_FILE)
parser.add_argument("--logfile", help="Log file", action="store", required=False)
args = parser.parse_args()

###########################
## Logging configuration ##
###########################
if args.logfile:
    logfile = args.logfile
else:
    logfile = args.outprefix + "_" + now("","long") + ".log"

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

#####################
### INITIAL SETUP ###
#####################
start_time = time.time()
regDB = args.regDB
build = args.build
out_prefix = args.outprefix
checkfile(regDB, logger=logger)

logger.info("#####  GREEN-DB QUERY version %s  #####", str(VERSION))
logger.info("# EFFECTIVE CONFIGURATION #")
logger.info("Python version: %s.%s", str(sys.version_info[0]), str(sys.version_info[1]))
if sys.version_info[0] < 3:
    logger.critical("Python version 3 is required")

#Detected configuration
logger.info("Genome build: %s", build)
logger.info("GREEN-DB: %s", regDB)

#Open connection to database
anno_db = create_connection(regDB)
if anno_db is None:
    logger.critical("Failed to connect to Regulatory regions db")
    sys.exit()

logger.info("# GREEN-DB QUERY START #")

### Query from region list without variants ###
if args.regIDs:
    region_IDs = readList(args.regIDs)
    logger.info("Region IDs input with %d IDs", len(region_IDs))
    out_tables = get_result_tables(anno_db, build, region_IDs, logger=logger)
    saveTables(TABLES,out_tables,logger)
    exitNow(start_time)

### Load var and regID from table ###
if args.table:
    checkfile(args.table, logger=logger) 
    logger.info("Variant table input: %s", args.table)
    region_IDs = set()
    variants = set()
    with open(args.table) as f:
        line = f.readline()
        while line:
            line = line.rstrip("\n")
            line = line.split("\t")
            
            variants.add(line[0])
            regions = set(line[1].split(","))
            if len(regions) > 0:
                region_IDs.update(regions)
            line = f.readline()
    logger.info("%d variants and %d regions read from file", len(variants), len(region_IDs))

### Load var and regID from VCF ###
if args.vcf:
    checkfile(args.vcf, logger=logger) 
    logger.info("VCF input: %s", args.vcf)

    logger.info("Reading variants and region IDs")
    region_IDs = set()
    variants = set()
    for variant in VCF(args.vcf):
        var_id = "_".join([str(x) for x in [variant.CHROM, variant.start, variant.REF, variant.ALT[0]]])
        nc_regid = variant.INFO.get("NC_regionID", False)
        nc_anno = variant.INFO.get("NC_ANNO", False)
        if nc_regid or nc_anno:
            variants.add(var_id)
        if nc_regid:
            region_IDs.update(nc_regid.split(","))
        elif nc_anno:
            nc_anno = nc_anno.split('|')
            region_IDs.update(nc_anno[1].split(","))
    
    logger.info("%d variants and %d regions read from file", len(variants), len(region_IDs))

### Query DB ###
step_time = time.time()
query_tables = get_result_tables(anno_db, build, region_IDs, var_id=True, logger=logger)
elapsed_time = time.time() - step_time
logger.info("Finished extracting regions from DB in %s sec", round(elapsed_time,3))

### Annotate var_id and prepare output tables ###
logger.info("Variant annotation started")
step_time = time.time()
interval_time = time.time()
n = 0
for v in variants:
    query_tables = fillVarID(v, query_tables)
    n = n + 1
    if n % 10 == 0:
        elapsed_time = time.time() - interval_time
        interval_time = time.time()
        logger.info("Annotated 10 variants in %s sec", round(elapsed_time,3))
elapsed_time = time.time() - step_time
logger.info("Finished variant annotation in %s sec", round(elapsed_time,3))

logger.info("Clean output tables")
for table in TABLES:
    if table in query_tables.keys():
        if not query_tables[table].empty:
            df = query_tables[table]
            df = df[df.var_id != "VAR_ID"]
            df['var_id'] = df['var_id'].str.replace('VAR_ID,', '')
            query_tables[table] = df

saveTables(TABLES,query_tables,logger)

exitNow(start_time)

#TODO Develop a plotting system