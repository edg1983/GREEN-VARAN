nextflow.enable.dsl=2

// Define allowed options
def allowed_builds  = ['GRCh37', 'GRCh38']
def known_scores  = params.annotations[params.build].scores.keySet() as String[]
def allowed_scores_values = known_scores + ['best', 'all']
def known_regions  = params.annotations[params.build].regions.keySet() as String[]
def allowed_regions_values = known_regions + ['best', 'all']
def resource_folder = "${params.resource_folder}/${params.build}"
def sql_folder = "${params.resource_folder}/SQlite"

// FUNCTIONS
def check_allowed_options(valuestocheck, allowedvalues) {
    for (x in valuestocheck) {
        if (!allowedvalues.contains(x)) {
            exit 1, "Invalid selection: ${x}. Valid options: ${allowedvalues.join(', ')}"
        }
    }
}

def print_dataset_items(dataset, main_folder, sources) {
    exists_string = ""
    dataset_file = sources[dataset].file
    dataset_path = file("${main_folder}/${dataset_file}")
    myfile = file(dataset_path)
    if (dataset_path.exists()) {
        exists_string = "*"
    }
    println "${exists_string}${dataset}: ${dataset_path}"
}

def define_selected_values(x, best, known_values) {
    switch(x) {
        case 'all':
            selected_values = known_values
            break
        case 'best':
            selected_values = best
            break
        default:
            selected_values = x.split(",")
        }
        check_allowed_options(selected_values, known_values)
    return selected_values
}

def check_annotation_file(dataset, annotations_def, folder) {
    dataset_file = annotations_def[dataset].file
    dataset_path = file("${folder}/${dataset_file}")
    if (!dataset_path.exists()) {
        log.info """\
            $dataset not found at $dataset_path
            """.stripIndent()
        if (params.nodownload) {
            log.error "Annotation file not found and download is disabled"
            exit 1 
        }
        return false
    } else {
        return true
    }
}

// Print help message when --help is used
params.help = false
if (params.help) {
    println """\
        GREEN-VARAN annotation download - PARAMETERS    
        ============================================
        --build GRCh37/GRCh38   :   Genome build
        --db                    :   Download GREEN-DB bed file
        --scores best/all/name  :   Download prediction scores
                                    best: ncER, FATHMM-MKL, ReMM
                                    all: all scores
                                    name download only the specified score (can be comma-separated list)
        --regions best/all/name :   Download functional regions
                                    best: TFBS, DNase, UCNE
                                    all: all regions
                                    name download only the specified regions (can be comma-separated list)
        --resource_folder       :   Specify a custom folder for the annotation files
        --list_data             :   Output the list of available scores / regions
        """
        .stripIndent()

    exit 0
}

// Print list of available scores/regions when --list_data is used
params.list_data = false
if (params.list_data) {
    println "Available datasets. Star indicate corresponding file is available locally"
    println "=== GREENDB ==="
    print_dataset_items('greendb_bed', resource_folder, params.annotations[params.build].database) 
    print_dataset_items('greendb', sql_folder, params.annotations.sqlite) 

    println "=== SCORES ==="
    for (x in known_scores) {
        print_dataset_items(x, resource_folder, params.annotations[params.build].scores) 
    }

    println "=== REGIONS ==="
    for (x in known_regions) {
        print_dataset_items(x, resource_folder, params.annotations[params.build].regions) 
    }

    println "=== GNOMAD AF ==="
    print_dataset_items('gnomAD', resource_folder, params.annotations[params.build].AF) 

    exit 0
}

log.info """\
    ==============================================
    GREEN-VARAN annotation - DOWNLOAD SOURCES    
    ==============================================
    build               : ${params.build}
    resource folder     : ${params.resource_folder}

    REQUESTED ANNOTATIONS:
        GREENDB         : ${params.db}
        Scores          : ${params.scores}
        Regions         : ${params.regions}
    ==============================================
    """
    .stripIndent()

// Check requested build is allowed
if (!allowed_builds.contains(params.build)) {
    exit 1, "Invalid genome build: ${params.operation}. Valid options: ${allowed_builds.join(', ')}"
}

// warning if no download is required
if (!params.scores && !params.regions && !params.db) {
    exit 1, "No download required, have you forgotten to set proper options?"
}

// inclusion statements
include { download_dataset as DOWNLOAD_DB } from './modules/utils' addParams( annotations : params.annotations[params.build].database, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_SQL } from './modules/utils' addParams( annotations : params.annotations.sqlite, resource_folder : sql_folder)
include { download_dataset as DOWNLOAD_SCORE } from './modules/utils' addParams( annotations : params.annotations[params.build].scores, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_REGION } from './modules/utils' addParams( annotations : params.annotations[params.build].regions, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_AF } from './modules/utils' addParams( annotations : params.annotations[params.build].AF, resource_folder : resource_folder)

//WORKFLOW
workflow {    
    //Download DB file if missing
    if (params.db) {
        if (!check_annotation_file('greendb_bed', params.annotations[params.build].database, resource_folder)) {
            DOWNLOAD_DB('greendb_bed')
        }
        if (!check_annotation_file('greendb', params.annotations.sqlite, sql_folder)) {
            DOWNLOAD_SQL('greendb')
        }
    }
    
    if (params.scores) {
        missing_scores = []
        selected_scores = define_selected_values(params.scores, ['ncER','FATHMM_MKLNC','ReMM'], known_scores)
        for (s in selected_scores) {
            if (!check_annotation_file(s, params.annotations[params.build].scores, resource_folder)) {
                missing_scores = missing_scores.plus(s)
            }
        }
        missing_scores_channel = Channel.fromList(missing_scores)
        DOWNLOAD_SCORE(missing_scores_channel)
    }
    
    if (params.regions) {
        missing_regions = []
        selected_regions = define_selected_values(params.regions, ['TFBS','DNase','UCNE'], known_regions)
        //log.info "Selected regions: $selected_regions"
        for (r in selected_regions) {
            if (!check_annotation_file(r, params.annotations[params.build].regions, resource_folder)) {
                missing_regions = missing_regions.plus(r)
            }
        }
        missing_regions_channel = Channel.fromList(missing_regions)
        DOWNLOAD_REGION(missing_regions_channel)
    }

    if (params.AF) {
        if (!check_annotation_file('gnomAD', params.annotations[params.build].AF, resource_folder)) {
            missing_AF_channel = Channel.of('gnomAD')
            DOWNLOAD_AF(missing_AF_channel)
        } 
    }
}

workflow.onComplete { 
	log.info ( workflow.success ? """\
    
    PIPELINE COMPLETED!    
    ============================================
    Annotated files downloaded to ${resource_folder}
    """.stripIndent() : "Oops .. something went wrong" )
}