nextflow.enable.dsl=2

// Define allowed options
def allowed_builds  = ['GRCh37', 'GRCh38']
def known_scores  = params.annotations[params.build].scores.keySet() as String[]
def allowed_scores_values = known_scores + ['best', 'all']
def known_regions  = params.annotations[params.build].regions.keySet() as String[]
def allowed_regions_values = known_regions + ['best', 'all']
def resource_folder = "${params.resource_folder}/${params.build}"

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

def define_selected_values(x, best, allowed_values) {
    switch(x) {
        case 'all':
            selected_values = params.known_scores
            break
        case 'best':
            selected_values = best
            break
        default:
            selected_values = x.split(",")
        }
        check_allowed_options(selected_values, allowed_values)
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
        GREEN-VARAN annotation pipeline - PARAMETERS    
        ============================================
        --input in.vcf.gz       :   Input VCF file(s), compressed and indexed
                                    You can input multiple files from a folder using quotes like
                                    --input 'mypath/*.vcf.gz'
        --build GRCh37/GRCh38   :   Genome build
        --out output_dir        :   Output directory
        --scores best/all/name  :   Annotate prediction scores
                                    best annotate ncER, FATHMM-MKL, ReMM
                                    all annotate all scores
                                    name annotate only the specified score (can be comma-separated list)
        --regions best/all/name :   Annotate functional regions
                                    best annotate TFBS, DNase, UCNE
                                    all annotate all regions
                                    name annotate only the specified regions (can be comma-separated list)
        --AF                    :   Annotate global AF from gnomAD genome v3
        --greenvaran_config     :   A json config file for GREEN-VARAN tool
        --greenvaran_dbschema   :   A json db schema file for GREEN-VARAN tool
        --resource_folder       :   Specify a custom folder for the annotation files
        --anno_toml             :   A custom annotation config file.
                                    This file is a toml file as specified for vcfanno.
                                    This will be added to other annotations defined with 
                                    scores and regions.
        --list_data             :   output the list of available scores / regions

        See https://github.com/brentp/vcfanno for guide on how to prepare toml files or 
        see an example in the config folder
        """
        .stripIndent()

    exit 1
}

// Print list of available scores/regions when --list_data is used
params.list_data = false
if (params.list_data) {
    println "Available datasets. Star indicate corresponding file is available locally"
    println "=== SCORES ==="
    for (x in known_scores) {
        print_dataset_items(x, resource_folder, params.annotations[params.build].scores) 
    }

    println "=== REGIONS ==="
    for (x in known_regions) {
        print_dataset_items(x, resource_folder, params.annotations[params.build].regions) 
    }

    println "=== GNOMAD AF ==="
    print_dataset_items('gnomAD', resource_folder, params.annotations[params.build].regions) 
    exit 0
}

log.info """\
    ==============================================
    GREEN-VARAN annotation - N F   P I P E L I N E    
    ==============================================
    input file          : ${params.input}
    build               : ${params.build}
    output              : ${params.out}
    greenvaran config   : ${params.greenvaran_config}
    greenvaran dbschema : ${params.greenvaran_dbschema}
    resource folder     : ${params.resource_folder}

    ACTIVE ANNOTATIONS:
        Scores          : ${params.scores}
        Regions         : ${params.regions}
        AF              : ${params.AF}
    ==============================================
    """
    .stripIndent()

// Check requested build is allowed
if (!allowed_builds.contains(params.build)) {
    exit 1, "Invalid genome build: ${params.operation}. Valid options: ${allowed_builds.join(', ')}"
}

// Check input files exist
checkPathParamList = [
    params.input, params.greenvaran_config, params.greenvaran_dbschema
]
for (param in checkPathParamList) { if (param) { file(param, checkIfExists: true) } }

if (params.anno_toml) {
    file(params.anno_config, checkIfExists: true)
}

//Make output dir and set prefix
outputdir = file(params.out)
outputdir.mkdirs()
outprefix = file(params.input).getSimpleName()

// inclusion statements
include { download_dataset as DOWNLOAD_DB } from './modules/utils' addParams( annotations : params.annotations[params.build].database, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_SCORE } from './modules/utils' addParams( annotations : params.annotations[params.build].scores, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_REGION } from './modules/utils' addParams( annotations : params.annotations[params.build].regions, resource_folder : resource_folder)
include { download_dataset as DOWNLOAD_AF } from './modules/utils' addParams( annotations : params.annotations[params.build].AF, resource_folder : resource_folder)

include { write_toml as WRITE_SCORE_TOML } from './modules/utils' addParams( annotations : params.annotations[params.build].scores, resource_folder : resource_folder)
include { write_toml as WRITE_REGION_TOML } from './modules/utils' addParams( annotations : params.annotations[params.build].regions, resource_folder : resource_folder)
include { write_toml as WRITE_AF_TOML } from './modules/utils' addParams( annotations : params.annotations[params.build].AF, resource_folder : resource_folder)

include { ANNOTATE } from './modules/annotate'
include { concat_toml as concat_scores_toml; concat_toml as concat_regions_toml; concat_toml} from './modules/utils'

//WORKFLOW
workflow {    
    def missing_data = [scores: [], regions: [], AF: []]
    def existing_data = [scores: [], regions: [], AF: []]

    //Read input VCF
    input_vcf = Channel
                .fromPath(params.input)
                .map { tuple(file("$it"), file("${it}.{csi,tbi}"))}

    //Download DB file if missing
    if (!check_annotation_file('greendb', params.annotations[params.build].database, resource_folder)) {
        DOWNLOAD_DB('greendb_bed')
    }
    
    //Check is annotations files exist, download them eventually and prepare toml files
    toml_files = Channel.empty()
    if (params.scores) {
        selected_scores = define_selected_values(params.scores, ['ncER','FATHMM_MKLNC','ReMM'], allowed_scores_values)
        //log.info "Selected scores: $selected_scores"
        for (s in selected_scores) {
            if (check_annotation_file(s, params.annotations[params.build].scores, resource_folder)) {
                existing_data.scores = existing_data.scores.plus(s)
            } else {
                missing_data.scores = missing_data.scores.plus(s)
            }
        }

        existing_scores_channel = Channel.fromList(existing_data.scores)
        if (missing_data.scores.size() > 0) {
            missing_scores_channel = Channel.fromList(missing_data.scores)
            DOWNLOAD_SCORE(missing_scores_channel)
            scores_channel = existing_scores_channel.mix(DOWNLOAD_SCORE.out)
        } else {
            scores_channel = existing_scores_channel
        }
        WRITE_SCORE_TOML(scores_channel, "max")
        toml_files = toml_files.concat(WRITE_SCORE_TOML.out)
    }
    
    if (params.regions) {
        selected_regions = define_selected_values(params.regions, ['TFSB','DNase','UCNE'], allowed_regions_values)
        //log.info "Selected regions: $selected_regions"
        for (s in selected_regions) {
            if (check_annotation_file(s, params.annotations[params.build].regions, resource_folder)) {
                existing_data.regions = existing_data.regions.plus(s)
            } else {
                missing_data.regions = missing_data.regions.plus(s)
            }
        }

        existing_regions_channel = Channel.fromList(existing_data.regions)
        if (missing_data.regions.size() > 0) {
            missing_regions_channel = Channel.fromList(missing_data.regions)
            DOWNLOAD_REGION(missing_regions_channel)
            regions_channel = existing_regions_channel.mix(DOWNLOAD_REGION.out)
        } else {
            regions_channel = existing_regions_channel
        }
        WRITE_REGION_TOML(regions_channel, "flag")
        toml_files = toml_files.concat(WRITE_REGION_TOML.out)
        //concat_regions_toml(WRITE_REGION_TOML.out.collect(), "regions")
    }

    if (params.AF) {
       if (check_annotation_file('gnomAD', params.annotations[params.build].AF, resource_folder)) {
            existing_data.AF = existing_data.AF.plus(s)
        } else {
            missing_data.AF = missing_data.AF.plus(s)
        }

        existing_AF_channel = Channel.fromList(existing_data.AF)
        if (missing_data.AF.size() > 0) {
            missing_AF_channel = Channel.fromList(missing_data.AF)
            DOWNLOAD_AF(missing_AF_channel)
            af_channel = existing_AF_channel.mix(DOWNLOAD_AF.out)
        } else {
            af_channel = existing_AF_channel
        }
        WRITE_AF_TOML(regions_channel, "max")
        toml_files = toml_files.concat(WRITE_AF_TOML.out)
        //concat_regions_toml(WRITE_REGION_TOML.out.collect(), "regions")
    }

    if (params.anno_toml) {
        anno_toml_channel = Channel.fromPath(params.anno_toml)
        toml_files = toml_files.concat(anno_toml_channel)
    }

    //When regions or scores are active perform vcfanno annotation
    if (params.regions || params.scores || params.AF) {
        //toml_files = Channel.fromPath("$outputdir/*.toml")
        concat_toml(toml_files.collect(), "annotations")
        ANNOTATE(input_vcf, concat_toml.out)
        annotated_vcf = ANNOTATE.out
    } else {
        annotated_vcf = input_vcf
    }
    
    //RUN GREEN-VARAN
    greendb = tuple(
        file("${resource_folder}/${params.annotations[params.build].database.greendb_bed.file}"),
        file("${resource_folder}/${params.annotations[params.build].database.greendb_bed.file}.csi")
    )
    green_varan(annotated_vcf,
        greendb,
        file(params.greenvaran_config), 
        file(params.greenvaran_dbschema), 
        params.greenvaran_options
    )
}

workflow.onComplete { 
	log.info ( workflow.success ? """\
    
    PIPELINE COMPLETED!    
    ============================================
    Annotated VCF file  --> ${outputdir}
    """.stripIndent() : "Oops .. something went wrong" )
}

process green_varan {
    label 'singlecore'
    publishDir "$outputdir", mode: 'move'
    
    input:
        tuple file(vcf_file), file(vcf_index)
        tuple file(db), file(db_index)
        file(config)
        file(dbschema)
        val(options)

    output:
        file "*.annotated.vcf.gz"
        file "*.log"

    script:
    """
    f="$vcf_file"
    prefix=\${f##*/}
    prefix=\${prefix%%.tmp.vcf.gz}

    greenvaran smallvars \
    -i $vcf_file \
    -o \${prefix}.annotated.vcf.gz \
    -c $config \
    --db $db \
    --dbschema $dbschema \
    $options
    """
}