nextflow.enable.dsl=2

params.resource_folder = "resources"
params.annotations = []
params.build = "GRCh38"

process download_dataset {
    label 'singlecore'
    publishDir "${params.resource_folder}", mode: 'move'
    
    input:
        val(dataset)

    output:
        tuple val(dataset), val("${params.resource_folder}/${dataset_file}"), emit: downloaded_file
        file "*.gz*"

    script:
    file(params.resource_folder).mkdirs() //this will create the directory if missing
    dataset_file = params.annotations[dataset].file
    webaddress = "${params.annotations[dataset].repository}/files/${dataset_file}"
    """
    wget $webaddress
    wget ${webaddress}.csi
    """
}

process write_toml {
    label 'singlecore'
    
    input:
        val(dataset)
        val(op)

    output:
        file "${dataset}.toml"

    script:
    id = "$dataset"
    dataset_file = params.annotations[id].file
    file = "${params.resource_folder}/${dataset_file}"
    
    if (params.annotations[id].columns && params.annotations[id].fields) {
        exit 1, "Only one of columns or fields can be set. Found both for ${id}"
    }
    
    if (params.annotations[id].columns) {
        input_fields = params.annotations[id].columns
        input_tag = "columns"
    } else if (params.annotations[id].fields) {
        single_fields = params.annotations[id].fields.split(',')
        input_fields = single_fields.collect { '"' + "$it" + '"' }.join(',')
        input_tag = "fields"
    } else {
        input_fields = "1"
        input_tag = "columns"
    }
    
    input_fields_string = "${input_tag}=[${input_fields}]"
    n_fields = input_fields.split(',').size()

    if (params.annotations[id].names) {
        single_names = params.annotations[id].names.split(',')
        if (single_names.size() != n_fields) {
            exit 1, "columns / fields and names must contain the same number of elements. Found $n_fields fields and ${single_names.size()} names for ${id}"
        }
        names = single_names.collect { '"' + "$it" + '"' }.join(',')
    } else {
        if (n_fields > 1) {
            names = (1..n_fields).collect { '"' + "${id}_$it" + '"' }.join(',')
        } else {
            names = '"' + "$id" + '"'
        }
    }

    if (params.annotations[id].ops) {
        single_ops = params.annotations[id].ops.split(',')
        if (single_ops.size() == 1) {
            ops = (1..n_fields).collect { '"' + "${params.annotations[id].ops}" + '"' }.join(',')
        } else if (single_ops.size() == n_fields) {
            ops = single_ops.collect { '"' + "$it" + '"' }.join(',')
        } else {
            exit 1, "ops must be a single value or must contain the same number of elements of columns/fields. Found $n_fields fields and ${single_ops.size()} ops for ${id}"
        }
    } else {
        if (n_fields > 1) {
            ops = (1..n_fields).collect { '"' + "$op" + '"' }.join(',')
        } else {
            ops = '"' + "$op" + '"'
        }
    }
    
    """
    echo "[[annotation]]" > ${dataset}.toml
    echo 'file="$file"' >> ${dataset}.toml
    echo '$input_fields_string' >> ${dataset}.toml
    echo 'names=[$names]' >> ${dataset}.toml
    echo 'ops=[$ops]' >> ${dataset}.toml
    """
}

process concat_toml {
    label 'singlecore'
    publishDir "${params.out}", mode: 'copy'
    
    input:
        file(block)
        val(outfile)

    output:
        file "${outfile}.toml"

    script:
    def toml_blocks = block.collect{ "$it" }.join(" ")
    """
    cat $toml_blocks > ${outfile}.toml
    """
}