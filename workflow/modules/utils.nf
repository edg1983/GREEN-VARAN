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
    if (params.annotations[id].column) {
        cols = params.annotations[id].column
    } else {
        cols = "1"
    }
    n_fields = cols.split(',').size()
    if (n_fields > 1) {
        names = (1..n_fields).collect { '"' + "${id}_$it" + '"' }.join(',') 
        ops = (1..n_fields).collect { '"' + "$op" + '"' }.join(',')
    } else {
        names = '"' + "$id" + '"'
        ops = '"' + "$op" + '"'
    }
    
    """
    echo "[[annotation]]" > ${dataset}.toml
    echo 'file="$file"' >> ${dataset}.toml
    echo 'columns=[$cols]' >> ${dataset}.toml
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