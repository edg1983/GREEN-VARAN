nextflow.enable.dsl=2

workflow ANNOTATE {
    take: 
        //a channel of toml files they will be concatenated for annotation
        vcf_file
        toml_files 
    
    main:
        annotate_vcf(vcf_file, toml_files)
        index_vcf(annotate_vcf.out)
    
    emit:
        index_vcf.out //a tuple: file(vcf), file(index)
}

process annotate_vcf {
    label 'multicore'
    
    input:
        tuple file(vcf_file), file(vcf_index)
        file(toml)

    output:
        file("*.tmp.vcf.gz")

    script:
    """
    f="$vcf_file"
    prefix=\${f##*/}
    prefix=\${prefix%%.gz}
    prefix=\${prefix%%.vcf}

    GOGC=2000 IRELATE_MAX_CHUNK=10000 IRELATE_MAX_GAP=1000 \
    vcfanno -p${params.ncpus} $toml $vcf_file | bgzip -c > \${prefix}.tmp.vcf.gz
    """
}

process index_vcf {
    label 'singlecore'
    
    input:
        file(vcf_file)

    output:
        tuple file("$vcf_file"), file("${vcf_file}.csi")

    script:
    """
    tabix -m12 --csi $vcf_file
    """
}