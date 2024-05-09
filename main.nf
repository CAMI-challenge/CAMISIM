#!/usr/bin/env nextflow

nextflow.enable.dsl=2

// Include the specific pipeline based on the parameter
if (params.pipeline == "metagenomic") {
    include { metagenomic } from "${projectDir}/workflow"
} else if (params.pipeline == "metatranscriptomic" ) {    
    include { metatranscriptomic } from "${projectDir}/pipelines/metatranscriptomic/metatranscriptomic"
}    

workflow {
    if (params.pipeline == "metagenomic") {
        metagenomic()
    } else if (params.pipeline == "metatranscriptomic") {
        metatranscriptomic()
    }
}