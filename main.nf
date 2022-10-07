
nextflow.enable.dsl = 2

include { CHECK_STRANDNESS } from './modules/check_strandness.nf'
include { FASTP } from './modules/fastp'
include { GENERATE_DECOY_TRANSCIPTROME ; SALMON_INDEX_REFERENCE ; SALMON_ALIGN_QUANT } from './modules/salmon.nf'

log.info """\
         RNAseq differential analysis using NextFlow 
         =============================
         outdir: ${params.outdir}
         """
         .stripIndent()
 
params.outdir = 'results'

workflow {
    
    read_pairs_ch = channel.fromFilePairs( params.reads, checkIfExists: true ) 
    CHECK_STRANDNESS( read_pairs_ch, params.reference_cdna_ensembl, params.reference_annotation_ensembl )
    FASTP( read_pairs_ch )
    GENERATE_DECOY_TRANSCIPTROME( params.reference_genome, params.reference_cdna )
    SALMON_INDEX_REFERENCE( GENERATE_DECOY_TRANSCIPTROME.out.decoy, GENERATE_DECOY_TRANSCIPTROME.out.gentrome )
    SALMON_ALIGN_QUANT( CHECK_STRANDNESS.out, FASTP.out.sample_trimmed, SALMON_INDEX_REFERENCE.out, params.reference_annotation )
}

