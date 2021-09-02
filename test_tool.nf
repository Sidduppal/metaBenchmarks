#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.reads = "/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz"

include{ FASTQC } from '/media/bigdrive1/sidd/autometa2_manuscript/metaBenchmarks/modules/nf-core/modules/fastqc/main.nf'

include {KRAKEN2_KRAKEN2} from '/media/bigdrive1/sidd/autometa2_manuscript/metaBenchmarks/modules/nf-core/modules/kraken2/kraken2/main.nf'

workflow sidd{
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")
    FASTQC(reads)
    KRAKEN2_KRAKEN2 (reads, db)
}

workflow chase{
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")
 
    Channel
        .fromPath(reads)
            .map { row ->
                    def meta = [:]
                    meta.id = row.simpleName
                    return [ meta, row ]
                }
            .set { ch_fasta }


    FASTQC(ch_fasta)
    //KRAKEN2_KRAKEN2 (ch_fasta, db)
}

workflow sidd2{
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")
 
    Channel
        .fromPath(reads)
            .map { read -> tuple(read.simpleName, read)}
            .set { ch_fasta }


    FASTQC(ch_fasta)
    //KRAKEN2_KRAKEN2 (ch_fasta, db)
}