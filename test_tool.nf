#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.reads = "/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz"

include{ FASTQC } from './modules/nf-core/modules/fastqc/main.nf'

include {KRAKEN2_KRAKEN2} from './modules/nf-core/modules/kraken2/kraken2/main.nf'

include {METABAT2} from './modules/local/metabat2.nf'

workflow sidd{
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")
    FASTQC(reads)
    KRAKEN2_KRAKEN2 (reads, db)
}

workflow chase{
    //reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    reads="/Users/sidd/Research/fun/nextflow/data/reads/*.fq.gz"
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

workflow metabat2 {
    bam=file("/Users/sidd/Research/autometa_v2/test_data/alignment.bam")
    assembly=file("/Users/sidd/Research/autometa_v2/test_data/78mbp_metagenome.fna")
    Channel
        .fromPath(assembly)
            .map { name ->
                    def meta = [:]
                    meta.id = name.simpleName
                    return [meta, name]
            }
        .set {ch_assembly}
    METABAT2(ch_assembly, bam)
}