#!/usr/bin/env nextflow
nextflow.enable.dsl=2

//params.reads = "/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz"

include{ FASTQC } from './modules/nf-core/modules/fastqc/main.nf'

include {KRAKEN2_KRAKEN2} from './modules/nf-core/modules/kraken2/kraken2/main.nf'

include {METABAT2} from './modules/local/metabat2.nf'

workflow fastqc{
    reads=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/forward_reads.fastq.gz")
    //reads="/Users/sidd/Research/fun/nextflow/data/reads/*.fq.gz"
 
    Channel
        .fromPath(reads)
            .map { read ->
                    def meta = [:]
                    meta.id = read.simpleName
                    return [ meta, read ]
                }
            .set { ch_fasta }


    FASTQC(ch_fasta)
}

workflow kraken{
    assembly=file("/media/bigdrive1/sidd/nextflow_trial/autometa_runs/78mbp_manual/interim/78mbp_metagenome.filtered.fna")
    db=file("/media/bigdrive1/Databases/kraken2/kraken2_db")

    Channel
        .fromPath(assembly)
            .map { fasta ->
                    def meta = [:]
                    meta.id = fasta.simpleName
                    meta.single_end = "True"
                    return [ meta, fasta ]
                }
            .set { ch_fasta }

    KRAKEN2_KRAKEN2 (ch_fasta, db)
}


workflow metabat2 {
    //bam=file("/Users/sidd/Research/autometa_v2/test_data/alignment.bam")
    //assembly=file("/Users/sidd/Research/autometa_v2/test_data/78mbp_metagenome.fna")
    bam=file("/media/bigdrive1/sidd/nextflow_trial/autometa_runs/78mbp_manual/interim/cov-alignmentsieg74wbx/alignment.bam")
    assembly=file("/media/bigdrive1/sidd/nextflow_trial/test_data/78/78mbp_metagenome.fna")
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