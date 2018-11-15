from os.path import join

rule generate_vcf_index:
    input: '{file}.vcf.gz'
    output: '{file}.vcf.gz.tbi'
    shell: 'tabix -p vcf {input}'

rule samtools_faidx:
    input: config['reference']['fasta']
    output: config['reference']['fasta'] + '.fai'
    wrapper: 'http://dohlee-bio.info:9193/samtools/faidx'

rule unzip_annotated_vcf:
    input: '{file}.snpeff_annotated.vcf.gz'
    output: '{file}.snpeff_annotated.vcf'
    threads: 1
    shell: 'gunzip {input}'
