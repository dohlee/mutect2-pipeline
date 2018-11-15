from os.path import join

rule mutect2_tumor_normal:
    input:
        # Required arguments.
        reference = config['reference']['fasta'],
        tumor = wxs_preprocessing(join(config['result_dir']['tumor'], '{tumor}.duplicates_marked.recalibrated.sorted.bam')),
        normal = wxs_preprocessing(join(config['result_dir']['normal'], '{normal}.duplicates_marked.recalibrated.sorted.bam')),
        tumor_idx = wxs_preprocessing(join(config['result_dir']['tumor'], '{tumor}.duplicates_marked.recalibrated.sorted.bam.bai')),
        normal_idx = wxs_preprocessing(join(config['result_dir']['tumor'], '{tumor}.duplicates_marked.recalibrated.sorted.bam.bai')),
        reference_dict = 'reference/Homo_sapiens_assembly38.dict',
        # Optional arguments. Omit unused files.
        germline_resource = config['resource']['known_variant_sites'],
        germline_resource_index = config['resource']['known_variant_sites'] + '.tbi',
        panel_of_normals = config['resource']['pon'],
    output:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.raw.mutect2.vcf.gz')
    params:
        extra = '--disable-read-filter MateOnSameContigOrNoMappedMateReadFilter ' \
                '--af-of-alleles-not-in-resource 0.0000025 ',
        # Optional parameters. Omit if unused.
        tumor_sample_name = '',
        normal_sample_name = '',
        java_options = ''
    resources: RAM = 16
    threads: config['threads']['gatk']['mutect2_tumor_normal']
    log: 'logs/gatk/mutect2-tumor-normal/{tumor}_vs_{normal}.log'
    wrapper:
        'http://dohlee-bio.info:9193/gatk/mutect2/tumor-normal'

rule get_pileup_summaries:
    input:
        # Required arguments.
        bam = wxs_preprocessing(join(config['result_dir']['sample'], '{sample}.duplicates_marked.recalibrated.sorted.bam')),
        common_biallelic_variants = config['resource']['common_biallelic_variants'],
        # Optional arguments. Omit unused files.
    output:
        join(config['result_dir']['sample'], '{sample}.pileupsummaries.table'),
    params:
        java_options = '-Xmx4g',
        extra = '',
    threads: config['threads']['gatk']['get_pileup_summaries']
    resources: RAM = 4
    log: 'logs/gatk/get-pileup-summaries/{sample}.log'
    wrapper:
        'http://dohlee-bio.info:9193/gatk/coverage/get-pileup-summaries'

rule calculate_contamination:
    input:
        # Required arguments.
        pileup_summaries = join(config['result_dir']['tumor'], '{tumor}.pileupsummaries.table'),
        normal_pileup_summaries = join(config['result_dir']['normal'], '{normal}.pileupsummaries.table'),
        # Optional arguments. Omit unused files.
    output:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.contamination.table')
    params:
        java_options = '-Xmx4g',
        extra = lambda wildcards, input: '--matched-normal ' + input.normal_pileup_summaries,
    threads: config['threads']['gatk']['calculate_contamination']
    resources: RAM = 4
    log: 'logs/gatk/calculate-contamination/{tumor}_vs_{normal}.log'
    wrapper:
        'http://dohlee-bio.info:9193/gatk/qc/calculate-contamination'

rule filter_mutect_calls:
    input:
        # Required input.
        vcf = join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.raw.mutect2.vcf.gz'),
        contamination_table = join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.contamination.table'),
    output:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.vcf.gz'),
    params:
        # Optional parameters. Omit if unused.
        extra = ''
    threads: config['threads']['gatk']['filter_mutect_calls']
    log: 'logs/gatk/filter-mutect-calls/{tumor}_vs_{normal}.log'
    wrapper:
        'http://dohlee-bio.info:9193/gatk/variant-filtering/filter-mutect-calls'

rule create_sequence_dictionary:
    input:
        config['reference']['fasta'],
    output:
        config['reference']['basename'] + '.dict'
    params:
        # Optional parameters. Omit if unused.
        java_options = '-Xmx4g'
    threads: config['threads']['gatk']['create_sequence_dictionary']
    resources: RAM = 4
    log: 'logs/gatk/create-sequence-dictionary/log.log'
    wrapper:
        'http://dohlee-bio.info:9193/gatk/reference/create_sequence_dictionary'

