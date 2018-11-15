from os.path import join

rule snpeff_gzipped:
    input:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.vcf.gz')
    output:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.snpeff_annotated.vcf.gz')
    params:
        # Required parameters.
        genome_version = config['snpeff']['genome_version'],
        # Optional parameters. Omit if unused.
        java_options = '-Xmx4g',
        # It true, there will be a significant speedup if there are a lot of samples.
        no_statistics = False,
        extra = ''
    threads: config['threads']['snpeff']
    log: 'logs/snpeff/ann/{tumor}_vs_{normal}.log'
    wrapper:
        'http://dohlee-bio.info:9193/snpeff/ann'

rule snpeff:
    input:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.vcf')
    output:
        join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.snpeff_annotated.vcf')
    params:
        # Required parameters.
        genome_version = config['snpeff']['genome_version'],
        # Optional parameters. Omit if unused.
        java_options = '-Xmx4g',
        # It true, there will be a significant speedup if there are a lot of samples.
        no_statistics = False,
        extra = ''
    threads: config['threads']['snpeff']
    log: 'logs/snpeff/ann/{tumor}_vs_{normal}.log'
    wrapper:
        'http://dohlee-bio.info:9193/snpeff/ann'

