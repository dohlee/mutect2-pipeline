from os.path import join

# WORKFLOW NOTIFICATIONS.
onstart: shell('')
onerror: shell('')
onsuccess: shell('')

# IMPORT SUBWORKFLOWS.
subworkflow wxs_preprocessing:
    workdir: 'wxs-preprocessing-pipeline'

configfile: 'config.yaml'
include: 'rules/gatk.smk'
include: 'rules/snpeff.smk'
include: 'rules/util.smk'

wildcard_constraints:
    sample = config['wildcard_constraints']['sample'],
    normal = config['wildcard_constraints']['normal'],
    tumor = config['wildcard_constraints']['tumor'],

# DEFINE YOUR SAMPLES HERE.
SAMPLES = []
TUMOR_SAMPLES, NORMAL_SAMPLES = [], []

rule all:
    input:
        expand(join(config['result_dir']['tumor'], '{tumor}_vs_{normal}.mutect2.snpeff_annotated.vcf'), zip, tumor=TUMOR_SAMPLES, normal=NORMAL_SAMPLES),
    
