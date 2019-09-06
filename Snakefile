configfile: 'config.yaml'

SCRIPT_DIR      = config.get('script_dir', 'scripts')

ALF_BIN         = config.get('alf_bin', 'alfsim')
ALF_CONF        = config.get('alf_conf_template', 'etc/alf_template.drw')
DISTMAT_DIR     = config.get('distmat_dir', 'distances')
GENOME_DIR      = config.get('genome_dir', 'genomes')
NO_GENES        = config.get('no_genes', [100])
NO_SPECIES      = config.get('no_species', 10)
PAMS            = config.get('pams', [10])
RESULT_DIR      = config.get('result_dir', 'results')
SIM_DATA_DIR    = config.get('simulation_data_dir', 'alf_data')
SI_K            = config.get('si_k_neighborhood', [10])
TP_RATE         = config.get('transposition_rate', 0.01)
TREE_DIR        = config.get('tree_dir', 'trees')
REPEATS         = config.get('no_repeats', 1)

rule all:
    input:
        expand('%s/s{no_species}_n{no_genes}_pam{pam}_k{si_k}/{no_repeats}.txt' %
                RESULT_DIR, no_species = NO_SPECIES, no_genes = NO_GENES, pam =
                PAMS, si_k=SI_K, no_repeats = range(REPEATS))


rule copy_alf_conf:
    input:
        ALF_CONF
    output:
        '%s/s{no_species}_n{no_genes}_pam{pam,[^/]+}/{no_repeat,\d+}.drw' %SIM_DATA_DIR
    shell:
        'OUT_DIR=$(dirname \"{output}\" | sed \'s/\//\\\\\\//\') && '
        'sed -e \"s/{{OUT_DIR}}/${{OUT_DIR}}\\/{wildcards.no_repeat}/\" '
        '-e \'s/{{PAM}}/{wildcards.pam}/\' '
        '-e \'s/{{NO_SPECIES}}/{wildcards.no_species}/\' '
        '-e \'s/{{NO_GENES}}/{wildcards.no_genes}/\' '
        '-e \'s/{{TP_RATE}}/%s/\' {input} > {output}' %TP_RATE


rule run_alf:
    input:
        '%s/s{no_species}_{alf_config}/{no_repeats}.drw' %SIM_DATA_DIR
    output:
        expand('%s/s{{no_species}}_{{alf_config}}/{{no_repeats,\d+}}/DB/{species}_aa.fa' %SIM_DATA_DIR,
                species = map(lambda x: 'SE%03i' %x, range(1, NO_SPECIES+1))),
        '%s/s{no_species}_{alf_config}/{no_repeats}/VP.tgz' %SIM_DATA_DIR,
        '%s/s{no_species}_{alf_config}/{no_repeats}/RealTree.nwk' %SIM_DATA_DIR
    log:
        'logs/alfsim_s{no_species}_{alf_config}_r{no_repeats}.log'
    shell:
        'IN_FILE="{input}"; rm -rf ${{IN_FILE%%.drw}};  %s {input} > {log}' %ALF_BIN


rule untar_VP:
    input:
        '%s/{alf_simul}/VP.tgz' %SIM_DATA_DIR,
    output:
        '%s/{alf_simul}/VP/HP.drw' %SIM_DATA_DIR
    shell:
        'tar -xzf {input} -C $(dirname {input})'


rule construct_gene_orders:
    input:
        fasta = expand('%s/s{{no_species}}_{{alf_config}}/DB/{species}_aa.fa' %
                SIM_DATA_DIR, species = map(lambda x: 'SE%03i' %x, range(1,
                    NO_SPECIES+1))),
        HP = '%s/s{no_species}_{alf_config}/VP/HP.drw' %SIM_DATA_DIR
    output:
        '%s/s{no_species}_{alf_config}.fa' %GENOME_DIR
    shell:
        '%s/vp_to_dists.py {input.HP} > {output}' %SCRIPT_DIR


rule construct_distance_matrices:
    input:
        '%s/{alf_simul}/{no_repeats}.fa' %GENOME_DIR
    output:
        '%s/{alf_simul}_k{si_k}/{no_repeats,\d+}.csv' %DISTMAT_DIR
    log:
        'logs/si_matrix_{alf_simul}_k{si_k}_r{no_repeats}.log'
    shell:
        '%s/si_matrix.py {input} {wildcards.si_k} > {output}' %SCRIPT_DIR


rule construct_tree:
    input:
        '%s/{alf_simul}.csv' %DISTMAT_DIR
    output:
        '%s/{alf_simul}.nwk' %TREE_DIR
    run:
        from dendropy import PhylogeneticDistanceMatrix as PDM

        D = PDM.from_csv(open(input[0]), delimiter='\t')
        T = D.nj_tree()
        T.write_to_path(output[0], schema='Newick')


rule compare_with_true_tree:
    input:
        inferred = '%s/{alf_simul}_k{si_k}/{no_repeats}.nwk' %TREE_DIR,
        real =  '%s/{alf_simul}/{no_repeats}/RealTree.nwk' %SIM_DATA_DIR
    output:
        '%s/{alf_simul}_k{si_k}/{no_repeats,\d+}.txt' %RESULT_DIR
    run:
        from dendropy import Tree, TaxonNamespace
        from dendropy.calculate import treecompare

        tns = TaxonNamespace()
        T0 = Tree.get(file=open(input[0]), schema='Newick',
                taxon_namespace=tns)
        T1 = Tree.get(file=open(input[1]), schema='Newick',
                taxon_namespace=tns)
        print(treecompare.symmetric_difference(T0, T1), file = open(output[0],
            'w'))

