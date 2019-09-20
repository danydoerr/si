configfile: 'config.yaml'

from itertools import product

SCRIPT_DIR      = config.get('script_dir', 'scripts')

ALF_BIN         = config.get('alf_bin', 'alfsim')
ALF_CONF        = config.get('alf_conf_template', 'etc/alf_template.drw')
DISTMAT_DIR     = config.get('distmat_dir', 'distances')
GENOME_DIR      = config.get('genome_dir', 'genomes')
NO_GENES        = config.get('no_genes', [100])
NO_SPECIES      = config.get('no_species', 10)
PAMS            = config.get('pams', [10])
REPEATS         = config.get('no_repeats', 1)
RESULT_DIR      = config.get('result_dir', 'results')
SIM_DATA_DIR    = config.get('simulation_data_dir', 'alf_data')
SI_K            = config.get('si_k_neighborhood', [10])
TP_RATE         = config.get('transposition_rate', 0.01)
TREE_DIR        = config.get('tree_dir', 'trees')

rule all:
    input:
        expand('%s/results_si_k{si_k}.csv' %RESULT_DIR, si_k=SI_K),
        '%s/results_dcj.csv' %RESULT_DIR,
        expand('%s/boxplot_s{no_species}_n{no_genes}.pdf' %RESULT_DIR,
                no_species = NO_SPECIES, no_genes = NO_GENES)
#        expand('%s/s{no_species}_n{no_genes}_pam{pam}/grappa/{no_repeats}.txt' %
#                RESULT_DIR, no_species = NO_SPECIES, no_genes = NO_GENES, pam =
#                PAMS, no_repeats = range(REPEATS)),


rule balanced_generate_tree:
    params:
        no_species = '{no_species}',
    output:
        temp('trees/true_tree_s{no_species,\d+}.nwk')
    shell:
        '%s/gen-balanced-tree.sh {params.no_species} 1 ' %SCRIPT_DIR +
        '> {output}'


rule rescale_tree:
    input:
        'trees/true_tree_s{no_species}.nwk'
    output:
        'trees/true_tree_s{no_species}_pam{pam,[0-9.]+}.nwk'
    shell:
        '%s/rescale.py -a{wildcards.pam} {input} > {output}' %SCRIPT_DIR

rule copy_alf_conf:
    input:
        conf = ALF_CONF,
        tree_file = 'trees/true_tree_s{no_species}_pam{pam}.nwk'
    output:
        '%s/s{no_species}_n{no_genes}_pam{pam,[^/]+}/{no_repeat,\d+}.drw' %SIM_DATA_DIR
    shell:
        'OUT_DIR=$(dirname \"{output}\" | sed \'s/\//\\\\\\//\') && '
        'TREE_FILE=$(echo \"{input.tree_file}\" | sed \'s/\//\\\\\\//\') && '
        'sed -e \"s/{{OUT_DIR}}/${{OUT_DIR}}\\/{wildcards.no_repeat}/\" '
        '-e \'s/{{PAM}}/{wildcards.pam}/\' '
        '-e \'s/{{NO_SPECIES}}/{wildcards.no_species}/\' '
        '-e \'s/{{NO_GENES}}/{wildcards.no_genes}/\' '
        '-e \"s/{{TREE_FILE}}/${{TREE_FILE}}/\" '
        '-e \'s/{{TP_RATE}}/%s/\' {input.conf} > {output}' %TP_RATE


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
        fasta = lambda wildcards: expand(
                '%s/s{no_species}_{alf_config}/DB/{species}_aa.fa' %
                SIM_DATA_DIR, no_species = (wildcards.no_species,),
                species = map(lambda x: 'SE%03i' %x, range(1,
                    int(wildcards.no_species)+1)),
                alf_config = (wildcards.alf_config, )),
        HP = '%s/s{no_species}_{alf_config}/VP/HP.drw' %SIM_DATA_DIR
    output:
        '%s/s{no_species,\d+}_{alf_config}.fa' %GENOME_DIR
    shell:
        '%s/vp_to_dists.py {input.HP} > {output}' %SCRIPT_DIR


rule construct_distance_matrices_si:
    input:
        '%s/{alf_simul}/{no_repeats}.fa' %GENOME_DIR
    output:
        '%s/{alf_simul}/si_k{si_k,\d+}/{no_repeats,\d+}.csv' %DISTMAT_DIR
    log:
        'logs/si_matrix_{alf_simul}_k{si_k}_r{no_repeats}.log'
    shell:
        '%s/si_matrix.py {input} {wildcards.si_k} > {output} 2> {log}' %SCRIPT_DIR


rule run_dcj:
    input:
        '%s/{alf_simul}/{no_repeats}.fa' %GENOME_DIR
    output:
        '%s/{alf_simul}/dcj/{no_repeats}.csv' %DISTMAT_DIR
    shell:
        '%s/compute_dcj.py {input} > {output}' %SCRIPT_DIR


rule run_grappa:
    input:
        '%s/{alf_simul}/{no_repeats}.fa' %GENOME_DIR
    output:
        '%s/{alf_simul}/grappa/{no_repeats}.out' %TREE_DIR
    log:
        'logs/grappa_{alf_simul}_r{no_repeats}.log'
    shell:
        '%s/run_grappa.sh {input} > {output} 2> {log}' %SCRIPT_DIR


rule extract_grappa_tree:
    input:
        '%s/{alf_simul}/grappa/{no_repeats}.out' %TREE_DIR
    output:
        '%s/{alf_simul}/grappa/{no_repeats}.nwk' %TREE_DIR
    shell:
        'grep -e \'^\s\?([0-9SE:,()]\+);$\' {input} | tail -n1 > {output}'


rule construct_tree_si:
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
        inferred = '%s/{alf_simul}/{subdir}/{no_repeats}.nwk' %TREE_DIR,
        real =  '%s/{alf_simul}/{no_repeats}/RealTree.nwk' %SIM_DATA_DIR
    output:
        '%s/{alf_simul}/{subdir,[^/]+}/{no_repeats,\d+}.txt' %RESULT_DIR
    run:
        from dendropy import Tree, TaxonNamespace
        from dendropy.calculate import treecompare

        tns = TaxonNamespace()
        T0 = Tree.get(file=open(input[0]), schema='Newick',
                taxon_namespace=tns)
        T1 = Tree.get(file=open(input[1]), schema='Newick',
                taxon_namespace=tns)
        #fp, fn = treecompare.false_positives_and_negatives(T0, T1)
        print(treecompare.symmetric_difference(T0, T1), file = open(output[0],
            'w'))


rule combine_into_results_table:
    input:
        expand('%s/s{no_species}_n{no_genes}_pam{pam}/{{subdir}}/{no_repeats}.txt' %
                RESULT_DIR, no_species = NO_SPECIES, no_genes = NO_GENES, pam =
                PAMS, si_k=SI_K, no_repeats = range(REPEATS))
    output:
        '%s/results_{subdir}.csv' %RESULT_DIR
    run:
        out = open(output[0], 'w')
        for x in product((RESULT_DIR,), (NO_SPECIES,), NO_GENES, PAMS,
                (wildcards.subdir,), range(REPEATS)):
            val = int(open('%s/s%s_n%s_pam%s/%s/%s.txt' %x).read())
            print('\t'.join(map(str, x[1:4] + x[5:] + (val, ))), file=out)
        out.close()


rule plot_performance:
    input:
        expand('%s/results_si_k{k}.csv' %RESULT_DIR, k = SI_K),
        '%s/results_dcj.csv' %RESULT_DIR,
    output:
        expand('%s/boxplot_s{no_species}_n{no_genes}.pdf' %RESULT_DIR,
                no_species = NO_SPECIES, no_genes = NO_GENES)
    run:
        import os
        if not os.environ.get('DISPLAY', None):
            import matplotlib; matplotlib.use('Agg')

        import matplotlib.pylab as plt
        import numpy as np

        for s, g in product(hasattr(NO_SPECIES, '__iter__') and NO_SPECIES or
                (NO_SPECIES, ), NO_GENES):
            plt.figure()
            title = '#species = %s, #genes = %s' %(s, g)
            max_dist = 0
            # plot SI distances for each choice of k
            for i, k in enumerate(SI_K):
                data = np.loadtxt('%s/results_si_k%s.csv' %(RESULT_DIR, k))
                res = data[(data[:, 0] == s) & (data[:, 1] == g), :]
                time = sorted(set(res[:, 2]))
                rf_dist = list(np.mean(res[res[:, 2] == t, -1]) for t in time)
                max_dist = max(max_dist, max(rf_dist))
                plt.plot(time, rf_dist, color='C%s' %i, label = 'k = %s'%k)

            # plot DCJ distances
            data = np.loadtxt('%s/results_dcj.csv' %RESULT_DIR)
            res = data[(data[:, 0] == s) & (data[:, 1] == g), :]
            time = sorted(set(res[:, 2]))
            rf_dist = list(np.mean(res[res[:, 2] == t, -1]) for t in time)
            max_dist = max(max_dist, max(rf_dist))
            plt.plot(time, rf_dist, color='C%s' %(i+1), label = 'DCJ')

            plt.ylim([0, max_dist])
            plt.title(title)
            plt.legend(loc='upper right')
            plt.xlabel('PAM')
            plt.ylabel('RF distance to true tree')
            plt.savefig('%s/boxplot_s%s_n%s.pdf' %(RESULT_DIR, s, g),
                    format='pdf')

