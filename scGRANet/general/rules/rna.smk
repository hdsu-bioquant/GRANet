if config['genome'] == 'hg19':
    scenic_data = {
      "motif_annotation": 'aux/scenic/rcistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
      "db_ranking_tss1" : 'aux/scenic/rcistarget/hg19-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2" : 'aux/scenic/rcistarget/hg19-tss-centered-5kb-7species.mc9nr.feather',
      "motif_annotationu": 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
      "db_ranking_tss1u" : 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2u" : 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.feather'
    }
elif config['genome'] == 'hg38':
    # will use the same files as hg19 because there is more info for hg19
    scenic_data = {
      "motif_annotation": 'aux/scenic/rcistarget/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
      "db_ranking_tss1" : 'aux/scenic/rcistarget/hg19-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2" : 'aux/scenic/rcistarget/hg19-tss-centered-5kb-7species.mc9nr.feather',
      "motif_annotationu": 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
      "db_ranking_tss1u" : 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2u" : 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg19/refseq_r45/mc9nr/gene_based/hg19-tss-centered-5kb-7species.mc9nr.feather'
    }
elif config['genome'] == 'mm9':
    scenic_data = {
      "motif_annotation": 'aux/scenic/rcistarget/motifs-v9-nr.mgi-m0.001-o0.0.tbl',
      "db_ranking_tss1" : 'aux/scenic/rcistarget/mm9-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2" : 'aux/scenic/rcistarget/mm9-tss-centered-5kb-7species.mc9nr.feather',
      "motif_annotationu": 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.mgi-m0.001-o0.0.tbl',
      "db_ranking_tss1u" : 'https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-10kb-7species.mc9nr.feather',
      "db_ranking_tss2u" : 'https://resources.aertslab.org/cistarget/databases/mus_musculus/mm9/refseq_r45/mc9nr/gene_based/mm9-tss-centered-5kb-7species.mc9nr.feather'
    }

#------------------------------------------------------------------------------#
#                                     SCENIC                                   #
#------------------------------------------------------------------------------#
# Run SCENIC over a scRNAseq matrix
rule scenic_step4_scenictoRDS:
    input:
        auc_mtx = join(DATAPATH, 'results/scrna/SCENIC/auc_mtx.csv')
    output:
        auc_mtxRDS = join(DATAPATH, 'results/scrna/SCENIC/auc_mtx.RDS')
    params:
        script = 'scripts/scRNA-seq/02_SCENIC_to_RDS.R',
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.auc_mtx} {output.auc_mtxRDS} 
        

        """



rule scenic_step3_aucell:
    input:
        regulons = join(DATAPATH, 'results/scrna/SCENIC/regulons.csv'),
        exprs    = join(DATAPATH, 'results/scrna/SCENIC/rnaseq_counts.tsv') 
    output:
        auc_mtx = join(DATAPATH, 'results/scrna/SCENIC/auc_mtx.csv')
    params:
        workdir = DATAPATH
    singularity:
        'docker://aertslab/pyscenic:0.10.4'
    shell:
        """
        
        pyscenic aucell {input.exprs} {input.regulons} \
        -o {output.auc_mtx} \
        --num_workers 6

        """
        
rule scenic_step2:
    input:
        adjacencies      = join(DATAPATH, 'results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv'),
        exprs            = join(DATAPATH, 'results/scrna/SCENIC/rnaseq_counts.tsv'), 
        motif_annotation = scenic_data['motif_annotation'],
        db_ranking_tss1  = scenic_data['db_ranking_tss1'],
        db_ranking_tss2  = scenic_data['db_ranking_tss2']
    output:
        regulons = join(DATAPATH, 'results/scrna/SCENIC/regulons.csv')
    params:
        workdir = DATAPATH,
        cores =  config['tfActity_int']['cores']
    singularity:
        'docker://aertslab/pyscenic:0.10.4'
    shell:
        """
        pyscenic ctx {input.adjacencies} {input.db_ranking_tss1} {input.db_ranking_tss2} \
        --annotations_fname {input.motif_annotation} \
        --expression_mtx_fname {input.exprs} \
        --mode "dask_multiprocessing" \
        --output {output.regulons} \
        --num_workers {params.cores}

        """

rule scenic_step1b_add_cor:
    input:
        adjacencies = join(DATAPATH, 'results/scrna/SCENIC/expr_mat.adjacencies.tsv'),
        exprs       = join(DATAPATH, 'results/scrna/SCENIC/rnaseq_counts.tsv') 
    output:
        adjacencies_withCor = join(DATAPATH, 'results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv')
    params:
        workdir = DATAPATH
    singularity:
        'docker://aertslab/pyscenic:0.10.4'
    shell:
        """
        pyscenic add_cor {input.adjacencies} {input.exprs} \
        --output {output.adjacencies_withCor} 

        """

rule scenic_step1_grn_arboreto:
    input:
        exprs = join(DATAPATH, 'results/scrna/SCENIC/rnaseq_counts.tsv') 
    output:
        adjacencies = join(DATAPATH, 'results/scrna/SCENIC/expr_mat.adjacencies.tsv')
    params:
        workdir = DATAPATH,
        TF_list = config['TF_list'],
        cores =  config['tfActity_int']['cores']
    singularity:
        'docker://aertslab/pyscenic:0.10.4'
    shell:
        """
        arboreto_with_multiprocessing.py -o {output.adjacencies} --num_workers {params.cores} {input.exprs} {params.TF_list}

        """

#------------------------------------------------------------------------------#
#                        Write matrix to use in SCENIC                         #
#------------------------------------------------------------------------------#
rule scenic_step0_seurat2scenic:
    input:
        seuratobj = seuratobj
    output:
        countsSCENIC = join(DATAPATH, 'results/scrna/SCENIC/rnaseq_counts.tsv') 
    params:
      script = 'scripts/scRNA-seq/01_seurat_to_SCENIC.R',
      hvg = config['GRNBoost']['HVG']
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        Rscript {params.script} {input.seuratobj} {output.countsSCENIC} \
        {params.hvg} 

        """

rule scenic_data:
    output:
        motif_annotation = scenic_data['motif_annotation'],
        db_ranking_tss1  = scenic_data['db_ranking_tss1'],
        db_ranking_tss2  = scenic_data['db_ranking_tss2']
    params:
        motif_annotationURL = scenic_data['motif_annotationu'],
        db_ranking_tss1URL  = scenic_data['db_ranking_tss1u'],
        db_ranking_tss2URL  = scenic_data['db_ranking_tss2u']
    shell:
        """
        cd aux/scenic/rcistarget/

        wget {params.motif_annotationURL}
        wget {params.db_ranking_tss1URL}
        wget {params.db_ranking_tss2URL}

        
        """
