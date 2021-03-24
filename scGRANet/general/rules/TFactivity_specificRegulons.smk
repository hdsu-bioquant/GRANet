#------------------------------------------------------------------------------#
#          Method 2: TF activity integrating scRNAseq and scATACseq data       #
#------------------------------------------------------------------------------#
rule tfm02_step4_RegulonStats:
    input:
        regulonAUC        = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulonAUC.RDS'),
        regulonAUC_devCor = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulonAUC_devCorrected.RDS'),
        regulonPySCENIC   = join(DATAPATH, 'results/scrna/SCENIC/auc_mtx.RDS')
    output:
        regulon_stats = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulon_stats.txt')
    params:
        script = 'scripts/tfActivity_int_method1/04_regulon_stats.R'
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.regulonAUC} {input.regulonAUC_devCor} \
        {input.regulonPySCENIC} {output.regulon_stats}
        
        """


rule tfm02_step3_RegulonActivity_NormMotifDeviation:
    input:
        archrproj  = join(DATAPATH, 'results/scatac/archr/ArchR02_MotifMatch/Save-ArchR-Project.rds'),
        exprs      = seuratobj,
        regulonAUC = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulonAUC.RDS')
    output:
        regulonAUC_devCor = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulonAUC_devCorrected.RDS')
    params:
        script = 'scripts/tfActivity_int_method2/03_RegulonActivity_NormMotifDeviation.R'
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.archrproj} {input.exprs} \
        {input.regulonAUC} {output.regulonAUC_devCor}
        
        """


rule tfm02_step2_AUCcell_intRegulonsByCellType:
    input:
        exprs    = seuratobj,
        #exprs    = join(DATAPATH, 'results/scrna/seurat/rna_seurat_int_annot_transfer.RDS'),
        regulons = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/tfRegulons_asDF.RDS')
    output:
        regulonAUC     = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/regulonAUC.RDS')
    params:
        script = 'scripts/tfActivity_int_method2/02_AUCcell_intRegulonsByCellType.R',
        nCores = config['tfActity_int']['cores'],
        seurat_annotCol = seurat_annotCol
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.exprs} {input.regulons} \
        {output.regulonAUC} {params.nCores} {params.seurat_annotCol}
        
        """

rule tfm02_step1_regulon_matchedMotifsByCellType:
    input:
        adjacencies = join(DATAPATH, 'results/scrna/SCENIC/expr_mat.adjacencies.cor.tsv'),
        archrproj   = join(DATAPATH, 'results/scatac/archr/ArchR02_MotifMatch/Save-ArchR-Project.rds')
    output:
        regulons = join(DATAPATH, 'results/integrated/TF_activity_method2_' + windowSize + '/tfRegulons_asDF.RDS')
    params:
        script  = 'scripts/tfActivity_int_method2/01_regulon_matchedMotifsByCellType.R',
        genome               = config['genome'],
        importance_threshold = config['tfActity_int']['regulons']['importance_threshold'],
        promoter_size        = config['tfActity_int']['regulons']['promoter_size'],
        min_regulon_size     = config['tfActity_int']['regulons']['min_regulon_size'],
        nCores               = config['tfActity_int']['cores'],
        cellWithPeak         = config['tfActity_int']['regulonsMethod2']['cellWithPeak']
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.adjacencies} {input.archrproj} \
        {output.regulons} \
        {params.genome} {params.importance_threshold} \
        {params.promoter_size} {params.min_regulon_size} \
        {params.nCores} {params.cellWithPeak}
        
        """
