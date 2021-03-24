
#------------------------------------------------------------------------------#
#               Create ArchR project and find motif positions                  #
#------------------------------------------------------------------------------#
rule atacMotifMatch:
    input:
        archrproj = join(DATAPATH, 'results/scatac/archr/ArchR01_transferLabels/Save-ArchR-Project.rds')
    output:
        archrprojout = join(DATAPATH, 'results/scatac/archr/ArchR02_MotifMatch/Save-ArchR-Project.rds')
    params:
        script = 'scripts/scATAC-seq/02_archR_motif_positions.R',
        nCores = config['tfActity_int']['cores']
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {input.archrproj} {output.archrprojout} {params.nCores}
        
        """
    


rule tranferRNA2atac:
    input:
        seuratobj  = seuratobj,
        arrowfiles = Arrowfiles.values()
    output:
        archrproj = join(DATAPATH, 'results/scatac/archr/ArchR01_transferLabels/Save-ArchR-Project.rds')
    params:
        script = 'scripts/scATAC-seq/01_create_archR_project.R',
        nCores = config['tfActity_int']['cores'],
        genome = config['genome'],
        seurat_annotCol = seurat_annotCol
    singularity:
        'docker://hdsu/r_scatac'
    shell:
        """
        
        Rscript {params.script} {output.archrproj} {input.seuratobj} \
        {params.seurat_annotCol} {params.genome} {params.nCores} \
        {input.arrowfiles}
        
        """
    
