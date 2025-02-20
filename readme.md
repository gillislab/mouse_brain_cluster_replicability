##  Cluster replicability in single-cell and single-nucleus atlases of the mouse brain

This repository contains the code used in the analysis for the study "Cluster Replicability in Single-Cell and Single-Nucleus Atlases of the Mouse Brain."

### Overview

This repository contains the code and analysis for "Cluster Replicability in Single-Cell and Single-Nucleus Atlases of the Mouse Brain." We assess the replicability of cell clusters in two large mouse brain atlases, each profiling over 4 million cells and 5000+ clusters. Using transcriptome-wide neighbor voting, we identify reciprocally matched clusters with consistent spatial localization and gene expression, conserved across species.

- Allen Institute Data: The AIT21_10Xv2, AIT21_10Xv3, and AIT21_multiome .h5ad files are not directly available online but the underlying data can be accessed through the Allen Brain Cell Atlas (https://github.com/AllenInstitute/abc_atlas_access).
- Other datasets: See the Data Files section below for a structured overview of additional datasets used in this project.

### System Requirements

- Several analyses require a **high-memory machine** (up to 2.6TB RAM).
- Code is a mix of **Python** and **R** notebooks.
- Python dependencies are listed in `project.env.yml`. The [GitHub version of `pyMN`](https://github.com/gillislab/pyMN) is required.

### Script order: 

- process_Zeng_AWS_files.ipynb: python notebook to combine the Allen AIT21.0 h5ad files (requires a high memory machine).  
- full_merge_and_HVG_write.ipynb: combines the single-cell and single-nucleus dataset and filters for highly variable genes (requires a high memory machine).  
- run_on_full_merged_ZengAWS.ipynb: calculate metaneighbor results (requires a high memory machine).  
- examine_region_overlap.ipynb: extracts the number of cells and clusters for each aligned major dissection region.  
- Figure_1_descriptive.R and Figure_1_markers.R: generates Figure 1 plots.  
- Macosko_all_marker_AUROCs.ipynb: calculates the AUROCs for the marker lists.  
- process_Zeng_AWS_files_all_genes_cpm.ipynb: creates CPM files for downstream tasks.  
- process_Macosko_all_genes_cpm.ipynb: creates CPM files for downstream tasks.  
- Marker_set_intersections.ipynb: calculates direct overlaps between the markers for the two datasets.  
- Macosko_all_marker_AUROCs.ipynb and Zeng_all_marker_AUROCs.ipynb: calculates AUROCs for all marker lists for all clusters (code duplication between these two, needs refactoring). These notebooks also contain commented-out code for marker subset creation to help with memory demands. Requires a high-memory machine.  
- split_for_metamarkers.ipynb: creates CPM files split into parts by genes due to memory constraints. Requires a high-memory machine. The number of splits can be adjusted if there are downstream problems with metamarkers.  
- Figure_1_markers.R: plots marker-based results and creates a cross-dataset matching based on markers (Zeng_markers_on_Macosko_one_to_one_mapping.csv).  
- make_centroids_mean_cpm_from_parts.R: creates the centroids (should be converted to Python).  
- Figure_1_centroid_expression.R: generates the gene-gene correlation plot and centroid correlation heatmap.  
- Figure_1_and_2_compare_other_mappings.R: compares to other integration efforts, generating an Euler diagram and bar graph.  
- Figure_2_metaneighbor_heatmaps.R: plots global and zoomed in heatmaps.  
- Figure_2_metaneighbor_plots.R: generates Supplement Figure 1 and other Figure 2 plots.  
- tophit_enrichment_tests.R: calculates higher-level enrichments of the reciprocal best hits (Supplement Tables 2 & 3).  
- centroid_correlations.bigcat.R: generates MERFISH cell calls using code from the scrattch.bigcat R package. Outputs are combined in Python using combine_MERFISH_calls_from_R.ipynb.  
- Figure_3_examine_MERFISH_calls_for_recips.R: processes Figure 3 plots, including the three clusters in Figure 5.  
- Table_for_high_confidence_filters.R: writes a table with the 612 high-confidence matches.  
- Table_annotation_for_tophits.R: annotates the top hit tables (2009 and 612 - Supplement Tables 1 & 8).  
- Figure_3_examine_MERFISH_recip_centroids.R: compares centroids in MERFISH coordinates.  
- Figure_4_coordinated_expression_in_recips_and_mean_exp.R: computes coordinated expression calculations.  
- Figure_4_coordinated_expression_GO.R: generates the GO plot in Figure 4 and Supplement Table 4.  
- make_BARseq_subsets.ipynb: creates the 109-gene subsets of the atlas data.  
- make_pretrained_using_full_merged.ipynb: creates MetaNeighbor pretrained models for both the highly variable gene set and the 109 BARseq gene panel. The highly variable gene set requires a large memory machine (>400GB).  
- run_on_barseq_pretrained.R: runs pretrained models against the BARseq clusters.  
- examine_barseq_pretrained_results.R: analyzes BARseq pretrained output and generates Supplement Table 5.  
- Figure_5_plot_BARseq_clusters.R: plots individual BARseq clusters.  
- meta_markers_in_parts.R: generates metamarkers for the two atlases.  
- meta_markers_in_reciprocals.R: calculates markers recurrent in both datasets based on Metaneighbor best reciprocal hits.  
- MetaMarkers_greedy_AUC_tests.ipynb: greedily tests markers for each cluster.  
- examine_greeedy_runs_meta_markers.R: generates Supplement Figure 1, counting the number of clusters where Metamarkers achieve a 0.95 best-vs-next AUROC threshold.  

### /data listing
Where possible, small data files (<100MB) have been included in the Git repository. Many of these files are supplements from other papers, and the directories are named accordingly. A tree listing of the data files used is:
```
├── [ 4.0K]  Allen_CCFv3
│   └── [  71K]  1-s2.0-S0092867420304025-mmc2.xlsx
├── [ 4.0K]  Carla_Winter_etal
│   └── [  23K]  Supp Table_5_He_AIBS_Mapping_Cluster_Summary.xlsx
├── [ 4.0K]  Hanqing_Liu_etal
│   └── [ 362K]  41586_2023_6805_MOESM7_ESM.csv
├── [ 4.0K]  Songpeng_Zu_etal
│   └── [ 4.0K]  SI Table 5 Transfer label scores for the integration of the snATAC-seq with the scRNA-seq data. Both cluster and subclass level scores are presented
│       ├── [ 167K]  Allen.scRNAseq.cl_subclassid_subclasslabel.csv
│       ├── [  16M]  sa2.all.cl2L4.transferLabelScore.csv
│       └── [ 1.3M]  sa2.all.subclass2L4.transferLabelScore.csv
├── [ 4.0K]  synGO
│   └── [ 4.0K]  release_20231201
│       └── [ 243K]  syngo_annotations.xlsx
└── [ 4.0K]  whole_mouse_brain
    ├── [ 4.0K]  macosko
    │   ├── [ 4.0K]  braincelldata.org
    │   │   └── [ 4.0K]  Mapping_matrices
    │   │       └── [  91K]  Puck_Num_01.mapping.MappedCellTypes.txt
    │   ├── [ 4.0K]  from_google_drive
    │   │   ├── [  95K]  CellType_MetaCluster.csv
    │   │   ├── [ 864K]  CellType_Metadata.tsv
    │   │   ├── [  41K]  Library_Metadata.tsv
    │   │   ├── [  31G]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.h5ad
    │   │   └── [ 155M]  Single_Nuc_Cluster_Avg_Expression.csv.gz
    │   └── [ 4.0K]  from_terra_non_uniform
    │       └── [ 913K]  CellType_Metadata.tsv
    ├── [ 4.0K]  processed
    │   ├── [ 4.0K]  barseq
    │   │   ├── [ 209M]  barseq_210630.rds
    │   │   ├── [ 2.7K]  barseq_gene_ids.csv
    │   │   ├── [  54M]  labels_20211208_withgaba_nonexc.csv
    │   │   └── [ 6.0K]  top_region_per_cluster.csv
    │   ├── [ 4.0K]  macosko
    │   │   ├── [ 155M]  merged_Zeng_AWS.Oct2023.pretrained_Macosko.csv.gz
    │   │   ├── [ 4.0K]  paper_supplements
    │   │   │   ├── [  17K]  41586_2023_6818_MOESM10_ESM.0_markers.txt
    │   │   │   └── [ 2.9M]  41586_2023_6818_MOESM10_ESM.xlsx
    │   │   └── [ 4.0K]  subsets
    │   │       ├── [ 4.8M]  Macosko_Mouse_Atlas_Single_Nuclei.pretrained_Macosko_109_genes.csv.gz
    │   │       ├── [  55G]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed_3309_markers_only.cpm.h5ad
    │   │       ├── [ 2.2G]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.barseq_109.h5ad
    │   │       ├── [  18G]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.cpm.part_0.h5ad
    │   │       ├── [  18G]  ........
    │   │       ├── [  18G]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.cpm.part_19.h5ad
    │   │       └── [ 598M]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.log2_cpm.mean_centroids.csv.gz
    │   ├── [ 5.6M]  Macosko_Mouse_Atlas_Single_Nuclei.Use_Backed.pretrained_barseq_128.csv.gz
    │   ├── [ 116G]  merged_Zeng_AWS.Oct2023.h5ad
    │   └── [ 4.0K]  zeng
    │       ├── [ 163M]  merged_Zeng_AWS.Oct2023.pretrained_Zeng.csv.gz
    │       └── [ 4.0K]  subsets
    │           ├── [  22G]  AIT21.0.merged.with_multiome_3820_markers_only.cpm.h5ad
    │           ├── [ 2.8G]  AIT21.0.merged.with_multiome.barseq_109.h5ad
    │           ├── [ 4.6G]  AIT21.0.merged.with_multiome.cpm.part_0.h5ad
    │           ├── [ 4.6G]  ........
    │           ├── [ 2.8G]  AIT21.0.merged.with_multiome.cpm.part_39.h5ad
    │           ├── [ 266G]  AIT21.0.merged.with_multiome.h5ad
    │           ├── [ 657M]  AIT21.0.merged.with_multiome.log2_cpm.mean_centroids.csv.gz
    │           └── [ 5.1M]  AIT21.0.merged.with_multiome.pretrained_Zeng_109_genes.csv.gz
    └── [ 4.0K]  zeng
        ├── [ 4.0K]  from_API
        │   ├── [ 123M]  ccf_coordinates_MERFISH-C57BL6J-638850.csv.gz
        │   └── [ 726K]  Cluster_colors_and_tree.abc_atlas_access.csv
        ├── [ 4.0K]  from_aws
        │   └── [ 4.0K]  AIT21.0
        │       ├── [  18K]  AIT21_annotation_freeze_081523.all_markers.txt
        │       ├── [ 155G]  AIT21_10Xv2.h5ad
        │       ├── [ 263G]  AIT21_10Xv3.h5ad
        │       ├── [ 206M]  AIT21_multiome.h5ad
        │       └── [ 1.8M]  AIT21_annotation_freeze_081523.tsv
        └── [ 4.0K]  MERFISH-C57BL6J-638850
            ├── [ 4.0K]  20230630
            │   └── [ 7.1G]  C57BL6J-638850-raw.h5ad
            └── [ 4.0K]  20230830
                ├── [ 538M]  cell_metadata.csv
                └── [  47K]  gene.csv
```