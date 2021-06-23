# ChesapeakeSulfurOxidation
Scripts associated with manuscript submitted to Environmental Microbiology 

These are some description of the scripts included. 
These are provided to show the methods used. 
The data required to run these are not included in this repo and are all in the public domain.
If you need assistance, please open an issue and I will respond. 


| Script Name  | General Purpose |
| ------------- | ------------- |
| 0_checkm_bins.sh  | Calculate Completeness/Contamination for Bins  |
| 2_gtdbtk_classify.sh  | Classify lineage with GTDB tk  |
| mapHELPER.sh + 4_readrecruitment.sh | Quantify coverage of bins |
| cb33_2017_qc.sh | QC shotgun libraries |
| assem_skel.sh | Assemble shotgun libraries |
| minimus_runner_CB33.sh | Deduplicate Contigs + Run Overlap Layout Consensus Assembly |
| cb33_binNquant.sh | Bin Genomes | 
| demux_pipeline.sh + check_headers.py + rev_comp_bcodes.py | Demultiplex amplicons |
| filter_pipeline.sh + FilterNTrim.R | QC amplicons |
| mergeRunandTaxa.sh + merge_chim_tax.R | Merge Amplicon samples, remove chimaeras, assign taxonomy |
| callOTUs_pipe.sh + fullDADApipe_PE.R | Full pipeline of amplicon analysis |
| make_tree_any_abund.sh | Generate 16S Amplicon Phylogeny | 
| dnld_hewson.sh | Download RNA seq data | 
| salmon_quant.sh | QC + Map RNA seq data |
| 0_DADA2_PostProcessing.ipynb | Pull sequences out of DADA outputs for making phylogeny |
| 1_Environmental_Data_Processing.ipynb | Aggregate environmental data |
| 2_Habitats.ipynb | Cluster Environmental data |
| 3_ParticleTransportAnalysis_and_MapFig.ipynb | Analysis of Particle data (not used in publication) + Map Figure Generation (used in publication) |
| 5_RDA_Analysis.ipynb | Cluster Microbial data |
| Predict_Hypoxia.ipynb | Used for recursive feature elimination |	 
| Pathway Analysis.ipynb + pathwayFunctions.py | Used to aggregate and analyze pathway data |
| pathwayAnalysis.py | Script version of Notebook of same name |








