####################
The collected files and scripts all relate to the manuscript as detailed below:


Title: Gene expression levels tune germline mutation rates through the compound effects of transcription-coupled repair and damage
Authors: Bo Xia1, and Itai Yanai1,2
Affiliations:
1 Institute for Computational Medicine, NYU Langone Health, New York, NY 10016, USA
2 Department of Biochemistry and Molecular Pharmacology, NYU Langone Health, New York, NY 10016, USA



####################
CONTENTS:


######
Folders:
\Data              : Relevant files related the scRNA-seq and germline SNV analysis.
\MATLAB_codes	   : Scripts
\SNV		   : germline variants extracted from the public database.


###
\DATA\
	GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat
	h_m_stage_expression.mat
	human_ensembl90_gene_features.txt
	human_genenames.tsv



###
\MATLAB_codes\
Main scripts:
	TS_REanalysis.m

Related functions:	
	F_compare_snv_asymmetry_v2_20200724.m
	F_compare_variants_20200724.m
	
Third party functions:
	\cbrewer
	cbrewer.m
	barwitherr.m	
	scatter_kde.m
	shadedErrorBar.m
	Violin.m
	violinplot.m
	


###
\SNV_Base_Frequency\
	human_1000G_genebody_12_result.12.stranded.tsv	
	human_ensembl90.genebody.singleFreq.tsv
	human_ensembl90.intron.singleFreq.tsv
	human_gene_body_result.12.stranded.tsv
	human_intron_12_result.12.stranded.tsv
	An_2018_SNV_control_result.12.stranded.tsv
	Jonsson_2017_Nature_SNV_all_result.12.stranded.tsv

	

####################
Contact:
Bo Xia:     Bo.Xia@nyulangone.org; xiabo90@gmail.com
Itai Yanai: Itai.Yanai@nyulangone.org
2021-02-01



