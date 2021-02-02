%% Reproduce Liu et al correlation analysis.
Liu = readtable('\Liu_and_Zhang_data\mutation_expressoin_data.xlsx');
[rho,pval] = partialcorr(Liu.devono_mutation_rate,Liu.testis_expression_level,...
                        [Liu.replicatoin_timing Liu.GC_content Liu.nucleosome_occupancy],...
                        'Type','Spearman')
[rho,pval] = corr(Liu.devono_mutation_rate,Liu.testis_expression_level,...
                        'Type','Spearman')



%% Liu and Zhang used raw gene expression datasets which includes ~70% empty barcodes. DONE!
%Inport Raw UMI matrices: Human TESE samples (Raw UMI files downloadable from NCBI-GEO: GSE125372)
%     fn1='\Analysis\2017-11-03_human_TESE\Original file\gene_expression.tsv'
%     fn2='\Analysis\2017-11-03_human_TESE\Original file\2018-08-17_human_TESE_techrep\gene_expression.tsv'
%     fn3='\Analysis\2018-10-30_human_TESE2\raw_data\human_TESE2_1st\gene_expression.tsv'
%     fn4='\Analysis\2018-10-30_human_TESE2\raw_data\human_TESE2_2nd\gene_expression.tsv'
    a1 = importdata(fn1);
    a2 = importdata(fn2);
    a3 = importdata(fn3);
    a4 = importdata(fn4);
    raw_exp = mean([a1.data a2.data a3.data a4.data]')';
    raw_names = a1.textdata(2:end,1);
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(raw_names,cellstr(human_gene_info.name),'stable');
    raw_names  = cellstr(human_gene_info.ensembl_id(xi,:));
% Liu and Zhang gene expression level
    fn_Liu='\Liu_and_Zhang_data\mutation_expressoin_data.xlsx'
    a5_Liu = importdata(fn_Liu);
    Liu_exp_denovo = a5_Liu.data.denovo_mutation(:,5);
    Liu_genenames_denovo = a5_Liu.textdata.denovo_mutation(2:end,1);
    Liu_exp_intron = a5_Liu.data.intron_SNP_density(:,5);
    Liu_genenames_intron = a5_Liu.textdata.intron_SNP_density(2:end,1);
% Check whether Liu and Zhang used raw matrices
    [~,xi,xj] = intersect(raw_names,Liu_genenames_denovo,'stable');
    [~,yi,yj] = intersect(raw_names,Liu_genenames_intron,'stable');
    if sum(raw_exp(xi)-Liu_exp_denovo(xj)) ==0 && sum(raw_exp(yi)-Liu_exp_intron(yj)) ==0 
        fprintf('Liu and Zhang wrongly used raw scRNA-seq matrices which includes:\n')
        fprintf('\t 5446 empty barcodes\n')
        fprintf('\t 460 non-germ cells(e.g. Leydig cells, PMC...)\n')
        fprintf('\t and 2094 male germ cells.\n')
        fprintf('\t SHOULD ONLY INCLUDE FILTERED MALE GERM CELLS!!!\n')
    end 
        

    
    
   

                    
%% human_variation_by_expressionlevel

if human_variation_by_expressionlevel
    %% LOAD data
    %######NEED TO ADD CODES OF HOW TO GENERATE THE FILES BELOW
    % human variants
load('YOUR_PATH\Data\h_m_stage_expression.mat');

    cd 'YOUR_PATH\MATLAB_codes\cbrewer'
    %cbrewer()
    OrRd=cbrewer('seq','OrRd',12);
    explevel_edgecm = [0.15 0.15 0.15;OrRd(4:11,:)];
    explevel_facecm = [0.5 0.5 0.5;OrRd(4:11,:)];
    genecluster_edgecm = [0.15 0.15 0.15;cbrewer('qual','Dark2',5)];
    genecluster_facecm = [0.5 0.5 0.5;cbrewer('qual','Set2',5)];
    cd 'YOUR_PATH\'
    
%% set up the expression level gene cluster
    
    h_mean_exp = log2(0.0001+mean(h_cluster_exp_meanUMI')');
    h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        h_expressed_mean = h_mean_exp(find(h_gene_cluster>0));
        [f,xi] = ksdensity(h_expressed_mean);
        figure;plot(xi,f,'LineWidth',1.5);set(gca,'Color','none');box off;
        ylabel('Kernel gene density');xlabel('Expression level(log2)');
        xlim([-12 10]);xticks([-10 -8 -6 -4 -2 0 2 4 6 8]);
        
    h_expressed_idx  = find(h_gene_cluster>0);
    h_explevel_cluster(1:length(h_gene_cluster),1) = 0;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-8) )) = 1;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-6 & h_expressed_mean>-8) )) = 2;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-4 & h_expressed_mean>-6) )) = 3;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=-2 & h_expressed_mean>-4) )) = 4;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=0 & h_expressed_mean>-2) )) = 5;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=2 & h_expressed_mean>0) )) = 6;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean<=4 & h_expressed_mean>2) )) = 7;
    h_explevel_cluster( h_expressed_idx( find(h_expressed_mean>4) )) = 8;
    histcounts(h_explevel_cluster)
    
    %define gene cluster composation in each expression level cluster
    for i=1:8
        temp = h_gene_cluster(find(h_explevel_cluster == i));
        percentage(i,:) = histcounts(temp)./length(temp);
    end
    figure;bar(percentage,'stacked');colormap(genecluster_facecm(2:6,:));
    set(gca,'Color','none');xlim([0 9]);
    


%% mutation assymetry by expression level
    load('YOUR_PATH\Data\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    human_geneID  = cellstr(human_gene_info.ensembl_id(xi,:));

%human intron: explevel cluster plot overall
    fn1='YOUR_PATH\SNV_Base_Frequency\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron_explevel'); 

% human_dbSNP_variants_in_total
%human gene body
    fn1='YOUR_PATH\SNV_Base_Frequency\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','no',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-genebody-explevel'); 
    
% An_ctl_&_Jonsson_all DNM rates across gene clusters
prot_idx = strmatch('protein_coding',cellstr(human_gene_info.type));
autosome={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
auto = ismember(cellstr(human_gene_info.chromosome(prot_idx,:)),autosome);
auto_prot_genes = cellstr(human_gene_info.ensembl_id(prot_idx(auto),:));
[~,~,j]=intersect(auto_prot_genes,human_geneID,'stable');
fn1="YOUR_PATH\SNV_Base_Frequency\Jonsson_2017_Nature_SNV_all_result.12.stranded.tsv"
fn2="YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv"
fn3="YOUR_PATH\SNV_Base_Frequency\An_2018_SNV_control_result.12.stranded.tsv"
    %by gene clusters
    An_Jonsson_all = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster(j),...
                                     human_geneID(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','yes',...
                                     'Include_meta_all','no',...
                                     'Filter_zero_mut_gene','no',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 
    An_Jonsson_all_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster(j),...
                                     human_geneID(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','yes',...
                                     'Include_meta_all','no',...
                                     'Filter_zero_mut_gene','yes',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 
end






%% Plot the frequencies of mutations per gene w/t or w/o zeroVar genes. DONE.
    load('YOUR_PATH\Data\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    human_geneID  = cellstr(human_gene_info.ensembl_id(xi,:));
% human_dbSNP_variants_in_total
%human gene body
    fn1='YOUR_PATH\SNV_Base_Frequency\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv'

    human_genebody = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','no',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-dbSNP-genecluster'); 
     human_genebody_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-dbSNP-genecluster'); 
%human gene body(1000GenomeSNV)
    fn1='YOUR_PATH\SNV_Base_Frequency\human_1000G_genebody_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody_1000G = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','only',...
                                     'Filter_zero_mut_gene','no',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-1000G'); 
     human_genebody_1000G_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','only',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-1000G'); 
% An_ctl_&_Jonsson_all DNM rates across gene clusters
prot_idx = strmatch('protein_coding',cellstr(human_gene_info.type));
autosome={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
auto = ismember(cellstr(human_gene_info.chromosome(prot_idx,:)),autosome);
auto_prot_genes = cellstr(human_gene_info.ensembl_id(prot_idx(auto),:));
[~,~,j]=intersect(auto_prot_genes,human_geneID,'stable');
fn1="YOUR_PATH\SNV_Base_Frequency\Jonsson_2017_Nature_SNV_all_result.12.stranded.tsv"
fn2="YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv"
fn3="YOUR_PATH\SNV_Base_Frequency\An_2018_SNV_control_result.12.stranded.tsv"
    %by gene clusters
    An_Jonsson_all = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster(j),...
                                     human_geneID(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','only',...
                                     'Include_meta_all','no',...
                                     'Filter_zero_mut_gene','no',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 
    An_Jonsson_all_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster(j),...
                                     human_geneID(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','only',...
                                     'Include_meta_all','no',...
                                     'Filter_zero_mut_gene','yes',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 

dbSNP_varpergene = log10(human_genebody.variants_level_structure.var_per_kb+1);
dbSNP_varpergene_filter0 = log10(human_genebody_filter_zero.variants_level_structure.var_per_kb+1);
G1000_varpergene = log10(human_genebody_1000G.variants_level_structure.var_per_kb+0.1);
G1000_varpergene_filter0 = log10(human_genebody_1000G_filter_zero.variants_level_structure.var_per_kb+0.1);
DNM_varpergene   = log10(An_Jonsson_all.variants_level_structure.var_per_kb+0.01);
DNM_varpergene_filter0   = log10(An_Jonsson_all_filter_zero.variants_level_structure.var_per_kb+0.01);

        figure;  
        subplot(2,3,1);histogram(dbSNP_varpergene,100); xlim([-0.5,3.5]);ylim([0 4000]);
        ylabel('Top: keep 0V genes');set(gca,'Color','none');
        subplot(2,3,4);histogram(dbSNP_varpergene_filter0,100); xlim([-0.5,3.5]);ylim([0 4000])
        ylabel('Bottom: remove 0V genes');xlabel('dbSNP150: log10(MutRate+1)');set(gca,'Color','none');
        subplot(2,3,2);histogram(G1000_varpergene,100); 
        xlim([-1.5,2.5]);ylim([0 4000]);set(gca,'Color','none');
        subplot(2,3,5);histogram(G1000_varpergene_filter0,100); 
        xlim([-1.5,2.5]);ylim([0 4000]);xlabel('1000Genome: log10(MutRate+0.1)');set(gca,'Color','none');
        subplot(2,3,3);histogram(DNM_varpergene,100); 
        xlim([-3,0.6]); ylim([0 5000]);set(gca,'Color','none');
        subplot(2,3,6);histogram(DNM_varpergene_filter0,100);
        xlim([-3,0.6]); ylim([0 5000]);xlabel('De novo mutations: log10(MutRate+0.01)');set(gca,'Color','none');
        suptitle('Variants per kb (log10) calculated from dbSNP150,1000G,DNM')
        % Gene number:  dbSNP: 19802; dbSNP_remove_0V: 19767; 
        %               1000G: 19802; 1000G_remove_0V: 19539; 
        %               DNM:   18905; DNM_remove_0V:   14131; 
        saveas(gcf,'YOUR_PATH\Results\Sup_fig_variants_per_gene_perkb_distribution.fig')


    
        
        
                           
%% Modelling DNM datasets by downsampling dbSNP150 variants. DONE.
load('YOUR_PATH\Data\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    human_geneID  = cellstr(human_gene_info.ensembl_id(xi,:));
% human_dbSNP_variants_in_total
%human gene body
    fn1='YOUR_PATH\SNV_Base_Frequency\human_gene_body_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv'
     human_genebody = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','only',...
                                     'Filter_zero_mut_gene','no',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-dbSNP-genecluster'); 
     human_genebody_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','only',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_edgecm,...
                                     'Title_prefix','Human-genebody-dbSNP-genecluster'); 

An_Jonsson_variant = An_Jonsson_all.variants_level_structure.var;
dbSNP_variant = human_genebody.variants_level_structure.var;
round(sum(dbSNP_variant)./sum(An_Jonsson_variant))   
fprintf('dbSNP150 has 1450 folds more genic variants than the DNM dataset. \n')


% 1. modelling mutation rate by downsampling.
    %First downgrade the variants to 1/145
mod_mutrate_gb = round(human_genebody.variants_level_structure.var/145);
sum(mod_mutrate_gb)
% 2. split by the indeces.
splited_gb = [];
for i=1:length(mod_mutrate_gb)
    splited_gb = [splited_gb ones([1 mod_mutrate_gb(i)])*i];
end

if sum(mod_mutrate_gb - histcounts(splited_gb,length(mod_mutrate_gb)) == 0)
    fprintf("Yes! Splitting the variants correctly.\n")
else
    fprintf("No! Check out inputs.\n")
end
% 3. random sample the array and then counting frequencies
downsampled_gb = randsample(splited_gb,round(length(splited_gb)/10));
downsamp_variant_gb = histcounts(downsampled_gb,length(mod_mutrate_gb))';

downsamp_struc_removezero = F_compare_variants_20200724(downsamp_variant_gb,...
                            human_genebody.variants_level_structure.glength,...
                            human_genebody.variants_level_structure.gname,...
                            human_genebody.variants_level_structure.gclu,...
                            human_genebody.variants_level_structure.gname,...
                            explevel_facecm,...
                            'Filter_zero_mut_gene','yes');
                           title('Genebody variant modeled rates (Filter Zero-mut genes)');
                           ylim([0 0.5])
downsamp_struc            = F_compare_variants_20200724(downsamp_variant_gb,...
                            human_genebody.variants_level_structure.glength,...
                            human_genebody.variants_level_structure.gname,...
                            human_genebody.variants_level_structure.gclu,...
                            human_genebody.variants_level_structure.gname,...
                            explevel_facecm,...
                            'Filter_zero_mut_gene','no');
                           title('Genebody variant modeled rates (W/o filtering Zero-mut genes)');
                           ylim([0 0.5])
                           
% 4. Downsampling to scatter the mutationrate with original rates.

    [~,i,j] = intersect(downsamp_struc.gname,human_genebody.variants_level_structure.gname,'stable');
    sum(i-j)  %--> equals 0, indicating matched gene coordinates.
    corr_downsample = [];pval =[];
    [corr_downsample(1) pval(1)] = corr(downsamp_struc.var_per_kb(i), ...
                    human_genebody.variants_level_structure.var_per_kb(j),'Type','spearman');
    [~,i,j] = intersect(downsamp_struc_removezero.gname,human_genebody.variants_level_structure.gname,'stable');
    [corr_downsample(2) pval(2)] = corr(downsamp_struc_removezero.var_per_kb(i),...
                    human_genebody.variants_level_structure.var_per_kb(j),'Type','spearman')

    
    temp = human_genebody.variants_level_structure.var_per_kb;
    figure;
    scatter_kde(log10(downsamp_struc.var_per_kb + 0.0001),log10(temp + 0.0001),'.');
    jet_custom = jet; colormap(jet_custom(50:240,:));
    xlabel('log10(downsampled-dbSNP mutrate + 0.0001)')
    ylabel('log10(dbSNP mutrate + 0.0001)')
    ylim([-4.5 3.5]);     xlim([-4.5 1.5]);      box off;     set(gca,'color','none');
    title('Rep.Vis. of downsampled mutrate vs original dbSNP mutrate')
    saveas(gcf,'YOUR_PATH\Results\Fig2A_downsample_scattering.fig')
    
% 5. random sampling, counting variant frequencies and calculate corrcoef for 100 times.
    gl = human_genebody.variants_level_structure.glength;
    [~,xi,xj] = intersect(human_geneID,human_genebody.variants_level_structure.gname,'stable');
    h_mean_exp = log2(0.0001+mean(h_cluster_exp_meanUMI(xi,:)')');
    corrcoef(human_genebody.variants_level_structure.var./gl, h_mean_exp)

    gl = An_Jonsson_all.variants_level_structure.glength;
    [~,xi,xj] = intersect(human_geneID,An_Jonsson_all.variants_level_structure.gname,'stable');
    h_mean_exp = log2(0.0001+mean(h_cluster_exp_meanUMI(xi,:)')');
    corrcoef(An_Jonsson_all.variants_level_structure.var./gl, h_mean_exp)

    [~,xi,xj] = intersect(human_genebody.variants_level_structure.gname,An_Jonsson_all.variants_level_structure.gname,'stable');
    corr(human_genebody.variants_level_structure.var_per_kb(xi),...
         An_Jonsson_all.variants_level_structure.var_per_kb(xj),'Type','spearman')
    a = find(An_Jonsson_all.variants_level_structure.var_per_kb >0);
    [~,xi,xj] = intersect(human_genebody.variants_level_structure.gname,An_Jonsson_all.variants_level_structure.gname(a),'stable');
    corr(human_genebody.variants_level_structure.var_per_kb(xi),...
         An_Jonsson_all.variants_level_structure.var_per_kb(xj),'Type','spearman')



    %Repeat 100 times of the downsampling experiments.
    gl = human_genebody.variants_level_structure.glength;
    dbSNP_gb_var = human_genebody.variants_level_structure.var./human_genebody.variants_level_structure.glength;
    original_vs_downsample = [];
    for i = 1:100
        downsampled_gb = randsample(splited_gb,round(length(splited_gb)/10));
        downsamp_variant_gb = histcounts(downsampled_gb,length(mod_mutrate_gb))';
        downsampled_varrate = downsamp_variant_gb./gl;
        %Spearman correlation between original variants and downsampled variants
        original_vs_downsample(i,1) = corr(dbSNP_gb_var, downsampled_varrate,'Type','spearman');
        %Spearman correlation between original variants and downsampled variants after removing zero-variant genes
        remove_0v = find(downsamp_variant_gb>0);
        original_vs_downsample(i,2) = corr(dbSNP_gb_var(remove_0v), downsampled_varrate(remove_0v),'Type','spearman');
    end
    %Plot the Spearman correlation coefficient rho between the original
    %dnSNP150 varant rates and the downsampled variant rates. 
        figure;
        boxplot(original_vs_downsample,'Symbol','','Notch','on');
        xlim([0.5 2.5]);ylim([0 0.55]);
        xticklabels({'Downsample w/t zeroVar','Downsample w/o zeroVar'})
        xlabel('dbSNP150 var dowmsampling rho')
        ylabel('rho: original vs dowmsampled ')
        box off;     set(gca,'color','none');
mean(original_vs_downsample)

        
        
        
 

    
%% Reproduce corr_coef analysis across all genes. Used intronic_SNP and DNM datasets. DONE. 
[~,i,j] = intersect(human_geneID,human_intron_explevel.variants_level_structure.gname,'stable');
mut_intron = human_intron_explevel.variants_level_structure.var_per_kb(j);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_intron,'stable');
rho=[];pval =[];
[rho(1,1),pval(1,1)] = corr(mut_intron(i),gexp(i),'Type','Spearman')
[rho(2,1),pval(2,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.intron_SNP_density(j,2),'Type','Spearman')
[rho(3,1),pval(3,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.intron_SNP_density(j,3),'Type','Spearman')
[rho(4,1),pval(4,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.intron_SNP_density(j,2:3),'Type','Spearman')
           
figure; bar(rho);set(gca,'color','none');box off; ylim([-0.12 0.06]);
xticklabels({'w/o Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});ylabel('Spearman rho, using intronic SNP')
saveas(gcf,'YOUR_PATH\Results\Fig_Spearman_rho_HumanIntronSNP.fig')


%Reproduce corr_coef analysis using human_DNM_filter_zeroVargene data
[~,i,j] = intersect(human_geneID,An_Jonsson_all_filter_zero.variants_level_structure.gname,'stable');
mut_intron = An_Jonsson_all_filter_zero.variants_level_structure.var_per_kb(j);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_denovo,'stable');
rho=[];pval =[];
[rho(1,1),pval(1,1)] = corr(mut_intron(i),gexp(i),'Type','Spearman')
[rho(2,1),pval(2,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.denovo_mutation(j,2),'Type','Spearman')
[rho(3,1),pval(3,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.denovo_mutation(j,3),'Type','Spearman')
[rho(4,1),pval(4,1)] = partialcorr(mut_intron(i),gexp(i),a5_Liu.data.denovo_mutation(j,2:3),'Type','Spearman')
                  
figure; bar(rho);set(gca,'color','none');box off;ylim([-0.12 0.06]);
xticklabels({'w/o Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});ylabel('Spearman rho, using human DNM')
saveas(gcf,'YOUR_PATH\Results\Fig_Spearman_rho_HumanDNM.fig')

   




%% GC content analysis. DONE. 
    load('YOUR_PATH\Data\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    human_geneID  = cellstr(human_gene_info.ensembl_id(xi,:));
    %human intron: explevel cluster plot overall
    fn1='YOUR_PATH\SNV_Base_Frequency\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron_explevel'); 
%(The output human_intron_explevel contrains the mutation rate by type.)
input_temp = human_intron_explevel;
GC_ratio = (input_temp.vglength(:,3)+ input_temp.vglength(:,4))./ input_temp.vglength(:,1);
temp = char(human_intron_192.tribase_names);
CpG_ratio = sum(human_intron_192.base_frequency(:,find(temp(:,2) == 'C' & temp(:,3) == 'G'))')' ./ human_intron_192.base_frequency(:,1);

    jet_custom = jet; jet_custom =jet_custom(50:240,:);
figure;
for i=1:12
    subplot(3,4,i);
    scatter_kde(GC_ratio,input_temp.vnum_perkb(:,i),'.');
	colormap(jet_custom);
    [GC_rho(i),GC_pval(i)] =corr(GC_ratio,input_temp.vnum_perkb(:,i),'Type','Spearman');
    [rho,pval] = partialcorr(GC_ratio,input_temp.vnum_perkb(:,i),CpG_ratio,'Type','Spearman')
    title(strcat(string(input_temp.vtype(i)),', rho:', num2str(GC_rho(i)),', pval:', num2str(GC_pval(i))))
    set(gca,'Color','none'); legend on;
end
suptitle('GC ratio-mutation rate correlation: intronic SNP')
saveas(gcf,'YOUR_PATH\Results\Sup_fig_GC ratio-mutation rate correlation-intronic SNP.fig')



%A>B/T>V type mutation rates does not correlate with GC ratio.)
input_temp = human_intron_explevel;
GC_ratio = (input_temp.vglength(:,3)+ input_temp.vglength(:,4))./ input_temp.vglength(:,1);
AT_mutrate_perkb = (sum(input_temp.vnum_pergene(:,1:6)')' ./ (input_temp.vglength(:,2) + input_temp.vglength(:,5))) *1000;
GC_mutrate_perkb = (sum(input_temp.vnum_pergene(:,7:12)')' ./ (input_temp.vglength(:,3) + input_temp.vglength(:,4))) *1000;
GC_ratio = (input_temp.vglength(:,3)+ input_temp.vglength(:,4))./ input_temp.vglength(:,1);
temp = char(human_intron_192.tribase_names);
CpG_ratio = sum(human_intron_192.base_frequency(:,find(temp(:,2) == 'C' & temp(:,3) == 'G'))')' ./ human_intron_192.base_frequency(:,1);
CpH_ratio = sum(human_intron_192.base_frequency(:,find(temp(:,2) == 'C' & temp(:,3) ~= 'G'))')' ./ human_intron_192.base_frequency(:,1);
    %Compute the correlation coefficients of mutrate to GC content
    [rho,pval] =corr(AT_mutrate_perkb,GC_ratio,'Type','Spearman')
    [rho,pval] =partialcorr(AT_mutrate_perkb,GC_ratio,CpG_ratio,'Type','Spearman')
    [rho,pval] =corr(GC_mutrate_perkb,GC_ratio,'Type','Spearman')
    [rho,pval] =partialcorr(GC_mutrate_perkb,GC_ratio,CpG_ratio,'Type','Spearman')

    [~,xi,xj] = intersect(human_geneID,input_temp.variants_level_structure.gname,'stable');
    [rho,pval] =corr(AT_mutrate_perkb(xj),h_mean_exp(xi),'Type','Spearman')
    [rho,pval] =partialcorr(AT_mutrate_perkb(xj),h_mean_exp(xi),CpG_ratio(xj),'Type','Spearman')
    [rho,pval] =corr(GC_mutrate_perkb(xj),h_mean_exp(xi),'Type','Spearman')
    [rho,pval] =partialcorr(GC_mutrate_perkb(xj),h_mean_exp(xi),CpG_ratio(xj),'Type','Spearman')

figure;
subplot(1,2,1); 
scatter_kde(GC_ratio,AT_mutrate_perkb,'.');colormap(jet_custom(50:240,:));
[rho,pval] =corr(GC_ratio,AT_mutrate_perkb,'Type','Spearman');
xlabel('Intronic GC ratio'); ylabel('Intronic A>B/T>V mut-rate')
title(strcat('A>B/T>V',', rho:', num2str(rho),', pval:', num2str(pval)));
set(gca,'Color','none');ylim([0 600])
subplot(1,2,2); 
scatter_kde(GC_ratio,GC_mutrate_perkb,'.');colormap(jet_custom(50:240,:));
[rho,pval] =corr(GC_ratio,GC_mutrate_perkb,'Type','Spearman');
xlabel('Intronic GC ratio'); ylabel('Intronic C>D/G>H mut-rate')
title(strcat('C>D/G>H',', rho:', num2str(rho),', pval:', num2str(pval)));
set(gca,'Color','none');ylim([0 600])
saveas(gcf,'YOUR_PATH\Results\GC ratio-AT-GC mutation types correlation-intronic SNP.fig')



%calculate CpG ratio to GC_rate: Very high!
input_temp = human_intron_explevel;
GC_ratio = (input_temp.vglength(:,3)+ input_temp.vglength(:,4))./ input_temp.vglength(:,1);
temp = char(human_intron_192.tribase_names);
CpG_ratio = sum(human_intron_192.base_frequency(:,find(temp(:,2) == 'C' & temp(:,3) == 'G'))')' ./ human_intron_192.base_frequency(:,1)*1000;
CpH_ratio = sum(human_intron_192.base_frequency(:,find(temp(:,2) == 'C' & temp(:,3) ~= 'G'))')' ./ human_intron_192.base_frequency(:,1)*1000;
[rho,pval] =corr(GC_ratio,CpG_ratio,'Type','Spearman')
[rho,pval] =corr(GC_ratio,log(CpG_ratio),'Type','Spearman')
figure; scatter_kde(GC_ratio,log10(CpG_ratio),'.');
    jet_custom = jet; jet_custom =jet_custom(50:240,:);
colormap(jet_custom);
set(gca,'Color','none');
xlabel('Intronic GC ratio'); ylabel('log10(Intronic CpG occurance per kb)')
title(strcat('GC ratio - CpG rate',', rho:', num2str(rho),', pval:', num2str(pval)))







%% Split expression level to two major groups: low-mid & high. DONE!
    fn_Liu='Liu_and_Zhang_data\mutation_expressoin_data.xlsx'
    a5_Liu = importdata(fn_Liu);
[~,i,j] = intersect(human_geneID,human_intron_explevel.variants_level_structure.gname,'stable');
mut_intron = human_intron_explevel.variants_level_structure.var_per_kb(j);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_intron,'stable');
gexp  = gexp(i);
mut_intron = mut_intron(i);
a5_Liu_data = a5_Liu.data.intron_SNP_density(j,:);


%Split genes to low-to-moderate_exp (70%) and high_exp (30%) genes
percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho =[];pval=[];
[rho(1,1),pval(1,1)] = corr(mut_intron(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,2),'Type','Spearman');
[rho(1,3),pval(1,3)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,3),'Type','Spearman');
[rho(1,4),pval(1,4)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,2:3),'Type','Spearman');
[rho(2,1),pval(2,1)] = corr(mut_intron(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,2),'Type','Spearman');
[rho(2,3),pval(2,3)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,3),'Type','Spearman');
[rho(2,4),pval(2,4)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,2:3),'Type','Spearman');
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using intronic SNP');ylim([-0.15 0.1])
title('Spearman rho: Split genes into low-mid & high expressed')




%Split with human_DNM datasets
[~,i,j] = intersect(human_geneID,An_Jonsson_all_filter_zero.variants_level_structure.gname,'stable');
mut_intron = An_Jonsson_all_filter_zero.variants_level_structure.var_per_kb(j);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_denovo,'stable');
gexp  = gexp(i);
mut_intron = mut_intron(i);
a5_Liu_data = a5_Liu.data.denovo_mutation(j,:);

percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho =[];pval=[];
[rho(1,1),pval(1,1)] = corr(mut_intron(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,2),'Type','Spearman');
[rho(1,3),pval(1,3)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,3),'Type','Spearman');
[rho(1,4),pval(1,4)] = partialcorr(mut_intron(i),gexp(i),a5_Liu_data(i,2:3),'Type','Spearman');
[rho(2,1),pval(2,1)] = corr(mut_intron(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,2),'Type','Spearman');
[rho(2,3),pval(2,3)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,3),'Type','Spearman');
[rho(2,4),pval(2,4)] = partialcorr(mut_intron(j),gexp(j),a5_Liu_data(j,2:3),'Type','Spearman');
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho(1:2,:));set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using DNM');ylim([-0.15 0.1])
title('Spearman rho: Split genes into low-mid & high expressed')







%% Coding strand and template strand mutation rates, and correlations with expression level
    load('YOUR_PATH\Data\h_m_stage_expression.mat');
    human_gene_info = struct2table(tdfread('YOUR_PATH\Data\human_genenames.tsv'));
    [~,~,xi] = intersect(human_genename,cellstr(human_gene_info.name),'stable');
    human_geneID  = cellstr(human_gene_info.ensembl_id(xi,:));
    %human intron: explevel cluster plot overall
    fn1='YOUR_PATH\SNV_Base_Frequency\human_intron_12_result.12.stranded.tsv'
    fn2='YOUR_PATH\SNV_Base_Frequency\human_ensembl90.intron.singleFreq.tsv'
     human_intron_explevel = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_explevel_cluster,...
                                     human_geneID,...
                                     'F_compare_variants','yes',...
                                     'Filter_zero_mut_gene','yes',...
                                     'FaceColor',explevel_facecm,...
                                     'Title_prefix','Human-intron_explevel'); 
%Correlation between SNP coding/template strand mutation rates with expression level
    fn_Liu='\Liu_and_Zhang_data\mutation_expressoin_data.xlsx'
    a5_Liu = importdata(fn_Liu);
[~,i,j] = intersect(human_geneID,human_intron_explevel.variants_level_structure.gname,'stable');
mut_intron = human_intron_explevel.vnum_pergene(j,:);
locus_length = human_intron_explevel.vglength(j,:);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_intron,'stable');
gexp  = gexp(i);
mut_intron = mut_intron(i,:);
locus_length = locus_length(i,:);
Ctl_profile = a5_Liu.data.intron_SNP_density(j,:);
CS_mutrate = sum(mut_intron(:,[1,3,5,7,9,11])')' ./ locus_length(:,1);
TS_mutrate = sum(mut_intron(:,[2,4,6,8,10,12])')' ./ locus_length(:,1);
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(CS_mutrate,gexp,'Type','Spearman')
[rho(2,1),pval(2,1)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,2),'Type','Spearman')
[rho(3,1),pval(3,1)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,3),'Type','Spearman')
[rho(4,1),pval(4,1)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,2:3),'Type','Spearman')
[rho(1,2),pval(1,2)] = corr(TS_mutrate,gexp,'Type','Spearman')
[rho(2,2),pval(2,2)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,2),'Type','Spearman')
[rho(3,2),pval(3,2)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,3),'Type','Spearman')
[rho(4,2),pval(4,2)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,2:3),'Type','Spearman')
figure; bar(rho');set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'CS','TS'});ylabel('Spearman rho, using intronic SNP');ylim([-0.3 0.25]);
title('Spearman rho: CS-TS with exp-level')

%SNP mutation rate on CS, and split genes by expression level
percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(CS_mutrate(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,2),'Type','Spearman')
[rho(1,3),pval(1,3)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,3),'Type','Spearman')
[rho(1,4),pval(1,4)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,2:3),'Type','Spearman')
[rho(2,1),pval(2,1)] = corr(CS_mutrate(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,2),'Type','Spearman')
[rho(2,3),pval(2,3)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,3),'Type','Spearman')
[rho(2,4),pval(2,4)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,2:3),'Type','Spearman')
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using intronic SNP')
title('Spearman rho: CS & Split genes');ylim([-0.3 0.22]);

%SNP mutation rate on TS, and split genes by expression level
percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(TS_mutrate(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,2),'Type','Spearman')
[rho(1,3),pval(1,3)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,3),'Type','Spearman')
[rho(1,4),pval(1,4)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,2:3),'Type','Spearman')
[rho(2,1),pval(2,1)] = corr(TS_mutrate(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,2),'Type','Spearman')
[rho(2,3),pval(2,3)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,3),'Type','Spearman')
[rho(2,4),pval(2,4)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,2:3),'Type','Spearman')
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using intronic SNP')
title('Spearman rho: TS & Split genes');	
                    
%Correlation between DNM coding/template strand mutation rates with expression level
% An_ctl_&_Jonsson_all DNM rates across gene clusters
prot_idx = strmatch('protein_coding',cellstr(human_gene_info.type));
autosome={'1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22'};
auto = ismember(cellstr(human_gene_info.chromosome(prot_idx,:)),autosome);
auto_prot_genes = cellstr(human_gene_info.ensembl_id(prot_idx(auto),:));
[~,~,j]=intersect(auto_prot_genes,human_geneID,'stable');
fn1="YOUR_PATH\SNV_Base_Frequency\Jonsson_2017_Nature_SNV_all_result.12.stranded.tsv"
fn2="YOUR_PATH\SNV_Base_Frequency\human_ensembl90.genebody.singleFreq.tsv"
fn3="YOUR_PATH\SNV_Base_Frequency\An_2018_SNV_control_result.12.stranded.tsv"
    %group by gene clusters
    An_Jonsson_all_filter_zero = F_compare_snv_asymmetry_v2_20200724(  fn1,   fn2,...
                                     h_gene_cluster(j),...
                                     human_geneID(j),...
                                     'Addup_SNV_file',fn3,...
                                     'Mode','single_gene',...
                                     'F_compare_variants','yes',...
                                     'Include_meta_all','no',...
                                     'Filter_zero_mut_gene','yes',...
                                     'Ylim_fold_factor',1.0,...
                                     'FaceColor',explevel_edgecm,... 
                                     'Title_prefix','An-Jonsson-all'); 

[~,i,j] = intersect(human_geneID,An_Jonsson_all_filter_zero.variants_level_structure.gname,'stable');
mut_intron = An_Jonsson_all_filter_zero.vnum_pergene(j,:);
locus_length = An_Jonsson_all_filter_zero.vglength(j,:);
geneID = human_geneID(i);
gexp  = h_mean_exp(i);
[~,i,j] = intersect(geneID,Liu_genenames_denovo,'stable');
gexp  = gexp(i);
mut_intron = mut_intron(i,:);
locus_length = locus_length(i,:);
Ctl_profile = a5_Liu.data.denovo_mutation(j,:);
CS_mutrate = sum(mut_intron(:,[1,3,5,7,9,11])')' ./ locus_length(:,1);
TS_mutrate = sum(mut_intron(:,[2,4,6,8,10,12])')' ./ locus_length(:,1);
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(CS_mutrate,gexp,'Type','Spearman')
[rho(1,2),pval(1,2)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,2),'Type','Spearman')
[rho(1,3),pval(1,3)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,3),'Type','Spearman')
[rho(1,4),pval(1,4)] = partialcorr(CS_mutrate,gexp,Ctl_profile(:,2:3),'Type','Spearman')
[rho(2,1),pval(2,1)] = corr(TS_mutrate,gexp,'Type','Spearman')
[rho(2,2),pval(2,2)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,2),'Type','Spearman')
[rho(2,3),pval(2,3)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,3),'Type','Spearman')
[rho(2,4),pval(2,4)] = partialcorr(TS_mutrate,gexp,Ctl_profile(:,2:3),'Type','Spearman')
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'CS','TS'});ylabel('Spearman rho, using DNM')
title('Spearman rho: DNM CS-TS with exp-level');ylim([-0.07 0.03]);                    
%SNP mutation rate on CS, and split genes by expression level
percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(CS_mutrate(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,2),'Type','Spearman')
[rho(1,3),pval(1,3)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,3),'Type','Spearman')
[rho(1,4),pval(1,4)] = partialcorr(CS_mutrate(i),gexp(i),Ctl_profile(i,2:3),'Type','Spearman')
[rho(2,1),pval(2,1)] = corr(CS_mutrate(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,2),'Type','Spearman')
[rho(2,3),pval(2,3)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,3),'Type','Spearman')
[rho(2,4),pval(2,4)] = partialcorr(CS_mutrate(j),gexp(j),Ctl_profile(j,2:3),'Type','Spearman')
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using DNM')
title('Spearman rho: CS & Split genes');ylim([-0.07 0.03]);

%SNP mutation rate on TS, and split genes by expression level
percent_cutoff = 70;
i = find(gexp<=prctile(gexp,percent_cutoff));
j = find(gexp>prctile(gexp,percent_cutoff));
rho = []; pval = [];
[rho(1,1),pval(1,1)] = corr(TS_mutrate(i),gexp(i),'Type','Spearman');
[rho(1,2),pval(1,2)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,2),'Type','Spearman')
[rho(1,3),pval(1,3)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,3),'Type','Spearman')
[rho(1,4),pval(1,4)] = partialcorr(TS_mutrate(i),gexp(i),Ctl_profile(i,2:3),'Type','Spearman')
[rho(2,1),pval(2,1)] = corr(TS_mutrate(j),gexp(j),'Type','Spearman');
[rho(2,2),pval(2,2)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,2),'Type','Spearman')
[rho(2,3),pval(2,3)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,3),'Type','Spearman')
[rho(2,4),pval(2,4)] = partialcorr(TS_mutrate(j),gexp(j),Ctl_profile(j,2:3),'Type','Spearman')
%Plot the corr coef by spliting the gene groups into low-mid & highs.
figure; bar(rho);set(gca,'color','none');box off;legend({'W/t Ctl','Ctl RT','Ctl GC','Ctl RT&GC'});
xticklabels({'low-mid exp','High exp'});ylabel('Spearman rho, using DNM')
title('Spearman rho: TS & Split genes'); ylim([-0.07 0.03]);





                   


%% Germ-soma expression 
    load('YOUR_PATH\Data\GTEx_Analysis_human_RNA-seq\GTEx_Analysis_2016-01-15_v7_RNASeQCv1.1.8_gene_median_tpm.mat')
    GTExAnalysis20160115v7RNASeQCv1_geneID ={}
    for i=1:length(GTExAnalysis20160115v7RNASeQCv1_gene(:,1))
        temp = strsplit(GTExAnalysis20160115v7RNASeQCv1_gene(i,1),'.');
        GTExAnalysis20160115v7RNASeQCv1_geneID(i,1) = cellstr(temp{1});
    end
    [input_geneID,xi,xj] = intersect(human_geneID,GTExAnalysis20160115v7RNASeQCv1_geneID,'stable');
    input_genename = human_genename(xi);
    sc_testis_exp = mean(h_cluster_exp_meanUMI(xi,:)')';
    sc_genecluster = h_gene_cluster(xi);
    %Set up GTEx non-testis protein-coding gene expression profiles
    [~,j]=setdiff(GTExAnalysis20160115v7RNASeQCv1_Sample,'Testis');
    GTEx_non_testis_exp = GTExAnalysis20160115v7RNASeQCv1(xj,j);
    k=strmatch('Testis',cellstr(GTExAnalysis20160115v7RNASeQCv1_Sample));
    GTEx_testis_exp = GTExAnalysis20160115v7RNASeQCv1(xj,k);
    
    % Generate gene index
    length(find(max(GTEx_non_testis_exp')<0.3))
    length(find(GTEx_testis_exp<0.3))
    max_gtex = max(GTEx_non_testis_exp')';
    cutoff = 0.3;
    clear GTEx_gene_cluster
    GTEx_gene_cluster(1:length(input_geneID)) = zeros;
    GTEx_gene_cluster(find(sc_genecluster==0 & max_gtex<cutoff)) = 1; %Germ-unexp; soma-unexp
    GTEx_gene_cluster(find(sc_genecluster==0 & max_gtex>=cutoff)) = 2; %Germ-unexp; soma-exp
    GTEx_gene_cluster(find(sc_genecluster>0 & max_gtex<cutoff)) = 3; %Germ-exp; soma-unexp
    GTEx_gene_cluster(find(sc_genecluster>0 & max_gtex>=cutoff)) = 4; %Germ-exp; soma-exp
    histcounts(GTEx_gene_cluster)
    input_genename(find(GTEx_gene_cluster == 3))
    figure;histogram(log2(sc_testis_exp(find(GTEx_gene_cluster == 3))))

    
    %For GO term analaysis:
    %Target gene set: uniquely expressed in the germ cells  (Germ-exp; soma-unexp)
        input_genename(find(GTEx_gene_cluster == 3))
    %Background gene set: all genes in the analysis
        input_genename
        
        