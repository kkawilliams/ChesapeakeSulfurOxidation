import os
import pandas as pd
import numpy as np
from sys import exit
os.chdir('/Volumes/KeithSSD/SulFox/scripts')
import seaborn as sns
import matplotlib as mpl
from matplotlib import pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from pathwayFunctions import *

base_path = "/Volumes/KeithSSD/SulFox/data"

# Load abundances:
bin_abundances = '/Volumes/KeithSSD/CBFunctions/data/outputs/bin_abund_n_metadata.p'
cbTr_abnd2, cb33_abnd2, metadata = metadata_and_abundances(bin_abundances)

# load lineage info & quality 
qualtax_file = '/Volumes/KeithSSD/CBFunctions/data/outputs/qual_tax_full.p'
qualtax_df = load_qualtax(qualtax_file)

# Locate prokka annotations 
annotation_folder = os.path.join(base_path, 'annotation')
bins_preselected = locate_annotated_files(annotation_folder)
bins_desired = [i for i in bins_preselected if i in set(qualtax_df.index[qualtax_df['Completeness'] >= 80])]

# 
report_files_path = os.path.join(base_path, "minpath_out")
target_pathways = os.path.join(base_path, "tab_seperated_files/PathwaysConsidered.csv")
suffix = '.report.txt'
mc_pathway_df = aggregate_and_filter_minpath_data(report_files_path, suffix, annotation_folder, target_pathways, bins_desired)

qt_idx = set(qualtax_df.index).intersection(set(mc_pathway_df.genome))
complete_list = qt_idx - set(qualtax_df.index[qualtax_df['Completeness'] < 80])
mc_pathway_df = mc_pathway_df[mc_pathway_df.locus_tag.notnull() & mc_pathway_df.genome.isin(complete_list)]

cog_presabs = categorical_pivot(mc_pathway_df, 'db_xref', 'genome')

kegg_fegenie_annots = '/Volumes/KeithSSD/CBFunctions/data/outputs/gene_annotations2.p'
prot_dir = '/Volumes/KeithSSD/CBFunctions/data/ORF_calls'
kf_pathway_df = load_kegg_database(kegg_fegenie_annots, qt_idx, prot_dir)

target_genomes = set(mc_pathway_df.genome).intersection(set(kf_pathway_df.genome))
pathway_df = merge_databases(mc_pathway_df, kf_pathway_df, target_genomes).drop_duplicates()
pathway_df['merger_col'] = pathway_df[['genome', 'contig', 'start', 'end']].apply(tuple, axis=1)
pathway_df['identifier'] = pathway_df['genome'] + "_" + pathway_df['contig'] + "_" + pathway_df['start'] + "_" + pathway_df['end']
annot_mat, gene_sets = parse_annotations_into_matrix(pathway_df)

seq_map = pull_genes(pathway_df, '/Volumes/KeithSSD/CBFunctions/data/all_bins','/Volumes/KeithSSD/SulFox/data/pulled_genes', '.fa')
seq_map = {"_".join(k):v for k, v in seq_map.items()}

salmon_df, salmon_dfs, normed_bins = load_salmon_data('/Volumes/KeithSSD/SulFox/data/salmon_out3', pathway_df)
salmon_df.drop(['SRR988008', 'SRR988009'], axis=1, inplace=True)
salmon_df = salmon_df.rename(columns=hewson_translator)
for d in [salmon_dfs, normed_bins]:
    for k in d.keys():
        subdf = d[k]
        d[k] = subdf.drop(['SRR988008', 'SRR988009'], axis=1).rename(columns=hewson_translator)        

copy_data_file = "/Volumes/KeithSSD/SulFox/data/gene_expression/gene_copy_number_data.txt"
copy_duplicates_file = "/Volumes/KeithSSD/SulFox/data/gene_expression/copy_number_duplicate_clusters.txt"
copy_num_df = load_cn_data(copy_data_file, copy_duplicates_file)

gene_id_lookup1 = pathway_df[['genome', 'GeneID']].dropna().groupby('genome').agg(set)['GeneID'].to_dict()
gene_id_lookup2 = pathway_df[['genome', 'merger_col']].dropna().groupby('genome').agg(set)['merger_col'].to_dict()
gene_id_lookup3 = {i:j for i, j in pathway_df[['identifier', 'GeneID']].dropna().values}

idtyp1_idtyp2 = {i:"_".join(k) for i, k in list(pathway_df[['GeneID', 'merger_col']].values)}
idtyp2_idtyp1 = {k:i for i, k in idtyp1_idtyp2.items()}
gene_id_lookup2 = {i:set(["_".join(x) for x in j]) for i, j in gene_id_lookup2.items()}

report_df = make_report(qualtax_df, complete_list, cbTr_abnd2, cb33_abnd2, salmon_df, copy_num_df, gene_id_lookup3, gene_sets)

# Convert the dataframe to an XlsxWriter Excel object.

add_genes_N = lambda x: " ".join(list(np.array(list(gene_cats.iloc[14:21, 1]))[annot_mat.loc[x, list(gene_cats.iloc[14:21, 1])]]))
add_genes_S = lambda x: len(list(np.array(list(gene_cats.iloc[22:-7, 1]))[annot_mat.loc[x, list(gene_cats.iloc[22:-7, 1])]]))
qualtax_df.loc[report_df.sort_values(by='mean rank', ascending=False).index[:40], ['Order', 'Family', 'Genus', 'Contamination']].join(pd.Series({i:add_genes_N(i) for i in report_df.sort_values(by='mean rank', ascending=False).index[:40]}, name='N')).join(pd.Series({i:add_genes_S(i) for i in report_df.sort_values(by='mean rank', ascending=False).index[:40]}, name='S'))

bins_of_focus = ['CBrest.bin.246', 'CBrest.bin.40', 'CB33.bin.111', 'CBrest.bin.375', 'CB33.bin.191',
                 'CBrest.bin.358', 'CB33.bin.222', 'CB33.bin.236', 'CB33.bin.196', 'CB33.bin.2']
qualtax_df.loc[bins_of_focus, ['Order', 'Family', 'Genus', 'Contamination']].join(pd.Series({i:add_genes_N(i) for i in bins_of_focus}, name='N')).join(pd.Series({i:add_genes_S(i) for i in bins_of_focus}, name='S'))
other_bins = [i for i in report_df.index if not i in bins_of_focus]

writer = pd.ExcelWriter("/Volumes/KeithSSD/SulFox/data/bin_prioritization_report.xlsx", engine='xlsxwriter')
report_df.loc[bins_of_focus, :].sort_values(by='mean rank', ascending=False).to_excel(writer, sheet_name='selected_bins')
report_df.loc[other_bins, :].sort_values(by='mean rank', ascending=False).to_excel(writer, sheet_name='other_bins')
writer.save()

special_gene_sets = set([annot for bof in bins_of_focus for annot in annot_mat.columns[annot_mat.loc[bof, :]]])
special_gene_sets = special_gene_sets - set(['', 'bcl_synth', 'taurine_degradation'])

housekeepers = ["proC", "recA", "rpoD", "rho", "glyA", "tpiA", "recF"]
writer = pd.ExcelWriter("/Volumes/KeithSSD/SulFox/data/priority_expression_report.xlsx", engine='xlsxwriter')

expression_summaries = {}
housekeepers_to_annotate = {}
for bof in bins_of_focus:
    special_genes = {}
    for sgs in special_gene_sets:
        for eid_ in gene_sets[sgs]:
            if eid_.startswith(bof+"_"):
                possible_name = "--".join(set(pathway_df.loc[pathway_df['identifier'] == eid_, 'gene'].dropna().values))
                possible_qid = "--".join(set(pathway_df.loc[pathway_df['identifier'] == eid_, 'query_id'].dropna().values))
                possible_prod = "--".join(set(pathway_df.loc[pathway_df['identifier'] == eid_, 'product'].dropna().values))
                special_genes[eid_] = {'bin':bof, 'enzyme':sgs, 'gene_type':'functional', 
                                       'gene': None, 'title': None, 'hmm':None}
                special_genes[eid_].update({'sequence':seq_map[eid_]})
                
                for field, possib in zip([possible_name, possible_qid, possible_prod], ['gene', 'hmm', 'title']):
                    special_genes[eid_].update({possib:field})
                if eid_ in set(salmon_dfs[bof].index):
                    special_genes[eid_].update(salmon_dfs[bof].loc[eid_, numeric_exp].to_dict())
    
    for h in housekeepers:
        hbool = (pathway_df.gene == h) & (pathway_df.genome == bof)
        if hbool.sum() > 0:
            hid_s = set(pathway_df.loc[hbool, 'identifier'])
            for eid_ in hid_s:
                special_genes[eid_] = {'bin':bof, 'enzyme':h, 'gene_type':'housekeeper', 
                                       'gene': None, 'title': None, 'hmm':None}
                special_genes[eid_].update(salmon_dfs[bof].loc[eid_, numeric_exp].to_dict())
                special_genes[eid_].update({'sequence':seq_map[eid_]})
                housekeepers_to_annotate[eid_] = seq_map[eid_]
    
    topranked = list(normed_bins[bof].sum(1).sort_values(ascending=False).index[:50])
    needed_cols = ['gene', 'product', 'query_id', 'identifier']
    catagg = lambda x: "--".join(set(x.dropna()))
    highlights_top = pathway_df.loc[pathway_df.identifier.isin(topranked), needed_cols].groupby('identifier').agg(catagg)
    highlights_top = highlights_top.drop([i for i in highlights_top.index if i in special_genes])
    for eid_ in highlights_top.index:
        if not eid_ in special_genes:
            special_genes[eid_] = {'bin':bof, 'gene_type':'ranked', 'hmm':None, 'enzyme': None,
                                   'title':highlights_top.loc[eid_, 'product'], 
                                   'gene':highlights_top.loc[eid_, 'gene']}
            special_genes[eid_].update(salmon_dfs[bof].loc[eid_, numeric_exp].to_dict())
            special_genes[eid_].update({'sequence':seq_map[eid_]})
        else:
            special_genes[eid_]['gene_type'] = special_genes[eid_]['gene_type'] + "+ranked"
    
    expression_summary = pd.DataFrame(special_genes).T
    bad_cols = np.array(numeric_exp)[expression_summary.loc[expression_summary.enzyme.isin(housekeepers), numeric_exp].sum() == 0]
    to_drop = highlights_top.index[(expression_summary.loc[list(highlights_top.index), bad_cols] > 0).sum(1) != 0]
    expression_summary.drop(to_drop, axis=0, inplace=True)
    expression_summary = expression_summary.loc[:, [i for i in expression_summary.columns if i != 'sequence']+['sequence']]
    expression_summary.to_excel(writer, sheet_name=f'{bof}')
    expression_summaries[bof] = expression_summary.copy()

writer.save()

gene_cats = pd.read_csv('/Volumes/KeithSSD/SulFox/data/genes_by_pathway.txt', header=None, sep="\t")
df = pd.DataFrame(index=bins_of_focus, columns=list(special_gene_sets | set(gene_cats[1]))).T.fillna(0)
set(df.index) - set(gene_cats[1])
set(gene_cats[1]) - set(df.index)

for b in df.columns:
    print(b, report_df.loc[b, 'taxonomic leaf name'])
    exp_sum_i = expression_summaries[b]
    valid_samples = np.array(numeric_exp)[(exp_sum_i.loc[exp_sum_i.enzyme.isin(housekeepers), numeric_exp] > 5).sum() > 0]
    
    for gene in df.index:
        if gene in housekeepers:
            hbool = (pathway_df.gene == gene) & (pathway_df.genome == b)
            if hbool.sum() > 0:
                df.loc[gene, b] += 1
        elif gene in set(annot_mat.columns):
            if gene in set(annot_mat.columns[annot_mat.loc[b, :]]):
                df.loc[gene, b] += 1
        else:
            raise ValueError("Missing gene from annot mat / housekeepers")
    
        if gene in set(exp_sum_i.enzyme):
            gslice = exp_sum_i.loc[exp_sum_i.enzyme == gene, valid_samples]
            gslice[gslice<5] = 0
            total_ge = gslice.sum().sum()
            if total_ge > (len(valid_samples)*5):
                df.loc[gene, b] += 1
 
df = df.loc[[i for i in gene_cats[1] if i in set(df.index)], :]
gene_names_and_cats = {i:"{} ({})".format(i, j.lower()) for i, j in gene_cats.set_index(1)[0].to_dict().items()}
bin_names_and_taxa = {i:"{} ({})".format(i, report_df.loc[i, 'taxonomic leaf name']) for i in df.columns}
df = df.rename(index=gene_names_and_cats, columns=bin_names_and_taxa)

# plot

myColors = ((1.0, 1.0, 1.0, 1.0), (0.0, 0.4, 0.0, 0.5), (0.0, 0.8, 0.0, 1.0))
cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))

sns.set(font="Calibri", font_scale=0.74)
plt.close('all')
fig = plt.figure(constrained_layout=True, figsize=(6*1.2, 8*1.2), dpi=300, num=1)
gs = fig.add_gridspec(8, 6)
ax = fig.add_subplot(gs[:-1, :])
ax2 = fig.add_subplot(gs[-1, 0])
#ax3 = fig.add_subplot(gs[-1, 0])
cm_fig = sns.heatmap(df, cbar_ax=ax2, ax=ax, cmap=cmap, linewidths=.2, cbar_kws={'drawedges':True})
ax.set_xticklabels(ax.get_xticklabels(),rotation=30, ha="right")
ax2.set_yticklabels(['Not Observed', 'Annotation Observed', '>1x expression observed'])
ax2.yaxis.set_ticks([ 0.3, 0.9, 1.6])
fig.savefig('/Volumes/KeithSSD/SulFox/figures/gene_repertoire.png', dpi=300)

cb33_cols = [i for i in set(copy_num_df.columns) & set(cb33_abnd2.columns) if not i.startswith('Zymo') and not i.startswith('SERC')]
cb33_tuples = [i.replace("CB33_", "").replace("June", "6").replace("Aug", "8").replace("July", "7").replace("_2017", "") for i in cb33_cols]
cb33_tuples = [i.replace("_R2", "").replace("R2", "").replace("_R7", "").replace("_R1", "").split("_") for i in cb33_tuples]
cb33_tuples = [[i[0], i[-1]] for i in cb33_tuples]
cb33_df = pd.DataFrame(index=cb33_cols, data=cb33_tuples)
cb33_df[2] = cb33_df[1].apply(lambda x: int(x.replace("m", "")))
cb33_df = cb33_df.sort_values(by=[0, 2])
cb33_df['label'] = cb33_df[0] + "-33-" + cb33_df[1]

cbRest_df = pd.DataFrame(index=cbRest_cols, data=cbRest_tuples)

cbRest_cols = [i for i in set(copy_num_df.columns) - set(cb33_abnd2.columns) if not 'Bay_Bridge' in i and not 'Negative' in i]
cbRest_tuples = [i.split("_") for i in cbRest_cols]
cbRest_df = pd.DataFrame(index=cbRest_cols, data=cbRest_tuples)
cbRest_df.loc[cbRest_df[0] == '4102017', 0] = '41017'
cbRest_df.loc[cbRest_df[0] == '6617', 0] = '60617'
cbRest_df.loc[cbRest_df[0] == '6517', 0] = '60517'
cbRest_df.loc[cbRest_df[1] == 'CB53-1', 1] = 'CB53'
cbRest_df[3] = cbRest_df[0].apply(lambda x: x[0])
cbRest_df[4] = cbRest_df[1].apply(lambda x: x.replace("CB", "").replace('C', ""))
cbRest_df = cbRest_df.sort_values(by=[2,3,4])
cbRest_df['label'] = cbRest_df[3] + "-" + cbRest_df[4] + "-" + cbRest_df[2]

cb33copies = copy_num_df[list(cb33_df.index)]
cbRestcopies = copy_num_df[list(cbRest_df.index)]

exp_cols = [i for i in salmon_df.columns if not i == 'merger_col']
exp_tuples = [i.split("_") for i in exp_cols]
exp_df = pd.DataFrame(index=exp_cols, data=exp_tuples)
exp_df[4] = exp_df[0].apply(lambda x: x.split("/")[0])
exp_df['label'] = exp_df[4] + "-43-" + exp_df[1]
exp_df = exp_df.sort_values(by=[1, 4])

import matplotlib.gridspec as gridspec

plt.close('all'); plt.clf();
fig2 = plt.figure(constrained_layout=True, figsize=(11, 8), dpi=300, num=2)
gs2 = gridspec.GridSpec(ncols=10, nrows=11, figure=fig2)
all_axes = {}
for b_i, bof in enumerate(bins_of_focus):
    spec_ax = {}
    spec_ax['name'] = fig2.add_subplot(gs2[b_i, 0], facecolor="white", frameon=True)
    _ = spec_ax['name'].set_xlim(-1, 1.)
    _ = spec_ax['name'].set_ylim(-1, 1.)
    _ = spec_ax['name'].axes.get_xaxis().set_ticks([])
    _ = spec_ax['name'].axes.get_yaxis().set_ticks([])
    _ = spec_ax['name'].text(0, .5, bof, ha="center", va="center")
    _ = spec_ax['name'].text(0, -.5, report_df.loc[bof, 'taxonomic leaf name'], ha="center", va="center")
    ##################################################
    spec_ax['cb33'] = fig2.add_subplot(gs2[b_i, 1:4], facecolor="white", frameon=True)
    cleaned_2 = set(pathway_df.loc[pathway_df['genome'] == bof, 'GeneID'].dropna())
    cleaned_4 = cleaned_2 & set(copy_num_df.index)
    cb33_slice = cb33copies.loc[cleaned_4, :].apply(np.median).values
    cb33_upper = cb33copies.loc[cleaned_4, :].apply(lambda x: np.percentile(x, 75)).values - cb33_slice
    cb33_lower = cb33_slice - cb33copies.loc[cleaned_4, :].apply(lambda x: np.percentile(x, 25)).values
    cb33_yerr  = np.vstack((cb33_upper, cb33_lower))
    _ = spec_ax['cb33'].errorbar(np.arange(len(cb33_df)), cb33_slice, yerr=cb33_yerr, 
            color='k', marker='.', capsize=1, capthick=1, ecolor='black')
    # _ = spec_ax['cb33'].plot(np.arange(len(cb33_df)), cb33_slice, color='k')
    if b_i == (len(bins_of_focus) - 1):
        _ = spec_ax['cb33'].set_xticks(np.arange(len(cb33_df)))
        _ = spec_ax['cb33'].set_xticklabels(list(cb33_df['label']), rotation=90, ha="center")
    else:
        _ = spec_ax['cb33'].axes.get_xaxis().set_ticks([])
    ##################################################
    spec_ax['cbrest'] = fig2.add_subplot(gs2[b_i, 4:8], facecolor="white", frameon=True)
    cbRest_slice = cbRestcopies.loc[cleaned_4, :].apply(np.median).values
    _ = spec_ax['cbrest'].plot(np.arange(len(cbRest_slice)), cbRest_slice, color='k')
    if b_i == (len(bins_of_focus) - 1):
        _ = spec_ax['cbrest'].set_xticks(np.arange(len(cbRest_slice)))
        _ = spec_ax['cbrest'].set_xticklabels(list(cbRest_df['label']), rotation=90, ha="center")
    else:
        _ = spec_ax['cbrest'].axes.get_xaxis().set_ticks([])
    ##################################################
    spec_ax['mrna'] = fig2.add_subplot(gs2[b_i, 8:10], facecolor="white", frameon=True)
    exp_slice = salmon_dfs[bof].loc[:, exp_df.index].apply(np.median).values
    _ = spec_ax['mrna'].plot(np.arange(len(exp_df)), exp_slice, color='k')
    if b_i == (len(bins_of_focus) - 1):
        _ = spec_ax['mrna'].set_xticks(np.arange(len(exp_slice)))
        _ = spec_ax['mrna'].set_xticklabels(list(exp_df['label']), rotation=90, ha="center")
    else:
        _ = spec_ax['mrna'].axes.get_xaxis().set_ticks([])
    ##################################################
    all_axes[bof] = spec_ax

all_axes['annotations'] = fig2.add_subplot(gs2[-1, :], facecolor="white", frameon=True)
_ = all_axes['annotations'].axes.get_xaxis().set_ticks([])
_ = all_axes['annotations'].axes.get_yaxis().set_ticks([])
fig2.savefig('/Volumes/KeithSSD/SulFox/figures/bin_location.png', dpi=300)

# Manually specify colorbar labelling after it's been generated



help(spec_ax['mrna'].errorbar)
#ax.figure.tight_layout()




# X - Y axis labels
ax.set_ylabel('FROM')
ax.set_xlabel('TO')

# Only y-axis labels need their rotation set, x-axis labels already have a rotation of 0
_, labels = plt.yticks()
plt.setp(labels, rotation=0)

add_legend
tight_layout
, , 

with open('/Volumes/KeithSSD/SulFox/data/housekeeping_genes.fa', 'w') as fh:
    for k, v in housekeepers_to_annotate.items():
        fh.write(">{}\n{}\n".format(k, v))

# read in kofam data for housekeepers
'/Volumes/KeithSSD/SulFox/data/housekeeping_genes.kofam.txt'


style_df =pd.DataFrame(index=cog_presabs.index, columns=['colors', 'markers', 'size'])
style_df['markers'] = 'o'
style_df['colors'] = '#000000'
style_df['size'] = 20
style_df['class'] = 'None'
packed_style = style_df[['markers', 'colors', 'size', 'class']].apply(tuple, axis=1).unique()
cog_fig = plot_presabs(cog_presabs, style_df, packed_style)
cog_fig.savefig('/Volumes/KeithSSD/SulFox/figures/TSNE_Fig.png')


    


#  

len(annotations['aerobic_resp'].intersection(annotations['rbcSL'])) / len(annotations['rbcSL'])
len(annotations['aerobic_resp'].intersection(annotations['rTCA'])) / len(annotations['rTCA'])
len(annotations['aerobic_resp'].intersection(annotations['3hp cycle'])) / len(annotations['3hp cycle'])
len(annotations['aerobic_resp'].intersection(annotations['3hp4hb cycle'])) / len(annotations['3hp4hb cycle'])





# new figure 

# replot figure 1

# 







self.functional_proteins = mergecomplex(self.functional_proteins, merge_cols)
def mergecomplex(x, cols):
    return x.groupby(cols).agg(lambda x: "/".join([i for i in set(x) if str(i) != 'nan']))
