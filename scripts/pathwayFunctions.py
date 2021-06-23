import os
import pickle
import time
import gc

import pandas as pd
import numpy as np
#import seaborn as sns
#import matplotlib.pyplot as plt
import scipy.stats as sst
import scipy.spatial.distance as ssd

#from matplotlib.lines import Line2D
#from matplotlib.patches import Patch
from collections import Counter, defaultdict
from sys import argv, exit
from Bio import SeqIO

NITROGEN_CYCLE = {'NO3 oxioreductase (nxR/narGH)': set(['K00370', 'K00371']),
                  'hydroxylamine dehydrogenase (hao)': set(['K10535']),
                  'ammonia oxidation (amoABC)': set(['K10944', 'K10945', 'K10946']),
                  'nitrate reductase (napAB/narIJ)': set(['K02567', 'K02568', 'K00374', 'K00373']),
                  'denitrification': set(['K00368', 'K15864', 'K00376', 'K02305', 'K04561']),
                  'dnra': set(['K00362', 'K00363', 'K03385', 'K15876'])}


SPECIFIC_KEGG_CODES = {'photosynthetic reaction center accesories (puhA/pufXC)': ['K13991', 'K13994', 'K13992'],
                       'light-harvesting complex 1 (pufAB)' : ['K08926' , 'K08927'],
                       'light-harvesting protein B-800-850 (pucAB)': ['K08930', 'K08939'],
                       'photosystem P840 reaction center (pscABCD)': ['K08940', 'K08941', 'K08942', 'K08943'],
                       'bacteriochlorophyll A protein (fmoA)': ['K08944'],
                       'chlorosome envelope proteins (csm)': ['K08945', 'K08946', 'K08947', 'K08948', 'K08949', 'K08950', 'K08951', 'K08952', 'K08953','K08954'],
                       'sulfur carrier complex (dsrEFH)': ['PF02635', 'PF04077', 'PF13686'],
                       'ABC transporters': ['K12536', 'K05648'], 
                       'succinate dehydrogenase (SDH)': ['K00234', 'K00235', 'K00236', 'K00237', 'K00239', 'K00240', 'K00241', 'K00242', 'K18859', 'K18860'],
                       'ATP citrate (pro-S)-lyase (ACLY)': ['K01648', 'K15230', 'K15231'],
                       '2-oxoglutarate dehydrogenase (sucB)': ['K00658'],
                       'malate dehydrogenase (mqo)': ['K00116'],
                       'methane/ammonia monooxygenase': ['K10944', 'K10945', 'K10946'],
                       'hydroxylamine reductase (hcp)': ['K05601'],
                       'nitrous-oxide reductase (nosZ)':['K00376'],
                       'assimilatory nitrate reductase (nasAB)': ['K00360', 'K00372'],
                       'nitrate reductase (NAD(P)H), (NR)': ['K10534'],
                       'ferredoxin-nitrate reductase (narB)': ['K00367'],
                       'NO3 reductase (narIJ)': ['K00374', 'K00373'],
                       'NO3 oxioreductase (nxR/narGH)': ['K00370', 'K00371'],
                       'periplasmic NO3 reductase (napAB)': ['K02567', 'K02568'],
                       'hydrazine oxidoreductase (hzo)': ['K20935'],
                       'nitrite reductase (nirKS)': ['K00368', 'K15864'],
                       'cytochrome c oxidase (COX1)': ['K00404', 'K02256', 'K02274', 'K02275', 'K02276', 'K02277', 'K15408', 'K15862'],
                       'cytochrome c nitrite reductase (nrfAH)': ['K03385', 'K15876'],
                       'nitrite reductase (NADH) (nirDB)': ['K00362', 'K00363'],
                       'nitric oxide reductase (norBC)': ['K04561', 'K02305', 'K15877'],
                       'hydroxylamine dehydrogenase (hao)': ['K10535'],
                       'thiosulfate/polysulfide reductase chain A (psrA)': ['K08352'],
                       'thiosulfate dehydrogenase (doxA)': ['K16936'],
                       'sulfur oxygenase/reductase (sor)': ['K16952'],
                       'sulfide dehydrogenase (fccB)': ['K17229'],
                       'sulfhydrogenase subunit gamma (sulfur reductase) (hydG)': ['K17995'],
                       'adenylyl-sulfate reductase (APR)': ['K05907'],
                       'anaerobic sulfite reductase (asrA)': ['K16950', 'K16951'],
                       'sulfur reductase molybdopterin subunit (sreA)': ['K17219'],
                       'sulfur dioxygenase (ETHE1)': ['K17725'],
                       'eukaryotic sulfide quinone oxidoreductase (SQOR)' : ['K22470'],
                       'ribulose-bisphosphate carboxylase small chain (rbcSL)':['K01601', 'K01602'],
                       'dissimilatory sulfite reductase D (DsrD)': ['PF08679'],
                       'alt sulfur carrier complex (dsrEFH)': ['K07235', 'K07236', 'K07237'],
                       'sulfite oxidase (SUOX)': ['K00387'], 
                       'PAPSS': ['K13811'],
                       'photosynthetic reaction center (pufLM)': ['K08928', 'K08929'],
                       'thiosulfate dehydrogenase (tsdA)':['K19713'],
                       'sulfate adenylyltransferase (sat/met3)': ['K00958'],
                       'adenylylsulfate reductase (aprAB)': ['K00394', 'K00395'],
                       'phosphosulfate reductase (cysH)': ['K00390'],
                       'sulfite dehydrogenase (sorA)': ['K05301'],
                       'sulfate adenylyltransferase cysN/D': ['K00955', 'K00956', 'K00957'],
                       'adenylylsulfate kinase cycC': ['K00860'],
                       'sulfite dehydrogenase (soeABC)': ['K21307', 'K21308', 'K21309'],
                       'NIR_SIR alternative dsrAB model': ['PF01077'],
                       'sulfide:quinone oxidoreductase (sqr)': ['K17218'],
                       'sulfur-oxidizing protein SoxYZ': ['K17226', 'K17227'],
                       'dissimilatory sulfite reductase (dsrAB)': ['K11180', 'K11181'],
                       'L-cysteine S-thiosulfotransferase (soxAX)': ['K17222', 'K17223'], 
                       'S-sulfosulfanyl-L-cysteine sulfohydrolase (soxB)': ['K17224'],
                       'sulfane dehydrogenase (soxCD)': [ 'K22622', 'K17225'],
                       'sulfite reductase (cysJI)': ['K00381', 'K00380'],
                       "sulfite reductase ferredoxin (sir)": ['K00392'],
                       'ferredoxin hydrogenase': ['K00532', 'K00533', 'K00534', 'K06441', 'K18016', 'K18017', 'K18023'],
                       'methane monooxygenase (mmoXYZBCD)': ['K16157', 'K16158', 'K16159', 'K16160', 'K16161', 'K16162'], 
                       'methanol dehydrogenase (xoxF/mxaF)': ['K14028', 'K23995']}

KOS_OF_INTEREST = [i for l in SPECIFIC_KEGG_CODES.values() for i in l ]

def merge_databases(mc_pathway_df ,kf_pathway_df, target_genomes):
    merge_cols = ['contig', 'end', 'start', 'genome']
    desired_mc_cols = merge_cols + ['strand', 'eC_number', 'gene', 'locus_tag', 'product', 'path_accession', 'path_name', 'metacyc_target']
    desired_kf_cols = merge_cols + ['ClusterID', 'GeneID', 'query_id', 'full_evalue', 'HMM', 'bitscore', 'annot_code']
    pathway_df = pd.merge(left=mc_pathway_df.loc[mc_pathway_df.genome.isin(target_genomes), desired_mc_cols], 
                          right=kf_pathway_df.loc[kf_pathway_df.genome.isin(target_genomes), desired_kf_cols], 
                          how='outer', on=merge_cols, sort=True)
    pathway_df['metacyc_target'] = pathway_df['metacyc_target'].fillna(False)
    return pathway_df


def metadata_and_abundances(bin_abundances):
    with open(bin_abundances, "rb" ) as fh:
        outpkg = pickle.load(fh)
    name_correction = outpkg['renamer']
    metadata = outpkg['meta']
    cb33_abnd = outpkg['cb33_abnd']
    cbTr_abnd = outpkg['cbTr_abnd']
    cb33_abnd2 = pd.read_csv('/Volumes/KeithSSD/CBFunctions/data/bin_data_old/CB33.abundances.txt', sep="\t", index_col=0)
    cbTr_abnd2 = pd.read_csv('/Volumes/KeithSSD/CBFunctions/data/bin_data_old/CBrest.abundance_table.txt', sep="\t", index_col=0)
    cb33_abnd2.rename(index={i:'CB33.'+i for i in cb33_abnd2.index}, columns=name_correction, inplace=True)
    cbTr_abnd2.rename(index={i:'CBrest.'+i for i in cbTr_abnd2.index}, columns=name_correction, inplace=True)   
    return [cbTr_abnd2, cb33_abnd2, metadata]



def plot_presabs(cog_presabs, style_df, packed_style):
    fashion_tsne = TSNE(random_state=123, n_jobs=12).fit_transform(cog_presabs.values)
    fashion_tsne = pd.DataFrame(index=cog_presabs.index, data=fashion_tsne, columns=['x', 'y'])

    colors = ['#4477DD', '#99DDFF', '#44BB99', '#BBCC33', '#AAAA00', '#EEDD88', '#EE8866', '#FFAABB', 
              '#DDDDDD', '#CC3311', '#AA4499']
    plt.close('all')
    plt.style.use('seaborn-paper')
    f, ax = plt.subplots(nrows=1, ncols=2, figsize=(6,3), dpi=120)
    ax[0].set_facecolor('w')
    ax[0].grid(color='k', linestyle='-', linewidth=0.5)
    for idx in fashion_tsne.index:
        ax[0].scatter(x=fashion_tsne['x'][idx], y=fashion_tsne['y'][idx], 
                      s=style_df['size'][idx], marker=style_df['markers'][idx], 
                      c=style_df['colors'][idx])
    #plt.xlim(-25, 25)
    #plt.ylim(-25, 25)
    ax[0].set_yticklabels([])
    ax[0].set_ylabel('TSNE Axis 2')
    ax[0].set_xticklabels([])
    ax[0].set_xlabel('TSNE Axis 1')
    legend_elements = [Line2D([0], [0], marker=m, color='w', label=l, markerfacecolor=c, markersize=s/2) for m, c, s, l in packed_style]
    ax[1].legend(handles=legend_elements, loc='center', fontsize=10, title="Sulfur Metabolic Genotype")
    ax[1].axis('tight')
    ax[1].set_yticklabels([])
    ax[1].set_xticklabels([])
    ax[1].axis('off')
    return f



def load_qualtax(file_name):
    select_cols = ['Completeness', 'Contamination', 'classification', 'N50 (scaffolds)', 
                   'Domain', 'Phylum', 'Class', 'Order', 'Family', 'Genus', 'Species']

    with open(file_name, "rb" ) as fh:
        qualtax = pickle.load(fh)

    #qualtax = qualtax[select_cols]
    return qualtax

def locate_annotated_files(annotation_folder):
    all_folders = os.listdir(annotation_folder)
    some_folders = [i for i in all_folders if not i.endswith('.ALT')]
    print("{} genome annotations identified".format(len(some_folders)))
    return some_folders 

def load_gene_annotations():
    # read in the annotated product 
    with open( '../data/outputs/gene_annotations2.p', "rb" ) as fh:
        clusterDF_ai = pickle.load(fh)

    print(clusterDF_ai.columsn)

pathway_lookup = "/Volumes/KeithSSD/SulFox/data/minPathFiles/updated.pathways.txt"
read_genome_src2 = lambda x: list(pd.read_csv(x, sep="\t").iloc[:, :2].apply(tuple, axis=1))
allPathNames = {j.strip():i.strip() for i, j in read_genome_src2(pathway_lookup)}

class Genome:
    def __init__(self, report_file, base_path):
        self.accession = report_file.replace(".report.txt", '')
        
        self.report_file = report_file

        self.report_path = os.path.join(base_path, report_file)

        self.name_alt = None; self.name_ = None

        self.pathway_set = None;

    def __str__(self):
        start_str = "\n"
        for k, v in self.__dict__.items():
            if k in ['pred_paths', 'true_paths', 'tp', 'fn', 'fp']:
                start_str = start_str + "{}\t{}, {}\n".format(k, type(v), len(self.pred_paths))
            else:
                start_str = start_str + "{}\t{}\n".format(k, v)

        return start_str

    def set_name(self, n):
        self.name_alt = n
        self.name_ = n.replace(" ", '_')
        return self
    
    def parse_predicted_pathways(self):
        expected_cols = ['path', 'any', 'naive', 'minpath', 'fam0', 'fam-found', 'name']

        all_recs, rec_cntr = {}, 0
        with open(self.report_path, 'r') as fh:
            for line_ in fh:
                rec_dict = {}
                cline_ = line_.split()
                for idx_, nam_ in zip(range(1, 15, 2), expected_cols):
                    try:
                        rec_dict[nam_] = cline_[idx_]
                    except IndexError:
                        print(cline_)
                        print(self.report_file)
                        raise IndexError('list index out of range')

                all_recs[rec_cntr] = rec_dict
                rec_cntr += 1

        pred_paths = pd.DataFrame(all_recs).T.drop(['any'], axis=1)
        assert list(pred_paths['path']) == list(pred_paths['name'])
        self.pathway_set = [Pathway(s) for s in set(pred_paths['name']) if s in allPathNames.keys()]
        print("\t{} pathways discovered: {}".format(self.accession, len(self.pathway_set)))
        return self

    def get_details(self, detail_path):
        self.details_path = detail_path + "/" + self.accession + '.details.txt'
        this_ec = None
        pathway_to_ec = {}
        with open(self.details_path, 'r') as fh:
            for line in fh:
                if line.startswith("path "):
                    this_path = line.split(" # ")[-1].strip()
                else:
                    this_ec = line.split()[-1].strip()
                
                path_set = pathway_to_ec.setdefault(this_path, set())
                if this_ec:
                    path_set.add(this_ec)
                    pathway_to_ec[this_path] = path_set

        captured = set()
        for s in self.pathway_set:
            if s.accession in pathway_to_ec.keys():
                s.ec_numbers.update(pathway_to_ec[s.accession])
                captured.add(s.accession)

        left_behind = set(pathway_to_ec.keys()) - captured
        del pathway_to_ec
        return self

    def find_gff(self, annotation_folder):
        gff_path_std = os.path.join(annotation_folder, self.accession)
        gff_path_alt = os.path.join(annotation_folder, self.accession + '.ALT')
        gff_folders = [gff_path_std, gff_path_alt]
        self.gffs = []
        for a_gff in gff_folders:
            if os.path.exists(a_gff):
                gff_files = [i for i in os.listdir(a_gff) if i.endswith('.gff')]
                if len(gff_files) > 0:
                    gff_file = os.path.join(a_gff, gff_files[0])
                    self.gffs.append(gff_file)

        #print("\t" + self.accession, 'has', len(self.gffs), 'gff files')
        return self

    def parse_gffs(self, target_pathways):
        target_df = pd.read_csv(target_pathways)
        self.functional_proteins = {}
        counter = 0
        for gff in self.gffs:
            for row in parse_gff(gff):
                self.functional_proteins[counter] = row
                counter += 1 

        self.functional_proteins = pd.DataFrame(self.functional_proteins).T
        self.ec_to_pathname = {}
        self.ec_to_pathnum = {}
        for s in self.pathway_set:
            for ec in s.ec_numbers:
                self.ec_to_pathname[ec] = s.name
                self.ec_to_pathnum[ec] = s.accession

        self.functional_proteins['path_accession'] = self.functional_proteins['eC_number'].map(self.ec_to_pathnum)
        self.functional_proteins['path_name'] = self.functional_proteins['eC_number'].map(self.ec_to_pathname)


        bool1 = self.functional_proteins['path_accession'].isin(target_df['Pathway Code'])
        bool2 = self.functional_proteins['path_name'].isin(target_df['Pathway Names'])
        self.functional_proteins['metacyc_target'] = pd.Series(index=bool1.index, data=[False]*len(bool1), dtype=bool)
        self.functional_proteins.loc[(bool1 | bool2), 'metacyc_target'] = True
        print(f"{self.functional_proteins.metacyc_target.sum()} pathway associated proteins (of interest)")
        self.functional_proteins['genome'] = pd.Series(index=bool1.index, data=[self.accession]*len(bool1), 
                                                       dtype=object)
        return self

def merge_identical_rows(mc_pathway_df):
    dt = np.dtype('U14')
    tag_names, tag_counts = np.unique(mc_pathway_df.locus_tag, return_counts=True)
    mc_unq = mc_pathway_df[mc_pathway_df.locus_tag.isin(tag_names[tag_counts == 1])]
    need_merging = tag_names[tag_counts > 1].astype(dt)
    locus_tags = mc_pathway_df.locus_tag.values.astype(dt)
    mc_index = mc_pathway_df.index.values
    start_ = time.perf_counter()
    dup_rows = []
    for i in need_merging:
        lt_bool = locus_tags == i
        lt_tags = mc_index[lt_bool]
        row_pair = mc_pathway_df.loc[lt_tags, :]
        null_total = row_pair.isnull().sum(1)
        # one has less information ?
        info_check = len(set(null_total)) == 2
        # the gene is assigned to different pathways ? 
        path_check = len(set(row_pair['path_name'].dropna())) < 2
        if info_check and path_check:
            dup_rows.append(mc_pathway_df.loc[null_total.sort_values().index[1], :].copy())
        else:
            halved_row = row_pair.iloc[0, :].copy()
            for col in ['product', 'eC_number', 'gene', 'db_xref', 'path_accession', 'path_name']:
                halved_row[col] = "/".join(list(row_pair[col].dropna()))
            halved_row['metacyc_target'] = row_pair['metacyc_target'].iloc[0] or row_pair['metacyc_target'].iloc[1]
            dup_rows.append(halved_row)
            
    end_ = time.perf_counter()
    print('Merging time: {:.2f} min'.format())
    superdf = pd.concat(dup_rows + [mc_unq], axis=1)
    return superdf



def parse_gff(f):
    fasta_f = f'{f}.temp.fasta'
    column_positions = {0:'contig', 3:'start', 4:'end', 6:'strand'}
    addl_columns = set(['locus_tag', 'gene', 'product', 'eC_number', 'db_xref'])
    
    in_gff_section = True
    all_records = []
    with open(fasta_f, 'w') as fhw:
        with open(f, 'r') as fh:
            for line in fh:
                if in_gff_section:
                    if (not line.startswith('##')) and ('hypothetical protein' not in line):
                        row_dict = {}
                        cells = line.strip().split("\t")
                        for a_idx, elem in enumerate(cells[:-1]):
                            c_idx = a_idx
                            if c_idx in column_positions.keys():
                                row_dict[column_positions[c_idx]] = elem
                        
                        more_cells = cells[-1].split(";")
                        for b_idx, elem in enumerate(more_cells):
                            c_idx = b_idx + a_idx + 1
                            colKey, colVal = elem.split("=")
                            if colKey in addl_columns:
                                row_dict[colKey] = colVal
                        
                        all_records.append(row_dict)
                    elif line.startswith('##FASTA'):
                        in_gff_section = False
                else:
                    _ = fhw.write(line)



    # break off fasta file

    # read fasta file 
    # iterate over list, check contig, add sequence 
    return all_records


def aggregate_and_filter_minpath_data(report_files_path, suffix, annotation_folder, target_pathways, desired):
    genomes_ = []
    for f in os.listdir(report_files_path):
        bin_name = os.path.basename(f).replace(suffix, '')
        if f.endswith(suffix) and bin_name in desired:
            g = Genome(f, report_files_path)
            g = g.parse_predicted_pathways()
            g = g.get_details(report_files_path)
            g = g.find_gff(annotation_folder)
            genomes_.append(g.parse_gffs(target_pathways))
    
    just_data = [g.functional_proteins.copy() for g in genomes_]
    pathway_df = pd.concat(just_data, ignore_index=True)
    assert len(pathway_df) == sum([len(j) for j in just_data])
    assert len(pathway_df.columns) == len(just_data[0].columns)
    return pathway_df


class Pathway:
    def __init__(self, name):
        self.accession = name
        self.name = allPathNames[name]
        self.ec_numbers = set()
        self.locus_tags = set()


def load_kegg_database(kegg_fegenie_annots, target_genomes, prot_dir):
    with open(kegg_fegenie_annots , "rb" ) as fh:
        kf_pathway_df = pickle.load(fh)    
    kf_pathway_df = kf_pathway_df[kf_pathway_df.query_id.notnull() | kf_pathway_df.HMM.notnull()]
    kf_pathway_df = kf_pathway_df[kf_pathway_df.BinID.isin(target_genomes)]
    kf_pathway_df['contig'] = kf_pathway_df['GeneID'].apply(lambda x: "_".join(x.split("_")[:-1]))
    
    start_end = {}
    for g in kf_pathway_df.BinID.unique():
        with open(f"{prot_dir}/{g}.fa-proteins.faa") as ph:
            for line in ph:
                if line.startswith(">"):
                    prot_, start_, end_ = [c.strip() for c in line[1:].split("#")[:3]]
                    assert not (prot_ in start_end.keys())
                    start_end[prot_] = (start_, end_)
    
    kf_pathway_df['start'] = kf_pathway_df['GeneID'].map(lambda x: start_end[x][0])
    kf_pathway_df['end'] = kf_pathway_df['GeneID'].map(lambda x: start_end[x][1])
    return kf_pathway_df.rename(columns={'BinID':'genome'})



def add_more_pathways(pathway_df):
    key_orgs = [i[:-3] for i in os.listdir("/Volumes/KeithSSD/SulFox/data/bins_of_interest") if i.endswith('.fa')]
    pathways_of_interest = sorted(list(pathway_df.loc[pathway_df.metacyc_target, 'path_name'].unique()))
    bins_of_interest = sorted(list(pathway_df.genome.unique()))

    kegg_codes_of_interest = list(SPECIFIC_KEGG_CODES.keys())
    to_remove = ['sulfur carrier complex (dsrEFH)', 'succinate dehydrogenase (SDH)', 'ATP citrate (pro-S)-lyase (ACLY)', 
                 '2-oxoglutarate dehydrogenase (sucB)', 'malate dehydrogenase (mqo)', 'NIR_SIR alternative dsrAB model', 
                 'methane monooxygenase (mmoXYZBCD)', 'methanol dehydrogenase (xoxF/mxaF)']
    kegg_codes_of_interest = [k for k in kegg_codes_of_interest if k not in to_remove]
    metabolism_df = pd.DataFrame(index=bins_of_interest, columns=pathways_of_interest+kegg_codes_of_interest).fillna(0)

    for b in bins_of_interest:
        bbool = pathway_df.genome == b
        for p in pathways_of_interest:
            pbool = pathway_df.path_name == p
            if len(pathway_df[bbool & pbool]) > 0:
                metabolism_df.loc[b, p] = 1
        for k in kegg_codes_of_interest:
            assoc_codes = SPECIFIC_KEGG_CODES[k]
            limit = ((len(assoc_codes)*(1/3))//1)+1
            kbool = pathway_df.query_id.isin(assoc_codes)
            if len(set(pathway_df.loc[(bbool & kbool), 'query_id'])) >= limit:
                metabolism_df.loc[b, k] = 1


def categorical_pivot(pathdf, col_col, row_col):
    df = pathdf.copy()
    df = df[df[col_col].notnull()].copy()
    df['Intercept'] = pd.Series(index=df.index, data=[1]*len(df), dtype=np.int64)
    flipDF = df.pivot_table(index=row_col, columns=col_col, values='Intercept', aggfunc='sum', fill_value=0)
    flipDF = (flipDF > 0).astype(int)
    return flipDF

def pull_genes(fulldf, sourcedir, outdir, extension='.fa'):
    ids = fulldf[['contig', 'end', 'start']].apply(tuple, axis=1)
    ids = pd.concat([ids, fulldf['genome']], axis=1)
    idsets = ids.groupby('genome').agg(set)[0]
    
    gene_hash = {}
    for b in idsets.index:
        targets_to_extract = idsets[b]
        gene_est = len(targets_to_extract)
        sourcefile = f"{sourcedir}/{b}{extension}"
        destfile = f"{outdir}/{b}.genes.fa"
        f = SeqIO.to_dict(SeqIO.parse(sourcefile, 'fasta'))
        if not os.path.exists(destfile):
            fhw = open(destfile, 'w')
        
        for t, e, s in targets_to_extract:
            header_name = f'{b}_{t}_{s}_{e}'
            subseq = str(f[t].seq)[int(s)-1:int(e)]
            if not os.path.exists(destfile):
                if len(subseq) > 100:
                    _ = fhw.write(">{}\n{}\n".format(header_name, subseq))
                else:
                    gene_est -= 1
            gene_hash[(b, t, s, e,)] = subseq
        print("Wrote {} ({:.1%}) sequences to {}".format(gene_est, gene_est/len(targets_to_extract), destfile))
    return gene_hash

def load_salmon_data(salmon_folder, pathway_df):
    # genes used to normalize
    housekeepers = ["proC", "recA", "rpoD", "rho", "glyA", "tpiA", "recF"]
    # 
    hk_per_bin = defaultdict(set)
    for h in housekeepers:
        hbool = (pathway_df.gene == h)
        gene_ids = pathway_df.loc[hbool, ['genome', 'contig', 'start', 'end']].apply(lambda x: "_".join(list(x)), axis=1).to_frame()
        gene_ids = (pathway_df.loc[hbool, ['genome']]).join(gene_ids)
        gene_ids = gene_ids.groupby('genome').agg(set)[0]
        for g in gene_ids.index:
            hk_per_bin[g].update(gene_ids[g])
    
    for g in hk_per_bin:
        assert len(hk_per_bin[g]) > 0
    
    salmon_files = defaultdict(list)
    for i in os.listdir(salmon_folder):
        file_i = os.path.join(salmon_folder, i, 'quant.sf')
        sample_i = i.split("_pass")[0]
        bin_i = i.split("pass_")[1].replace(".quant", "")
        df = pd.read_csv(file_i, sep="\t", index_col=0).loc[:, ['TPM']].rename(columns={'TPM':sample_i})
        salmon_files[bin_i].append(df)
    
    salmon_dfs = {b:pd.concat(v, axis=1, verify_integrity=True) for b, v in salmon_files.items()}
    normed_bins = {bin_i:a_df.div(a_df.loc[hk_per_bin[bin_i], :].T.mean(1)).fillna(0).replace(np.inf, 0) for bin_i, a_df in salmon_dfs.items()}
    salmon_df = pd.concat(normed_bins.values(), verify_integrity=True)
    salmon_df['merger_col'] = pd.Series({i:tuple(i.split("_")) for i in salmon_df.index})
    return salmon_df, salmon_dfs, normed_bins


def parse_annotations_into_matrix(pathway_df):
    sarNO = ~pathway_df['product'].astype(str).apply(lambda x: 'arcosine' in x)
    
    annotations = {'psrA': pathway_df['product'].astype(str).str.lower().apply(lambda x: 'polysulfide reductase' in x) | 
                           (pathway_df['gene'].astype(str).str.startswith('psrA') & ~pathway_df['product'].astype(str).str.lower().str.contains('fatty')),
                   'phs': pathway_df.gene.astype(str).str.contains('phs'),
                   'fccA': pathway_df.gene.astype(str).str.startswith('fccA') & ~pathway_df['product'].astype(str).str.contains('umarate'),
                   'fccB': pathway_df.gene.astype(str).str.startswith('fccB') & ~pathway_df['product'].astype(str).str.contains('umarate'),
                   'hydG': pathway_df['product'].astype(str).apply(lambda x: 'ulfhydrogenase' in x) | (pathway_df.query_id.astype(str) == 'K17995'),
                   'aprA': pathway_df.gene.astype(str).str.startswith('aprA') | (pathway_df.query_id.astype(str) == 'K00394'),
                   'aprB': pathway_df.gene.astype(str).str.startswith('aprB') | (pathway_df.query_id.astype(str) == 'K00395'),
                   'asrA': pathway_df.gene.astype(str).str.startswith('asrA') | (pathway_df.query_id.astype(str) == 'K16950'),
                   'asrB': pathway_df.gene.astype(str).str.startswith('asrB') | (pathway_df.query_id.astype(str) == 'K16951'),
                   'dsrE/tusD': pathway_df.query_id.isin(['K07235']) | pathway_df.gene.astype(str).str.startswith("dsrE") | pathway_df['product'].astype(str).str.lower().str.contains('tusd') | pathway_df.gene.astype(str).str.startswith("tusD"),
                   'dsrF/tusC': pathway_df.query_id.isin(['K07236']) | pathway_df.gene.astype(str).str.startswith("dsrF") | pathway_df['product'].astype(str).str.lower().str.contains('tusc') | pathway_df.gene.astype(str).str.startswith("tusC"),
                   'dsrH/tusB': pathway_df.query_id.isin(['K07237']) | pathway_df.gene.astype(str).str.startswith("dsrH") | pathway_df['product'].astype(str).str.lower().str.contains('tusb') | pathway_df.gene.astype(str).str.startswith("tusB"),
                   'tsdA': (pathway_df.query_id.astype(str) == 'K19713'),
                   'sat/cysN/cysD/cysC/cysNC/cysH': pathway_df.query_id.isin(['K00958', 'K00860', 'K00955', 'K00956', 'K00957','K00390']) &  ~pathway_df['product'].astype(str).apply(lambda x: 'alanine' in x) & ~pathway_df['product'].astype(str).apply(lambda x: 'elongation' in x),
                   'soeA': pathway_df.gene.astype(str).str.startswith('soe') & (pathway_df.query_id.astype(str) == 'K21307'),
                   'soeB': pathway_df.gene.astype(str).str.startswith('soe') & (pathway_df.query_id.astype(str) == 'K21308'),
                   'soeC': pathway_df.gene.astype(str).str.startswith('soe') & (pathway_df.query_id.astype(str) == 'K21309'),
                   'sqr': (pathway_df.query_id.astype(str) == 'K17218'),
                   'soxL': pathway_df.gene.astype(str).str.startswith('soxL'),
                   'msrB': pathway_df.gene.astype(str).str.startswith('msrB'),
                   'msrA': pathway_df.gene.astype(str).str.startswith('msrA'),
                   'dsrC': (pathway_df.gene.astype(str).str.startswith('dsrC') | pathway_df.gene.astype(str).str.startswith('tusE') | (pathway_df['product'] == 'sulfite reductase subunit gamma') | pathway_df['product'].astype(str).str.lower().str.contains('tuse')),
                   'soxY': ((pathway_df.query_id.astype(str) == 'K17226') | pathway_df['gene'].astype(str).str.startswith('soxY') | pathway_df['product'].astype(str).str.lower().str.contains('soxy')) & sarNO,
                   'soxZ': ((pathway_df.query_id.astype(str) == 'K17227') | pathway_df['gene'].astype(str).str.startswith('soxZ') | pathway_df['product'].astype(str).str.lower().str.contains('soxz')) & sarNO,
                   "dsrA": (pathway_df.query_id.astype(str) == 'K11181') | pathway_df['gene'].astype(str).str.startswith('dsrB') | pathway_df['gene'].astype(str).str.startswith('dsvB'),
                   "dsrB": (pathway_df.query_id.astype(str) == 'K11180') | pathway_df['gene'].astype(str).str.startswith('dsrA') | pathway_df['gene'].astype(str).str.startswith('dsvA'),
                   "dsrL": (pathway_df['gene'].astype(str).str.startswith('dsrL') | pathway_df['product'].astype(str).str.lower().str.contains('dsrL')),
                   'soxA': ((pathway_df.query_id.astype(str) == 'K17222') | pathway_df['gene'].astype(str).str.startswith('soxA') | pathway_df['product'].astype(str).str.lower().str.contains('soxa')) & sarNO,
                   'soxX': ((pathway_df.query_id.astype(str) == 'K17223') | pathway_df['gene'].astype(str).str.startswith('soxX') | pathway_df['product'].astype(str).str.lower().str.contains('soxx')) & sarNO,
                   'soxB': ((pathway_df.query_id.astype(str) == 'K17224') | pathway_df['gene'].astype(str).str.startswith('soxB') | pathway_df['product'].astype(str).str.lower().str.contains('soxb')) & sarNO,
                   'soxD': ((pathway_df.query_id.astype(str) == 'K22622') | pathway_df['gene'].astype(str).str.startswith('soxD') | pathway_df['product'].astype(str).str.lower().str.contains('soxd')) & sarNO,
                   'soxC': ((pathway_df.query_id.astype(str) == 'K17225') | pathway_df['gene'].astype(str).str.startswith('soxC') | pathway_df['product'].astype(str).str.lower().str.contains('soxc')) & sarNO,
                   'rdhA/thtR/sseB': ((pathway_df.query_id.astype(str) == 'K01011') | pathway_df['gene'].astype(str).str.startswith('rdhA')),
                   'sir/cysJI': (pathway_df['gene'].astype(str).str.startswith('sir_') | (pathway_df['gene'].astype(str) == 'sir') | pathway_df['gene'].astype(str).str.startswith('cysJ') | pathway_df['gene'].astype(str).str.startswith('cysI')) | pathway_df.query_id.isin(['K00381', 'K00380', 'K00392']),
                   'rbcS': (pathway_df.query_id.astype(str) == 'K01602') | (pathway_df.gene.astype(str).apply(lambda x: 'ibulose' in x) & pathway_df.gene.astype(str).apply(lambda x: 'isphosphate' in x) & pathway_df.gene.astype(str).apply(lambda x: 'small' in x.lower())),
                   'rbcL': (pathway_df.query_id.astype(str) == 'K01601') | (pathway_df.gene.astype(str).apply(lambda x: 'ibulose' in x) & pathway_df.gene.astype(str).apply(lambda x: 'isphosphate' in x) & pathway_df.gene.astype(str).apply(lambda x: 'large' in x.lower()))}
    
    complexes = {'fccAB': ("fccA", "fccB",), #> 1
                 'aprAB': ("aprA", "aprB",), #> 1
                 'msrAB': ("msrA", "msrB",), #> 1
                 'asrAB': ("asrA", "asrB",), #> 1
                 'dsrEFH/tusBCD': ("dsrE/tusD", "dsrF/tusC", "dsrH/tusB",), # > 1
                 'soeABC': ("soeA", "soeB", "soeC",), # > 1
                 'soxYZ': ("soxY", "soxZ",),
                 'dsrAB': ("dsrA", "dsrB",),
                 'soxAX': ("soxA", "soxX",),
                 'soxCD': ("soxD", "soxC",),
                 'rbcSL': ("rbcS", "rbcL",)}
    
    quinol_oxidases_names = ['appB', 'appC', 'cioB', 'cydA', 'cydB', 'cydX', 'cyoA', 'cyoB', 'cyoC', 'cyoD', 'qoxC']
    quinol_oxidases_ids = ['K00426', 'K00425', 'K02299']
    high_o2_oxidases = ['K02276', 'K02298', 'K02297', 'K02300']

    full_pathways = {'taurine_degradation': (pathway_df.path_name.astype(str).str.contains('taurine degradation'), 1, 'locus_tag'),
                     'coxLMS/cutLMS': (pathway_df['product'].astype(str).str.contains('arbon monoxide dehydrogenase'), 1, 'locus_tag'),
                     'appBC/cydABX/qoxC':   (pathway_df.query_id.isin(quinol_oxidases_ids) | 
                                             pathway_df.gene.astype(str).str.startswith('appB') |  # anoxic conditions
                                             pathway_df.gene.astype(str).str.startswith('appC') |  # anoxic conditions
                                             pathway_df.gene.astype(str).str.startswith('cioB') | # microaerobic
                                             pathway_df.gene.astype(str).str.startswith('cydA') | # microaerobic
                                             pathway_df.gene.astype(str).str.startswith('cydB') | # microaerobic
                                             pathway_df.gene.astype(str).str.startswith('cydX') | # microaerobic
                                             pathway_df.gene.astype(str).str.startswith('qoxC')  & 
                                             (pathway_df.query_id != "K02033"), 1, 'locus_tag'),
                     'accABCD/pccAB': (pathway_df.gene.astype(str).str.startswith('accA') | 
                                       pathway_df.gene.astype(str).str.startswith('accB') | 
                                       pathway_df.gene.astype(str).str.startswith('accC') | 
                                       pathway_df.gene.astype(str).str.startswith('accD') | 
                                       pathway_df.gene.astype(str).str.startswith('pccA') | 
                                       pathway_df.gene.astype(str).str.startswith('pccB') , 1, 'locus_tag'),
                     "napAB+narIJ": ((pathway_df.query_id.isin(NITROGEN_CYCLE['nitrate reductase (napAB/narIJ)']) & pathway_df.query_id.notnull()), 1, 'query_id'),
                     'amoABC': ((pathway_df.query_id.isin(NITROGEN_CYCLE['ammonia oxidation (amoABC)']) & pathway_df.query_id.notnull()), 1, 'query_id'),
                     'nirK/nirS+norBC': (pathway_df.query_id.isin({'K15864', 'K04561', 'K00368', 'K02305',}) |
                                         pathway_df.gene.astype(str).str.startswith('norB') | 
                                         pathway_df.gene.astype(str).str.startswith('norC') |
                                         pathway_df.gene.astype(str).str.startswith('nirK') |
                                         pathway_df.gene.astype(str).str.startswith('nirS'), 1, 'query_id'),
                     'nosZ': ((pathway_df.query_id == 'K00376') | pathway_df.gene.astype(str).str.startswith('nosZ'), 1, 'query_id'),
                     'nirABD': (pathway_df.query_id.isin({'K00363', 'K00362', 'K00366'}) |
                                      pathway_df.gene.astype(str).str.startswith('nirA') |
                                      pathway_df.gene.astype(str).str.startswith('nirB') |
                                      pathway_df.gene.astype(str).str.startswith('nirD'), 1, 'query_id'),
                     'nrfAH': (pathway_df.query_id.isin({'K15876', 'K03385'}) | pathway_df.gene.astype(str).str.startswith('nrfH') |
                               pathway_df.gene.astype(str).str.startswith('nrfA'), 0, 'query_id'),
                     'nifKHDW/vnfA': (pathway_df['product'].astype(str).str.lower().str.contains('nitrogenase'), 1, 'query_id'),
                     'coxABC/ctaCDE/cyoABCD': ((pathway_df['product'].astype(str).str.lower().str.contains('cytochrome c oxidase') & 
                                               (pathway_df['gene'].astype(str).str.startswith('cta') | pathway_df['gene'].astype(str).str.startswith('cox'))) |
                                              (pathway_df.gene.astype(str).str.startswith('cyoA') | 
                                               pathway_df.gene.astype(str).str.startswith('cyoB') | 
                                               pathway_df.gene.astype(str).str.startswith('cyoC') | 
                                               pathway_df.gene.astype(str).str.startswith('cyoD'))  | 
                                               pathway_df.query_id.isin(high_o2_oxidases), 1, 'locus_tag'),
                     'katGBE/sodABC/dfx':(pathway_df.gene.astype(str).str.startswith('grxC') | 
                                      pathway_df.gene.astype(str).str.startswith('dfx') | 
                                      pathway_df['product'].astype(str).str.lower().str.contains('alkylhydroperoxidase') | 
                                      pathway_df['product'].astype(str).str.lower().str.contains('oxide dismutase') | 
                                      pathway_df['product'].astype(str).str.lower().str.contains('oxid dismutase') | 
                                      pathway_df['product'].astype(str).str.lower().str.contains('catalase'), 1, 'locus_tag'),
                     'psaABCDEFIJKL':(pathway_df.gene.astype(str).str.startswith('psa'), 9, 'locus_tag'),
                     'petABCDGLMN':(pathway_df.gene.astype(str).str.startswith('pet'), 10, 'locus_tag'),
                     'pufABCLMX':(pathway_df.gene.astype(str).str.startswith('puf'), 1, 'locus_tag'),
                     'psbA-Z':(pathway_df.gene.astype(str).str.startswith('psb'), 3, 'locus_tag'),
                     'hoxS/hndABC': (pathway_df.gene.astype(str).str.startswith('hox') | 
                                     pathway_df.gene.astype(str).str.startswith('hnd'), 1, 'locus_tag'),
                     'hupLBS/hoxZKG': (pathway_df.gene.astype(str).str.startswith('hup'), 1, 'locus_tag'),
                     'cpeABCDERSTYZ': (pathway_df.gene.astype(str).str.startswith('cpe'), 1, 'locus_tag'),
                     'fdnG/fdhF/fdhAB/fdoGHI/fdsAB': (pathway_df['product'].astype(str).str.lower().str.contains('formate dehydrogenase') & 
                                        ~pathway_df['product'].astype(str).str.lower().str.contains('mitochondrial'), 1, 'locus_tag'),
                     'hao': (pathway_df.gene.astype(str).str.startswith('hao'), 0, 'locus_tag')}
    
    in_complex = set()
    is_complex = set()
    genomes_with = {}
    genes_within = {}
    for k, v in complexes.items():
        is_complex.add(k)
        in_complex.update([i for i in v])
        
        if len(v) == 2:
            annotations[k] = annotations[v[0]] | annotations[v[1]]
        
            have_complex = np.unique(list(pathway_df[annotations[v[0]]]['genome'].unique()) + \
                                      list(pathway_df[annotations[v[1]]]['genome'].unique()), return_counts=True)
        elif len(v) == 3:
            annotations[k] = annotations[v[0]] | annotations[v[1]] | annotations[v[2]]
            have_complex = np.unique(list(pathway_df[annotations[v[0]]]['genome'].unique()) + \
                                      list(pathway_df[annotations[v[1]]]['genome'].unique()) + \
                                      list(pathway_df[annotations[v[2]]]['genome'].unique()), return_counts=True)
        
        genes_within[k] = set(pathway_df.loc[annotations[k], 'identifier'])
        genomes_with[k] = set(have_complex[0][have_complex[1] > 1])
    
    for k, v in annotations.items():
        if (k not in in_complex) and (k not in is_complex):
            genomes_with[k] = set(pathway_df[v]['genome'])
            genes_within[k] = set(pathway_df.loc[annotations[k], 'identifier'])
    
    for fp, (bool_i, cutoff_i, query_col) in full_pathways.items():
        agg_df = pathway_df.loc[bool_i, ['genome', query_col]].groupby('genome').agg(lambda x: len(set(x.values)))[query_col]
        genomes_with[fp] = set(agg_df.index[agg_df > cutoff_i])
        genes_within[fp] = set(pathway_df.loc[bool_i & pathway_df.genome.isin(genomes_with[fp]), 'identifier'])
    
    new_annotations = {}
    for key, val in genomes_with.items():
        new_annotations[key] = {k:True if k in val else False for k in pathway_df.genome.unique()}
    
    annot_mat = pd.DataFrame(new_annotations)
    return annot_mat, genes_within 

hewson_translator = {'SRR1011421': '9/21/11_Oxic_17m_1',
                     'SRR1011422': '9/21/11_Oxic_17m_2',
                     'SRR980293': '5/17/10_Oxic_13m_1',
                     'SRR984833': '5/17/10_Oxic_13m_2',
                     'SRR984835': '5/17/10_Oxic_3m_1',
                     'SRR984837': '5/17/10_Oxic_3m_2',
                     'SRR984839': '6/7/10_Anox_16.5m_1',
                     'SRR984891': '6/7/10_Anox_16.5m_2',
                     'SRR985022': '7/11/10_Anox_17.5m_1',
                     'SRR985023': '7/11/10_Anox_17.5m_2',
                     'SRR985024': '7/11/10_Oxic_3m_1',
                     'SRR985025': '7/11/10_Oxic_3m_2',
                     'SRR985026': '8/5/10_Anox_22m_1',
                     'SRR985027': '8/5/10_Anox_22m_2',
                     'SRR985028': '8/30/10_Oxic_3m_1',
                     'SRR985029': '8/30/10_Oxic_3m_2',
                     'SRR985030': '8/30/10_Anox_20m_1',
                     'SRR985031': '8/30/10_Anox_20m_2',
                     'SRR985301': '10/18/10_Oxic_13m_1',
                     'SRR985423': '10/18/10_Oxic_13m_2',
                     'SRR988005': '7/08/11_Oxic_3m_1',
                     'SRR988006': '7/08/11_Oxic_3m_2',
                     'SRR988007': '5/24/11_Anox_18m_1',
                     'SRR988008': '4/18/11_Oxic_13.5m_1',
                     'SRR988009': '4/18/11_Oxic_13.5m_2',
                     'SRR988014': '5/24/11_Anox_18m_2',
                     'SRR988024': '6/14/11_Oxic_3m_2',
                     'SRR988025': '6/14/11_Anox_18m_1',
                     'SRR988026': '6/14/11_Anox_18m_2'}


def load_cn_data(data_file, duplicate_file=None):
    cov_df = pd.read_csv(data_file, sep="\t", index_col=0)
    included_rows = set(cov_df.index)
    add_back = []
    if duplicate_file:
        print("Adding back duplicates")
        dup_df = pd.read_csv(duplicate_file, sep="\t")
        for idx in range(len(dup_df)):
            retained_, removed_ = list(dup_df.iloc[idx, [0, 1]].values)
            retained_row = cov_df.loc[retained_, :].copy()
            retained_row.name = removed_
            if not removed_ in included_rows:
                add_back.append(retained_row)
            else:
                assert sum(cov_df.loc[removed_, :].values == retained_row.values) == len(retained_row)
    
    if len(add_back):
        assert len(add_back) == len(dup_df)
        removed_set = set(dup_df.iloc[:, 1])
        assert len(removed_set) == len(dup_df)
        assert len(removed_set.intersection(set(cov_df.index))) == 0
        add_df = pd.concat(add_back, axis=1).T
        cov_df_plus = pd.concat((cov_df, add_df,))
        assert cov_df_plus.shape[1] == cov_df.shape[1]
        assert len(cov_df_plus) == (len(removed_set) + len(included_rows))
    
        if len(set(hewson_translator.keys).intersection(cov_df_plus.columns)):
            cov_df_plus = cov_df_plus.rename(columns=hewson_translator)
        return cov_df_plus
    else:
        if len(set(hewson_translator.keys()).intersection(cov_df.columns)):
            cov_df = cov_df.rename(columns=hewson_translator)
        return cov_df

        
def make_report(qualtax_df, complete_list, cbTr_abnd2, cb33_abnd2, salmon_df, copy_num_df, gene_id_lookup3, gene_sets):
    report_df = qualtax_df.loc[complete_list, ['Completeness', 'Contamination', 'N50 (scaffolds)']]
    addl_columns = ['taxonomic leaf rank', 'taxonomic leaf name', 'max copy number', 'median nonzero expression', 
                    'genes with sum >5x rpm expression', 'genes with mean >5x rpm expression', 
                    'median copy number (target genes)', 'sum expression (target genes)',
                    '# of enzymes max >1x (exp)', 'enzymes max >1x (exp)',
                    '# of enzymes max >10x (exp)', 'enzymes max >10x (exp)',
                    '# of enzymes max >1 rpm (copy)', 'enzymes max >1 rpm (copy)',
                    '# of enzymes max >10 rpm (copy)', 'enzymes max >10 rpm (copy)']
    
    for c in addl_columns:
        if c.startswith("enzymes"):
            report_df[c] = pd.Series(index=report_df.index, dtype=str, data=['']*len(report_df))
        else:
            report_df[c] = pd.Series(index=report_df.index, dtype=str)
    
    cbTrexp = [i for i in cbTr_abnd2.columns if not 'control' in i.lower()]
    cb33exp = [i for i in cb33_abnd2.columns if not 'ctrl' in i.lower()]
    
    numeric_exp = [i for i in salmon_df.columns if not i in ['merger_col', 'gene_name']]
    
    exp_keys = set(salmon_df.index)
    copy_num_keys = set(copy_num_df.index)
    exp_keys_translated = set([i for i in exp_keys if i in gene_id_lookup3 and gene_id_lookup3[i] in copy_num_keys])
    
    for bof in complete_list:
        done_flag = False
        for c in ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Domain']:
            cname = str(qualtax_df.loc[bof, c])
            if (cname.lower() in ['nan', '']):
                pass            
            elif not done_flag:
                done_flag = True
                report_df.loc[bof, 'taxonomic leaf rank'] = c
                report_df.loc[bof, 'taxonomic leaf name'] = cname
        
        print(bof, report_df.loc[bof, 'taxonomic leaf name'])
        if bof.startswith('CBrest'):
            slice_data = cbTr_abnd2.loc[bof, cbTrexp].copy()
            report_df.loc[bof, 'max copy number'] = np.max(slice_data)
        elif bof.startswith('CB33'):
            slice_data = cb33_abnd2.loc[bof, cb33exp].copy()
            report_df.loc[bof, 'max copy number'] = np.max(slice_data)
        
        all_exp_found = exp_keys.intersection(gene_id_lookup2[bof])
        report_df.loc[bof, 'genes with sum >5x rpm expression'] = (salmon_df.loc[all_exp_found, numeric_exp].sum(1) > 5).sum()
        report_df.loc[bof, 'genes with mean >5x rpm expression'] = (salmon_df.loc[all_exp_found, numeric_exp].mean(1) > 5).sum()
        report_df.loc[bof, 'median nonzero expression'] = np.median([i for i in salmon_df.loc[all_exp_found, numeric_exp].mean(1).values if i > 0])
        
        copy_list = []
        expression_list = []
        for k, pre_keys in gene_sets.items():
            local_keys = [i for i in pre_keys if i.startswith(bof+"_")]
            expindx = set(local_keys) & set(salmon_df.index)
            if len(expindx):
                expvalue = np.max(salmon_df.loc[expindx, numeric_exp].sum())
                expression_list.append(expvalue)
                if expvalue >= 1:
                    report_df.loc[bof, 'enzymes max >1x (exp)'] = " ".join([report_df.loc[bof, 'enzymes max >1x (exp)'], f"{k}"])
                if expvalue >= 10:
                    report_df.loc[bof, 'enzymes max >10x (exp)'] = " ".join([report_df.loc[bof, 'enzymes max >10x (exp)'], f"{k}"])
            
            cnindx = set([gene_id_lookup3[lk] for lk in set(local_keys) & exp_keys_translated]) 
            if len(cnindx):
                cnvalue = np.max(copy_num_df.loc[cnindx, :].sum())
                copy_list.append(cnvalue)
                if cnvalue >= 1:
                    report_df.loc[bof, 'enzymes max >1 rpm (copy)'] = " ".join([report_df.loc[bof, 'enzymes max >1 rpm (copy)'], f"{k}"])
                if cnvalue >= 10:
                    report_df.loc[bof, 'enzymes max >10 rpm (copy)'] = " ".join([report_df.loc[bof, 'enzymes max >10 rpm (copy)'], f"{k}"])
        
        if len(copy_list):
            report_df.loc[bof, 'median copy number (target genes)'] = np.median(np.array(copy_list))
        
        if len(expression_list):
            report_df.loc[bof, 'sum expression (target genes)'] = np.sum(np.array(expression_list))
        
        report_df.loc[bof, '# of enzymes max >1x (exp)'] = len(report_df.loc[bof, 'enzymes max >1x (exp)'].split())
        report_df.loc[bof, '# of enzymes max >10x (exp)'] = len(report_df.loc[bof, 'enzymes max >10x (exp)'].split())
        report_df.loc[bof, '# of enzymes max >1 rpm (copy)'] = len(report_df.loc[bof, 'enzymes max >1 rpm (copy)'].split())
        report_df.loc[bof, '# of enzymes max >10 rpm (copy)'] = len(report_df.loc[bof, 'enzymes max >10 rpm (copy)'].split())
        
    cols_to_rank = ['median nonzero expression', 'max copy number', 'median copy number (target genes)', 'sum expression (target genes)']
    
    rank_cols = [i+" Rank" for i in cols_to_rank]
    
    for cr, rc in zip(cols_to_rank, rank_cols):
        report_df[rc] = report_df[cr].rank()
    
    report_df['mean rank'] = report_df[rank_cols].mean(1)
    return report_df





