import os
import time

import random

import pickle
import pandas as pd
import networkx as nx
import numpy as np

import h5py

class GRPreprocessor():
    

output_folder = 'gogo_hierarchical_30'
s = time.time()
print('...loading data')
"""
Set up the necessary paths
"""
root = '../../../ski_data/{}/original_node_tables/'.format(output_folder)
table_path = '../../../ski_data/{}/contracted_tables/'.format(output_folder)
graph_path = '../../../ski_data/{}/contracted_graphs/'.format(output_folder)
record_path = '../../../ski_data/{}/contraction_records/'.format(output_folder)
"""
Load master_edge_table, GR50, affinity data, pathway data, disease scores, and kinase information
"""
master_edge_table = pd.read_csv('../../../ski_data/master_edge_table.csv')
gr50 = pd.read_csv('../../../ski_data/gr50.csv') # only used create the instance
ba = pd.read_csv('../../../ski_data/affinity.csv')
raw_pathways = pd.read_csv('../../../ski_data/gogo_hierarchical_30.csv')
raw_pathways.rename(columns={'cluster_id':'pathways'}, inplace=True)
kinases = []
with open('../../../ski_data/kinases.lst','r') as f:
    for item in f.readlines():
        kinases.append(item.replace('\n',''))
"""
create instance paths based on the combination of different cell line and drug name
"""
combination = ['HCC1395','fascaplysin']
cell_line = combination[0]
drug = combination[1]
instance_id = '_'.join(combination)
ge_folder = '../../../ski_data/gene_expressions/'
ge = pd.read_csv(ge_folder + cell_line + '.csv')
try:
    os.mkdir(table_path + instance_id)
except:
    print('Folder exists')
"""
Generate instance node table (original_node_table)
"""
node_table = ge
affinity = ba.loc[ba['Drug_Name']==drug][['ENSEMBLID','Affinity', 'Class']]
affinity.rename(columns={'ENSEMBLID':'ensembl'}, inplace=True)
z_scores = pd.read_csv('../../../ski_data/disease_scores/{}.csv'.format(cell_line))
node_table_aff = node_table.merge(affinity, how='left') # combine affinity data into node_table
node_table_na = node_table_aff.fillna(0) # fill NaN affinities as 0
node_table_z = node_table_na.merge(z_scores, how='left', on='ensembl') # combine disease scores into node_table
node_table_final = node_table_z.merge(raw_pathways) # combine pathway information into node_table
node_table_final.fillna({'disgenet_score':0, 'disease_score':0}, inplace=True)
node_table_final['kinaseness'] = [1.0]*len(node_table_final) # combine kinase information into node_table
node_table_final.at[node_table_final['ensembl'].isin(kinases), 'kinaseness'] = 2.0
# node_table_final.to_csv(root + instance_id + '.csv', index=False) # save instance
"""
Start the graph contraction
"""
print('...start contraction')
new_edge_table = master_edge_table
# add gene expressions information to edge_table
ged = node_table_final[['ensembl', 'gene_exp']].to_dict('list')
ged1 = dict(zip(ged['ensembl'], ged['gene_exp']))
new_edge_table.at[:,'ge1'] = [ged1[x] for x in new_edge_table['protein1'].tolist()]
new_edge_table.at[:,'ge2'] = [ged1[x] for x in new_edge_table['protein2'].tolist()]
# add pathway information to edge_table
pathways = node_table_final[['ensembl', 'pathways']].to_dict('list')
pathways1 = dict(zip(pathways['ensembl'], pathways['pathways']))
new_edge_table.at[:,'inter'] = new_edge_table.apply(lambda row: pathways1[row['protein1']]==pathways1[row['protein2']], axis=1)
# find edges that are NOT contractable (gene expression not the same | either node is a kinase | pathways not the same)
keepers = new_edge_table[(new_edge_table['protein1'].isin(kinases)) | 
                              (new_edge_table['protein2'].isin(kinases)) | 
                              (new_edge_table['ge1'] != new_edge_table['ge2'])|
                              (new_edge_table['inter'] == False)]
print('Percentage of edges to keep {}'.format(len(keepers)/(len(master_edge_table))))
# Create a temp graph
tmp_G = nx.from_pandas_edgelist(new_edge_table, source='protein1', target='protein2',
                                edge_attr = new_edge_table.columns.tolist()[2])
# Remove the UNcontractable edges
keepers_tuple = [(x,y) for x,y in zip(keepers['protein1'].tolist(), keepers['protein2'].tolist())]
tmp_G.remove_edges_from(keepers_tuple)
# Find connected components that have more than one node
#Gs = list(nx.connected_component_subgraphs(tmp_G))
Gs = [tmp_G.subgraph(c).copy() for c in nx.connected_components(tmp_G)]
connected_Gs = [x for x in Gs if len(x.nodes) > 1]
print('Number of connected subgraphs {}'.format(len(connected_Gs)))
# Actual contraction part
new_node_table = node_table_final.copy()
d = {i:i for i in new_node_table['ensembl'].tolist()}
ids = 0
for G in connected_Gs:
    if ids % (int(len(connected_Gs)/10)) == 0 and ids >0:
        print('Currently at {} connected graphs'.format(ids))
    nodes = list(G.nodes)
    dest = random.choice(nodes)
    tmp_row = new_node_table[new_node_table['ensembl'] == dest]
    idx = tmp_row.index[0]
    node_attrs = new_node_table[new_node_table['ensembl'].isin(nodes)]
    medians = node_attrs.iloc[:,-4:-2].replace(0, np.nan).median() # can be max as well
#    maxs = node_attrs.iloc[:,-4:-1].replace(0, np.nan).max()
    medians.fillna(0, inplace=True)
    tmp_row.at[idx, ['disgenet_score', 'disease_score']] = medians.to_numpy()
    tmp_row.at[idx, 'ensembl'] = 'v{}'.format(ids)

    for item in nodes:
        d[item] = 'v{}'.format(ids)
    new_node_table.loc[new_node_table['ensembl'].isin(nodes),'kinaseness'] = 0.0
    new_node_table.loc[new_node_table['ensembl'].isin(nodes),'ensembl'] = 'v{}'.format(ids)
    new_edge_table.loc[new_edge_table['protein1'].isin(nodes),'protein1'] = 'v{}'.format(ids)
    new_edge_table.loc[new_edge_table['protein2'].isin(nodes),'protein2'] = 'v{}'.format(ids)
    new_edge_table = new_edge_table[new_edge_table['protein1'] != new_edge_table['protein2']]
    ids += 1
    
# Replace edge score with median of all contracted scores
print('...post processing')
new_node_table.drop_duplicates('ensembl', inplace=True)
grouped = new_edge_table.groupby(['protein1','protein2'], as_index=False)
new_edge_table1 = grouped.agg({'combined_score':'median'})
# Remove scenarios where A-B and B-A are the same
new_edge_table1.drop_duplicates(['protein1', 'protein2'], inplace=True)
new_edge_table1['check_string'] = new_edge_table1.apply(lambda row: ''.join(sorted([str(row['protein1']), str(row['protein2'])])), axis=1)
new_edge_table1.drop_duplicates('check_string', inplace=True)
new_edge_table1.drop('check_string', inplace=True, axis=1)
# Replace NaN with 0 for disease scores
new_node_table.fillna({'disgenet_score':0, 'disease_score':0}, inplace=True)
# Create new graph
print('...creating new graphs')
GG = nx.from_pandas_edgelist(new_edge_table, source='protein1', target='protein2',
                            edge_attr = new_edge_table.columns.tolist()[2])
"""
Graph statistics
"""
Nodes = len(GG.nodes)
Edges = len(GG.edges)
print('Number of nodes: {}'.format(Nodes))
print('Number of edges: {}'.format(Edges))
total_degree = sum(dict(GG.degree).values())
avg_degree = total_degree/Nodes
print('Average node degree: {}'.format(avg_degree))
density = 2*Edges/(Nodes*(Nodes-1))
print('Graph density: {}'.format(density))
"""
Assign node attributes
"""
print('...assigning attributes')
attr_d = new_node_table.to_dict(orient='index')
keys = list(attr_d.keys())
for k in keys:
   tmp_id = attr_d[k]['ensembl']
   attr_d[tmp_id] = attr_d.pop(k)
nx.set_node_attributes(GG, attr_d)
"""
Save stuff
"""
nx.write_gexf(GG, graph_path + instance_id + '.gexf')
new_node_table.to_csv(table_path + instance_id + '/' + 'contracted_node_table_{}.csv'.format(instance_id), index=False)
new_edge_table.to_csv(table_path + instance_id + '/' + 'contracted_edge_table_{}.csv'.format(instance_id), index=False)
d1 = {'original':[], 'new':[]}
for k,v in d.items():
    d1['original'].append(k)
    d1['new'].append(v)
contraction_record = pd.DataFrame(d1)
contraction_record.to_csv(record_path + instance_id + '.csv', index=False)
"""
Generate data
"""
# node_ids = new_node_table['ensembl'].tolist()
# G_node_ids = list(GG.nodes)
# new_index = [node_ids.index(x) for x in G_node_ids]
# new_node_table = new_node_table.reindex(index=new_index)

# feature = new_node_table.loc[:,['gene_exp', 'Affinity','disease_score', 'disgenet_score', 'kinaseness']]
# feature_matrix = feature.to_numpy()
# with h5py.File('./test.h5', 'w') as f:
#     f.create_dataset('X', data=feature_matrix, compression='gzip', compression_opts=9)

print(time.time() - s)