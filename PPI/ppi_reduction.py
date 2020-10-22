import argparse

import time

import random

import pandas as pd
import networkx as nx
import numpy as np

def get_neighbors(df, protein):
    """
    Find neighbors based from the edge list for a given node
    df: DataFrame -- master edge table
    protein: string -- protein name
    """
    nei_df1 = df[df['protein1'] == protein]
    neighbors1 = nei_df1['protein2'].tolist()
    nei_df2 = df[df['protein2'] == protein]
    neighbors2 = nei_df2['protein1'].tolist()
    neighbors = neighbors1 + neighbors2
    return neighbors

def main(opt):

    output_name = opt.o
    s = time.time()
    print('...loading data')
    """
    Load master_edge_table and kinase information
    """
    master_edge_table = pd.read_csv('./master.edges')
    raw_pathways = pd.read_csv('./gogo_bpo.groups')
    raw_pathways.rename(columns={'cluster_id':'pathways'}, inplace=True)
    kinases = []
    with open('../../../ski_data/kinases.lst','r') as f:
        for item in f.readlines():
            kinases.append(item.replace('\n',''))
    """
    Read the original node table
    """
    node_table_final = pd.read_csv(opt.i)
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
    nx.write_gexf(GG,  '{}.gexf'.format(output_name))
    new_node_table.to_csv('reduced_{}.nodes'.format(output_name), index=False)
    new_edge_table.to_csv('reduced_{}.edges'.format(output_name), index=False)
    d1 = {'original':[], 'new':[]}
    for k,v in d.items():
        d1['original'].append(k)
        d1['new'].append(v)
    contraction_record = pd.DataFrame(d1)
    contraction_record.to_csv('reduced_{}.records'.format(output_name), index=False)
    
    print(time.time() - s)
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--i', type=str, help='Input node table filename')
    parser.add_argument('--o', type=str, help='Output filename')
    opt = parser.parse_args()
    main(opt)