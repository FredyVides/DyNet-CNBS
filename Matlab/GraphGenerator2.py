# -*- coding: utf-8 -*-
"""
Created on Fri Sep 22 09:10:48 2023

@author: fredy.vides
"""
def GraphGenerator2(AdjA,AdjB,subset_sizes):
    import itertools
    import matplotlib.pyplot as plt
    import networkx as nx
    from numpy import zeros
    
    #subset_sizes = [18, 66]
    subset_color = [
        "limegreen",
        "darkorange",
        "limegreen",
        "darkorange",
        "limegreen",
        "darkorange",
        "limegreen",
        "darkorange",
    ]
    
    
    def multilayered_graph(*subset_sizes):
        extents = nx.utils.pairwise(itertools.accumulate((0,) + subset_sizes))
        layers = [range(start, end) for start, end in extents]
        G = nx.Graph()
        for i, layer in enumerate(layers):
            G.add_nodes_from(layer, layer=i)
        for layer1, layer2 in nx.utils.pairwise(layers):
            G.add_edges_from(itertools.product(layer1, layer2))
        return G
    
    
    G = multilayered_graph(*subset_sizes)
    A0 = zeros((sum(subset_sizes),sum(subset_sizes)))
    A0[subset_sizes[0]:sum(subset_sizes[:-1]),:subset_sizes[0]] = AdjA
    A0 = A0 + A0.T
    G0 = nx.from_numpy_array(A0, create_using = nx.MultiGraph())
    G.edges = G0.edges
    A0[sum(subset_sizes[:-1]):,subset_sizes[0]:sum(subset_sizes[:-1])] = AdjB
    A0 = A0 + A0.T
    G0 = nx.from_numpy_array(A0, create_using = nx.MultiGraph())
    G.edges = G0.edges
    
    color = [subset_color[data["layer"]] for v, data in G.nodes(data=True)]
    pos = nx.multipartite_layout(G, subset_key="layer")
    plt.figure(figsize=(8, 8))
    nx.draw_networkx(G, pos, node_color=color, with_labels=False)
    plt.show()