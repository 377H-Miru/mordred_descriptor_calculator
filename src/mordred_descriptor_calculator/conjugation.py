import networkx as nx
import numpy as np

def calc_conjugation_features(mol):
    """
    Calculate conjugation descriptors for an RDKit Mol object.
    Returns a dictionary of conjugation features.
    """
    res = {
        "Conjugation_Count": 0,
        "Conjugation_MaxAtomCount": 0,
        "Conjugation_MaxLength": 0,
        "Conjugation_BLA": np.nan,
        "Conjugation_GraphEnergy": 0.0
    }
    if mol is None:
        return res
    try:
        conjugated_bonds = [b for b in mol.GetBonds() if b.GetIsConjugated()]
        if not conjugated_bonds:
            return res
        
        has_3d = (mol.GetNumConformers() > 0)
        conf = mol.GetConformer() if has_3d else None
        
        graph = nx.Graph()
        for b in conjugated_bonds:
            u, v = b.GetBeginAtomIdx(), b.GetEndAtomIdx()
            if has_3d and conf is not None:
                w = conf.GetAtomPosition(u).Distance(conf.GetAtomPosition(v))
            else:
                w = 1.0
            graph.add_edge(u, v, weight=w)
            
        systems = list(nx.connected_components(graph))
        if not systems:
            return res
        
        largest = max(systems, key=len)
        sub = graph.subgraph(largest)
        res["Conjugation_Count"] = len(systems)
        res["Conjugation_MaxAtomCount"] = len(largest)
        if sub.number_of_nodes() > 1:
            try:
                res["Conjugation_MaxLength"] = nx.diameter(sub)
            except Exception:
                pass
        
        if has_3d:
            lengths = [d['weight'] for u, v, d in sub.edges(data=True)]
            if lengths:
                res["Conjugation_BLA"] = float(np.std(lengths))
                
        adj = nx.to_numpy_array(sub, weight=None)
        res["Conjugation_GraphEnergy"] = float(np.sum(np.abs(np.linalg.eigvalsh(adj))))
        return res
    except Exception:
        return res
