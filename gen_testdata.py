#!/usr/bin/env python3
"""
Generate test datasets for SFD validation.

Creates multiple datasets with diverse simplet structures:
- Dataset A: Small (20 vertices) - sparse, tree-like
- Dataset B: Medium (50 vertices) - mixed, with triangles and tetrahedra
- Dataset C: Large (200 vertices) - dense clusters + sparse bridges
- Dataset D: Medium (80 vertices) - scale-free-like (hub structure)

All datasets are designed so that exact enumeration is feasible (< 1 minute)
while being large enough to test V2 sampling quality.
"""

import os
import random

PROJECT_DIR = os.path.dirname(os.path.abspath(__file__))
TESTDATA_DIR = os.path.join(PROJECT_DIR, "testdata")

def save_dataset(name, simplices):
    """Save simplices in DBLP format (nverts.txt + simplices.txt)"""
    dirpath = os.path.join(TESTDATA_DIR, name)
    os.makedirs(dirpath, exist_ok=True)
    
    with open(os.path.join(dirpath, "nverts.txt"), "w") as f:
        for s in simplices:
            f.write(f"{len(s)}\n")
    
    with open(os.path.join(dirpath, "simplices.txt"), "w") as f:
        for s in simplices:
            for v in s:
                f.write(f"{v}\n")
    
    # Count structures
    edges = set()
    tris = set()
    tets = set()
    for s in simplices:
        s = sorted(s)
        if len(s) == 2:
            edges.add(tuple(s))
        elif len(s) == 3:
            tris.add(tuple(s))
            for i in range(3):
                for j in range(i+1, 3):
                    edges.add((s[i], s[j]))
        elif len(s) == 4:
            tets.add(tuple(s))
            for i in range(4):
                for j in range(i+1, 4):
                    edges.add((s[i], s[j]))
                    for k in range(j+1, 4):
                        tris.add((s[i], s[j], s[k]))
    
    vertices = set()
    for s in simplices:
        vertices.update(s)
    
    print(f"  {name}: {len(vertices)} vertices, {len(edges)} edges, "
          f"{len(tris)} triangles, {len(tets)} tetrahedra, "
          f"{len(simplices)} input simplices")
    
    return dirpath


def dataset_a():
    """
    Dataset A: 20 vertices, sparse tree-like + some triangles
    Expected dominant types: 1 (edge), 2 (path), 5 (star), 6 (path-4)
    """
    simplices = []
    
    # Tree backbone: 0-1-2-3-4-5-6-7-8-9
    for i in range(9):
        simplices.append([i, i+1])
    
    # Star from vertex 5: 5-10, 5-11, 5-12
    for v in [10, 11, 12]:
        simplices.append([5, v])
    
    # A few triangles (filled)
    simplices.append([0, 1, 2])    # filled triangle
    simplices.append([10, 11, 12]) # empty triangle (edges but no 3-simplex... wait, we add as 3-simplex)
    
    # Some extra edges for paths
    simplices.append([3, 13])
    simplices.append([7, 14])
    simplices.append([14, 15])
    simplices.append([13, 16])
    simplices.append([16, 17])
    simplices.append([17, 18])
    simplices.append([18, 19])
    
    # Empty triangle (just edges, no filled)
    simplices.append([8, 9])  # already exists
    simplices.append([8, 14])
    # 9-14 not added -> path through 8
    
    return save_dataset("dataset_a", simplices)


def dataset_b():
    """
    Dataset B: 50 vertices, mixed structures with triangles and tetrahedra
    Expected: all 18 types should appear
    """
    simplices = []
    
    # Cluster 1: Dense clique (K5 = vertices 0-4) -> tetrahedra
    for i in range(5):
        for j in range(i+1, 5):
            simplices.append([i, j])
    # Add some filled triangles
    simplices.append([0, 1, 2])
    simplices.append([0, 1, 3])
    simplices.append([0, 2, 3])
    simplices.append([1, 2, 3])
    # Tetrahedron
    simplices.append([0, 1, 2, 3])
    # More triangles
    simplices.append([2, 3, 4])
    simplices.append([1, 3, 4])
    
    # Cluster 2: Triangle mesh (vertices 10-19)
    for i in range(10, 19):
        simplices.append([i, i+1])
    simplices.append([10, 12])
    simplices.append([12, 14])
    simplices.append([14, 16])
    simplices.append([16, 18])
    # Filled triangles
    simplices.append([10, 11, 12])
    simplices.append([12, 13, 14])
    simplices.append([14, 15, 16])
    # Empty triangles exist via edges: 10-11-12 has edges + filled
    # 11-12-13 has edges 11-12, 12-13, but no 11-13 -> path, not triangle
    
    # Cluster 3: Star hub (vertex 25, connected to 26-34)
    for v in range(26, 35):
        simplices.append([25, v])
    # Some edges among leaves
    simplices.append([26, 27])
    simplices.append([28, 29])
    simplices.append([30, 31])
    simplices.append([32, 33])
    # Triangle in star
    simplices.append([25, 26, 27])
    
    # Cluster 4: 4-cycles and paths (vertices 35-44)
    # 4-cycle: 35-36-37-38-35
    simplices.append([35, 36])
    simplices.append([36, 37])
    simplices.append([37, 38])
    simplices.append([38, 35])
    # Another 4-cycle: 39-40-41-42-39
    simplices.append([39, 40])
    simplices.append([40, 41])
    simplices.append([41, 42])
    simplices.append([42, 39])
    # Bridges
    simplices.append([38, 39])
    # Path: 43-44
    simplices.append([43, 44])
    simplices.append([42, 43])
    
    # Bridges between clusters
    simplices.append([4, 10])   # Cluster 1 - 2
    simplices.append([19, 25])  # Cluster 2 - 3
    simplices.append([34, 35])  # Cluster 3 - 4
    simplices.append([44, 0])   # Cluster 4 - 1 (cycle)
    
    # Extra vertices for diversity
    simplices.append([45, 46])
    simplices.append([46, 47])
    simplices.append([47, 48])
    simplices.append([48, 49])
    simplices.append([45, 25])  # Connect to hub
    simplices.append([49, 10])  # Connect to cluster 2
    
    # K4 without tetrahedron (Type 17 candidate)
    # All 4 triangles filled but no 3-simplex
    simplices.append([45, 46, 47])
    simplices.append([45, 46, 48])
    simplices.append([45, 47, 48])
    simplices.append([46, 47, 48])
    simplices.append([45, 48])  # ensure edge
    
    return save_dataset("dataset_b", simplices)


def dataset_c():
    """
    Dataset C: 200 vertices, dense clusters connected by sparse bridges
    Stress tests V2's ability to traverse heterogeneous structure
    """
    random.seed(42)
    simplices = []
    
    # 5 dense clusters of 20 vertices each (0-19, 20-39, 40-59, 60-79, 80-99)
    for cluster in range(5):
        base = cluster * 20
        vertices = list(range(base, base + 20))
        
        # Random edges within cluster (density ~30%)
        for i in range(len(vertices)):
            for j in range(i+1, len(vertices)):
                if random.random() < 0.3:
                    simplices.append([vertices[i], vertices[j]])
        
        # Some filled triangles
        for _ in range(5):
            triple = random.sample(vertices, 3)
            triple.sort()
            simplices.append(triple)
            # Also add edges
            simplices.append([triple[0], triple[1]])
            simplices.append([triple[1], triple[2]])
            simplices.append([triple[0], triple[2]])
        
        # 1-2 tetrahedra per cluster
        if cluster < 3:
            quad = random.sample(vertices[:10], 4)
            quad.sort()
            simplices.append(quad)
            # Add all sub-simplices
            for a in range(4):
                for b in range(a+1, 4):
                    simplices.append([quad[a], quad[b]])
                    for c in range(b+1, 4):
                        simplices.append([quad[a], quad[b], quad[c]])
    
    # Sparse bridge region (vertices 100-149): long paths
    for i in range(100, 149):
        simplices.append([i, i+1])
    
    # Connect clusters to bridge
    simplices.append([19, 100])
    simplices.append([39, 110])
    simplices.append([59, 120])
    simplices.append([79, 130])
    simplices.append([99, 140])
    
    # Star hub (vertex 150, connected to 151-170)
    for v in range(151, 171):
        simplices.append([150, v])
    simplices.append([150, 149])  # connect to bridge
    simplices.append([150, 0])    # connect to cluster 0
    
    # Additional sparse region (171-199)
    for i in range(171, 199):
        if random.random() < 0.4:
            simplices.append([i, i+1])
    # Connect to rest
    simplices.append([171, 170])
    simplices.append([199, 80])
    
    return save_dataset("dataset_c", simplices)


def dataset_d():
    """
    Dataset D: 80 vertices, scale-free-like with hub structure
    One big hub + medium hubs + many leaves
    Similar to co-authorship networks
    """
    random.seed(123)
    simplices = []
    
    # Main hub: vertex 0, connected to 1-20
    for v in range(1, 21):
        simplices.append([0, v])
    
    # Medium hubs: vertices 1, 5, 10, connected to some others
    for v in range(21, 31):
        simplices.append([1, v])
    for v in range(31, 39):
        simplices.append([5, v])
    for v in range(39, 46):
        simplices.append([10, v])
    
    # Small hubs
    for v in range(46, 51):
        simplices.append([21, v])
    for v in range(51, 55):
        simplices.append([31, v])
    
    # Leaf chains
    for i in range(55, 65):
        simplices.append([i, i+1])
    simplices.append([55, 46])
    
    # Some triangles among hubs
    simplices.append([0, 1, 5])
    simplices.append([0, 5, 10])
    simplices.append([0, 1, 10])
    simplices.append([1, 5, 10])  # filled
    
    # Tetrahedron
    simplices.append([0, 1, 5, 10])
    
    # More triangles
    simplices.append([1, 21, 22])
    simplices.append([5, 31, 32])
    simplices.append([21, 46, 47])
    
    # Extra edges among leaves for 4-cycles
    simplices.append([2, 3])
    simplices.append([3, 4])
    simplices.append([4, 2])   # triangle 2-3-4
    simplices.append([6, 7])
    simplices.append([7, 8])
    
    # Isolated cluster (66-79) with own structure
    for i in range(66, 75):
        simplices.append([i, i+1])
    simplices.append([66, 75])
    simplices.append([75, 76])
    simplices.append([76, 77])
    simplices.append([77, 78])
    simplices.append([78, 79])
    simplices.append([79, 66])  # cycle
    
    # Filled triangles in isolated cluster
    simplices.append([66, 67, 68])
    simplices.append([70, 71, 72])
    
    # Connect isolated cluster
    simplices.append([66, 15])
    simplices.append([75, 39])
    
    return save_dataset("dataset_d", simplices)


# ============================================================
# MAIN
# ============================================================

if __name__ == "__main__":
    print("Generating test datasets...")
    print()
    
    os.makedirs(TESTDATA_DIR, exist_ok=True)
    
    # Also keep the original tiny testdata
    save_dataset("tiny", [
        [0,1], [1,2], [2,3], [3,4], [0,4],
        [0,1,2], [2,3,4], [0,3,4]
    ])
    print()
    
    path_a = dataset_a()
    path_b = dataset_b()
    path_c = dataset_c()
    path_d = dataset_d()
    
    print()
    print("All datasets created in testdata/")
    print()
    print("Usage:")
    print("  ./exact_sfd testdata/dataset_b/nverts.txt testdata/dataset_b/simplices.txt results/exact.res")
    print("  ./unified_sfd_v2 testdata/dataset_b/nverts.txt testdata/dataset_b/simplices.txt results/v2.res 0")
    print()
    print("Or use run_20times.py with DATASET environment variable:")
    print("  DATASET=dataset_b python3 run_20times.py")