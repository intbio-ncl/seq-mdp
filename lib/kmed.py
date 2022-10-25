import json

import numpy as np
from sklearn_extra.cluster import KMedoids

METHOD: str = "single"

mat = np.load("PF00155_mat.npy")
result = KMedoids(method='pam', n_clusters=100, metric="precomputed", max_iter=300).fit(mat)

with open("PF00155_headings.json") as f:
    headings = json.load(f)
    print([headings[idx] for idx in result.medoid_indices_])
