"""Municipal adjacency graph and the BYM2 scaling factor.

Queen contiguity on the IBGE municipal mesh, restricted to the 5,565 municipalities of the
analysis table (data/analysis_muni.csv), then welded into a single connected component by
joining each isolated piece to its nearest centroid. An ICAR prior is improper on every
disconnected component, so this is required rather than cosmetic.

The BYM2 parameterisation (Riebler et al. 2016) needs the ICAR component scaled so that
its generalised variance is 1; the scaling factor is the geometric mean of the diagonal
of the generalised inverse of the graph Laplacian Q = D - W, computed here by eigen
decomposition of the dense 5,565 x 5,565 matrix.

Outputs (v2/data/):
    adjacency_pairs.csv      undirected neighbour pairs (i < j), 0-indexed in the
                             order of municipalities.csv
    municipalities.csv       CODMUNRES, longitude, latitude, in model order
    bym2_scaling.json        scaling factor and graph summary
    adjacency_report.txt

Run (after 03_build_dataset.py, whose municipality order this file follows):
    python v2/src/02_adjacency.py
"""

import json
import time
from pathlib import Path

import geopandas as gpd
import numpy as np
import pandas as pd
from libpysal.weights import Queen
from shapely.geometry import Point

HERE = Path(__file__).resolve().parent
V2 = HERE.parent
A1 = V2.parent

# Nazária (PI), emancipated from Teresina in 2008, is absent from the mesh; coordinates are
# the median of its eight CNES establishments.
FALLBACK_POINTS = {"220672": (-5.3578, -42.8193)}


def component_labels(n, neighbours):
    labels = np.full(n, -1)
    current = 0
    for start in range(n):
        if labels[start] != -1:
            continue
        stack = [start]
        labels[start] = current
        while stack:
            node = stack.pop()
            for nb in neighbours[node]:
                if labels[nb] == -1:
                    labels[nb] = current
                    stack.append(nb)
        current += 1
    return labels


def main():
    t0 = time.time()
    # Model order is the sorted CODMUNRES of the analysis table written by 03_build_dataset.py.
    modelled = pd.read_csv(V2 / "data" / "analysis_muni.csv", dtype={"CODMUNRES": str}).CODMUNRES.values
    if not (np.sort(modelled) == modelled).all():
        raise SystemExit("analysis table is not sorted by CODMUNRES")

    munis = gpd.read_file(A1 / "data" / "ref" / "br_municipios.geojson")
    munis["cod6"] = munis["id"].astype(str).str[:6]
    munis = munis[munis.cod6.isin(modelled)].copy()
    invalid = int((~munis.geometry.is_valid).sum())
    if invalid:
        munis["geometry"] = munis.geometry.make_valid()
    munis = munis.set_index("cod6").reindex(modelled).reset_index()
    # The published mesh predates a few municipal emancipations. A municipality without a
    # polygon gets a point from the median coordinates of its CNES establishments and is
    # attached to its three nearest neighbours by centroid distance.
    no_polygon = munis.loc[munis.geometry.isna(), "cod6"].tolist()
    for code in no_polygon:
        if code not in FALLBACK_POINTS:
            raise SystemExit(f"municipality {code} has no polygon and no fallback point")
        lat, lon = FALLBACK_POINTS[code]
        munis.loc[munis.cod6 == code, "geometry"] = Point(lon, lat).buffer(1e-4)
    munis = gpd.GeoDataFrame(munis, geometry="geometry", crs="EPSG:4326")

    w = Queen.from_dataframe(munis, use_index=False)
    neighbours = {i: set(v) for i, v in w.neighbors.items()}
    islands = [i for i, v in neighbours.items() if not v]

    proj = munis.to_crs(5880)
    cent = proj.geometry.centroid
    xy = np.column_stack([cent.x.values, cent.y.values])
    joins = []
    for code in no_polygon:
        i = int(np.where(munis.cod6.values == code)[0][0])
        d = np.hypot(xy[:, 0] - xy[i, 0], xy[:, 1] - xy[i, 1])
        d[i] = np.inf
        for j in np.argsort(d)[:3]:
            neighbours[i].add(int(j)); neighbours[int(j)].add(i)
            joins.append((i, int(j), float(d[j]) / 1000.0))
    while True:
        labels = component_labels(len(xy), neighbours)
        if labels.max() == 0:
            break
        # Join the smallest component to its nearest municipality in any other component.
        sizes = np.bincount(labels)
        small = int(np.argmin(sizes))
        target = labels == small
        other = ~target
        best = None
        for i in np.where(target)[0]:
            d = np.hypot(xy[other, 0] - xy[i, 0], xy[other, 1] - xy[i, 1])
            k = int(np.argmin(d))
            j = int(np.where(other)[0][k])
            if best is None or d[k] < best[0]:
                best = (float(d[k]), i, j)
        _, i, j = best
        neighbours[i].add(j)
        neighbours[j].add(i)
        joins.append((i, j, best[0] / 1000.0))

    pairs = sorted({(min(i, j), max(i, j)) for i, v in neighbours.items() for j in v})
    n = len(munis)
    pd.DataFrame(pairs, columns=["i", "j"]).to_csv(V2 / "data" / "adjacency_pairs.csv", index=False)
    wgs = cent.to_crs(4326)
    pd.DataFrame({"CODMUNRES": munis.cod6.values, "lon": wgs.x.values, "lat": wgs.y.values}) \
        .to_csv(V2 / "data" / "municipalities.csv", index=False)

    # BYM2 scaling factor: geometric mean of diag(Q^+), Q = D - W (Riebler et al. 2016).
    W = np.zeros((n, n))
    for i, j in pairs:
        W[i, j] = W[j, i] = 1.0
    Q = np.diag(W.sum(1)) - W
    vals, vecs = np.linalg.eigh(Q)
    keep = vals > 1e-8 * vals.max()
    diag_pinv = np.einsum("ij,j,ij->i", vecs[:, keep], 1.0 / vals[keep], vecs[:, keep])
    scaling = float(np.exp(np.mean(np.log(diag_pinv))))
    degrees = W.sum(1)
    summary = {
        "n": int(n), "pairs": len(pairs), "mean_degree": float(degrees.mean()),
        "min_degree": int(degrees.min()), "max_degree": int(degrees.max()),
        "queen_islands": len(islands), "components_welded": len(joins),
        "zero_eigenvalues": int((~keep).sum()), "scaling_factor": scaling,
        "joins_km": [round(k, 1) for _, _, k in joins],
    }
    (V2 / "data" / "bym2_scaling.json").write_text(json.dumps(summary, indent=2))
    report = "\n".join([
        "Municipal adjacency (queen contiguity, welded)", "=" * 50,
        f"municipalities       : {n}",
        f"invalid geometries   : {invalid} (repaired)",
        f"undirected pairs     : {len(pairs):,}",
        f"mean degree          : {degrees.mean():.2f} (min {int(degrees.min())}, max {int(degrees.max())})",
        f"without polygon      : {len(no_polygon)} {no_polygon} (point from CNES, 3 nearest neighbours)",
        f"queen islands        : {len(islands)}",
        f"joins added          : {len(joins)}",
        *[f"    joined {munis.cod6.iloc[i]} <-> {munis.cod6.iloc[j]}  ({km:.0f} km)" for i, j, km in joins],
        f"zero eigenvalues of Q: {int((~keep).sum())} (1 = connected)",
        f"BYM2 scaling factor  : {scaling:.4f}",
        f"elapsed              : {time.time() - t0:.0f}s",
    ])
    (V2 / "data" / "adjacency_report.txt").write_text(report + "\n")
    print(report)


if __name__ == "__main__":
    main()
