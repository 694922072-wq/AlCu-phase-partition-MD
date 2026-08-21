from __future__ import annotations

import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from ovito.io import import_file
from ovito.modifiers import AtomicStrainModifier


ROOT = Path(r"D:\leng\AlCu_MAIN_PPP_Zhou")
REPORT = ROOT / "reports" / "interface_resolved_4models"
MECH = REPORT / "01_mechanical_response"
REGION = REPORT / "02_region_definition_refined"
STATS = REPORT / "03_interface_D2min_deltaPE_shear"
OUT = REPORT / "04_void_cluster"
DEBUG = REPORT / "10_debug_batch2"
OUT.mkdir(parents=True, exist_ok=True)
DEBUG.mkdir(parents=True, exist_ok=True)

DUMPS = {
    "M3_SYM": ROOT / "dumps" / "M3_SYM_tensile_x25_xy_z_clean_light.lammpstrj",
    "M4_SYM_RATIO": ROOT / "dumps" / "M4_SYM_RATIO_tensile_x25_xy_z_clean_light.lammpstrj",
    "M4_SYM": ROOT / "dumps" / "M4_SYM_tensile_x25_xy_z_clean_light.lammpstrj",
    "M4_LIT": ROOT / "dumps" / "M4_LIT_tensile_x25_clean_light.lammpstrj",
}
MODEL_ORDER = ["M3_SYM", "M4_SYM_RATIO", "M4_SYM", "M4_LIT"]
CUTOFFS = [3.5, 4.0]
WIDTHS = [10, 5]
SCOPES = ["model_specific", "global"]
CLUSTER_TYPES = ["damage", "strict_damage"]


class UnionFind:
    def __init__(self, n: int) -> None:
        self.parent = np.arange(n, dtype=np.int64)
        self.rank = np.zeros(n, dtype=np.int8)

    def find(self, x: int) -> int:
        p = int(self.parent[x])
        if p != x:
            self.parent[x] = self.find(p)
        return int(self.parent[x])

    def union(self, a: int, b: int) -> None:
        ra, rb = self.find(a), self.find(b)
        if ra == rb:
            return
        if self.rank[ra] < self.rank[rb]:
            self.parent[ra] = rb
        elif self.rank[ra] > self.rank[rb]:
            self.parent[rb] = ra
        else:
            self.parent[rb] = ra
            self.rank[ra] += 1


def particle_array(data, name: str) -> np.ndarray:
    return np.asarray(data.particles[name])


def cell_info(data) -> tuple[np.ndarray, np.ndarray]:
    mat = np.asarray(data.cell.matrix, dtype=float)
    lengths = np.array([mat[0, 0], mat[1, 1], mat[2, 2]], dtype=float)
    origin = np.array([mat[0, 3], mat[1, 3], mat[2, 3]], dtype=float)
    return lengths, origin


def pbc_cluster_indices(pos: np.ndarray, lengths: np.ndarray, origin: np.ndarray, cutoff: float) -> list[np.ndarray]:
    n = len(pos)
    if n == 0:
        return []
    if n == 1:
        return [np.array([0], dtype=np.int64)]
    rel = (pos - origin) % lengths
    ncells = np.maximum(np.floor(lengths / cutoff).astype(int), 1)
    cell_size = lengths / ncells
    cidx = np.floor(rel / cell_size).astype(int)
    cidx = np.minimum(cidx, ncells - 1)
    cells: dict[tuple[int, int, int], list[int]] = defaultdict(list)
    for i, c in enumerate(cidx):
        cells[(int(c[0]), int(c[1]), int(c[2]))].append(i)
    uf = UnionFind(n)
    c2 = cutoff * cutoff
    neighbor_offsets = [(dx, dy, dz) for dx in (-1, 0, 1) for dy in (-1, 0, 1) for dz in (-1, 0, 1)]
    for key, ids in cells.items():
        for off in neighbor_offsets:
            nb = ((key[0] + off[0]) % ncells[0], (key[1] + off[1]) % ncells[1], (key[2] + off[2]) % ncells[2])
            if nb not in cells:
                continue
            for i in ids:
                for j in cells[nb]:
                    if j <= i:
                        continue
                    d = rel[j] - rel[i]
                    d -= lengths * np.round(d / lengths)
                    if float(np.dot(d, d)) <= c2:
                        uf.union(i, j)
    groups: dict[int, list[int]] = defaultdict(list)
    for i in range(n):
        groups[uf.find(i)].append(i)
    return [np.array(v, dtype=np.int64) for v in groups.values()]


def load_map(width: int) -> pd.DataFrame:
    gz = REGION / f"atom_region_id_map_refined_4models_{width}A.csv.gz"
    df = pd.read_csv(gz)
    df["id"] = df["id"].astype(np.int64)
    df["region_id"] = df["region_id"].astype(int)
    return df


def region_priority(name: str, region_type: str) -> int:
    if any(s in name for s in ["near_IMC", "near_Al2Cu", "Al_center_core", "Al_matrix_core"]):
        return 1
    if "interface" in name or "interface" in region_type:
        return 2
    if "interior_excl_interfaces" in name:
        return 3
    return 4


def primary_region_map(map_df: pd.DataFrame, bounds: pd.DataFrame) -> pd.DataFrame:
    merged = map_df.merge(bounds[["model", "region_id", "region_type", "z_min", "z_max"]], on=["model", "region_id"], how="left")
    merged["priority"] = [region_priority(n, t) for n, t in zip(merged["region_name"], merged["region_type"])]
    merged = merged.sort_values(["model", "id", "priority", "region_id"])
    return merged.drop_duplicates(["model", "id"], keep="first")[["model", "id", "region_id", "region_name", "region_type"]].copy()


def nearest_interface(bounds_model: pd.DataFrame, z_initial: float) -> tuple[str, float, str, float]:
    interfaces = bounds_model[bounds_model["region_type"].astype(str).str.contains("interface", case=False, na=False)].copy()
    if interfaces.empty:
        return "unavailable", np.nan, "unavailable", np.nan
    interfaces["z_center"] = (interfaces["z_min"] + interfaces["z_max"]) / 2.0
    interfaces["distance"] = (interfaces["z_center"] - z_initial).abs()
    row = interfaces.sort_values("distance").iloc[0]
    bounds = bounds_model.copy()
    boundary_rows = []
    for r in bounds.itertuples(index=False):
        boundary_rows.append((f"{r.region_name}:z_min", abs(float(r.z_min) - z_initial)))
        boundary_rows.append((f"{r.region_name}:z_max", abs(float(r.z_max) - z_initial)))
    bname, bdist = min(boundary_rows, key=lambda x: x[1])
    return str(row["region_name"]), float(row["distance"]), bname, float(bdist)


def cluster_records(
    model: str,
    event_row: pd.Series,
    scope: str,
    cluster_type: str,
    cutoff: float,
    ids_sel: np.ndarray,
    pos_sel: np.ndarray,
    z0_by_id: pd.Series,
    primary: pd.DataFrame,
    bounds_model: pd.DataFrame,
    total_atoms: int,
    lengths: np.ndarray,
    origin: np.ndarray,
) -> tuple[list[dict], list[dict], list[dict]]:
    if len(ids_sel) == 0:
        return [], [], []
    clusters = pbc_cluster_indices(pos_sel, lengths, origin, cutoff)
    preg = primary.set_index("id")
    records = []
    dist_rows = []
    atom_id_rows = []
    for cid, local_idx in enumerate(sorted(clusters, key=lambda a: len(a), reverse=True), start=1):
        ids = ids_sel[local_idx]
        pos = pos_sel[local_idx]
        size = int(len(ids))
        center = pos.mean(axis=0)
        spans = pos.max(axis=0) - pos.min(axis=0)
        rg = float(np.sqrt(((pos - center) ** 2).sum(axis=1).mean())) if size else np.nan
        z0_vals = z0_by_id.reindex(ids).to_numpy(dtype=float)
        z0_mean = float(np.nanmean(z0_vals)) if len(z0_vals) else np.nan
        reg_names = preg.reindex(ids)["region_name"].fillna("unmapped").to_list()
        counts = Counter(reg_names)
        total_members = sum(counts.values())
        dom_name, dom_count = counts.most_common(1)[0]
        second_name, second_count = counts.most_common(2)[1] if len(counts) > 1 else ("none", 0)
        nearest_if, nearest_if_dist, nearest_boundary, nearest_boundary_dist = nearest_interface(bounds_model, z0_mean)
        dist_json = json.dumps({k: int(v) for k, v in counts.items()}, ensure_ascii=False)
        records.append(
            {
                "model": model,
                "event": event_row["event_name"],
                "strain_percent": event_row["actual_strain_percent"],
                "frame_index": event_row["frame_index"],
                "threshold_scope": scope,
                "cluster_type": cluster_type,
                "cutoff_A": cutoff,
                "cluster_id": cid,
                "cluster_size_atoms": size,
                "cluster_fraction_total_atoms": size / total_atoms,
                "cluster_center_x": float(center[0]),
                "cluster_center_y": float(center[1]),
                "cluster_center_z_current": float(center[2]),
                "cluster_center_z_initial_mean": z0_mean,
                "cluster_span_x": float(spans[0]),
                "cluster_span_y": float(spans[1]),
                "cluster_span_z": float(spans[2]),
                "cluster_radius_gyration": rg,
                "dominant_region_name": dom_name,
                "dominant_region_fraction": dom_count / total_members if total_members else np.nan,
                "second_region_name": second_name,
                "second_region_fraction": second_count / total_members if total_members else np.nan,
                "nearest_initial_interface": nearest_if,
                "distance_to_nearest_initial_interface_A": nearest_if_dist,
                "nearest_region_boundary": nearest_boundary,
                "distance_to_nearest_region_boundary_A": nearest_boundary_dist,
                "region_distribution_json": dist_json,
                "interpretation_note": "candidate damage cluster, not confirmed crack-initiation plane",
            }
        )
        for name, count in counts.items():
            dist_rows.append(
                {
                    "model": model,
                    "event": event_row["event_name"],
                    "threshold_scope": scope,
                    "cluster_type": cluster_type,
                    "cutoff_A": cutoff,
                    "cluster_id": cid,
                    "region_name": name,
                    "atom_count_in_cluster_region": int(count),
                    "fraction_in_cluster": count / total_members if total_members else np.nan,
                }
            )
        if cid <= 3 and cutoff == 3.5:
            atom_id_rows.append(
                {
                    "model": model,
                    "event": event_row["event_name"],
                    "threshold_scope": scope,
                    "cluster_type": cluster_type,
                    "cutoff_A": cutoff,
                    "cluster_id": cid,
                    "atom_ids": " ".join(map(str, ids.tolist())),
                }
            )
    return records, dist_rows, atom_id_rows


def frame_metrics(pipe, frame: int, pe0_by_id: pd.Series | None = None) -> tuple[pd.DataFrame, np.ndarray, np.ndarray]:
    data = pipe.compute(frame)
    ids = particle_array(data, "Particle Identifier").astype(np.int64)
    pos = particle_array(data, "Position").astype(float)
    d2 = particle_array(data, "Nonaffine Squared Displacement").astype(float)
    pe = particle_array(data, "c_peatom").astype(float)
    if "Shear Strain" in list(data.particles.keys()):
        shear = particle_array(data, "Shear Strain").astype(float)
    else:
        shear = np.full_like(d2, np.nan)
    df = pd.DataFrame({"id": ids, "D2min": d2, "PE": pe, "shear": shear})
    if pe0_by_id is None:
        df["PE0"] = df["PE"]
    else:
        df["PE0"] = df["id"].map(pe0_by_id).to_numpy(dtype=float)
    df["delta_PE"] = df["PE"] - df["PE0"]
    lengths, origin = cell_info(data)
    return df, pos, np.vstack([lengths, origin])


def main() -> int:
    selected = pd.read_csv(MECH / "selected_keyframes_4models_batch11.csv")
    thresholds_model = pd.read_csv(STATS / "thresholds_model_specific_batch11.csv")
    thresholds_global = pd.read_csv(STATS / "thresholds_global_batch11.csv")
    bounds = {w: pd.read_csv(REGION / f"region_bounds_refined_4models_{w}A.csv") for w in WIDTHS}
    maps = {w: load_map(w) for w in WIDTHS}
    primary = {w: primary_region_map(maps[w], bounds[w]) for w in WIDTHS}

    all_records = {10: [], 5: []}
    all_dist = {10: [], 5: []}
    sensitivity = []
    atom_id_rows = []

    for model in MODEL_ORDER:
        print(f"[cluster] loading {model}")
        pipe = import_file(str(DUMPS[model]), multiple_frames=True)
        pipe.modifiers.append(AtomicStrainModifier(cutoff=4.0, reference_frame=0, output_nonaffine_squared_displacements=True))
        frame0_df, pos0, box0 = frame_metrics(pipe, 0)
        pe0_by_id = frame0_df.set_index("id")["PE"]
        z0_by_id = pd.Series(pos0[:, 2], index=frame0_df["id"].to_numpy())
        model_frames = selected[selected["model"] == model].sort_values("event_order")
        for event_row in model_frames.to_dict("records"):
            frame = int(event_row["frame_index"])
            print(f"[cluster] {model} {event_row['event_name']} frame={frame}")
            fdf, pos, box = frame_metrics(pipe, frame, pe0_by_id)
            lengths, origin = box[0], box[1]
            total_atoms = len(fdf)
            for width in WIDTHS:
                bounds_model = bounds[width][bounds[width]["model"] == model].copy()
                primary_model = primary[width][primary[width]["model"] == model].copy()
                for scope in SCOPES:
                    if scope == "model_specific":
                        th = thresholds_model[(thresholds_model["model"] == model) & (thresholds_model["band_width_A"] == width)].iloc[0]
                    else:
                        th = thresholds_global[thresholds_global["band_width_A"] == width].iloc[0]
                    high_d2 = fdf["D2min"].to_numpy() >= float(th["D2min_threshold_p95"])
                    high_delta = fdf["delta_PE"].to_numpy() >= float(th["deltaPE_threshold_p95"])
                    high_shear = fdf["shear"].to_numpy() >= float(th["shear_threshold_p95"])
                    masks = {
                        "damage": high_d2 & high_delta,
                        "strict_damage": high_d2 & high_delta & high_shear,
                    }
                    ids_all = fdf["id"].to_numpy()
                    for cluster_type, mask in masks.items():
                        ids_sel = ids_all[mask]
                        pos_sel = pos[mask]
                        for cutoff in CUTOFFS:
                            recs, dists, atom_rows = cluster_records(
                                model,
                                pd.Series(event_row),
                                scope,
                                cluster_type,
                                cutoff,
                                ids_sel,
                                pos_sel,
                                z0_by_id,
                                primary_model,
                                bounds_model,
                                total_atoms,
                                lengths,
                                origin,
                            )
                            all_records[width].extend(recs)
                            all_dist[width].extend(dists)
                            if width == 10:
                                atom_id_rows.extend(atom_rows)
                            if recs:
                                largest = max(recs, key=lambda r: r["cluster_size_atoms"])
                                nclusters = len(recs)
                                selected_count = int(mask.sum())
                                largest_size = largest["cluster_size_atoms"]
                            else:
                                nclusters = 0
                                selected_count = int(mask.sum())
                                largest_size = 0
                            sensitivity.append(
                                {
                                    "model": model,
                                    "event": event_row["event_name"],
                                    "strain_percent": event_row["actual_strain_percent"],
                                    "band_width_A": width,
                                    "threshold_scope": scope,
                                    "cluster_type": cluster_type,
                                    "cutoff_A": cutoff,
                                    "selected_atom_count": selected_count,
                                    "cluster_count": nclusters,
                                    "largest_cluster_size_atoms": largest_size,
                                    "largest_cluster_fraction_total_atoms": largest_size / total_atoms,
                                }
                            )

    for width in WIDTHS:
        df = pd.DataFrame(all_records[width])
        dist = pd.DataFrame(all_dist[width])
        if df.empty:
            df = pd.DataFrame()
            largest = pd.DataFrame()
        else:
            largest = df.sort_values("cluster_size_atoms", ascending=False).groupby(["model", "event", "threshold_scope", "cluster_type", "cutoff_A"], as_index=False).head(1)
            largest = largest.sort_values(["model", "threshold_scope", "cluster_type", "cutoff_A", "strain_percent"])
        df.to_csv(OUT / f"damage_clusters_all_4models_{width}A.csv", index=False, encoding="utf-8-sig")
        largest.to_csv(OUT / f"largest_damage_cluster_by_event_{width}A.csv", index=False, encoding="utf-8-sig")
        dist.to_csv(OUT / f"damage_cluster_region_distribution_{width}A.csv", index=False, encoding="utf-8-sig")
    pd.DataFrame(sensitivity).to_csv(OUT / "damage_cluster_threshold_sensitivity.csv", index=False, encoding="utf-8-sig")
    if atom_id_rows:
        with gzip.open(OUT / "largest_cluster_atom_ids_top3.csv.gz", "wt", encoding="utf-8") as f:
            pd.DataFrame(atom_id_rows).to_csv(f, index=False)
    (DEBUG / "damage_cluster_method_notes.md").write_text(
        "Damage clusters were computed from selected high-D2min/high-deltaPE atoms using a pure-Python periodic cell-list connected-component algorithm. "
        "Strict damage additionally requires high shear. Cluster outputs are candidate damage clusters, not cracks.\n",
        encoding="utf-8",
    )
    print("Batch 2 damage cluster analysis complete.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
