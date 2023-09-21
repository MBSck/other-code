import shutil
from pathlib import Path
from typing import Set

import numpy as np
from tqdm import tqdm
from spectral_cube import SpectralCube


def check_adjacent(indices_set: Set, other_indices_set: Set) -> bool:
    """Checks if any of the pixels in the sets are adjacent."""
    adjacencies = []
    for count, index in enumerate(indices_set):
        adjacencies.append([])
        for other_index in other_indices_set:
            adjacencies[count].append(abs(index[0]-other_index[0]) < 1
                                      and abs(index[1]-other_index[1]) < 1)
    return np.any(adjacencies)


def set_adjacent(cluster_ids: np.ndarray,
                 cluster_indices: np.ndarray) -> None:
    """Checks if any indices are adjacent"""
    clusters_to_merge = {cluster_id: [] for cluster_id in cluster_ids}
    for cluster_id, indices in zip(cluster_ids, cluster_indices):
        indices_set = set(zip(*indices))
        for other_cluster_id, other_indices in zip(cluster_ids, cluster_indices):
            other_indices_set = set(zip(*other_indices))
            if check_adjacent(indices_set, other_indices_set):
                clusters_to_merge[cluster_id].append(other_cluster_id)
    new_cluster_ids, new_cluster_indices = [], []
    all_cluster_ids = cluster_ids.copy().tolist()
    for cluster_id, other_cluster_ids in clusters_to_merge.items():
        if cluster_id in all_cluster_ids:
            tmp_indices = cluster_indices[np.where(cluster_ids == cluster_id)]
            new_cluster_ids.append(cluster_id)
            all_cluster_ids.remove(cluster_id)
            for other_cluster_id in other_cluster_ids:
                if other_cluster_id in all_cluster_ids:
                    all_cluster_ids.remove(other_cluster_id)
        # TODO: Include merging of lists here.
        else:
            continue

    return cluster_ids, cluster_indices


def get_clusters(cluster: np.ndarray):
    """Reads the clusters' information from a (.fits)-file."""
    if cluster[cluster != -1].size == 0:
        return None, None
    cluster_ids = np.unique(cluster[cluster != -1])
    cluster_indices = [np.where(cluster == cluster_id)
                       for cluster_id in cluster_ids]
    cluster_ids, cluster_indices = set_adjacent(cluster_ids, cluster_indices)
    return cluster_ids, cluster_indices


def check_intersect(indices: np.ndarray,
                    shifted_indices: np.ndarray) -> bool:
    """Checks if there is an intersect between the clusters."""
    indices = set(zip(*indices))
    shifted_indices = set(zip(*shifted_indices))
    return indices.intersection(shifted_indices)


def write_non_overlap_clusters(combined: np.ndarray,
                               cluster_indices: np.ndarray,
                               cluster_ids: np.ndarray) -> None:
    """Writes the non-overlapping clusters in a (.fits)-file
    from a list of clusters.
    """
    if cluster_indices is None:
        return
    if cluster_indices[0].size == 0:
        return
    for indices, cluster_id in zip(cluster_indices, cluster_ids):
        combined[indices] = cluster_id


def combine_runs(combined_file: np.ndarray,
                 cluster_file: Path, cluster_shifted_path: Path) -> None:
    """Creates a new combines (.fits)-file with the real clusters.

    Removes overlapping clusters retaining the larger ones and
    inserts all non-overlapping clusters.
    """
    cluster = SpectralCube.read(cluster_file)
    cluster_shifted = SpectralCube.read(cluster_shifted_path)
    combined = SpectralCube.read(combined_file)

    for slice_index in tqdm(range(combined.shape[0]), "Remove overlap Clusters"):
        combined_slice = combined[slice_index].value
        original_ids, original_indices = get_clusters(cluster[slice_index].value)
        shifted_ids, shifted_indices = get_clusters(cluster_shifted[slice_index].value)
        if original_indices is None or shifted_indices is None:
            continue
        original_ids_copy = original_ids.copy().tolist()
        shifted_ids_copy = shifted_ids.copy().tolist()

        original_cluster_indices = None
        shifted_cluster_indices = None

        for original_index, original_id in enumerate(original_ids):
            original_cluster_indices = original_indices[original_index]
            if original_cluster_indices[0].size == 0:
                continue

            for shifted_index, shifted_id in enumerate(shifted_ids):
                shifted_cluster_indices = shifted_indices[shifted_index]
                if shifted_cluster_indices[0].size == 0:
                    continue
                if check_intersect(original_cluster_indices, shifted_cluster_indices):
                    if original_cluster_indices[0].size > shifted_cluster_indices[0].size:
                        combined_slice[original_cluster_indices] = original_id
                    else:
                        combined_slice[shifted_cluster_indices] = shifted_id

                    if original_id in original_ids_copy:
                        original_ids_copy.remove(original_id)
                    if shifted_id in shifted_ids_copy:
                        shifted_ids_copy.remove(shifted_id)

        write_non_overlap_clusters(
                combined_slice, original_cluster_indices, original_ids_copy)
        write_non_overlap_clusters(
                combined_slice, shifted_cluster_indices, shifted_ids_copy)


if __name__ == "__main__":
    data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/cube")
    cluster_path = data_dir / "clusters_run_01.fits"
    shifted_cluster_path = data_dir / "clusters_run_02.fits"
    combined_clusters_path = data_dir / "final_clusters_cube.fits"
    combined_clusters_path_copy = data_dir / "final_clusters_cube_copy.fits"
    # shutil.copyfile(combined_clusters_path, combined_clusters_path_copy)
    combine_runs(combined_clusters_path_copy, cluster_path, shifted_cluster_path)
