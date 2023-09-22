import shutil
from typing import List
from pathlib import Path

import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube
from tqdm import tqdm


def remove_edge_clusters(cluster: Path, output_dir: Path,
                         edges_to_remove: List[str] = ["left", "right"]) -> None:
    """Creates a (.fits)-file without clusters at the edges
    of the cubes.
    Removes the clusters that (in any channel) contact the left
    or right edges of the frame.
    """
    if not output_dir.exists():
        output_dir.mkdir()
    output_file = output_dir / f"no_edge_{cluster.name}"
    shutil.copyfile(cluster, output_file)

    left_edge_clusters, right_edge_clusters = [], []
    with fits.open(cluster) as hdul:
        data = hdul[0].data
    if "left" in edges_to_remove:
        left_edge_clusters = data[:, :, 0][data[:, :, 0] != -1]
    if "right" in edges_to_remove:
        right_edge_clusters = data[:, :, -1][data[:, :, -1] != -1]
    edge_clusters = np.unique(np.concatenate((right_edge_clusters, left_edge_clusters)))
    data[np.isin(data, edge_clusters)] = -1
    with fits.open(output_file, "update") as hdul_out:
        hdul_out[0].data = data


def get_clusters(cluster: np.ndarray):
    """Reads the clusters' information from a (.fits)-file."""
    if cluster[cluster != -1].size == 0:
        return None, None
    cluster_ids = np.unique(cluster[cluster != -1])
    cluster_indices = [np.where(cluster == cluster_id)
                       for cluster_id in cluster_ids]
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

    for slice_index in tqdm(range(combined.shape[0]), "Remove overlapping clusters"):
        combined_slice = combined[slice_index].value
        original_ids, original_indices = get_clusters(cluster[slice_index].value)
        shifted_ids, shifted_indices = get_clusters(cluster_shifted[slice_index].value)
        if original_indices is None or shifted_indices is None:
            continue

        original_ids_set, shifted_ids_set = set(original_ids), set(shifted_ids)
        original_cluster_indices, shifted_cluster_indices = None, None

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

                    if original_id in original_ids_set:
                        original_ids_set.remove(original_id)
                    if shifted_id in shifted_ids_set:
                        shifted_ids_set.remove(shifted_id)

        write_non_overlap_clusters(
                combined_slice, original_cluster_indices, original_ids_set)
        write_non_overlap_clusters(
                combined_slice, shifted_cluster_indices, shifted_ids_set)


if __name__ == "__main__":
    first_run = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/cube/sub_cubes")
    first_run_paths = sorted(list(first_run.glob("*.fits")), key=lambda x: x.stem)
    for index, cube in tqdm(enumerate(first_run_paths), "First run"):
        if index == 0:
            edges_to_remove = ["right"]
        elif index == len(first_run_paths)-1:
            edges_to_remove = ["left"]
        else:
            edges_to_remove = ["left", "right"]
        remove_edge_clusters(cube, first_run / "no_edge", edges_to_remove)

    second_run = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/cube/sub_cubes/second_run")
    second_run_paths = sorted(list(second_run.glob("*.fits")), key=lambda x: x.stem)
    for cube in tqdm(second_run.glob("*.fits"), "Second run"):
        if index == 0:
            edges_to_remove = ["right"]
        elif index == len(first_run_paths)-1:
            edges_to_remove = ["left"]
        else:
            edges_to_remove = ["left", "right"]
        remove_edge_clusters(cube, second_run / "no_edge", edges_to_remove)

