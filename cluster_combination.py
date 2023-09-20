import shutil
from pathlib import Path
from typing import Dict

import numpy as np
from astropy.wcs import WCS
from tqdm import tqdm
from spectral_cube import SpectralCube


def get_clusters(cluster: Path) -> Dict:
    """Reads the clusters' information from a (.fits)-file."""
    cluster_ids = np.unique(cluster[cluster != -1])
    cluster_indices = [np.where(cluster == cluster_id)
                       for cluser_id in tqdm(cluster_ids, "Check")]
    return cluster_ids, cluster_indices


# def check_intersect(coordinate_range: CoordinateRange,
                    # shifted_coordinate: np.ndarray) -> bool:
    # """Checks if there is an intersect between the clusters."""
    # latitude, longitude = shifted_coordinate[0], shifted_coordinate[1]
    # in_latitude = np.logical_and(coordinate_range.lat_min <= latitude, 
                                 # coordinate_range.lat_max >= latitude)
    # in_longitude = np.logical_and(coordinate_range.lon_min <= longitude,
                                  # coordinate_range.lon_max >= longitude)
    # return np.any(np.logical_and(in_latitude, in_longitude))


def get_indices(coordinate: np.ndarray, w: WCS) -> np.ndarray:
    """Gets indices from world coordinates."""
    indices = w.all_world2pix(*coordinate, 0, 0)[:2]
    return [np.array(list(map(int, indices[0]))),
            np.array(list(map(int, indices[1])))]


def add_non_overlapping_clusters(
        output: np.ndarray, w: WCS,
        cluster_coordinates: np.ndarray, cluster_ids: np.ndarray) -> None:
    """Sets the clusters in a (.fits)-file from a list of clusters."""
    for cluster_index, cluster_id in enumerate(cluster_ids):
        coordinate = cluster_coordinates[cluster_index]
        if coordinate[0].size == 0:
            continue
        indices = get_indices(coordinate, w)
        output[indices] = cluster_id


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
        original_cluster_ids, original_indices = get_clusters(
            cluster[slice_index].value)
        shifted_cluster_ids, shifted_indices = get_clusters(
            cluster_shifted[slice_index].value)
        breakpoint()
        original_cluster_ids_copy = original_cluster_ids.copy().tolist()
        shifted_cluster_ids_copy = shifted_cluster_ids.copy().tolist()
        for original_cluster_index, original_cluster_id in enumerate(original_cluster_ids):
            original_coordinate = original_slice_coordinates[original_cluster_index]
            if original_coordinate[0].size == 0:
                continue
            # coordinate_range = CoordinateRange(original_coordinate)
            # for shifted_cluster_index, shifted_cluster_id in enumerate(shifted_cluster_ids):
                # shifted_coordinate = shifted_slice_coordinates[shifted_cluster_index]
                # if shifted_coordinate[0].size == 0:
                    # continue
                # if check_intersect(coordinate_range, shifted_coordinate):
                    # indices_first = get_indices(original_coordinate, w)
                    # indices_second = get_indices(shifted_coordinate, w)
                    # if original_coordinate[0].size < shifted_coordinate[0].size:
                        # combined[slice_index][indices_first] = -1
                        # combined[slice_index][indices_second] = shifted_cluster_id
                    # else:
                        # combined[slice_index][indices_first] = original_cluster_id
                        # combined[slice_index][indices_second] = -1
                    # if original_cluster_id in original_cluster_ids_copy:
                        # original_cluster_ids_copy.remove(original_cluster_id)
                    # if shifted_cluster_id in shifted_cluster_ids_copy:
                        # shifted_cluster_ids_copy.remove(shifted_cluster_id)

        add_non_overlapping_clusters(combined, w, slice_index,
                                     original_slice_coordinates, original_cluster_ids_copy)
        add_non_overlapping_clusters(combined, w, slice_index,
                                     shifted_slice_coordinates, shifted_cluster_ids_copy)


if __name__ == "__main__":
    data_dir = Path("/data/beegfs/astro-storage/groups/matisse/scheuck/data/cube")
    cluster_path = data_dir / "clusters_run_01.fits"
    shifted_cluster_path = data_dir / "clusters_run_02.fits"
    combined_clusters_path = data_dir / "final_clusters_cube.fits"
    combined_clusters_path_copy = data_dir / "final_clusters_cube_copy.fits"
    shutil.copyfile(combined_clusters_path, combined_clusters_path_copy)
    combine_runs(combined_clusters_path_copy, cluster_path, shifted_cluster_path)
