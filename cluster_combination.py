import shutil
from typing import Set, List, Optional
from pathlib import Path

import numpy as np
from astropy.io import fits
from spectral_cube import SpectralCube
from tqdm import tqdm


def remove_edge_clusters(cluster: Path, output_dir: Optional[Path] = None,
                         pixels_from_edge: Optional[int] = 2,
                         edges_to_remove: Optional[List[str]] = ["left", "right"]) -> None:
    """Creates a (.fits)-file without clusters at the edges
    of the cubes.
    Removes the clusters that (in any channel) contact the left
    or right edges of the frame.

    Parameters
    ----------
    cluster : pathlib.Path
        The (.fits)-file to be without edge clusters.
    output_dir: pathlih.Path
        The directory to contain the "no_edge" (.fits)-files.
    pixels_to_edge : int, optional
        The amount of pixels from the edge(s) to be defined "in contact" with it.
    edges_to_remove : list of str, optional
        A list containing either "left" or "right" or both to indicate which edge
        should be checked for clusters that contact it.
    """
    if output_dir is None:
        output_dir = cluster.parent / "no_edge"
    if not output_dir.exists():
        output_dir.mkdir()

    output_file = output_dir / f"no_edge_{cluster.name}"
    shutil.copyfile(cluster, output_file)

    left_edge_clusters, right_edge_clusters = [], []
    with fits.open(cluster) as hdul:
        data = hdul[0].data
    if "left" in edges_to_remove:
        left_edge = pixels_from_edge-1
        left_edge_clusters = data[:, :, left_edge][data[:, :, left_edge] != -1]
    if "right" in edges_to_remove:
        right_edge = -pixels_from_edge
        right_edge_clusters = data[:, :, right_edge][data[:, :, right_edge] != -1]
    edge_clusters = np.unique(np.concatenate((right_edge_clusters, left_edge_clusters)))
    data[np.isin(data, edge_clusters)] = -1
    with fits.open(output_file, "update") as hdul_out:
        hdul_out[0].data = data


def run_edge_removal(run_dir: Path,
                     glob_tag: Optional[str] = "cluster_*.fits",
                     **kwargs) -> None:
    """Runs the edge removal for a run.

    Will create a folder "no_edge" that will contain the individual
    edge cluster free (.fits)-files and one combined file that contains the
    information of all clusters.

    Parameters
    ----------
    run_dir : pathlib.Path
        A directory containing the (.fits)-files corresponding to one run.
    glob_tag : str, optional
        The tag that is used to glob the files to be edge removed.
    output_dir: pathlih.Path
        The directory to contain the "no_edge" (.fits)-files.
    pixels_to_edge : int, optional
        The amount of pixels from the edge(s) to be defined "in contact" with it.
    """
    run_paths = sorted(list(run_dir.glob(glob_tag)), key=lambda x: x.stem)
    for index, cube in tqdm(enumerate(run_paths), "Removing Edges"):
        if index == 0:
            edges_to_remove = ["right"]
        elif index == len(run_paths)-1:
            edges_to_remove = ["left"]
        else:
            edges_to_remove = ["left", "right"]
        remove_edge_clusters(cube, edges_to_remove=edges_to_remove, **kwargs)


def get_clusters(cluster: np.ndarray):
    """Reads the clusters' ids and their indices from a 2D numpy array."""
    if cluster[cluster != -1].size == 0:
        return None, None
    cluster_ids = np.unique(cluster[cluster != -1])
    cluster_indices = [np.where(cluster == cluster_id)
                       for cluster_id in cluster_ids]
    return cluster_ids, cluster_indices


def check_intersect(indices: np.ndarray,
                    shifted_indices: np.ndarray) -> bool:
    """Checks if there is an intersect between the clusters' indices."""
    indices = set(zip(*indices))
    shifted_indices = set(zip(*shifted_indices))
    return indices.intersection(shifted_indices)


def write_non_overlap_clusters(combined_slice: np.ndarray,
                               indices: np.ndarray,
                               cluster_ids: np.ndarray,
                               cluster_id_set: Set) -> np.ndarray:
    """Writes the non-overlapping clusters in a (.fits)-file
    from a list of clusters.

    Parameters
    ----------
    combined_slice : numpy.ndarray
        A 2D array to contain the wanted information from both runs.
    indices : numpy.ndarray
        The indices of the specified run.
    cluster_ids : numpy.ndarray
        The cluster ids of the specified run.
    cluster_ids_set : set
        The cluster ids of the specified run that are to be written to the combined run.

    Returns
    -------
    all_indices : List[numpy.ndarray]
        A list of the indices that are written to the combined slice (for debugging).
    """
    if indices is None:
        return
    all_indices = []
    for cluster_id in cluster_id_set:
        index = np.where(cluster_ids == cluster_id)[0][0]
        if indices[index][0].size == 0:
            continue
        combined_slice[indices[index]] = cluster_id
        all_indices.append(indices[index])
    return all_indices


def combine_runs(cluster_file: Path,
                 cluster_shifted_path: Path, output_path: Path) -> None:
    """Creates a new combines (.fits)-file with the real clusters.

    Removes overlapping clusters retaining the larger ones and
    inserts all non-overlapping clusters.

    Parameters
    ----------
    cluster_file : pathlib.Path
        The file of the first run.
    cluster_shifted_file : pathlib.Path
        The file of the second/shifted run.
    output_file : pathlib.Path
        The path and name for the output file which will be created.
    """
    cluster = SpectralCube.read(cluster_file)
    cluster_shifted = SpectralCube.read(cluster_shifted_path)

    combined_array = []
    for slice_index in tqdm(range(combined.shape[0]), "Remove overlapping clusters"):
        appended, combined_slice = False, -np.ones(cluster[slice_index].shape)
        original_ids, original_indices = get_clusters(cluster[slice_index].value)
        shifted_ids, shifted_indices = get_clusters(cluster_shifted[slice_index].value)
        if original_indices is None or shifted_indices is None:
            if not appended:
                combined_array.append(combined_slice)
                appended = True
            continue

        original_ids_set, shifted_ids_set = set(original_ids), set(shifted_ids)
        original_cluster_indices, shifted_cluster_indices = None, None

        for original_index, original_id in enumerate(original_ids):
            original_cluster_indices = original_indices[original_index]
            if original_cluster_indices[0].size == 0:
                if not appended:
                    combined_array.append(combined_slice)
                    appended = True
                continue

            for shifted_index, shifted_id in enumerate(shifted_ids):
                shifted_cluster_indices = shifted_indices[shifted_index]
                if shifted_cluster_indices[0].size == 0:
                    if not appended:
                        combined_array.append(combined_slice)
                        appended = True
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

        non_overlap_original = write_non_overlap_clusters(
                combined_slice, original_indices, original_ids, original_ids_set)
        non_overlap_shifted = write_non_overlap_clusters(
                combined_slice, shifted_indices, shifted_ids, shifted_ids_set)
        if not appended:
            combined_array.append(combined_slice)
            appended = True
    combined = SpectralCube(data=np.array(combined_array), wcs=cluster.wcs)
    combined.write(output_path, overwrite=True)


if __name__ == "__main__":
    data_dir = Path("/data/beegfs/astro-storage/groups/beuther/syed/HISA_study/THOR_survey/hisa_extraction/scimes/cluster_extraction/")
    second_run = data_dir / "second_run"

    # run_edge_removal(data_dir)
    run_edge_removal(second_run)

    # first_cluster = data_dir / "all_clusters_run_01.fits"
    # second_cluster = data_dir / "all_clusters_run_02.fits"
    # combine_runs(data_dir / "test.fits", first_cluster, second_cluster)

