from libcpp.vector cimport vector
from libc.stdint cimport uint64_t, uint32_t
import os
import numpy as np
cimport numpy as np

def zwatershed(np.ndarray[np.float32_t, ndim=4] affs, thresholds, np.ndarray[uint32_t, ndim=3] gt = None):
    '''
    Compute segmentations from an affinity graph for several thresholds.

    Passed volumes need to be converted into contiguous memory arrays. This will
    be done for you if needed, but you can save memory by making sure your
    volumes are already C_CONTIGUOUS.

    Parameters
    ----------

        affs: numpy array, float32, 4 dimensional

            The affinities as an array with affs[channel][z][y][x].

        thresholds: list of int

            The thresholds to compute segmentations for. For each threshold, one
            segmentation is returned.

        gt  : numpy array, uint32, 3 dimensional (optional)

            An optional ground-truth segmentation as an array with gt[z][y][x].
            If given, metrics

    Returns
    -------

        segmentations

            List of segmentations (numpy arrays, uint64, 3 dimensional).

        (segmentations, metrics)

            Tuple of segmentations and metrics, if ground-truth volume was
            given.
    '''

    cdef vector[uint64_t*]            segmentation_data
    cdef np.ndarray[uint64_t, ndim=3] segmentation
    cdef uint32_t*                    gt_data = NULL

    # the C++ part assumes contiguous memory, make sure we have it (and do 
    # nothing, if we do)
    if not affs.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous affinity arrray (avoid this by passing C_CONTIGUOUS arrays)")
        affs = np.ascontiguousarray(affs)
    if gt is not None and not gt.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous ground-truth arrray (avoid this by passing C_CONTIGUOUS arrays)")
        gt = np.ascontiguousarray(gt)

    print("Preparing segmentation volumes...")

    segmentations = []
    volume_shape = (affs.shape[1], affs.shape[2], affs.shape[3])
    thresholds.sort()
    for i in range(len(thresholds)):
        segmentation = np.zeros(volume_shape, dtype=np.uint64)
        segmentations.append(segmentation)
        segmentation_data.push_back(&segmentation[0,0,0])

    if gt is not None:
        gt_data = &gt[0,0,0]

    print("Processing thresholds")

    metrics = process_thresholds(
        thresholds,
        affs.shape[1], affs.shape[2], affs.shape[3],
        &affs[0, 0, 0, 0],
        segmentation_data,
        gt_data)

    if gt is not None:
        stats = {
            'V_Rand'      : 0,
            'V_Rand_split': [],
            'V_Rand_merge': [],
            'V_Info'      : 0,
            'V_Info_split': [],
            'V_Info_merge': []
        }
        for metric in metrics:
            rand_f_score = 2.0/(1.0/metric.rand_split + 1.0/metric.rand_merge)
            voi_score = 2.0/(1.0/metric.voi_split + 1.0/metric.voi_merge)
            stats['V_Rand'] = max(stats['V_Rand'], rand_f_score)
            stats['V_Rand_split'].append(metric.rand_split)
            stats['V_Rand_merge'].append(metric.rand_merge)
            stats['V_Info'] = max(stats['V_Info'], voi_score)
            stats['V_Info_split'].append(metric.voi_split)
            stats['V_Info_merge'].append(metric.voi_merge)
        return (segmentations, stats)
    return segmentations

cdef extern from "zwatershed.h":

    struct Metrics:
        double voi_split
        double voi_merge
        double rand_split
        double rand_merge

    struct ZwatershedState:
        pass

    vector[Metrics] process_thresholds(
            vector[size_t] thresholds,
            size_t width, size_t height, size_t depth,
            np.float32_t* affs,
            vector[uint64_t*]& segmentation_data,
            uint32_t* gt_data)

    ZwatershedState get_initial_state(
            size_t width, size_t height, size_t depth,
            np.float32_t* affs)

####################
# PREVIOUS METHODS #
####################

def zwshed_with_stats_arb(np.ndarray[uint64_t, ndim=3] gt, np.ndarray[uint64_t, ndim=1] node1,
                          np.ndarray[uint64_t, ndim=1] node2, np.ndarray[float, ndim=1] edgeWeight, threshes,
                          save_threshes,
                          int h5, seg_save_path="NULL/"):

    raise RuntimeError("This method is currently not adjusted to the new c-order C++ backend")

    # NOTE:
    # This method and the downstream *_arb versions of zwatershed need to be 
    # rewritten/modified:
    # * use c-order
    # * return stats if gt is not None (make last argument)
    # * remove h5 support (delegated to user), remove relevant args
    # * rename to 'zwatershed_arb'
    # * write documentation

def zwatershed_basic_h5(np.ndarray[np.float32_t, ndim=4] affs, seg_save_path="NULL/"):

    makedirs(seg_save_path)

    # the C++ part assumes contiguous memory, make sure we have it (and do 
    # nothing, if we do)
    if not affs.flags['C_CONTIGUOUS']:
        print("Creating memory-contiguous affinity arrray (avoid this by passing C_CONTIGUOUS arrays)")
        affs = np.ascontiguousarray(affs)

    state = get_initial_state(
        affs.shape[1], affs.shape[2], affs.shape[3],
        &affs[0,0,0,0])

    # NOTE:
    # Similar to zwatershed() above, this method needs to be adjusted for the 
    # c-order backend. In particular:
    # * handle affinities and ground-truth the same way
    # * return a region graph representation to store it in h5 file

    ## get initial seg,rg
    #affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    #dims = affs.shape
    #seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    ## TODO: call process_thresholds, but return region graph after initial 
    ## watershed
    #cdef ZwatershedResult result = zwshed_initial(seg_empty, affs)
    #counts = result.counts
    #rg = result.rg
    #f = h5py.File(seg_save_path + 'basic.h5', 'w')
    #f["seg"] = np.array(result.seg, dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
    #f["counts"]=counts
    #rg_edges = np.array(rg.edges, dtype=np.uint64)
    #f["rg_edges"]=rg_edges.reshape(len(rg_edges)/2,2)
    #f["rg_weights"]=np.array(rg.weights, dtype=np.float32)
    #f.close()

def makedirs(seg_save_path):
    if not seg_save_path.endswith("/"):
        seg_save_path += "/"
    if not os.path.exists(seg_save_path):
        os.makedirs(seg_save_path)
