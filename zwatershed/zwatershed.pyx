from libcpp.list cimport list
from libcpp.map cimport map
from libcpp.vector cimport vector
from libcpp.string cimport string
from libcpp.pair cimport pair
from libc.stdint cimport uint64_t, uint32_t
import numpy as np
import os
cimport numpy as np
import h5py

def zwatershed_unified(np.ndarray[np.float32_t, ndim=4] affs, threshes, np.ndarray[uint32_t, ndim=3] gt = None):

    cdef vector[uint64_t*]            segmentation_data
    cdef np.ndarray[uint64_t, ndim=3] segmentation
    cdef uint32_t*                    gt_data

    print("Preparing segmentation volumes...")

    segmentations = []
    volume_shape = (affs.shape[1], affs.shape[2], affs.shape[3])
    threshes.sort()
    for i in range(len(threshes)):
        segmentation = np.zeros(volume_shape, dtype=np.uint64, order='F')
        segmentations.append(segmentation)
        segmentation_data.push_back(&segmentation[0,0,0])

    print("Processing thresholds")

    if gt is not None:

        gt_data = &gt[0,0,0]
        metrics = process_thresholds(
            threshes,
            affs.shape[1], affs.shape[2], affs.shape[3],
            &affs[0, 0, 0, 0],
            segmentation_data,
            gt_data)
        return (segmentations, metrics)

    process_thresholds(
        threshes,
        affs.shape[1], affs.shape[2], affs.shape[3],
        &affs[0, 0, 0, 0],
        segmentation_data)
    return segmentations

#-------------- interface methods --------------------------------------------------------------
def zwatershed_and_metrics(gt, affs, threshes, save_threshes):
    threshes.sort()
    return zwshed_with_stats(gt, affs, threshes, save_threshes, h5=0)

def zwatershed_and_metrics_h5(gt, affs, threshes, save_threshes, seg_save_path):
    threshes.sort()
    return zwshed_with_stats(gt, affs, threshes, save_threshes, h5=1, seg_save_path=seg_save_path)

def zwatershed(affs, threshes):
    threshes.sort()
    return zwshed_no_stats(affs, threshes, threshes, h5=0)

def zwatershed_h5(affs, threshes, seg_save_path):
    threshes.sort()
    zwshed_no_stats(affs, threshes, threshes, h5=1, seg_save_path=seg_save_path)

def zwatershed_and_metrics_arb(gt, node1, node2, edgeWeight, threshes, save_threshes):
    threshes.sort()
    return zwshed_with_stats_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, h5=0)

def zwatershed_and_metrics_h5_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, seg_save_path):
    threshes.sort()
    return zwshed_with_stats_arb(gt, node1, node2, edgeWeight, threshes, save_threshes, h5=1,
                                 seg_save_path=seg_save_path)

def zwatershed_arb(seg_shape, node1, node2, edgeWeight, save_threshes):
    save_threshes.sort()
    return zwshed_no_stats_arb(seg_shape, node1, node2, edgeWeight, save_threshes, save_threshes, h5=0)

def zwatershed_h5_arb(seg_shape, node1, node2, edgeWeight, save_threshes, seg_save_path):
    save_threshes.sort()
    return zwshed_no_stats_arb(seg_shape, node1, node2, edgeWeight, save_threshes, save_threshes, h5=1,
                               seg_save_path=seg_save_path)
                               
def zwatershed_basic_h5(affs, seg_save_path):
    zwatershed_basic_h5(affs, seg_save_path=seg_save_path)

#-------------- helper methods --------------------------------------------------------------
def zwatershed_basic_h5(np.ndarray[np.float32_t, ndim=4] affs, seg_save_path="NULL/"):
    makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    cdef ZwatershedResult result = zwshed_initial(seg_empty, affs)
    counts = result.counts
    rg = result.rg
    f = h5py.File(seg_save_path + 'basic.h5', 'w')
    f["seg"] = np.array(result.seg, dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
    f["counts"]=counts
    rg_edges = np.array(rg.edges, dtype=np.uint64)
    f["rg_edges"]=rg_edges.reshape(len(rg_edges)/2,2)
    f["rg_weights"]=np.array(rg.weights, dtype=np.float32)
    f.close()

def zwshed_with_stats(np.ndarray[uint64_t, ndim=3] gt, np.ndarray[np.float32_t, ndim=4] affs, threshes, save_threshes,
                      int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    gt = np.array(gt, order='F')
    cdef ZwatershedResult result = zwshed_initial(gt, affs)
    cdef np.ndarray[uint64_t, ndim=1] counts_out = result.counts

    counts_len = len(result.counts)
    dims = affs.shape

    # get segs, stats
    segs, splits, merges, info_splits, info_merges = [], [], [], [], []
    for i in range(len(threshes)):
        if(len(result.rg.edges) > 0):
            result = merge_with_stats(dims[0], dims[1], dims[2], &gt[0, 0, 0], result.rg, result.seg, &counts_out[0], counts_len, threshes[i])
        counts_out = np.array(result.counts, dtype='uint64')
        counts_len = len(counts_out)
        if threshes[i] in save_threshes:
            seg = np.array(result.seg, dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
        splits = splits + [result.stats[0]]
        merges = merges + [result.stats[1]]
        info_splits = info_splits + [result.stats[2]]
        info_merges = info_merges + [result.stats[3]]
    max_f_score = 2 / (1 / splits[0] + 1 / merges[0])
    max_v_info = 2 / (1 / info_splits[0] + 1 / info_merges[0])
    for j in range(len(splits)):
        f_score = 2 / (1 / splits[j] + 1 / merges[j])
        if f_score > max_f_score:
            max_f_score = f_score
        info_score = 2 / (1 / info_splits[j] + 1 / info_merges[j])
        if info_score > max_v_info:
            max_v_info = info_score
        
    returnMap = {'V_Rand': max_f_score, 'V_Rand_split': splits, 'V_Rand_merge': merges, 'V_Info_split':info_splits, 'V_Info_merge':info_merges, 'V_Info':max_v_info}
    if h5:
        return returnMap
    else:
        return segs, returnMap

def zwshed_no_stats(np.ndarray[np.float32_t, ndim=4] affs, threshes, save_threshes, int h5,
                    seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    affs = np.asfortranarray(np.transpose(affs, (1, 2, 3, 0)))
    dims = affs.shape
    seg_empty = np.empty((dims[0], dims[1], dims[2]), dtype='uint64')
    cdef ZwatershedResult result = zwshed_initial(seg_empty, affs)
    cdef np.ndarray[uint64_t, ndim=1] counts_out = np.array(result.counts, dtype=np.uintp)
    segs = []

    # get segs, stats
    for i in range(len(threshes)):
        if(len(result.rg.edges) > 0):
            result = merge_no_stats(dims[0], dims[1], dims[2], result.rg, result.seg, &counts_out[0], len(result.counts), threshes[i])
        counts_out = np.array(result.counts, dtype='uint64')
        counts_len = len(counts_out)
        if threshes[i] in save_threshes:
            seg = np.array(result.seg, dtype='uint64').reshape((dims[2], dims[1], dims[0])).transpose(2, 1, 0)
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
    if not h5:
        return segs

def zwshed_initial(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[np.float32_t, ndim=4] affs):
    dims = affs.shape
    return zwshed_initial_c(dims[0], dims[1], dims[2], &affs[0, 0, 0, 0])

def makedirs(seg_save_path):
    if not seg_save_path.endswith("/"):
        seg_save_path += "/"
    if not os.path.exists(seg_save_path):
        os.makedirs(seg_save_path)

# arb methods ---------------
def zwshed_with_stats_arb(np.ndarray[uint64_t, ndim=3] gt, np.ndarray[uint64_t, ndim=1] node1,
                          np.ndarray[uint64_t, ndim=1] node2, np.ndarray[float, ndim=1] edgeWeight, threshes,
                          save_threshes,
                          int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    gt = np.array(gt, order='C')
    n_edge = node1.size
    map = zwshed_initial_arb(gt, node1, node2, n_edge, edgeWeight)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']

    counts_len = len(map['counts'])
    dims = gt.shape

    # get segs, stats
    segs, splits, merges = [], [], []
    for i in range(len(threshes)):
        # TODO: use RegionGraph
        #if(len(rgn_graph) > 0):
            #map = merge_with_stats_arb(dims[0], dims[1], dims[2], &gt[0, 0, 0], &rgn_graph[0, 0],
                                   #rgn_graph.shape[0], &seg_in[0], &counts_out[0], counts_len, threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[0], dims[1], dims[2]))
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
        splits = splits + [map['stats'][0]]
        merges = merges + [map['stats'][1]]
    max_f_score = 2 / (1 / splits[0] + 1 / merges[0])
    for j in range(len(splits)):
        f_score = 2 / (1 / splits[j] + 1 / merges[j])
        if f_score > max_f_score:
            max_f_score = f_score
    returnMap = {'V_Rand': max_f_score, 'V_Rand_split': splits, 'V_Rand_merge': merges}
    if h5:
        return returnMap
    else:
        return segs, returnMap

def zwshed_no_stats_arb(dims, np.ndarray[uint64_t, ndim=1] node1,
                        np.ndarray[uint64_t, ndim=1] node2, np.ndarray[float, ndim=1] edgeWeight, threshes,
                        save_threshes,
                        int h5, seg_save_path="NULL/"):
    if h5:
        makedirs(seg_save_path)

    # get initial seg,rg
    n_edge = node1.size
    seg_empty = np.zeros(dims,dtype='uint64')
    map = zwshed_initial_arb(seg_empty, node1, node2, n_edge, edgeWeight)
    cdef np.ndarray[uint64_t, ndim=1] seg_in = map['seg']
    cdef np.ndarray[uint64_t, ndim=1] counts_out = map['counts']
    cdef np.ndarray[np.float32_t, ndim=2] rgn_graph = map['rg']
    counts_len = len(map['counts'])

    # get segs, stats
    segs = []
    for i in range(len(threshes)):
        # TODO: use RegionGraph
        #if(len(rgn_graph) > 0):
            #map = merge_no_stats_arb(dims[0], dims[1], dims[2], &rgn_graph[0, 0],
                                 #rgn_graph.shape[0], &seg_in[0], &counts_out[0], counts_len, threshes[i])
        seg = np.array(map['seg'], dtype='uint64').reshape((dims[0], dims[1], dims[2]))
        graph = np.array(map['rg'], dtype='float32')
        counts_out = np.array(map['counts'], dtype='uint64')
        counts_len = len(counts_out)
        seg_in = np.array(map['seg'], dtype='uint64')
        rgn_graph = graph.reshape(len(graph) / 3, 3)
        if threshes[i] in save_threshes:
            if h5:
                f = h5py.File(seg_save_path + 'seg_' + str(threshes[i]) + '.h5', 'w')
                f["main"] = seg
                f.close()
            else:
                segs.append(seg)
    if not h5:
        return segs

def zwshed_initial_arb(np.ndarray[uint64_t, ndim=3] seg, np.ndarray[uint64_t, ndim=1] node1,
                       np.ndarray[uint64_t, ndim=1] node2, size_t n_edge, np.ndarray[float, ndim=1] edgeWeight):
    cdef np.ndarray[uint64_t, ndim=1] counts = np.empty(1, dtype='uint64')
    dims = seg.shape
    map = zwshed_initial_c_arb(dims[0], dims[1], dims[2], &node1[0], &node2[0], &edgeWeight[0], n_edge)
    graph = np.array(map['rg'], dtype='float32')
    return {'rg': graph.reshape(len(graph) / 3, 3), 'seg': np.array(map['seg'], dtype='uint64'),
            'counts': np.array(map['counts'], dtype='uint64')}



#-------------- c++ methods --------------------------------------------------------------
cdef extern from "zwatershed.h":

    struct Metrics:
        double voi_split
        double voi_merge
        double rand_split
        double rand_merge

    vector[Metrics] process_thresholds(vector[size_t] thresholds, size_t width, size_t height, size_t depth, np.float32_t* affs, vector[uint64_t*]& segmentation_data, uint32_t* gt_data)
    vector[Metrics] process_thresholds(vector[size_t] thresholds, size_t width, size_t height, size_t depth, np.float32_t* affs, vector[uint64_t*]& segmentation_data)

    #struct RegionGraphEdge {

        #uint64_t id1;
        #uint64_t id2;
    #};

    #struct RegionGraph {

        #std::vector<RegionGraphEdge> edges;
        #std::vector<float>           weights;
    #};

    #struct ZwatershedResult {

        #std::vector<RegionGraph> rg;
        #std::vector<uint64_t> seg;
        #std::vector<size_t> counts;
        #std::vector<double> stats;
    #};

    struct RegionGraphEdge:
        uint64_t id1
        uint64_t id2

    struct RegionGraph:
        vector[RegionGraphEdge] edges
        vector[float]           weights

    struct ZwatershedResult:
        RegionGraph rg
        vector[uint64_t] seg
        vector[size_t] counts
        vector[double] stats

    ZwatershedResult zwshed_initial_c(size_t dimX, size_t dimY, size_t dimZ, np.float32_t*affs)
    ZwatershedResult merge_with_stats(size_t dx, size_t dy, size_t dz, np.uint64_t*gt,
                                                 RegionGraph& rgn_graph, const vector[uint64_t]& seg,
                                                 uint64_t*counts, size_t counts_len, size_t thresh)
    ZwatershedResult merge_no_stats(size_t dx, size_t dy, size_t dz, RegionGraph& rgn_graph, const vector[uint64_t]& seg,
                                               uint64_t*counts, size_t counts_len, size_t thresh)
    map[string, list[float]] zwshed_initial_c_arb(size_t dimX, size_t dimY, size_t dimZ, uint64_t*node1,
                                                  uint64_t*node2, float*edgeWeight, size_t n_edge)
    map[string, vector[double]] merge_with_stats_arb(size_t dx, size_t dy, size_t dz, np.uint64_t*gt,
                                                     RegionGraph& rgn_graph, uint64_t*seg,
                                                     uint64_t*counts, size_t counts_len, size_t thresh)
    map[string, vector[double]] merge_no_stats_arb(size_t dx, size_t dy, size_t dz,
                                                   RegionGraph& rgn_graph, uint64_t*seg,
                                                   uint64_t*counts, size_t counts_len, size_t thresh)
