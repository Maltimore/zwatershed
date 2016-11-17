#ifndef ITERATIVE_REGION_MERGING_H__
#define ITERATIVE_REGION_MERGING_H__

#include <vector>
#include <map>
#include <queue>
#include <limits>

#include "types.hpp"

template <typename NodeIdType, typename AffinityType, typename ScoreType = AffinityType>
class IterativeRegionMerging {

public:

	typedef region_graph<NodeIdType, AffinityType>    RegionGraphType;
	typedef RegionGraphEdge<NodeIdType, AffinityType> EdgeType;
	typedef std::size_t                               EdgeIdType;

	static const EdgeIdType NoEdge = std::numeric_limits<EdgeIdType>::max();

	/**
	 * Create a region merging for the given initial RAG.
	 */
	IterativeRegionMerging(region_graph_ptr<NodeIdType, AffinityType> initialRegionGraph) :
		_regionGraphPointer(initialRegionGraph),
		_regionGraph(*initialRegionGraph),
		_edgeQueue(EdgeCompare(_edgeScores)),
		_mergedUntil(0) {

		NodeIdType maxNodeId = 0;

		// prepare reverse lookup for nodes to edges
		for (EdgeIdType e = 0; e < _regionGraph.size(); e++) {

			NodeIdType u = _regionGraph[e].id1;
			NodeIdType v = _regionGraph[e].id2;

			_incidentEdges[u].push_back(e);
			_incidentEdges[v].push_back(e);

			maxNodeId = std::max(maxNodeId, std::max(u, v));
		}

		_affiliatedEdges.resize(_regionGraph.size());
		_edgeScores.resize(_regionGraph.size());
		_nextNodeId = maxNodeId + 1;
	}

	/**
	 * Merge a RAG with the given edge scoring function until the given threshold.
	 */
	template <typename EdgeScoringFunction>
	void mergeUntil(
			const EdgeScoringFunction& edgeScoringFunction,
			ScoreType threshold) {

		// compute scores of each edge not scored so far
		if (_mergedUntil == 0)
			for (EdgeIdType e = 0; e < _regionGraph.size(); e++)
				scoreEdge(e, edgeScoringFunction);

		// while there are still unhandled edges
		while (_edgeQueue.size() > 0) {

			// get the next cheapest edge to merge
			EdgeIdType next = _edgeQueue.top();
			_edgeQueue.pop();

			// stop, if the threshold got exceeded
			ScoreType score = _edgeScores[next];
			if (score >= threshold)
				break;

			NodeIdType u = _regionGraph[next].id1;
			NodeIdType v = _regionGraph[next].id2;

			// skip if incident regions already got merged
			if (!isRoot(u) || !isRoot(v))
				continue;

			mergeRegions(u, v, edgeScoringFunction);
		}

		_mergedUntil = threshold;
	}

	/**
	 * Get the segmentation corresponding to the current merge level.
	 *
	 * The provided segmentation has to hold the initial segmentation, or any 
	 * segmentation created by previous calls to extractSegmentation(). In other 
	 * words, it has to hold IDs that have been seen before.
	 */
	template <typename SegmentationVolume>
	void extractSegmentation(SegmentationVolume& segmentation) {

		for (std::size_t i = 0; i < segmentation.num_elements(); i++)
			segmentation.data()[i] = getRoot(segmentation.data()[i]);
	}

private:

	/**
	 * Compare two edges based on their score. To be used in the priority queue.
	 */
	class EdgeCompare {

	public:

		EdgeCompare(const std::vector<ScoreType>& edgeScores) :
			_edgeScores(edgeScores) {}

		bool operator()(const EdgeIdType a, const EdgeIdType b) {

			return _edgeScores[a] > _edgeScores[b];
		}

	private:

		const std::vector<ScoreType>& _edgeScores;
	};

	/**
	 * Merge regions a and b.
	 */
	template <typename EdgeScoringFunction>
	void mergeRegions(
			NodeIdType a,
			NodeIdType b,
			const EdgeScoringFunction& edgeScoringFunction) {

		// create a new node c = a + b
		NodeIdType c = _nextNodeId;
		_nextNodeId++;

		// set parents
		_rootPaths[a] = c;
		_rootPaths[b] = c;

		std::vector<EdgeIdType>& cEdges = _incidentEdges[c];

		// connect c to neighbors of a and b, and update affiliatedEdges

		// for each child region
		for (NodeIdType child : { a, b } ) {

			// for all neighbors of child
			for (EdgeIdType neighborEdge : _incidentEdges[child]) {

				NodeIdType neighbor = (
						_regionGraph[neighborEdge].id1 == child ?
						_regionGraph[neighborEdge].id2 :
						_regionGraph[neighborEdge].id1);

				// don't consider already merged regions
				if (!isRoot(neighbor))
					continue;

				EdgeIdType newEdge = findEdge(c, neighbor, cEdges);

				if (newEdge == NoEdge) {

					// add the edge from c -> neighbor, with temporary affinity of 0
					newEdge = addEdge(c, neighbor, 0);
					cEdges.push_back(newEdge);
				}

				// add affiliated edges to new edge
				if (_affiliatedEdges[neighborEdge].size() == 0) // initial edge
					_affiliatedEdges[newEdge].push_back(_regionGraph[neighborEdge]);
				else
					std::copy(
							_affiliatedEdges[neighborEdge].begin(),
							_affiliatedEdges[neighborEdge].end(),
							std::back_inserter(_affiliatedEdges[newEdge]));

				// clear affiliated edges of merged region edges -- they are not 
				// needed anymore
				_affiliatedEdges[neighborEdge].clear();
			}
		}

		// score all new edges == incident edges to c
		// it is expected that the scoring function updated the affinities, if it 
		// needs them
		for (EdgeIdType e : cEdges)
			scoreEdge(e, edgeScoringFunction);
	}

	/**
	 * Create a new edge between u and v.
	 */
	EdgeIdType addEdge(NodeIdType u, NodeIdType v, AffinityType affinity) {

		_regionGraph.push_back(EdgeType(u, v, affinity));
		EdgeIdType newEdge = _regionGraph.size() - 1;

		_affiliatedEdges.push_back(std::vector<EdgeType>());
		_edgeScores.push_back(0);

		return newEdge;
	}

	/**
	 * Find the edge connecting u and v. Returns NoEdge, if there is none.
	 */
	inline EdgeIdType findEdge(NodeIdType u, NodeIdType v) {

		return findEdge(u, v, _regionGraph);
	}

	/**
	 * Same as findEdge(u, v), but restricted to edges in pool.
	 */
	inline EdgeIdType findEdge(NodeIdType u, NodeIdType v, const std::vector<EdgeIdType>& pool) {

		NodeIdType min = std::min(u, v);
		NodeIdType max = std::max(u, v);

		for (EdgeIdType e : pool) {

			if (std::min(_regionGraph[e].id1, _regionGraph[e].id2) == min &&
				std::max(_regionGraph[e].id1, _regionGraph[e].id2) == max)
				return e;
		}

		return NoEdge;
	}

	/**
	 * Score edge i.
	 */
	template <typename EdgeScoringFunction>
	void scoreEdge(EdgeIdType e, const EdgeScoringFunction& edgeScoringFunction) {

		const EdgeType& edge = _regionGraph[e];
		std::vector<EdgeType>& affiliatedEdges = _affiliatedEdges[e];

		_edgeScores[e] = edgeScoringFunction(edge, affiliatedEdges);
		_edgeQueue.push(e);
	}

	inline bool isRoot(NodeIdType id) {

		// if there is no root path, it is a root
		return (_rootPaths.count(id) != 0);
	}

	/**
	 * Get the root node of a merge-tree.
	 */
	NodeIdType getRoot(NodeIdType id) {

		NodeIdType root = id;

		while (!isRoot(root))
			root = _rootPaths.at(root);

		// not compressed, yet
		if (_rootPaths.at(id) != root)
			while (id != root) {

				NodeIdType next = _rootPaths.at(id);
				_rootPaths[id] = root;
				id = next;
			}

		return root;
	}

	// shared pointer just to protect lifetime of region graph
	region_graph_ptr<NodeIdType, AffinityType> _regionGraphPointer;
	RegionGraphType& _regionGraph;

	// reverse look-up from nodes to edges
	std::map<NodeIdType, std::vector<EdgeIdType>> _incidentEdges;

	// for every new edge between regions u and v, the edges of the initial RAG 
	// between any child of u and any child of v
	//
	// initial edges will have this empty
	//
	// congruent to _regionGraph
	std::vector<std::vector<EdgeType>> _affiliatedEdges;

	// the score of each edge
	//
	// congruent to _regionGraph
	std::vector<ScoreType> _edgeScores;

	// sorted list of edges indices, cheapest edge first
	std::priority_queue<EdgeIdType, std::vector<EdgeIdType>, EdgeCompare> _edgeQueue;

	// paths from nodes to the roots of the merge-tree they are part of
	//
	// root nodes are not in the map
	//
	// paths will be compressed when read
	std::map<NodeIdType, NodeIdType> _rootPaths;

	// current state of merging
	ScoreType _mergedUntil;

	NodeIdType _nextNodeId;
};

#endif // ITERATIVE_REGION_MERGING_H__

