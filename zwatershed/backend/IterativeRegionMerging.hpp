#ifndef ITERATIVE_REGION_MERGING_H__
#define ITERATIVE_REGION_MERGING_H__

#include <vector>
#include <map>
#include <queue>

#include "RegionGraph.hpp"

template <typename NodeIdType, typename AffinityType, typename ScoreType = AffinityType>
class IterativeRegionMerging {

public:

	typedef RegionGraph<NodeIdType>              RegionGraphType;
	typedef typename RegionGraphType::EdgeType   EdgeType;
	typedef typename RegionGraphType::EdgeIdType EdgeIdType;

	/**
	 * Create a region merging for the given initial RAG.
	 */
	IterativeRegionMerging(
			std::shared_ptr<RegionGraphType> initialRegionGraph) :
		_regionGraphPointer(initialRegionGraph),
		_regionGraph(*initialRegionGraph),
		_affiliatedEdges(*initialRegionGraph),
		_edgeScores(*initialRegionGraph),
		_edgeQueue(EdgeCompare(_edgeScores)),
		_mergedUntil(0) {}

	/**
	 * Merge a RAG with the given edge scoring function until the given threshold.
	 */
	template <typename EdgeScoringFunction>
	void mergeUntil(
			const EdgeScoringFunction& edgeScoringFunction,
			ScoreType threshold) {

		if (threshold <= _mergedUntil) {

			std::cout << "already merged until " << threshold << ", skipping" << std::endl;
			return;
		}

		// compute scores of each edge not scored so far
		if (_mergedUntil == 0)
			for (EdgeIdType e = 0; e < _regionGraph.edges().size(); e++)
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

			NodeIdType u = _regionGraph.edge(next).u;
			NodeIdType v = _regionGraph.edge(next).v;

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

		EdgeCompare(const typename RegionGraphType::template EdgeMap<ScoreType>& edgeScores) :
			_edgeScores(edgeScores) {}

		bool operator()(const EdgeIdType a, const EdgeIdType b) {

			return _edgeScores[a] > _edgeScores[b];
		}

	private:

		const typename RegionGraphType::template EdgeMap<ScoreType>& _edgeScores;
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
		NodeIdType c = _regionGraph.addNode();

		// set parents
		_rootPaths[a] = c;
		_rootPaths[b] = c;

		// connect c to neighbors of a and b, and update affiliatedEdges

		// for each child region
		for (NodeIdType child : { a, b } ) {

			// for all neighbors of child
			for (EdgeIdType neighborEdge : _regionGraph.incEdges(child)) {

				NodeIdType neighbor = _regionGraph.getOpposite(child, neighborEdge);

				// don't consider already merged regions
				if (!isRoot(neighbor))
					continue;

				// do we already have an edge to this neighbor?
				EdgeIdType newEdge = _regionGraph.findEdge(c, neighbor);

				// if not, add the edge c
				if (newEdge == RegionGraphType::NoEdge)
					newEdge = _regionGraph.addEdge(c, neighbor);

				// TODO: delegate affiliated edge book-keeping to scoring 
				// function

				// add affiliated edges to new edge
				if (_affiliatedEdges[neighborEdge].size() == 0)
					// neighborEdge is an initial edge, add it to the new edge 
					// affiliated edge list
					_affiliatedEdges[newEdge].push_back(neighborEdge);
				else
					// neighborEdge is a compound edge, copy its affiliated 
					// edges to the new affiliated edge list
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
		for (EdgeIdType e : _regionGraph.incEdges(c))
			scoreEdge(e, edgeScoringFunction);
	}

	/**
	 * Score edge i.
	 */
	template <typename EdgeScoringFunction>
	ScoreType scoreEdge(EdgeIdType e, const EdgeScoringFunction& edgeScoringFunction) {

		std::vector<EdgeIdType>& affiliatedEdges = _affiliatedEdges[e];
		ScoreType score = edgeScoringFunction(e, affiliatedEdges);

		_edgeScores[e] = score;
		_edgeQueue.push(e);

		return score;
	}

	inline bool isRoot(NodeIdType id) {

		// if there is no root path, it is a root
		return (_rootPaths.count(id) == 0);
	}

	/**
	 * Get the root node of a merge-tree.
	 */
	NodeIdType getRoot(NodeIdType id) {

		// early way out
		if (isRoot(id))
			return id;

		// walk up to root

		NodeIdType root = _rootPaths.at(id);
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
	std::shared_ptr<RegionGraphType> _regionGraphPointer;
	RegionGraphType& _regionGraph;

	// for every new edge between regions u and v, the edges of the initial RAG 
	// between any child of u and any child of v
	//
	// initial edges will have this empty
	//
	// congruent to _regionGraph
	typename RegionGraphType::template EdgeMap<std::vector<EdgeIdType>> _affiliatedEdges;

	// the score of each edge
	//
	// congruent to _regionGraph
	typename RegionGraphType::template EdgeMap<ScoreType> _edgeScores;

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

