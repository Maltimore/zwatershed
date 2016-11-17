#ifndef REGION_GRAPH_H__
#define REGION_GRAPH_H__

#include <vector>
#include <limits>

template <typename ID>
struct RegionGraphEdge {

	typedef ID NodeIdType;

	NodeIdType u;
	NodeIdType v;

	RegionGraphEdge() : u(0), v(0) {}
	RegionGraphEdge(NodeIdType u_, NodeIdType v_) : u(u_), v(v_) {}
};

// forward declaration
template <typename ID>
class RegionGraph;

template<typename ID>
class RegionGraphEdgeMapBase {

protected:

	RegionGraphEdgeMapBase(RegionGraph<ID>& regionGraph) :
		_regionGraph(regionGraph) {

		_regionGraph.registerEdgeMap(this);
	}

	virtual ~RegionGraphEdgeMapBase() {

		_regionGraph.deregisterEdgeMap(this);
	}

private:

	friend RegionGraph<ID>;

	virtual void onNewEdge(std::size_t id) = 0;

	RegionGraph<ID>& _regionGraph;
};

template<typename ID, typename T, typename Container>
class RegionGraphEdgeMap : public RegionGraphEdgeMapBase<ID> {

public:

	RegionGraphEdgeMap(RegionGraph<ID>& regionGraph) :
		RegionGraphEdgeMapBase<ID>(regionGraph),
		_values(regionGraph.edges().size()) {}

	inline const T& operator[](std::size_t i) const { return _values[i]; }
	inline T& operator[](std::size_t i) { return _values[i]; }

private:

	void onNewEdge(std::size_t id) {

		_values.emplace_back();
	}

	Container _values;
};

template <typename ID>
class RegionGraph {

public:

	typedef ID                          NodeIdType;
	typedef std::size_t                 EdgeIdType;

	typedef RegionGraphEdge<NodeIdType> EdgeType;

	template <typename T, typename Container = std::vector<T>>
	using EdgeMap = RegionGraphEdgeMap<ID, T, Container>;

	static const EdgeIdType NoEdge = std::numeric_limits<EdgeIdType>::max();

	RegionGraph(ID numNodes = 0) :
		_numNodes(numNodes),
		_incEdges(numNodes) {}

	ID addNode() {

		NodeIdType id = _numNodes;
		_numNodes++;

		_incEdges.emplace_back();

		return id;
	}

	std::size_t addEdge(NodeIdType u, NodeIdType v) {

		std::size_t id = _edges.size();
		_edges.push_back(EdgeType(u, v));

		_incEdges[u].push_back(id);
		_incEdges[v].push_back(id);

		for (RegionGraphEdgeMapBase<ID>* map : _edgeMaps)
			map->onNewEdge(id);

		return id;
	}

	inline const EdgeType& edge(EdgeIdType e) const { return _edges[e]; }

	inline const std::vector<EdgeType>& edges() const { return _edges; }

	inline const std::vector<EdgeIdType>& incEdges(ID node) const { return _incEdges[node]; }

	inline NodeIdType getOpposite(NodeIdType n, EdgeIdType e) const {

		return (_edges[e].u == n ? _edges[e].v : _edges[e].u);
	}

	/**
	 * Find the edge connecting u and v. Returns NoEdge, if there is none.
	 */
	inline EdgeIdType findEdge(NodeIdType u, NodeIdType v) {

		return findEdge(u, v, (_incEdges[u].size() < _incEdges[v].size() ? _incEdges[u] : _incEdges[v]));
	}

	/**
	 * Same as findEdge(u, v), but restricted to edges in pool.
	 */
	inline EdgeIdType findEdge(NodeIdType u, NodeIdType v, const std::vector<EdgeIdType>& pool) {

		NodeIdType min = std::min(u, v);
		NodeIdType max = std::max(u, v);

		for (EdgeIdType e : pool)
			if (std::min(_edges[e].u, _edges[e].v) == min &&
				std::max(_edges[e].u, _edges[e].v) == max)
				return e;

		return NoEdge;
	}

private:

	friend RegionGraphEdgeMapBase<ID>;

	void registerEdgeMap(RegionGraphEdgeMapBase<ID>* edgeMap) {

		_edgeMaps.push_back(edgeMap);
	}

	void deregisterEdgeMap(RegionGraphEdgeMapBase<ID>* edgeMap) {

		auto it = std::find(_edgeMaps.begin(), _edgeMaps.end(), edgeMap);
		if (it != _edgeMaps.end())
			_edgeMaps.erase(it);
	}

	ID _numNodes;

	std::vector<EdgeType> _edges;

	std::vector<std::vector<EdgeIdType>> _incEdges;

	std::vector<RegionGraphEdgeMapBase<ID>*> _edgeMaps;
};

#endif // REGION_GRAPH_H__

