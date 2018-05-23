#ifndef GRAPH_H
// ReSharper disable once CppInconsistentNaming
#define _CRT_SECURE_NO_WARNINGS
#define GRAPH_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <set>
#include <fstream>
#include <stack>
#include <queue>

class Dsu {
	std::vector<unsigned int> parent_;
	std::vector<unsigned int> rank_;
public:
	Dsu(const unsigned int size) {
		parent_.resize(size);
		rank_.resize(size);
		for (unsigned int i = 0; i < size; ++i)
			parent_[i] = i;

	}
	unsigned int find(const unsigned int x) {
		if (x == parent_[x])
			return x;
		return parent_[x] = find(parent_[x]);
	}
	void unite(unsigned int x, unsigned int y) {
		x = find(x);
		y = find(y);
		if (x != y) {
			if (rank_[x] < rank_[y])
				std::swap(x, y);
			parent_[y] = x;
			if (rank_[x] == rank_[y])
				++rank_[x];
		}
	}
};


class IGraph {
public:
	virtual void readGraph(const std::string &file_name) = 0;

	virtual void addEdge(int from, int to, int weight) = 0;

	virtual void removeEdge(int from, int to) = 0;

	virtual int changeEdge(int from, int to, int new_weight) = 0;

	virtual void writeGraph(const std::string &file_name) = 0;

	virtual void reset() = 0;

	virtual void reset(unsigned int vertices_count, unsigned int is_directional, unsigned int is_weighted) = 0;

	virtual ~IGraph() = default;

	virtual void feelGraph(IGraph * g) = 0;

};



class AdjacencyListGraph;

class AdjacencyMatrixGraph;

class EdgeListGraph;


class AdjacencyMatrixGraph : public IGraph {

	std::vector<std::vector<int>> graph_;
	unsigned int vertices_count_ = 0;
	unsigned int is_directional_ = 0;
	unsigned int is_weighted_ = 0;

public:


	AdjacencyMatrixGraph(const unsigned int vectices_count = 0, const unsigned int is_directional = 0, const unsigned int is_weighted = 0) {
		AdjacencyMatrixGraph::reset(vectices_count, is_directional, is_weighted);
	}

	void feelGraph(IGraph *graph) override {
		graph->reset(vertices_count_, is_directional_, is_weighted_);
		for (unsigned int u = 1; u <= vertices_count_; ++u)
			for (unsigned int v = 1; v <= vertices_count_; ++v)
				if (graph_[u][v])
					graph->addEdge(u, v, graph_[u][v]);
	}

	void readGraph(const std::string &file_name) override {
		reset();
		std::ifstream stream;
		stream.open(file_name);
		{
			char type;
			stream >> type >> vertices_count_ >> is_directional_ >> is_weighted_;
			graph_.resize(vertices_count_ + 1);
			for (unsigned int i = 1; i <= vertices_count_; ++i)
				graph_[i].resize(vertices_count_ + 1);
			for (unsigned int from = 1; from <= vertices_count_; ++from)
				for (unsigned int to = 1; to <= vertices_count_; ++to)
					stream >> graph_[from][to];
		}
		stream.close();
	}


	void writeGraph(const std::string &file_name) override {
		std::ofstream stream;
		stream.open(file_name);
		{
			stream << "C " << vertices_count_ << std::endl << is_directional_ << " " << is_weighted_ << std::endl;
			for (unsigned int from = 1; from <= vertices_count_; ++from) {
				for (unsigned int to = 1; to <= vertices_count_; ++to) {
					stream << graph_[from][to];
					if (to == vertices_count_) stream << std::endl;
					else stream << " ";
				}
			}
		}
		stream.close();
	}

	void addEdge(const int from, const int to, const int weight) override {
		// TODO: Correct way?
		// isWeighted = isWeighted || weight > 1 ? 1 : 0;
		// or
		// if(!isWeighted) weight = 1;
		if (!is_directional_) graph_[to][from] = weight;
		graph_[from][to] = weight;
	}

	void removeEdge(const int from, const int to) override {
		if (!is_directional_) graph_[to][from] = 0;
		graph_[from][to] = 0;
	}

	int changeEdge(const int from, const int to, const int new_weight) override {
		// TODO: Correct way?
		// isWeighted = isWeighted || newWeight > 1 ? 1 : 0;
		// or
		// if(!isWeighted) newWeight = min(newWeight, 1);
		if (!is_directional_) graph_[to][from] = new_weight;
		const int old_weight = graph_[from][to];
		graph_[from][to] = new_weight;
		return old_weight;
	}

	void reset() override {
		for (unsigned int i = 1; i <= vertices_count_; ++i) graph_[i].clear();
		graph_.clear();
		vertices_count_ = 0;
	}

	void reset(const unsigned int vertices_count, const unsigned int is_directional, const unsigned int is_weighted) override {
		reset();
		vertices_count_ = vertices_count;
		is_directional_ = is_directional;
		is_weighted_ = is_weighted;
		graph_.resize(vertices_count + 1);
		for (unsigned int i = 1; i <= vertices_count; ++i) graph_[i].resize(vertices_count + 1);
	}
};

class AdjacencyListGraph : public IGraph {

	std::vector<std::map<int, int>> graph_;
	unsigned int vertices_count_ = 0;
	unsigned int is_directional_ = 0;
	unsigned int is_weighted_ = 0;

	static std::vector<int> extractIntSequence(char *str) {
		std::vector<int> sequence;
		int number = 0;
		bool empty = true;
		while (str != nullptr && *str != '\0' && *str != '\n') {
			empty = false;
			if (isdigit(*str))
				number = number * 10 + (*str - '0');
			else {
				sequence.push_back(number);
				number = 0;
			}
			++str;
		}
		if (!empty)
			sequence.push_back(number);
		return sequence;
	}
	static std::vector<int> extractIntSequence(std::string str) {
		std::vector<int> sequence;
		int number = 0;
		for (char i : str) {
			if (isdigit(i)) {
				number = number * 10 + (i - '0');
			}
			else {
				sequence.push_back(number);
				number = 0;
			}
		}
		if (!str.empty())
			sequence.push_back(number);
		return sequence;
	}




public:
	AdjacencyListGraph(const unsigned int vectices_count = 0, const unsigned int is_directional = 0, const unsigned int is_weighted = 0) {
		AdjacencyListGraph::reset(vectices_count, is_directional, is_weighted);
	}

	void readGraph(const std::string &file_name) override {
		reset();
		std::ifstream stream;
		stream.open(file_name);
		{
			char type;
			stream >> type >> vertices_count_ >> is_directional_ >> is_weighted_;
			graph_.resize(vertices_count_ + 1);
			int weight = 1;
			int to = 0;
			std::string line;
			for (unsigned int from = 0; from <= vertices_count_; ++from) {
				std::getline(stream, line);
				std::vector<int> numbers = extractIntSequence(line);
				for (unsigned int idx = 0; idx < numbers.size(); ++idx) {
					to = numbers[idx];
					if (is_weighted_) weight = numbers[++idx];
					graph_[from].emplace(to, weight);
					if(!is_directional_)
						graph_[to].emplace(from, weight);

				}
			}
		}
		stream.close();
	}

	void writeGraph(const std::string &file_name) override {
		std::ofstream stream;
		stream.open(file_name);
		{
			stream << "L " << vertices_count_ << std::endl << is_directional_ << " " << is_weighted_ << std::endl;
			for (unsigned int i = 1; i <= vertices_count_; ++i) {
				bool is_first = true;
				for (const auto to : graph_[i]) {
					if (!is_first) stream << " ";
					stream << to.first;
					if (is_weighted_) stream << " " << to.second;
					is_first = false;
				}
				stream << std::endl;
			}
		}
		stream.close();
	}

	void addEdge(const int from, const int to, const int weight) override {
		graph_[from][to] = weight;
		if (!is_directional_) graph_[to][from] = weight;
	}

	void removeEdge(const int from, const int to) override {
		graph_[from].erase(to);
		if (!is_directional_) graph_[to].erase(from);
	}

	int changeEdge(const int from, const int to, const int new_weight) override {
		const int old_weight = graph_[from][to];
		graph_[from][to] = new_weight;
		if (!is_directional_) graph_[to][from] = new_weight;
		return old_weight;
	}

	void reset() override {
		for (unsigned int i = 1; i <= vertices_count_; ++i) graph_[i].clear();
		graph_.clear();
		vertices_count_ = 0;
	}

	void reset(const unsigned int vertices_count, const unsigned int is_directional, const unsigned int is_weighted) override {
		reset();
		vertices_count_ = vertices_count;
		is_directional_ = is_directional;
		is_weighted_ = is_weighted;
		graph_.resize(vertices_count + 1);
	}


	void feelGraph(IGraph *graph) override {
		graph->reset(vertices_count_, is_directional_, is_weighted_);
		for (unsigned int u = 1; u <= vertices_count_; ++u)
			for (const auto v : graph_[u]) {
				graph->addEdge(u, v.first, v.second);
				if (!is_directional_)
					graph->addEdge(v.first, u, v.second);
			}
	}

	AdjacencyListGraph* getSpaingTreePrima() {
		AdjacencyListGraph* result = new AdjacencyListGraph(vertices_count_, 0, is_weighted_);

		std::set<std::pair<int, int>> prior;
		std::vector<int> dists(vertices_count_ + 1, std::numeric_limits<int>::max());
		std::vector<int> ends(vertices_count_ + 1, -1);
		dists[1] = 0;
		prior.emplace(0, 1);
		while (!prior.empty()) {
			if (prior.empty())
				exit(1); // graph has few cc
			const int v = prior.begin()->second;
			prior.erase(prior.begin());
			if (ends[v] != -1)
				result->addEdge(v, ends[v], dists[v]);

			for (auto &edge : graph_[v]) {
				int to = edge.first, w = edge.second;
				if (w < dists[to])
				{
					const std::pair<int, int> item = std::make_pair(dists[to], to);
					if (prior.count(item))
						prior.erase(prior.find(item));
					prior.emplace(w, to);
					dists[to] = w;
					ends[to] = v;
				}
			}

		}

		return result;
	}

	std::vector<int> getVerticesDegrees() const {
		std::vector<int> degrees(vertices_count_ + 1);
		for (unsigned int i = 1; i <= vertices_count_; ++i) {
			for (const auto v : graph_[i]) {
				if (is_directional_)
					++degrees[i];
				++degrees[v.first];
			}
		}
		return degrees;
	}

	std::vector<int> getComponentsSizes() {
		Dsu dsu(vertices_count_ + 1);
		for (unsigned int u = 1; u <= vertices_count_; ++u) {
			for (const auto v : graph_[u]) {
				dsu.unite(u, v.first);
			}
		}
		std::vector<int> sizes(vertices_count_ + 1);
		for (unsigned int i = 1; i <= vertices_count_; ++i)
			++sizes[dsu.find(i)];
		return sizes;
	}

	void bridgeDfs(int v, int & timer, std::vector<bool> & visited, std::vector<int> & tin, std::vector<int> & fup, std::set<std::pair<int, int>> & bridges, const int parent = -1) {
		visited[v] = true;
		tin[v] = fup[v] = timer++;
		for (const auto edge : graph_[v]) {
			auto to = edge.first;
			if (to == parent)
				continue;
			if (visited[to])
				fup[v] = std::min(fup[v], tin[to]);
			else {
				bridgeDfs(to, timer, visited, tin, fup, bridges, v);
				fup[v] = std::min(fup[v], fup[to]);
				if (fup[to] > tin[v])
					bridges.emplace(v, to);
			}
		}
	}

	bool isBridge(const int u, const int v) {
		std::vector<bool> visited(vertices_count_ + 1);
		std::vector<int> tin(vertices_count_ + 1);
		std::vector<int> fup(vertices_count_ + 1);
		std::set<std::pair<int, int>> bridges;
		int timer = 0;
		bridgeDfs(v, timer, visited, tin, fup, bridges);
		for (const auto bridge : bridges) {
			if ((bridge.first == u && bridge.second == v) || (bridge.first == v && bridge.second == u))
				return true;
		}
		return false;
	}

	std::vector<int> getEuleranTourFleri(int u) {
		std::vector<int> res;
		res.push_back(u);
		std::vector<int> degrees = getVerticesDegrees();
		while (true) {
			const int deg = degrees[u];
			if (!deg)
				break;
			bool moved = false;
			for (const auto edge : graph_[u]) {
				const int v = edge.first;
				if (deg == 1 || !isBridge(u, v)) {
					moved = true;
					removeEdge(u, v);
					--degrees[u], --degrees[v];
					u = v;
					res.push_back(u);
					break;
				}
			}
			if (!moved) {
				const int v = graph_[u].begin()->first;
				removeEdge(u, v);
				--degrees[u], --degrees[v];
				u = v;
				res.push_back(u);
			}

		}

		return res;
	}

	std::vector<int> getEuleranTourEffective(int u) {
		std::vector<int> res;
		std::stack<int> stack;
		stack.push(u);
		while (!stack.empty()) {
			u = stack.top();
			for (auto& edge : graph_[u]) {
				int v = edge.first;
				stack.push(v);
				removeEdge(u, v);
				break;
			}
			if (u == stack.top()) {
				stack.pop();
				res.push_back(u);
			}
		}
		return res;
	}

	int checkBipart(std::vector<char> &marks) const {
		std::queue<int> qu;
		if (vertices_count_ < 2)
			return 0;
//		if(marks.size() != vertices_count_ + 1)
//			marks.resize(vertices_count_ + 1);
		for (unsigned int i = 1; i <= vertices_count_; ++i) {
			marks[i] = 'u';
		}
		for (unsigned int i = 1; i <= vertices_count_; ++i) {
			if (marks[i] == 'u') {
				qu.push(i);
				marks[i] = 'a';
				while (!qu.empty()) {
					const int v = qu.front(); qu.pop();
					for (const auto edge : graph_[v]) {
						const int to = edge.first;
						if (marks[to] == 'u') {
							marks[to] = marks[v] == 'a' ? 'b' : 'a';
							qu.push(to);
						}
						else if (marks[to] == marks[v]) {
							return 0;
						}
					}
				}
			}
		}
		return 1;
	}

	std::vector<bool> getAnyMatching(std::vector<int> & left_part, std::vector<int> & matching) const {
		std::vector<bool> used(vertices_count_ + 1);
		for (auto u : left_part) {
			for (const auto edge : graph_[u]) {
				if (matching[edge.first] == -1) {
					matching[edge.first] = u;
					used[u] = true;
					break;
				}
			}
		}
		return used;
	}

	bool tryKuhn(int v, std::vector<int> & colors, std::vector<int> & matching, int color) const {
		if (colors[v] == color) return false;
		colors[v] = color;
		for (const auto edge : graph_[v]) {
			const int to = edge.first;
			if (matching[to] == -1 || tryKuhn(matching[to], colors, matching, color)) {
				matching[to] = v;
				return true;
			}
		}
		return false;
	}

	std::vector<std::pair<int, int> > getMaximumMatchingBipart() const {
		std::vector<std::pair<int, int> > res;
		std::vector<char> marks(vertices_count_ + 1);
		if (!checkBipart(marks)) {
			return res;
		}
		std::vector<int> left_part;
		for (unsigned int i = 1; i <= vertices_count_; ++i) {
			if (marks[i] == marks[1])
				left_part.push_back(i);
		}

		std::vector<int> colors(vertices_count_ + 1);
		std::vector<int> matching(vertices_count_ + 1, -1);

		std::vector<bool> used = getAnyMatching(left_part, matching);
		int color = 0;
		for (auto u : left_part) {
			if (!used[u])
				tryKuhn(u, colors, matching, ++color);
		}

		for (unsigned int v = 1; v <= vertices_count_; ++v) {
			if (matching[v] != -1) {
				res.emplace_back(matching[v], v);
			}
		}

		return res;
	}

	struct Edge {
		int from, to, capacity, flow;
		Edge(const int from, const int to, const int capacity = 0, const int flow = 0) :from(from), to(to), capacity(capacity), flow(flow) {}
	};



	int dfsDinitz(int v, int min_capacity, int sink, std::vector<int> &distances, std::vector<int> & undeleted_ids, std::vector<std::vector<int>> & g_as_edge_ids, std::vector<Edge> & edges) const {
		if (min_capacity == 0 || v == sink)
			return min_capacity;
		for (; undeleted_ids[v] < static_cast<int>(g_as_edge_ids[v].size()); ++undeleted_ids[v]) {
			const int id = g_as_edge_ids[v][undeleted_ids[v]];
			const int to = edges[id].to;
			if (distances[to] != distances[v] + 1)  continue;
			const int pushed = dfsDinitz(to, std::min(min_capacity, edges[id].capacity - edges[id].flow), sink, distances, undeleted_ids, g_as_edge_ids, edges);
			if (pushed) {
				edges[id].flow += pushed;
				edges[id ^ 1].flow -= pushed;
				return pushed;
			}
		}
		return 0;
	}

	bool bfsDinitz(const int source, const int sink, std::vector<int> &distances, std::vector<std::vector<int>> & g_as_edge_ids, std::vector<Edge> & edges) const {
		for (unsigned int i = 0; i <= vertices_count_; ++i)
			distances[i] = std::numeric_limits<int>::max();
		distances[source] = 0;
		std::queue<int> que;
		que.push(source);
		while (!que.empty() && distances[sink] == std::numeric_limits<int>::max()) {
			const int v = que.front(); que.pop();
			for (const auto edge_id : g_as_edge_ids[v]) {
				int to = edges[edge_id].to;
				if (distances[to] == std::numeric_limits<int>::max() && edges[edge_id].flow < edges[edge_id].capacity) {
					que.push(to);
					distances[to] = distances[v] + 1;
				}
			}
		}
		return distances[sink] != std::numeric_limits<int>::max();
	}

	AdjacencyListGraph* flowDinitz(const int source, const int sink) const {
		AdjacencyListGraph* result = new AdjacencyListGraph(vertices_count_, is_directional_, is_weighted_);

		std::vector<std::vector<int>> g_as_edge_ids(vertices_count_ + 1);
		std::vector<Edge> edges;

		for (unsigned int u = 1; u <= vertices_count_; ++u) {
			for (auto edge : graph_[u]) {
				g_as_edge_ids[u].push_back(static_cast<int>(edges.size()));
				edges.emplace_back(u, edge.first, edge.second);
				g_as_edge_ids[edge.first].push_back(static_cast<int>(edges.size()));
				edges.emplace_back(edge.first, u);
			}
		}
		std::vector<int> undeleted_ids(vertices_count_ + 1);
		std::vector<int> distances(vertices_count_ + 1, std::numeric_limits<int>::max());
		int flow = 0;
		while (bfsDinitz(source, sink, distances, g_as_edge_ids, edges)) {
			for (unsigned int i = 0; i <= vertices_count_; ++i)
				undeleted_ids[i] = 0;
			while (const int pushed = dfsDinitz(source, std::numeric_limits<int>::max(), sink, distances, undeleted_ids, g_as_edge_ids, edges)) {
				flow += pushed;
			}
		}
		
		for(unsigned int i = 0; i < edges.size(); i+=2) {
			result->addEdge(edges[i].from, edges[i].to, edges[i].flow);
		}

		return result;
	}

	int dfsFordFulkerson(int v, int min_capacity, int sink, int color, std::vector<int> &colors, std::vector<std::vector<int>> & g_as_edge_ids, std::vector<Edge> & edges) const {
		if (v == sink)
			return min_capacity;
		colors[v] = color;
		for(auto edge_id: g_as_edge_ids[v]) {
			const int to = edges[edge_id].to;
			if(colors[to] != color && edges[edge_id].flow < edges[edge_id].capacity) {
				const int delta = dfsFordFulkerson(to, std::min(min_capacity, edges[edge_id].capacity - edges[edge_id].flow), sink, color, colors, g_as_edge_ids, edges);
				if(delta > 0) {
					edges[edge_id].flow += delta;
					edges[edge_id^1].flow -= delta;
					return delta;
				}
			}
		}
		return 0;
	}

	AdjacencyListGraph* flowFordFulkerson(const int source, const int sink) const {
		AdjacencyListGraph* result = new AdjacencyListGraph(vertices_count_, is_directional_, is_weighted_);

		std::vector<std::vector<int>> g_as_edge_ids(vertices_count_ + 1);
		std::vector<Edge> edges;

		for (unsigned int u = 1; u <= vertices_count_; ++u) {
			for (auto edge : graph_[u]) {
				g_as_edge_ids[u].push_back(static_cast<int>(edges.size()));
				edges.emplace_back(u, edge.first, edge.second);
				g_as_edge_ids[edge.first].push_back(static_cast<int>(edges.size()));
				edges.emplace_back(edge.first, u);
			}
		}

		int flow = 0;
		std::vector<int> colors(vertices_count_ + 1, 0);
		int color = 0;
		while(const int pushed = dfsFordFulkerson(source, std::numeric_limits<int>::max(), sink, ++color, colors, g_as_edge_ids, edges)) {
			flow += pushed;
		}

		for (unsigned int i = 0; i < edges.size(); i += 2) {
			result->addEdge(edges[i].from, edges[i].to, edges[i].flow);
		}

		return result;
	}
};

class EdgeListGraph : public IGraph {

	std::map<std::pair<int, int>, int> graph_;
	unsigned int vertices_count_ = 0;
	unsigned int is_directional_ = 0;
	unsigned int is_weighted_ = 0;


	void rearange(int &from, int &to) const {
		const int tmp = std::max(from, to);
		from = std::min(from, to);
		to = tmp;
	}
public:

	EdgeListGraph(const unsigned int vectices_count = 0, const unsigned int is_directional = 0, const unsigned int is_weighted = 0) {
		EdgeListGraph::reset(vectices_count, is_directional, is_weighted);

	}

	void readGraph(const std::string &file_name) override {
		reset();
		std::ifstream stream;
		stream.open(file_name);
		{
			char type;
			unsigned int edges_count;
			stream >> type >> vertices_count_ >> edges_count >> is_directional_ >> is_weighted_;
			int from, to, weight = 1;
			for (unsigned int i = 1; i <= edges_count; ++i) {
				stream >> from >> to;
				if (is_weighted_) stream >> weight;
				if (!is_directional_)
					rearange(from, to);
				graph_.emplace(std::make_pair(from, to), weight);
			}
		}
		stream.close();
	}

	void writeGraph(const std::string &file_name) override {
		std::ofstream stream;
		stream.open(file_name);
		{
			stream << "E " << vertices_count_ << " " << graph_.size() << std::endl << is_directional_ << " " << is_weighted_ << std::endl;
			for (const auto edge : graph_) {
				int from = edge.first.first;
				int to = edge.first.second;
				const int weight = edge.second;
				stream << from << " " << to;
				if (is_weighted_) stream << " " << weight;
				stream << std::endl;
			}
		}
		stream.close();
	}

	void addEdge(int from, int to, int weight) override {
		if (!is_directional_)
			rearange(from, to);
		graph_.emplace(std::make_pair(from, to), weight);

	}

	void removeEdge(int from, int to) override {
		if (!is_directional_)
			rearange(from, to);
		const auto edge = graph_.find(std::make_pair(from, to));
		if (edge != graph_.end())
			graph_.erase(edge);
	}

	int changeEdge(int from, int to, const int new_weight) override {
		if (!is_directional_)
			rearange(from, to);
		const auto edge = graph_.find(std::make_pair(from, to));
		const int old_weight = edge->second;
		edge->second = new_weight;
		return old_weight;
	}

	void reset() override {
		graph_.clear();
	}

	void reset(unsigned int vertices_count, unsigned int is_directional, unsigned int is_weighted) override {
		reset();
		vertices_count_ = vertices_count;
		is_directional_ = is_directional;
		is_weighted_ = is_weighted;
	}

	void feelGraph(IGraph *graph) override {
		graph->reset(vertices_count_, is_directional_, is_weighted_);
		for (const auto item : graph_) {
			const int u = item.first.first;
			const int v = item.first.second;
			const int w = item.second;
			graph->addEdge(u, v, w);
		}
	}

	AdjacencyListGraph* getSpaingTreeKruscal() {
		AdjacencyListGraph * result = new AdjacencyListGraph(vertices_count_, 0, is_weighted_);

		typedef std::pair<std::pair<int, int>, int> Edge;
		std::vector<Edge> edges(graph_.begin(), graph_.end());
		std::sort(edges.begin(), edges.end(), [](const Edge & current, const Edge & other) {
			return current.second < other.second;
		});

		Dsu dsu(vertices_count_ + 1);
		for (auto& edge : edges) {
			const int from = edge.first.first;
			const int to = edge.first.second;
			if (dsu.find(to) != dsu.find(from)) {
				dsu.unite(from, to);
				result->addEdge(from, to, edge.second);
			}
		}

		return result;
	}


	AdjacencyListGraph* getSpaingTreeBoruvka() {
		AdjacencyListGraph * result = new AdjacencyListGraph(vertices_count_, 0, is_weighted_);

		typedef std::pair<std::pair<int, int>, int> Edge;
		std::vector<Edge> edges(graph_.begin(), graph_.end());
		std::sort(edges.begin(), edges.end(), [](const Edge & current, const Edge & other) {
			return current.second < other.second;
		});
		Dsu dsu(vertices_count_ + 1);
		std::vector<int> steps(vertices_count_ + 1);
		std::vector<std::pair<int, std::pair<int, int>>> min_vals(vertices_count_ + 1);
		int step = 0;
		bool state_changed = true;
		while (state_changed) {
			state_changed = false;
			++step;

			for (auto& edge : edges) {
				const int from = edge.first.first;
				const int to = edge.first.second;
				const int weight = edge.second;
				const unsigned int parent = dsu.find(from);
				const unsigned int to_parent = dsu.find(to);
				if (to_parent != parent) {
					if (steps[parent] != step) {
						steps[parent] = step;
						min_vals[parent] = std::make_pair(weight, std::make_pair(from, to));
					}
					if (weight < min_vals[parent].first) {
						min_vals[parent] = std::make_pair(weight, std::make_pair(from, to));
					}
					if(!is_directional_) {
						if (steps[to_parent] != step) {
							steps[to_parent] = step;
							min_vals[to_parent] = std::make_pair(weight, std::make_pair(from, to));
						}
						if (weight < min_vals[to_parent].first) {
							min_vals[to_parent] = std::make_pair(weight, std::make_pair(from, to));
						}
					}
				}
			}

			for (unsigned int i = 1; i <= vertices_count_; ++i) {
				if (steps[i] == step) {
					const int from = min_vals[i].second.first;
					const int to = min_vals[i].second.second;
					const int weight = min_vals[i].first;
					if(dsu.find(from) != dsu.find(to)) {
						result->addEdge(from, to, weight);
						dsu.unite(from, to);
						state_changed = true;
					}
				}
			}

		}

		return result;
	}




};

class Graph : public IGraph {
	IGraph *graph = nullptr;

	static char extractType(const std::string &file_name) {
		char type;
		std::ifstream stream;
		stream.open(file_name);
		{
			stream >> type;
		}
		stream.close();
		return type;
	}


public:

	Graph(IGraph * graph) : graph(graph) {

	}

	Graph(const unsigned int vertices_count = 0) {
		graph = new AdjacencyListGraph(vertices_count);
	}

	void transformToAdjList() {
		IGraph * g = new AdjacencyListGraph();
		graph->feelGraph(g);
		graph->reset();
		delete graph;
		graph = g;
	}

	void transformToAdjMatrix() {
		IGraph * g = new AdjacencyMatrixGraph();
		graph->feelGraph(g);
		graph->reset();
		delete graph;
		graph = g;
	}

	void transformToListOfEdges() {
		IGraph * g = new EdgeListGraph();
		graph->feelGraph(g);
		graph->reset();
		delete graph;
		graph = g;
	}

	void writeGraph(const std::string &file_name) override {
		graph->writeGraph(file_name);
	}

	void readGraph(const std::string &file_name) override {
		reset();
		switch (extractType(file_name)) {
		case 'L':
			graph = new AdjacencyListGraph();
			break;
		case 'C':
			graph = new AdjacencyMatrixGraph();
			break;
		case 'E':
			graph = new EdgeListGraph();
			break;
		default:
			break;
		}
		graph->readGraph(file_name);
	}

	void addEdge(const int from, const int to, const int weight) override {
		graph->addEdge(from, to, weight);
	}

	void removeEdge(const int from, const int to) override {
		graph->removeEdge(from, to);
	}

	int changeEdge(const int from, const int to, const int new_weight) override {
		return graph->changeEdge(from, to, new_weight);
	}

	void reset() override {
		if (graph)
			graph->reset();
	}

	void reset(const unsigned int vertices_count, const unsigned int is_directional, const unsigned int is_weighted) override {
		graph->reset(vertices_count, is_directional, is_weighted);
	}

	void feelGraph(IGraph *graph) override {
		graph->feelGraph(graph);
	}

	Graph getSpaingTreePrima() const {
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return Graph(static_cast<IGraph *>(g.getSpaingTreePrima()));
	}

	Graph getSpaingTreeKruscal() const {
		EdgeListGraph g = EdgeListGraph();
		graph->feelGraph(&g);
		return Graph(static_cast<IGraph*>(g.getSpaingTreeKruscal()));
	}

	Graph getSpaingTreeBoruvka() const {
		EdgeListGraph g = EdgeListGraph();
		graph->feelGraph(&g);
		return Graph(static_cast<IGraph*>(g.getSpaingTreeBoruvka()));
	}

	int checkEuler(bool &circle_exist) const {
		circle_exist = false;
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);

		std::vector<int> degrees = g.getVerticesDegrees();
		std::vector<int> components_sizes = g.getComponentsSizes();

		int big_components = 0;
		int odd_number = 0;
		int result = 1;
		for (unsigned int i = 1; i < degrees.size(); ++i) {
			if (degrees[i] % 2) {
				result = i;
				++odd_number;
            		}
			big_components += components_sizes[i] > 1;
		}
		if (odd_number > 2 || big_components > 1)
			return 0;
		circle_exist = odd_number == 0;


		return result;
	}

	std::vector<int> getEuleranTourFleri() const {
		bool has_cycle = false;
		const int v = checkEuler(has_cycle);
		if (!v)
			return std::vector<int>();
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return g.getEuleranTourFleri(v);
	}


	std::vector<int> getEuleranTourEffective() const {
		bool has_cycle = false;
		const int v = checkEuler(has_cycle);
		if (!v)
			return std::vector<int>();
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return g.getEuleranTourEffective(v);
	}

	int checkBipart(std::vector<char> &marks) const {
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return g.checkBipart(marks);
	}

	std::vector<std::pair<int, int> > getMaximumMatchingBipart() const {
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return g.getMaximumMatchingBipart();
	}


	Graph flowFordFulkerson(int source, int sink) const {
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return Graph(static_cast<IGraph *>(g.flowFordFulkerson(source, sink)));
	}

	Graph flowDinitz(int source, int sink) const {
		AdjacencyListGraph g = AdjacencyListGraph();
		graph->feelGraph(&g);
		return Graph(static_cast<IGraph *>(g.flowDinitz(source, sink)));
	}
};


#endif //GRAPH_H
