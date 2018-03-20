#ifndef GRAPH_H
// ReSharper disable once CppInconsistentNaming
#define _CRT_SECURE_NO_WARNINGS
#define GRAPH_H

#include <string>
#include <vector>
#include <algorithm>
#include <map>
#include <tuple>
#include <set>
#include <fstream>

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
	unsigned int find(unsigned int x) {
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
		for (unsigned int u = 0; u < vertices_count_; ++u)
			for (unsigned int v = 0; v < vertices_count_; ++v)
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
			graph_.resize(vertices_count_);
			for (unsigned int i = 0; i < vertices_count_; ++i)
				graph_[i].resize(vertices_count_);
			for (unsigned int from = 0; from < vertices_count_; ++from)
				for (unsigned int to = 0; to < vertices_count_; ++to)
					stream >> graph_[from][to];
		}
		stream.close();
	}


	void writeGraph(const std::string &file_name) override {
		std::ofstream stream;
		stream.open(file_name);
		{
			stream << "C " << vertices_count_ << std::endl << is_directional_ << " " << is_weighted_ << std::endl;
			for (unsigned int from = 0; from < vertices_count_; ++from) {
				for (unsigned int to = 0; to < vertices_count_; ++to) {
					stream << graph_[from][to];
					if (to + 1 == vertices_count_) stream << std::endl;
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
		for (unsigned int i = 0; i < vertices_count_; ++i) graph_[i].clear();
		graph_.clear();
		vertices_count_ = 0;
	}

	void reset(const unsigned int vertices_count, const unsigned int is_directional, const unsigned int is_weighted) override {
		reset();
		vertices_count_ = vertices_count;
		is_directional_ = is_directional;
		is_weighted_ = is_weighted;
		graph_.resize(vertices_count);
		for (unsigned int i = 0; i < vertices_count; ++i) graph_[i].resize(vertices_count);
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
			graph_.resize(vertices_count_);
			int weight = 1;
			int to = 0;
			std::string line;
			for (unsigned int from = 0; from < vertices_count_; ++from) {
				std::getline(stream, line);
				std::vector<int> numbers = extractIntSequence(line);
				for (unsigned int idx = 0; idx < numbers.size(); ++idx) {
					to = numbers[idx] - 1;
					if (is_weighted_) weight = numbers[++idx];
					graph_[from][to] = weight;
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
			for (unsigned int i = 0; i < vertices_count_; ++i) {
				bool is_first = true;
				for (const auto to : graph_[i]) {
					if (!is_first) stream << " ";
					stream << to.first + 1;
					if (is_weighted_) stream << to.second;
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
		for (unsigned int i = 0; i < vertices_count_; ++i) graph_[i].clear();
		graph_.clear();
		vertices_count_ = 0;
	}

	void reset(const unsigned int vertices_count, const unsigned int is_directional, const unsigned int is_weighted) override {
		reset();
		vertices_count_ = vertices_count;
		is_directional_ = is_directional;
		is_weighted_ = is_weighted;
		graph_.resize(vertices_count);
	}


	void feelGraph(IGraph *graph) override {
		graph->reset(vertices_count_, is_directional_, is_weighted_);
		for (unsigned int u = 0; u < vertices_count_; ++u)
			for (const auto v : graph_[u]) {
				graph->addEdge(u, v.first, v.second);
				if (!is_directional_)
					graph->addEdge(v.first, u, v.second);
			}
	}


	AdjacencyListGraph* getSpaingTreePrima() {
		AdjacencyListGraph* result = new AdjacencyListGraph(vertices_count_, 0, is_weighted_);

		std::set<std::pair<int, int>> prior;
		std::vector<int> dists(vertices_count_, std::numeric_limits<int>::max());
		std::vector<int> ends(vertices_count_, -1);
		dists[0] = 0;
		prior.emplace(0, 0);
		for (unsigned int i = 0; i < vertices_count_; ++i) {
			if (prior.empty())
				exit(1); // graph has few cc
			const int v = prior.begin()->second;
			prior.erase(prior.begin());
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

		for (unsigned int i = 0; i < vertices_count_; ++i)
			if (ends[i] != -1)
				result->addEdge(i, ends[i], dists[i]);
		return result;
	}

};

class EdgeListGraph : public IGraph {

	std::vector<std::tuple<int, int, int>> graph_;
	unsigned int vertices_count_ = 0;
	unsigned int is_directional_ = 0;
	unsigned int is_weighted_ = 0;
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
			graph_.reserve(edges_count);
			for (unsigned int i = 0; i < edges_count; ++i) {
				stream >> from >> to;
				if (is_weighted_) stream >> weight;
				graph_.emplace_back(from - 1, to - 1, weight);
			}		
		}
		stream.close();
	}

	void writeGraph(const std::string &file_name) override {
		std::ofstream stream;
		stream.open(file_name);
		{
			stream << "E " << vertices_count_ << " " << graph_.size() << std::endl << is_directional_ << " " << is_weighted_ << std::endl;
			int from, to, weight;
			for (const auto edge : graph_) {
				std::tie(from, to, weight) = edge;
				++from, ++to;
				stream << from << " " << to;
				if (is_weighted_) stream << " " << weight;
				stream << std::endl;
			}
		}
		stream.close();
	}

	void addEdge(int from, int to, int weight) override {
		graph_.emplace_back(from, to, weight);

	}

	void removeEdge(int from, int to) override {
		const auto edge = std::find_if(graph_.begin(), graph_.end(), [from, to](std::tuple<int, int, int> &edge) {
			return (std::get<0>(edge) == from) && (std::get<1>(edge) == to);
		});
		graph_.erase(edge);
	}

	int changeEdge(int from, int to, int new_weight) override {
		auto edge = std::find_if(graph_.begin(), graph_.end(), [from, to](std::tuple<int, int, int> &edge) {
			return (std::get<0>(edge) == from) && (std::get<1>(edge) == to);
		});
		std::swap(std::get<2>(*edge), new_weight);
		return new_weight;
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
		int u, v, w;
		for (auto &item : graph_) {
			std::tie(u, v, w) = item;
			graph->addEdge(u, v, w);
		}
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
		graph->addEdge(from - 1, to - 1, weight);
	}

	void removeEdge(const int from, const int to) override {
		graph->removeEdge(from - 1, to - 1);
	}

	int changeEdge(const int from, const int to, const int new_weight) override {
		return graph->changeEdge(from - 1, to - 1, new_weight);
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

	Graph getSpaingTreePrima() {
		transformToAdjList();
		AdjacencyListGraph * g = dynamic_cast<AdjacencyListGraph*>(graph);
		return Graph(static_cast<IGraph *>(g->getSpaingTreePrima()));
	}
	//    Graph getSpaingTreeKruscal() {
	//
	//    }
	//    Graph getSpaingTreeBoruvka() {
	//
	//    }


};


#endif //GRAPH_H