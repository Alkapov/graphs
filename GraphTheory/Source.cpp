#include <iostream>

#include "Graph.h"

#if _MSC_VER
const std::string ADJACENCY_MATRIX = "samples/adjacency_matrix.txt";
const std::string ADJACENCY_LIST = "samples/adjacency_list.txt";
const std::string EDGES_LIST = "samples/edges_list.txt";
const std::string TEMP_FILE = "samples/temp.txt";
#else
const std::string ADJACENCY_MATRIX = "../samples/adjacency_matrix.txt";
const std::string ADJACENCY_LIST = "../samples/adjacency_list.txt";
const std::string EDGES_LIST = "../samples/edges_list.txt";
const std::string TEMP_FILE = "../samples/temp.txt";
#endif

int main() {
	Graph graph;
	graph.readGraph(EDGES_LIST);
	bool check;
	int s = graph.checkEuler(check);
	std::cout << check << " " << s << std::endl;
	//Graph tree = graph.getSpaingTreeBoruvka();
	//			
	//tree.transformToListOfEdges();

	////    graph.addEdge(1, 1, 1000);
	////    printf("%d", graph.changeEdge(1, 2, 12200));
	//tree.writeGraph(TEMP_FILE);

	return 0;
}