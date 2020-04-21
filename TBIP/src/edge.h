#ifndef EDGE_H_
#define EDGE_H_
class Edge{
public:
	int edgeId;
	int sourceNode;
	int desNode;
	double dis;
	Edge(int edgeId, int sourceNode, int desNode, double dis);
};

#endif /* EDGE_H_ */
