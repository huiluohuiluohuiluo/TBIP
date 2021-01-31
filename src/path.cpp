#include <iostream>

#include "path.h"
using namespace std;

Path::Path(double dis, list<int> edges){
	this->dis = dis;
	this->edges = edges;
}
