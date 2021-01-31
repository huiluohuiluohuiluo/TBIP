#ifndef PATHSTR_H_
#define PATHSTR_H_
#include <vector>
#include <set>
using namespace std;

//class Path_Str {
//
//public:
//	double dis;
//	set<int> edges;
//	int endEdge;
//
//	Path_Str(double dis, set<int> edges, int endEdge);
//
//	bool operator <(const struct Path_Str & path) const {
//		set<int> edges1 = this->edges;
//		set<int> edges2 = path.edges;
//
//		set<int>::iterator iter1 = edges1.begin();
//		set<int>::iterator iter2 = edges2.begin();
//		while (iter1 != edges1.end() && iter2 != edges2.end()) {
//			if (*iter1 == *iter2) {
//				++iter1;
//				++iter2;
//			} else {
////				cout << "not equal" << endl;
//				return *iter1 < *iter2;
//			}
//		}
////		cout << "equal" << endl;
//		return false;
//	}
//};

struct Path_Str {
	double dis;
	set<int> edges;
	int endEdge;
};

#endif /* PATHSTR_H_ */
