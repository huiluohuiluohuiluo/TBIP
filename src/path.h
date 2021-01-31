#ifndef PATH_H_
#define PATH_H_
#include <vector>
#include <list>
using namespace std;

class Path {

public:
	double dis;
	list<int> edges;

	Path(double dis, list<int> edges);


	bool operator <(const struct Path & path) const {
		list<int> edges1 = this->edges;
		list<int> edges2 = path.edges;

		list<int>::iterator iter1 = edges1.begin();
		list<int>::iterator iter2 = edges2.begin();
		while (iter1 != edges1.end() && iter2 != edges2.end()) {
			if (*iter1 == *iter2) {
				++iter1;
				++iter2;
			} else {
//				cout << "not equal" << endl;
				return *iter1 < *iter2;
			}
		}
//		cout << "equal" << endl;
		return false;
	}
};

#endif /* PATH_H_ */
