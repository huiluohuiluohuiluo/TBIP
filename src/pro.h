#ifndef PRO_H_
#define PRO_H_
class Pro{
	int edgeId_head;
	int edgeId_end;
	double pr;

public:
	Pro(int edgeId_head, int edgeId_end, double pr);
};

Pro::Pro(int edgeId_head, int edgeId_end, double pr){
	this->edgeId_head = edgeId_head;
	this->edgeId_end = edgeId_end;
	this->pr = pr;

}

#endif /* PRO_H_ */
