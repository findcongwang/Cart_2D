#ifndef ADT_H
#define ADT_H

class Polyhedron;
class pt4D;

class ADTNode{
  public:	
	Polyhedron* poly;
	pt4D *region;
	int depth; //maybe short_int
	ADTNode* L_branch;
	ADTNode* R_branch;
	
	ADTNode(Polyhedron* poly_in, int depth_in, pt4D *region_in);
};

#endif