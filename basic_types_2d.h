#ifndef BASIC_TYPES_H
#define BASIC_TYPES_H

#include <iostream>
#include <vector>

enum vert_cond {OUT, IN, ON_BOUND, COND_COUNT};
class ADTNode;

class pt4D{	//4D point class to represent boxes
  public:
	double x_min;
	double x_max;
	double y_min;
	double y_max;
	pt4D(double x1, double x2, double y1, double y2);
};

class vert_precise{
  public:
	double x;
	double y;
	int ID;
	vert_precise(double x_in, double y_in, int ID_in);
};

class edge{
  public:
	int pt1; //nodeID 1
	int pt2; //nodeID 2
	bool walked;	// NO: 0, YES: 1
	
	//edge connectives
	ADTNode *node1;
	ADTNode *node2;
	edge* child1;
	edge* child2;
	edge(int in1, int in2);
};

class n_gon{	//n-sided polygon defined by a vector of edges
  public:
	int n;
	std::vector<edge*> *edges;
	n_gon(int n);
};

class Refinement{
  public:
  	std::vector<std::vector<pt4D*>*>  *table;	//table of pt4D to mark refinement regions
  	Refinement(std::vector<std::vector<pt4D*>*> *refinement_table);
};

class Minsize_Refinement{
  public:
  	pt4D *box;
  	double min_size;
  	Minsize_Refinement(pt4D *box_in, double min_size_in);
  	~Minsize_Refinement();
};

class Adaptive_variable{
  public:
  	int max_refs;	//Maxmium level of refinements allowed
  	double alpha;	//Angle parameter in which determines if a cell need to be divided (radians)
  	Adaptive_variable(int max_refs_in, double alpha_in);
};

class size_box{
  public:
  	double min;
  	double max;
  	pt4D *box;	//affected box
  	size_box(double min_in, double max_in, pt4D *box_in);
};

class body_edge;

class vert_intersection: public vert_precise{
  public:
	bool on_cell_bound;	//flag for if the vert is valid to use for in/out determination
	vert_cond status;	// OUT: 0, IN: 1, ON_BOUND: 2
	edge* host;
	vert_intersection(double x_in, double y_in, edge* host_in, vert_cond status_in, bool on_cellbound, int ID_in);
};

class body_edge: public edge{
  public:
	bool is_cart;	// NO: 0, YES: 1
	std::vector<vert_intersection*> *list_of_intersecton_points;
	body_edge(int in1, int in2, std::vector<vert_intersection*> *list_in, bool flag);
};

std::ostream& operator<< (std::ostream &out, pt4D &pt);
std::ostream& operator<< (std::ostream &out, vert_precise &vert);
std::ostream& operator<< (std::ostream &out, edge &edge);
std::ostream& operator<< (std::ostream &out, vert_intersection &vert);
std::ostream& operator<< (std::ostream &out, body_edge &edge);

#endif
