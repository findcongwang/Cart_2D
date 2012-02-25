#ifndef POLYHEDRON_H
#define POLYHEDRON_H

#include <iostream>
#include <vector>
#include "basic_types_2d.h"
#include "ADT_2d.h"

class Polyhedron{
  public:
  	bool generic;
  	bool is_cut;
	pt4D *bounding_box;
	Polyhedron();
};

class Graph: public Polyhedron{
  public:
	int* list_of_verts;	//vertID to the vector of vert_precise
	int* list_of_edges;	//edgeID to the vector of edges
	int* num_entries;	//pointer to int
	int** num_indices;	//list of pointers to int
	//virtual void info() = 0;
};

class Face: public Graph{	//this is a new class
	Polyhedron *parent;	//the polyhedron that this face belongs to
	Face(Polyhedron *parent, int *vert_lst, int *edge_lst );
};

class CartCell: public Graph{
  public:
	double length;
	int pos;	//vertID to the vector of vert_ints //pos of b_l
	std::vector<Polyhedron*> *intersect_candidates;
	bool painted;
	bool in_geo;
  	edge *top;
	edge *bot;
	edge *left;
	edge *right;
	vert_precise *t_l;
	vert_precise *t_r;
	vert_precise *b_l;
	vert_precise *b_r;
	CartCell(pt4D *bound, double l, int coord, bool init);
	~CartCell();
};

class Segment: public Graph{	//triangle in 3d
  public:
	Segment(pt4D *bound, int *vert_lst, int *edge_lst);
	~Segment();
};

class CutCell: public Graph{
  public:
	Face **list_of_faces;
	
	CutCell(int *vert_lst, int *edge_lst, int* entries, int** indices);
	CutCell(int *vert_lst, int *edge_lst, int* entries, int** indices, Face** face_list);
	~CutCell();
};

class cell_piece: public CutCell{
  public:
	std::vector<vert_intersection*> *list_of_intersecton_points;
	std::vector<Polyhedron*> *intersect_segments;
	std::vector<int> *old_vert_bound;
	std::vector<int> *old_edge_bound;
	
	cell_piece(int *vert_lst, int *edge_lst, int** indices,
			   std::vector<vert_intersection*> *list_of_intersecton_points,
			   std::vector<Polyhedron*> *intersect_segments);
};

//make destructor for cell_piece

std::ostream& operator<< (std::ostream &out, CartCell &cart_cell);
std::ostream& operator<< (std::ostream &out, Segment &segment);
std::ostream& operator<< (std::ostream &out, cell_piece &piece);



#endif
