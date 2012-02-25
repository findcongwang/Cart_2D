/*
 * Polyhedron.cc
 *
 *  Created on: Sep 7, 2011
 *      Author: c233wang
 */

#include "Polyhedron_2d.h"

Polyhedron::Polyhedron(){
	generic = true;
	is_cut = false;
}

CartCell::CartCell(pt4D *bound, double l, int coord, bool init_in){
	//a Cartesian cell is described by its low-valued vertex and length
	bounding_box = bound;
	length = l;
	pos = coord;	//vert_ID of a vertex in the vert_int vector
	painted = false;
	in_geo = init_in;
	intersect_candidates = NULL;	//initializes to null
	num_entries = NULL;
	num_indices = NULL;
	generic = false;
	top = NULL;
	bot = NULL;
	left = NULL;
	right = NULL;
	t_l = NULL;
	t_r = NULL;
	b_l = NULL;
	b_r = NULL;
}

CartCell::~CartCell(){
	//delete used memory
	delete intersect_candidates;
}


Segment::Segment(pt4D *bound, int *vert_lst, int *edge_lst){
	//described by indexes to vertices and edges involved
	bounding_box = bound;
	list_of_verts = vert_lst;
	list_of_edges = edge_lst;
}


Segment::~Segment(){
	//delete used memory
	delete[] list_of_verts;
	delete[] list_of_edges;
}

Face::Face(Polyhedron *poly, int *vert_lst, int *edge_lst ){
	parent = poly;
	list_of_verts = vert_lst;
	list_of_edges = edge_lst;
}

CutCell::CutCell(int *vert_lst, int *edge_lst, int* entries, int** indices){
	//described by indexes to vertices and edges involved
	//bounding_box no longer needed
	num_entries = entries;
	num_indices = indices;
	list_of_verts = vert_lst;
	list_of_edges = edge_lst;
	list_of_faces = NULL;
	is_cut = true;
	generic = false;
}

CutCell::CutCell(int *vert_lst, int *edge_lst, int* entries, int** indices, Face** face_list){
	//described by indexes to vertices and edges involved
	//bounding_box no longer needed
	num_entries = entries;
	num_indices = indices;
	list_of_verts = vert_lst;
	list_of_edges = edge_lst;
	list_of_faces = face_list;
	is_cut = true;
	generic = false;
}

CutCell::~CutCell(){
	//delete used memory
	delete[] list_of_verts;
	delete[] list_of_edges;
	
	if(list_of_faces != NULL){
		delete[] list_of_faces;
	}
}

cell_piece::cell_piece(int *vert_lst, int *edge_lst, int** indices,
			   		   std::vector<vert_intersection*> *list_pts,
			   		   std::vector<Polyhedron*> *list_segs)
		   :CutCell(vert_lst, edge_lst, NULL, indices){
	list_of_intersecton_points = list_pts;
	intersect_segments = list_segs;

	old_vert_bound = NULL;
	old_edge_bound = NULL;
	//make destructor for cell_piece
}

std::ostream& operator<< (std::ostream &out, CartCell &cart_cell){
	//output cell fields
	out << "<Cell> nodeID: " << cart_cell.pos << std::endl <<
			" Length: " << cart_cell.length <<
			" Box: " << *cart_cell.bounding_box << std::endl;
	return out;
}

std::ostream& operator<< (std::ostream &out, Segment &segment){
	//output IDs of vertices and edges
	out << *segment.bounding_box << std::endl;
	for(int i = 0; i < 2; i++){
		out <<
			"Vert " << i << ": " << segment.list_of_verts[i] << 
			" Edge " << i << ": " << segment.list_of_edges[i] << std::endl;
	}
	return out;
}

std::ostream& operator<< (std::ostream &out, cell_piece &piece){
	for(int i = 0; i < *(piece.num_indices[0]); i++){
		out << piece.list_of_verts[i] << std::endl;
		out << piece.list_of_edges[i] << std::endl;
	}
	for(int j = 0; j < piece.list_of_intersecton_points->size(); j++){
		out << piece.list_of_intersecton_points->at(j) << std::endl;
	}
	for(int j = 0; j < piece.intersect_segments->size(); j++){
		out << piece.intersect_segments->at(j) << std::endl;
	}
}


/*
std::ostream& operator<< (std::ostream &out, Triangle &triangle){
	//output IDs of vertices and edges
	for(int i = 0; i < 3; i++){
		out << *triangle.bounding_box << 
			"Vert " << i << ": " << triangle.list_of_verts[i] << 
			" Edge " << i << ": " << triangle.list_of_edges[i] << std::endl;
	}
	return out;
}

/*
std::ostream& operator<< (std::ostream &out, BodyCell &body_cell){
	//output IDs of vertices and edges
	for(int i = 0; i < body_cell.num_verts; i++){
			out << body_cell.list_of_verts[i] << " " <<
					body_cell.list_of_edges[i] << std::endl;
	}

	for(int i = 0; i < body_cell.num_edges; i++){
			out << body_cell.list_of_verts[i] << " " <<
					body_cell.list_of_edges[i] << std::endl;
	}

	return out;
}
*/