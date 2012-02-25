#include "basic_types_2d.h"

pt4D::pt4D(double x1, double x2, double y1, double y2){
	x_min = x1;
	x_max = x2;
	y_min = y1;
	y_max = y2;
}

n_gon::n_gon(int n_in){
	n = n_in;
}

vert_precise::vert_precise(double x_in, double y_in, int ID_in){
	x = x_in;
	y = y_in;
	ID = ID_in;
}

edge::edge(int in1, int in2){
	pt1 = in1;
	pt2 = in2;
	walked = false;
	node1 = NULL;
	node2 = NULL;
	child1 = NULL;
	child2 = NULL;
}

Refinement::Refinement(std::vector<std::vector<pt4D*>*> *refinement_table){
	table = refinement_table;
}

Minsize_Refinement::Minsize_Refinement(pt4D *box_in, double min_size_in){
	box = box_in;
	min_size = min_size_in;
}

Minsize_Refinement::~Minsize_Refinement(){
	delete box;
}

Adaptive_variable::Adaptive_variable(int max_refs_in, double alpha_in){
	max_refs = max_refs_in;
	alpha = alpha_in;
}

size_box::size_box(double min_in, double max_in, pt4D *box_in){
	min = min_in;
	max = max_in;
	box = box_in;
}

vert_intersection::vert_intersection(double x_in, double y_in, edge* host_in, vert_cond status_in, bool on_cellbound, int ID_in)
				  :vert_precise(x_in, y_in, ID_in){
	on_cell_bound = on_cellbound;
	host = host_in;
	status = status_in;
}

body_edge::body_edge(int in1, int in2, std::vector<vert_intersection*> *list_in, bool flag)
		  :edge(in1, in2){
	list_of_intersecton_points = list_in;
	is_cart = flag;
}

std::ostream& operator<< (std::ostream &out, pt4D &pt){
	//prints "coordinates in 4D: (x_min, ..., z_max)"
	out << "(" << pt.x_min << ", " << pt.x_max << ", "
	    << pt.y_min << ", " << pt.y_max << ")" << std::endl;
}

std::ostream& operator<< (std::ostream &out, vert_precise &vert){
	//prints "[x, y]"
	out << "[" << vert.x << ", " << vert.y << "]";
	return out;
}

std::ostream& operator<< (std::ostream &out, edge &edge){
	//prints "EDGE(#-#)"
	out << "EDGE(" << edge.pt1 << "-" << edge.pt2 << ")" << std::endl;
	return out;
}

std::ostream& operator<< (std::ostream &out, vert_intersection &vert){
	out << "##vert_intersection##" << std::endl;
	out << "[" << vert.x << ", " << vert.y << "]" << std::endl;
	out << "on_cell_bound: " << vert.on_cell_bound << std::endl;
	out << "status: " << vert.status << std::endl;
	if(vert.host != NULL){
		out << "host: " << *(vert.host) << std::endl;
	}
	else{
		out << "host: NULL" << std:: endl;
	}
}

std::ostream& operator<< (std::ostream &out, body_edge &edge){
	out << "BODY_EDGE(" << edge.pt1 << "-" << edge.pt2 << ")" << std::endl;
	out << "is_cart: " << edge.is_cart << std::endl;
	for(int j = 0; j < edge.list_of_intersecton_points->size(); j++){
		out << edge.list_of_intersecton_points->at(j) << std::endl;
	}
	
}