/*
 *	Cart.cpp
 *
 *  Created on: Sep 7, 2011
 *      Author: c233wang (Cong Wang)
 *
 *	For usage, please refer to the README file
 *	To make the code more readable, I recommend collapsing all function blocks first
 *
 */
#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <string>
#include <queue>
#include <cstdio>
#include <cassert>
#include <cmath>
#include <climits>

enum gene_type {MANUAL, MIN_SIZE, ADAPTIVE, TYPE_COUNT};
double *GBOX_X_MIN = new double;
double *GBOX_X_MAX = new double;
double *GBOX_Y_MIN = new double;
double *GBOX_Y_MAX = new double;
double *INIT_SIZE = new double;
double *INIT_MIN = new double;
double *EPSILON = new double;
bool *INTERNAL = new bool;
gene_type *PROGRAM_GENE_TYPE = new gene_type;

#include "basic_types_2d.h"
#include "ADT_2d.h"
#include "Polyhedron_2d.h"

using namespace std;

//variables & constants:
double const PI = 4 * atan(1);
const int DIM = 2;
vector<vert_precise*> precise_verts;	//stores vertex coordinates
vector<vert_precise*> precise_verts_backup;
vector<edge*> precise_edges;	//edges in precise_verts
vector<int> print_vertIDs;	//verts of leaves
vector<edge*> print_edges;	//edges of leaves
vector<edge*> print_boundary;
vector<n_gon*> print_cutcells; //cutcells to print
int geo_num = 0;
vector<Segment*> segments;	//input triangles
vector<ADTNode*> non_cut_cell_nodes;
vector<ADTNode*> non_cut_cell_nodes_backup;
vector<ADTNode*> cut_cell_nodes;	//add cutcell pointers to here as we go
vector<ADTNode*> cut_cell_nodes_backup;
ADTNode *segments_ADT = NULL;
ADTNode *final_mesh = NULL;
const int boundary_cap = 3;	// 2 is simply closing off
int CARTCELL_NUM_ENTRIES = 1;
int CARTCELL_NUM_INDICES = 4;
vert_precise *mesh_coord;

//functions:

//returns the base 2 log of num
double log_base2(double num){
	return ( log(num)/log(2.0) );
}

//returns the length of the smallest child edge from an edge
double deepest_length(edge *eg){
	
	if( (eg->child1 == NULL) && (eg->child2 == NULL) ){
		//base case
		CartCell *c1 = static_cast<CartCell*>(eg->node1->poly);
		CartCell *c2 = static_cast<CartCell*>(eg->node2->poly);
		
		if(eg->node1 == NULL){
			return c2->length;
		}
		if(eg->node1 == NULL){
			return c1->length;
		}
		
		return min(c1->length, c2->length);
	}
	
	if( (eg->child1 != NULL) && (eg->child2 != NULL) ){
		//both has cells
		return min( deepest_length(eg->child1), deepest_length(eg->child2) );
	}
	
	if(eg->child1 != NULL){
		return deepest_length(eg->child1);
	}
	
	return deepest_length(eg->child2);
}

//returns the level difference of an edge and a cell length
int level_difference(edge *eg, double length){
	if(eg == NULL){
		return 0;
	}
	double other_length = deepest_length(eg);
	return log_base2(length/other_length);
}

//returns the depth of an edge
int edge_depth(edge *eg){
	if(eg == NULL || (eg->child1 == NULL && eg->child2 == NULL)){
		return 0;
	}
	return max(1 + edge_depth(eg->child1), 1 + edge_depth(eg->child2));
}

//returns the node pointer equivalent to self in a edge 
ADTNode *find_self(edge *eg, ADTNode *self){
	if(eg->node1 == self){
		return eg->node1;
	}
	return eg->node2;
}

//returns the other node pointer in a edge 
ADTNode *find_other(edge *eg, ADTNode *self){
	if(eg->node1 == self){
		return eg->node2;
	}
	return eg->node1;
}

//returns the address of the pointer equivalent to self (to find and set all copies of self to null)
ADTNode **find_other_ref(edge *eg, ADTNode *self){
	if(eg->node1 == self){
		return &eg->node2;
	}
	return &eg->node1;
}

//return the node pointing the neighbouring edges of a cell and one of its edge
ADTNode *find_other(edge *eg, CartCell *self_cell){
	if( (eg->node1 != NULL) && (static_cast<CartCell*>(eg->node1->poly) == self_cell) ){
		return eg->node2;
	}
	return eg->node1;
}

//refines a node (cell) by dividing it in 4, the is_dummy flags if the cell is from the initial mesh
void refine(ADTNode *node, bool is_dummy){
	
	//determine depth of the new subcells
	int new_depth;
	if(is_dummy){
		new_depth = 0;
	}
	else{
		new_depth = node->depth + 1;
	}

	CartCell *cell = static_cast<CartCell*>(node->poly);
	pt4D *current_bound = cell->bounding_box;
	double med = cell->length / 2;
	
	////cout << "node parameters initalized" << endl;
	
	////cout << "dividing in 2" << endl;
	
	//<--- divide cell in 2
	pt4D *l_bound = new pt4D(current_bound->x_min, current_bound->x_max - med,
							 current_bound->y_min, current_bound->y_max);
	pt4D *r_bound = new pt4D(current_bound->x_min + med, current_bound->x_max,
							 current_bound->y_min, current_bound->y_max);
	Polyhedron *rect1 = new Polyhedron();
	Polyhedron *rect2 = new Polyhedron();
	rect1->bounding_box = l_bound;
	rect2->bounding_box = r_bound;
	node->L_branch = new ADTNode(rect1, 0, l_bound);
	node->R_branch = new ADTNode(rect2, 0, r_bound);
	// end of divide cell in 2 --->
	
	//cout << "dividing in 4" << endl;
	
	pt4D *b_c1 = new pt4D(l_bound->x_min, l_bound->x_max,
						  l_bound->y_min, l_bound->y_max - med);
	pt4D *b_c2 = new pt4D(l_bound->x_min, l_bound->x_max,
						  l_bound->y_min + med, l_bound->y_max);
	pt4D *b_c3 = new pt4D(r_bound->x_min, r_bound->x_max,
						  r_bound->y_min, r_bound->y_max - med);
	pt4D *b_c4 = new pt4D(r_bound->x_min, r_bound->x_max,
						  r_bound->y_min + med, r_bound->y_max);
			
	vert_precise *v2 = new vert_precise(precise_verts.at(cell->pos)->x,
										precise_verts.at(cell->pos)->y + med, precise_verts.size());
	precise_verts.push_back(v2);							
	
	vert_precise *v3 = new vert_precise(precise_verts.at(cell->pos)->x + med,
										precise_verts.at(cell->pos)->y, precise_verts.size());
	precise_verts.push_back(v3);									
										
	vert_precise *v4 = new vert_precise(precise_verts.at(cell->pos)->x + med,
										precise_verts.at(cell->pos)->y + med, precise_verts.size());
	precise_verts.push_back(v4);									
	
	vert_precise *v_top = new vert_precise(precise_verts.at(cell->pos)->x + med,
										   precise_verts.at(cell->pos)->y + cell->length, precise_verts.size());
	precise_verts.push_back(v_top);
	
	vert_precise *v_right = new vert_precise(precise_verts.at(cell->pos)->x + cell->length,
										     precise_verts.at(cell->pos)->y + med, precise_verts.size());
	precise_verts.push_back(v_right);	
	
	//cout << "new verts were made" << endl;
	
	//cout << *b_c1 << endl;
	//cout << med << endl;
	//cout << cell->b_l->ID << endl;
	
	//CartCell(pt4D *bound, double l, int coord);
	CartCell *c1 = new CartCell(b_c1, med, cell->b_l->ID, *INTERNAL);
	CartCell *c2 = new CartCell(b_c2, med, v2->ID, *INTERNAL);	
	CartCell *c3 = new CartCell(b_c3, med, v3->ID, *INTERNAL);
	CartCell *c4 = new CartCell(b_c4, med, v4->ID, *INTERNAL);
	
	//new cells in this formation
	//. - . - .
	//| 2 | 4 |
	//. _ . _ .
	//| 1 | 3 |
	//. - . - .
	
	//cout << "new cartcells were made" << endl;
	
	//ADTNode(Polyhedron* poly_in, int depth_in, pt4D *region_in);
	ADTNode *n1 = new ADTNode(c1, new_depth, b_c1);
	ADTNode *n2 = new ADTNode(c2, new_depth, b_c2);
	ADTNode *n3 = new ADTNode(c3, new_depth, b_c3);
	ADTNode *n4 = new ADTNode(c4, new_depth, b_c4);

	//cout << "new nodes were made" << endl;
	
	node->L_branch->L_branch = n1;
	node->L_branch->R_branch = n2;
	node->R_branch->L_branch = n3;
	node->R_branch->R_branch = n4;
	
	//cout << "initializing edge data" << endl;
	//linking up all edge relations
	
	edge *top_left;
	edge *top_right;
	edge *left_up;
	edge *mid_up;
	edge *right_up;
	edge *mid_left;
	edge *mid_right;
	edge *left_down;
	edge *mid_down;
	edge *right_down;
	edge *bot_left;
	edge *bot_right;
	
	//cout << *cell->bounding_box;
	//cout << "child1 = null?: " << (cell->top->child1 == NULL) << endl;
	//cout << "child2 = null?: " << (cell->top->child2 == NULL) << endl;
	//cout << "node1 = null?: " << (cell->top->node1 == NULL) << endl;
	//cout << "node2 = null?: " << (cell->top->node2 == NULL) << endl;
	
	if(cell->top->child1 != NULL){
		
		top_left = cell->top->child1;
		if(top_left->node1->region->y_min < node->region->y_max){
			top_left->node1 = n2;
		}
		else{
			top_left->node2 = n2;
		}
	}
	else{
		top_left = new edge(cell->t_l->ID, v_top->ID);
		cell->top->child1 = top_left;
		
		//cout << "node1 = node?: " << (cell->top->node1 == node) << endl;
		//cout << "node2 = node?: " << (cell->top->node2 == node) << endl;
		
		if(cell->top->node1 != NULL && cell->top->node1->region->y_min < node->region->y_max){
			top_left->node1 = cell->top->node2;
		}
		else{
			top_left->node1 = cell->top->node1;
		}
		top_left->node2 = n2;
	}
	
	//cout << "done child1" << endl;
	
	if(cell->top->child2 != NULL){
		top_right = cell->top->child2;
		if(top_right->node1->region->y_min < node->region->y_max){
			top_right->node1 = n4;
		}
		else{
			top_right->node2 = n4;
		}
	}
	else{
		top_right = new edge(v_top->ID, cell->t_r->ID);
		cell->top->child2 = top_right;
		
		if(cell->top->node1 != NULL && cell->top->node1->region->y_min < node->region->y_max){
			top_right->node1 = cell->top->node2;
		}
		else{
			top_right->node1 = cell->top->node1;
		}
		top_right->node2 = n4;
	}
	
	//cout << "processed cell->top" << endl;

	//cout << "child1 = null?: " << (cell->bot->child1 == NULL) << endl;
	//cout << "child2 = null?: " << (cell->bot->child2 == NULL) << endl;
	//cout << "node1 = null?: " << (cell->bot->node1 == NULL) << endl;
	//cout << "node2 = null?: " << (cell->bot->node2 == NULL) << endl;
	
	if(cell->bot->child1 != NULL){
		bot_left = cell->bot->child1;
		if(bot_left->node1->region->y_max > node->region->y_min){
			bot_left->node1 = n1;
		}
		else{
			bot_left->node2 = n1;
		}
	}
	else{
		bot_left = new edge(cell->b_l->ID, v3->ID);
		cell->bot->child1 = bot_left;
		
		//cout << "node1 = node?: " << (cell->bot->node1 == node) << endl;
		//cout << "node2 = node?: " << (cell->bot->node2 == node) << endl;
		
		if(cell->bot->node1 != NULL && cell->bot->node1->region->y_max > node->region->y_min){
			bot_left->node1 = cell->bot->node2;
		}
		else{
			bot_left->node1 = cell->bot->node1;
		}
		bot_left->node2 = n1;
	}
	
	//cout << "done child1" << endl;
	
	if(cell->bot->child2 != NULL){
		bot_right = cell->bot->child2;
		if(bot_right->node1->region->y_max > node->region->y_min){
			bot_right->node1 = n3;
		}
		else{
			bot_right->node2 = n3;
		}
	}
	else{
		bot_right = new edge(v3->ID, cell->b_r->ID);
		cell->bot->child2 = bot_right;
		
		if(cell->bot->node1 != NULL && cell->bot->node1->region->y_max > node->region->y_min){
			bot_right->node1 = cell->bot->node2;
		}
		else{
			bot_right->node1 = cell->bot->node1;
		}
		bot_right->node2 = n3;
	}
	
	//cout << "processed cell->bot" << endl;
		
	if(cell->left->child1 != NULL){
		
		left_up = cell->left->child1;
		if(left_up->node1->region->x_max > node->region->x_min){
			left_up->node1 = n2;
		}
		else{
			left_up->node2 = n2;
		}
	}
	else{		
		left_up = new edge(v2->ID, cell->t_l->ID);
		cell->left->child1 = left_up;
		
		if(cell->left->node1 != NULL && cell->left->node1->region->x_max > node->region->x_min){
			left_up->node1 = cell->left->node2;
		}
		else{
			left_up->node1 = cell->left->node1;
		}
		left_up->node2 = n2;
	}
	
	//cout << "done child1" << endl;
	
	if(cell->left->child2 != NULL){
		left_down = cell->left->child2;
		if(left_down->node1->region->x_max > node->region->x_min){
			left_down->node1 = n1;
		}
		else{
			left_down->node2 = n1;
		}
	}
	else{
		left_down = new edge(cell->b_l->ID, v2->ID);
		cell->left->child2 = left_down;
		
		if(cell->left->node1 != NULL && cell->left->node1->region->x_max > node->region->x_min){
			left_down->node1 = cell->left->node2;
		}
		else{
			left_down->node1 = cell->left->node1;
		}
		left_down->node2 = n1;
	}
	
	////cout << "n1 left " << *(static_cast<CartCell*>(left_down->node1->poly)->bounding_box) << endl;
	////cout << "n2 left " << *(static_cast<CartCell*>(left_up->node1->poly)->bounding_box) << endl;
	
	//cout << "processed cell->left" << endl;
	
	if(cell->right->child1 != NULL){
		right_up = cell->right->child1;
		if(right_up->node1->region->x_min < node->region->x_max){
			right_up->node1 = n4;
		}
		else{
			right_up->node2 = n4;
		}
	}
	else{
		right_up = new edge(v_right->ID, cell->t_r->ID);
		cell->right->child1 = right_up;
		
		if(cell->right->node1 != NULL && cell->right->node1->region->x_min < node->region->x_max){
			right_up->node1 = cell->right->node2;
		}
		else{
			right_up->node1 = cell->right->node1;
		}
		right_up->node2 = n4;
	}
	
	//cout << "done child1" << endl;
	
	if(cell->right->child2 != NULL){
		right_down = cell->right->child2;
		if(right_down->node1->region->x_min < node->region->x_max){
			right_down->node1 = n3;
		}
		else{
			right_down->node2 = n3;
		}
	}
	else{
		right_down = new edge(cell->b_r->ID, v_right->ID);
		cell->right->child2 = right_down;
		
		if(cell->right->node1 != NULL && cell->right->node1->region->x_min < node->region->x_max){
			right_down->node1 = cell->right->node2;
		}
		else{
			right_down->node1 = cell->right->node1;
		}
		right_down->node2 = n3;
	}
	
	//cout << "processed cell->right" << endl;
	
	mid_up = new edge(v4->ID, v_top->ID);
	mid_left = new edge(v2->ID, v4->ID);
	mid_right = new edge(v4->ID, v_right->ID);
	mid_down = new edge(v3->ID, v4->ID);
	
	mid_up->node1 = n2;
	mid_up->node2 = n4;
	mid_left->node1 = n2;
	mid_left->node2 = n1;
	mid_right->node1 = n4;
	mid_right->node2 = n3;
	mid_down->node1 = n1;
	mid_down->node2 = n3;
	
	c1->top = mid_left;
	c1->bot = bot_left;
	c1->left = left_down;
	c1->right = mid_down;
	c2->top = top_left;
	c2->bot = mid_left;
	c2->left = left_up;
	c2->right = mid_up;
	c3->top = mid_right;
	c3->bot = bot_right;
	c3->left = mid_down;
	c3->right = right_down;
	c4->top = top_right;
	c4->bot = mid_right;
	c4->left = mid_up;
	c4->right = right_up;
	
	/*
	//cout << "n1 left ";
	if(find_other(c1->left, n1) != NULL){
		//cout << *(static_cast<CartCell*>(find_other(c1->left, n1)->poly)->bounding_box);
		//cout << *(static_cast<CartCell*>(find_self(c1->left, n1)->poly)->bounding_box);
	}
	else{
		//cout << "NULL" << endl;
	}
	
	
	//cout << "n2 left ";
	if(find_other(c2->left, n2) != NULL){
		//cout << *(static_cast<CartCell*>(find_other(c2->left, n2)->poly)->bounding_box);
		//cout << *(static_cast<CartCell*>(find_self(c2->left, n2)->poly)->bounding_box);
	}
	else{
		//cout << "NULL" << endl;
	}
	
	if(cell->left->node1 != NULL){
		//cout << "node1: " << *cell->left->node1->region;
	}
	if(cell->left->node2 != NULL){
		//cout << "node2: " << *cell->left->node2->region;
	}*/
	
	//linking
	c1->t_l = v2;
	c1->t_r = v4;
	c1->b_l = cell->b_l;
	c1->b_r = v3;
	c2->t_l = cell->t_l;
	c2->t_r = v_top;
	c2->b_l = v2;
	c2->b_r = v4;
	c3->t_l = v4;
	c3->t_r = v_right;
	c3->b_l = v3;
	c3->b_r = cell->b_r;
	c4->t_l = v_top;
	c4->t_r = cell->t_r;
	c4->b_l = v4;
	c4->b_r = v_right;

	return;
	
}

//helper to new_generate_dummy_mesh
void new_generate_dummy_mesh_helper(ADTNode* node, double size, double min_size){
	
	//cout << "refining the node" << endl;
	refine(node, true);
	
	double new_size = size / 2.0;
	
	if(new_size > min_size){
		//cout << "refining children of the node" << endl;
		new_generate_dummy_mesh_helper(node->L_branch->L_branch, new_size, min_size);
		new_generate_dummy_mesh_helper(node->L_branch->R_branch, new_size, min_size);
		new_generate_dummy_mesh_helper(node->R_branch->L_branch, new_size, min_size);
		new_generate_dummy_mesh_helper(node->R_branch->R_branch, new_size, min_size);
	}
	return;
}

//the initial mesh generation function, needs to be cropped by new_crop_mesh
ADTNode* new_generate_dummy_mesh(double size, double min_size, int pos_ID){
	//generates dummy square mesh with all depth 0
	
	vert_precise *t_l = new vert_precise(precise_verts.at(pos_ID)->x,
										 precise_verts.at(pos_ID)->y + size,
										 precise_verts.size());
	precise_verts.push_back(t_l);
	vert_precise *t_r = new vert_precise(precise_verts.at(pos_ID)->x + size,
										 precise_verts.at(pos_ID)->y + size,
										 precise_verts.size());
	precise_verts.push_back(t_r);
	vert_precise *b_l = precise_verts.at(pos_ID);
	vert_precise *b_r = new vert_precise(precise_verts.at(pos_ID)->x + size,
										 precise_verts.at(pos_ID)->y,
										 precise_verts.size());
	precise_verts.push_back(b_r);
	
	pt4D *sBox = new pt4D(b_l->x, b_l->x + size, b_l->y, b_l->y + size);
	CartCell *cell = new CartCell(sBox, size, pos_ID, *INTERNAL);
	ADTNode* init_mesh = new ADTNode(cell, 0, sBox);
	
	//setup big box of *INIT_SIZE and pass down for recursive division
	edge *top = new edge(t_l->ID, t_r->ID);
	edge *bot = new edge(b_l->ID, b_r->ID );
	edge *left = new edge(b_l->ID, t_l->ID);
	edge *right = new edge(b_r->ID, t_r->ID);
	
	cell->t_l = t_l;
	cell->t_r = t_r;
	cell->b_l = b_l;
	cell->b_r = b_r;
	cell->top = top;
	cell->bot = bot;
	cell->left = left;
	cell->right = right;
	
	top->node1 = init_mesh;
	bot->node1 = init_mesh;
	left->node1 = init_mesh;
	right->node1 = init_mesh;
	
	//cout << "root node is ok" << endl;
	
	new_generate_dummy_mesh_helper(init_mesh, size, min_size);
	
	return init_mesh;
}

//<--- line intersection algorithm
static bool inters_helper(double xi, double yi, double xj, double yj, double xk, double yk){
	return ((xi <= xk || xj <= xk) && (xk <= xi || xk <= xj) &&
		    (yi <= yk || yj <= yk) && (yk <= yi || yk <= yj));
}
static char compute_dir(double xi, double yi, double xj, double yj, double xk, double yk){
	double a = (xk - xi) * (yj - yi);
	double b = (xj - xi) * (yk - yi);
	return a < b ? -1 : a > b ? 1 : 0;
}
bool do_lines_intersect(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	// Do line segments (x1, y1)--(x2, y2) and (x3, y3)--(x4, y4) intersect?
	char d1 = compute_dir(x3, y3, x4, y4, x1, y1);
	char d2 = compute_dir(x3, y3, x4, y4, x2, y2);
	char d3 = compute_dir(x1, y1, x2, y2, x3, y3);
	char d4 = compute_dir(x1, y1, x2, y2, x4, y4);
	return (((d1 > 0 && d2 < 0) || (d1 < 0 && d2 > 0)) &&
			((d3 > 0 && d4 < 0) || (d3 < 0 && d4 > 0))) ||
			(d1 == 0 && inters_helper(x3, y3, x4, y4, x1, y1)) ||
			(d2 == 0 && inters_helper(x3, y3, x4, y4, x2, y2)) ||
			(d3 == 0 && inters_helper(x1, y1, x2, y2, x3, y3)) ||
			(d4 == 0 && inters_helper(x1, y1, x2, y2, x4, y4));
}
bool do_lines_intersect(Segment *seg, edge *eg){
	return do_lines_intersect(precise_verts.at(seg->list_of_verts[0])->x, 
							  precise_verts.at(seg->list_of_verts[0])->y,
							  precise_verts.at(seg->list_of_verts[1])->x,
							  precise_verts.at(seg->list_of_verts[1])->y,
							  precise_verts.at(eg->pt1)->x, precise_verts.at(eg->pt1)->y,
							  precise_verts.at(eg->pt2)->x, precise_verts.at(eg->pt2)->y);
}
bool do_lines_intersect(Segment *seg1, Segment *seg2){
	return do_lines_intersect(precise_verts.at(seg1->list_of_verts[0])->x, 
							  precise_verts.at(seg1->list_of_verts[0])->y,
							  precise_verts.at(seg1->list_of_verts[1])->x,
							  precise_verts.at(seg1->list_of_verts[1])->y,
							  precise_verts.at(seg2->list_of_verts[0])->x, 
							  precise_verts.at(seg2->list_of_verts[0])->y,
							  precise_verts.at(seg2->list_of_verts[1])->x,
							  precise_verts.at(seg2->list_of_verts[1])->y);
}
//end of line intersection algorithm --->

//returns the shortest distance of a point pt0 to the line pt1-pt2
double pt_line_dist(vert_precise *pt0, vert_precise *pt1, vert_precise *pt2){
	return ( fabs((pt2->x - pt1->x) * (pt1->y -  pt0->y) - (pt1->x - pt0->x) * (pt2->y - pt1->y)) /
				sqrt( pow((pt2->x - pt1->x),2) + pow((pt2->y - pt1->y),2) ) );
}

//returns the midpoint of a line a-b
vert_precise *find_midpoint(vert_precise *a, vert_precise *b){
	double new_x = (a->x + b->x) / 2;
	double new_y = (a->y + b->y) / 2;
	
	vert_precise *midpoint = new vert_precise(new_x, new_y, precise_verts.size());
	return midpoint;
}

//returns the vector of vertIDs of the new boundary after picking the furtherest
vector<int> *max_dist_vert(vector<int> *vec, vert_precise *pt1, vert_precise *pt2, int cap){
	//returns ID of the furtherest point in vec from pt1<->pt2
	double *dist_array;
	dist_array = new double[vec->size()];
	vector<int> *result = new vector<int>;
	
	/*
	if(vec->size() == 2){
		vert_precise *mid_point = find_midpoint(pt1, pt2);
		precise_verts.push_back(mid_point);
		result = mid_point->ID;
	}
	*/
	
	//map the verts to their distances
	for(int i = 1; i < (vec->size() - 1); i++){
		dist_array[i] = pt_line_dist(precise_verts.at(vec->at(i)), pt1, pt2);
	}
	
	result->push_back(vec->front());
	
	//check if wing tail case (both intersection points on same edge)
	if((pt1->x == pt2->x) || (pt1->y == pt2->y)){
		if(cap == 2){
			cap++;
		}
	}
	
	//keep the IDs of cap number of highest distances
	for(int i = 0; i < (cap - 2); i++){
		double max = -1.0;
		int  ID = -1;
		
		for(int j = 1; j < (vec->size() - 1); j++){
			if(dist_array[j] > max){
				max = dist_array[j];
				ID = j;
			}
		}
		
		result->push_back(vec->at(ID));
		dist_array[ID] = -1;
	}
	
	result->push_back(vec->back());
	
	delete [] dist_array;
	
	return result;
}

//returns the angle between lines (x1,y1)-(x2,y2) and (x3,y3)-(x4,y4)
double angle_between(double x1, double y1, double x2, double y2, double x3, double y3, double x4, double y4){
	 double dot_product = (((x1-x2)*(x3-x4)+(y1-y2)*(y3-y4))/(sqrt(pow((x1-x2),2)+pow((y1-y2),2))*sqrt(pow((x3-x4),2)+pow((y3-y4),2))));
	 return acos(dot_product); 
}
//differnt inputs for the checking angles
double angle_between(Segment *seg1, Segment *seg2){
	return angle_between(precise_verts.at(seg1->list_of_verts[0])->x,
						 precise_verts.at(seg1->list_of_verts[0])->y,
						 precise_verts.at(seg1->list_of_verts[1])->x,
						 precise_verts.at(seg1->list_of_verts[1])->y,
						 precise_verts.at(seg2->list_of_verts[0])->x,
						 precise_verts.at(seg2->list_of_verts[0])->y,
						 precise_verts.at(seg2->list_of_verts[1])->x,
						 precise_verts.at(seg2->list_of_verts[1])->y);
}
double angle_between(edge *e1, edge *e2){
	return angle_between(precise_verts.at(e1->pt1)->x,
						 precise_verts.at(e1->pt1)->y,
						 precise_verts.at(e1->pt2)->x,
						 precise_verts.at(e1->pt2)->y,
						 precise_verts.at(e2->pt1)->x,
						 precise_verts.at(e2->pt1)->y,
						 precise_verts.at(e2->pt2)->x,
						 precise_verts.at(e2->pt2)->y);
}
double angle_between(double x1, double y1, double x2, double y2, edge *eg){
	return angle_between(x1, y1, x2, y2,
						 precise_verts.at(eg->pt1)->x,
						 precise_verts.at(eg->pt1)->y,
						 precise_verts.at(eg->pt2)->x,
						 precise_verts.at(eg->pt2)->y);
}
double angle_between(double x1, double y1, double x2, double y2, vert_precise *a, vert_precise *b){
	return angle_between(x1, y1, x2, y2, a->x, a->y, b->x, b->y);
}

//checks if a slope of an edge eg is negative
bool is_slope_neg(edge *eg){
	if( ((precise_verts.at(eg->pt1)->x - precise_verts.at(eg->pt2)->x) != 0) &&
		((precise_verts.at(eg->pt1)->y - precise_verts.at(eg->pt2)->y)/
		 (precise_verts.at(eg->pt1)->x - precise_verts.at(eg->pt2)->x)) < 0){
		 return true;
	}
	return false;
}

//<--- remove dup algorithm (sort & unique)
bool is_equal(vert_precise* point, vert_precise* ref_point){
	//compare by value, x1 == x2, y1 == y2
	return (fabs(point->x - ref_point->x) < *EPSILON
		 && fabs(point->y - ref_point->y) < *EPSILON);
}
bool less_than(vert_precise* a, vert_precise* b){
	if((a->x) < (b->x)){
		return true;
	}
	else if((a->x) == (b->x)){
		return ((a->y) < (b->y));
	}
	return false;
}
bool more_than(vert_precise* a, vert_precise* b){
	if((a->x) > (b->x)){
		return true;
	}
	else if((a->x) == (b->x)){
		return ((a->y) > (b->y));
	}
	return false;
}
void removeDuplicates(vector<vert_intersection*> *vec){
   sort(vec->begin(), vec->end(), less_than);
   vec->erase(unique(vec->begin(), vec->end(), is_equal), vec->end());
}
void removeDuplicates(vector<vert_precise*> *vec){
   sort(vec->begin(), vec->end(), less_than);
   vec->erase(unique(vec->begin(), vec->end(), is_equal), vec->end());
}
//end of remove dup algorithm --->

//returns the number of nodes in an ADT
int num_of_nodes(ADTNode *tree){
	//count number of nodes in a tree
	if(tree == NULL){
		return 0;
	}
	else{
		return (1 + num_of_nodes(tree->L_branch) + num_of_nodes(tree->R_branch));
	}
}

//frees an ADT
void free_ADTNode(ADTNode* node){
	if(node == NULL){
		return;
	}
	free_ADTNode(node->L_branch);
	free_ADTNode(node->R_branch);
	delete node;
}

//prints an ADT
void print_ADT(ADTNode* tree){
	if(tree == NULL){
		return;
	}
	//cout << *(tree->region) << *(tree->poly->bounding_box);
	//cout << tree->depth << endl << endl;
	print_ADT(tree->L_branch);
	print_ADT(tree->R_branch);
}

//prints an edge
void edge_stat(edge* eg){
	if(eg == NULL){
		return;
	}
	cout << "(" << precise_verts.at(eg->pt1)->x << ", " << precise_verts.at(eg->pt1)->y << ")" << " -- "
		 << "(" << precise_verts.at(eg->pt2)->x << ", " << precise_verts.at(eg->pt2)->y << ")" << endl;
}

//prints a segment
void edge_stat(Segment* seg){
	if(seg == NULL){
		//cout << "seg is NULL" << endl;
		return;
	}
	cout << "(" << precise_verts.at(seg->list_of_verts[0])->x << ", " << precise_verts.at(seg->list_of_verts[0])->y << ")" << " -- "
		 << "(" << precise_verts.at(seg->list_of_verts[1])->x << ", " << precise_verts.at(seg->list_of_verts[1])->y << ")" << endl;
}

//comparsion functions, for sorting by x corrdinate or y corrdinate
bool comp_in_x(Segment *s1, Segment *s2){
	//compare segments by the highest x coordinate
	double v1 = max(precise_verts.at(s1->list_of_verts[0])->x, 
					precise_verts.at(s1->list_of_verts[1])->x);
					
	double v2 = max(precise_verts.at(s2->list_of_verts[0])->x, 
					precise_verts.at(s2->list_of_verts[1])->x);
	return (v1 < v2);
}
bool comp_in_y(Segment *s1, Segment *s2){
	//compare segments by the highest y coordinate
	double v1 = max(precise_verts.at(s1->list_of_verts[0])->y, 
					precise_verts.at(s1->list_of_verts[1])->y);
					
	double v2 = max(precise_verts.at(s2->list_of_verts[0])->y, 
					precise_verts.at(s2->list_of_verts[1])->y);
	return (v1 < v2);
}

//overlap functions, strict means all comparsions are strict
bool overlap(pt4D *A, pt4D *B){
	return ( (A->x_min < B->x_max) && (A->x_max > B->x_min)
		  && (A->y_min < B->y_max) && (A->y_max > B->y_min) );
}
bool overlap_nonstrict(pt4D *A, pt4D *B){
	return ( (A->x_min <= B->x_max) && (A->x_max >= B->x_min)
		  && (A->y_min <= B->y_max) && (A->y_max >= B->y_min) );
}
bool overlap_x_or_y(pt4D *A, pt4D *B){
	return ( (A->x_min <= B->x_max) && (A->x_max >= B->x_min)
		  || (A->y_min <= B->y_max) && (A->y_max >= B->y_min) );
}
bool contained(pt4D *container, pt4D *obj){
	return ( (container->x_min <= obj->x_min) && (container->x_max >= obj->x_max)
		  && (container->y_min <= obj->y_min) && (container->y_max >= obj->y_max) );
}
bool contained(vert_precise *pt, pt4D *container){
	return ( (container->x_min <= pt->x) && (container->x_max >= pt->x)
		  && (container->y_min <= pt->y) && (container->y_max >= pt->y) );
}
bool contained_strict(vert_precise *pt, pt4D *container){
	return ( (container->x_min < pt->x) && (container->x_max > pt->x)
		  && (container->y_min < pt->y) && (container->y_max > pt->y) );
}

//prints a cellpiece or CutCell class
void cell_piece_info(CutCell *piece){

	int part_entries = 1;
	int part_sum = 0;
	
	if(piece->num_entries != NULL){
		part_entries = *(piece->num_entries);
		//cout << "cutcell, not cellpiece" << endl;
	}
	
	//cout << "part_entries: " << part_entries << endl;
	for(int i = 1; i <= part_entries; i++){
		int part_elenum = *(piece->num_indices[i-1]);
				
		for(int temp_index = 0; temp_index < part_elenum; temp_index++){
			//cout << *(precise_verts.at(piece->list_of_verts[part_sum + temp_index])) << endl;
			//edge_stat(precise_edges.at(piece->list_of_edges[part_sum + temp_index]));
			//cout << endl;
		}
		part_sum += part_elenum;
	}
		
}
void cell_piece_print(CutCell *piece){

	int part_entries = 1;
	
	if(piece->num_entries != NULL){
		part_entries = *(piece->num_entries);
	}
		
	//cout << "part_entries: " << part_entries << endl;
	for(int i = 1; i <= part_entries; i++){
		int part_elenum = *(piece->num_indices[i-1]);
		
		for(int temp_index = 0; temp_index < part_elenum; temp_index++){
			//edge_stat(precise_edges.at(piece->list_of_edges[temp_index]));
			print_edges.push_back(precise_edges.at(piece->list_of_edges[temp_index]));
		}
	}
		
}

//removes the referrence of a point pt from a cell_piece piece's intersection list
void remove_pt(vert_intersection *pt, cell_piece *piece){
	vector<vert_intersection*> *temp = piece->list_of_intersecton_points;
	for(int i = 0, cap = temp->size(); i < cap; ++i, cap = temp->size()){
		if( is_equal(pt, temp->at(i)) ){
			temp->at(i) = NULL;
			temp->erase(temp->begin() + i); 
			i--;
			//cout << "vert erased from interv_list" << endl;
			//break;
		}
	}
}

//checks of a point is on the line
bool is_pt_on_line(vert_precise *pt, Segment *seg){

	//cout << *pt << endl;
	//edge_stat(seg);

	if( is_equal(pt, precise_verts.at(seg->list_of_verts[0])) 
	 || is_equal(pt, precise_verts.at(seg->list_of_verts[1]))){
		return true;
	}
	
	double x1 = precise_verts.at(seg->list_of_verts[0])->x;
	double y1 = precise_verts.at(seg->list_of_verts[0])->y;
	double x2 = precise_verts.at(seg->list_of_verts[1])->x;
	double y2 = precise_verts.at(seg->list_of_verts[1])->y;
	
	pt4D *temp = seg->bounding_box;
	if(contained(pt, temp)){
		if( (temp->x_max - temp->x_min) <= *EPSILON ){
			//vert line
			return true;
		}
		const double M = (y2 - y1) / (x2 - x1);	//slope
		const double B = -(M * x1) + y1;	//y-intercept
		//check if pt is on the line
		
		//cout << "variance: " << fabs(pt->y - (M * pt->x + B)) << endl;
		
		return ( fabs(pt->y - (M * pt->x + B) ) <= *EPSILON );
	}
	//if not even contained, not on line
	return false;
}
bool is_pt_on_line(vert_precise *pt, edge *eg){	
	if( is_equal(pt, precise_verts.at(eg->pt1)) 
	 || is_equal(pt, precise_verts.at(eg->pt2))){
		return true;
	}
	double x1 = precise_verts.at(eg->pt1)->x;
	double y1 = precise_verts.at(eg->pt1)->y;
	double x2 = precise_verts.at(eg->pt2)->x;
	double y2 = precise_verts.at(eg->pt2)->y;	
		
	double x_min = min(x1, x2);
	double x_max = max(x1, x2);
	double y_min = min(y1, y2);
	double y_max = max(y1, y2);
		 
	pt4D *temp = new pt4D(x_min, x_max, y_min, y_max);
	if(contained(pt, temp)){
		if( (temp->x_max - temp->x_min) <= *EPSILON ){
			//vert line
			delete temp;
			return true;
		}
		const double M = (y2 - y1) / (x2 - x1);	//slope
		const double B = -(M * x1) + y1;	//y-intercept
		//check if pt is on the line
		delete temp;
		
		//cout << "variance: " << fabs(pt->y - (M * pt->x + B)) << endl;
		
		return ( fabs(pt->y - (M * pt->x + B) ) <= *EPSILON );
	}
	//if not even contained, not on line
	delete temp;
	return false;
}

//terminating conditional functions passed to the walk function 
bool cond_onbound(vert_precise* point, void* dummy){
	return (static_cast<vert_intersection*>(point)->status == ON_BOUND);
}
bool is_in (vert_precise* point, void* ref_edge){
	if(ref_edge == NULL){
		return false;
	}
	edge *eg = static_cast<edge*>(ref_edge);
	return is_pt_on_line(point, eg);
}

//checks if a vertice vert is inside a cell
vert_cond is_inside_helper(double x, double y, double l, double r, double b, double t){
	if((x > l) && (x < r) && (y > b) && (y < t)){
		//cout << "vert_cond is IN" << endl;
		return IN;
	}
	else if((x < l) || (x > r) || (y < b) || (y > t)){
		//cout << "vert_cond is OUT" << endl;
		return OUT;
	}
	else{
		return ON_BOUND;
	}
}
vert_cond is_inside(CartCell *cell, vert_precise* vert){
	vert_precise* temp = precise_verts.at(cell->pos);
	
	//cout << "vert_cond is_inside is called" << endl;
	//cout << vert->x << endl;
	//cout << vert->y << endl;
	//cout << "cell dimensions:" << endl;	
	//cout << temp->x << endl;
	//cout << temp->x + cell->length << endl;
	//cout << temp->y << endl;
	//cout << temp->y + cell->length << endl;
	//cout << "^ the coordinates passed ^ " << endl;
	//cout << "result: " << is_inside_helper(vert->x, vert->y, 
	//						temp->x, temp->x + cell->length,
	//						temp->y, temp->y + cell->length) << endl;
	
	return is_inside_helper(vert->x, vert->y, 
							temp->x, temp->x + cell->length,
							temp->y, temp->y + cell->length);
}

//returns edge formed by the single segment that cuts the cell, NULL if case is trivial (does not cut the cell)
body_edge *trivial_check(cell_piece *piece, Segment *seg){
	int count = 0;
	vector<vert_intersection*> *temp = piece->list_of_intersecton_points;
	vector<vert_intersection*> bag;
	
	for(int i = 0; i < temp->size(); i++){
		if(is_pt_on_line(temp->at(i) , seg)){
			count++;
			bag.push_back(temp->at(i));
			if(count == 2){
				body_edge *new_edge = new body_edge(bag.at(0)->ID, bag.at(1)->ID, NULL, false);
				return new_edge;
			}
		}
	}
	return NULL;
}

//removes all trivial intersections in a cell_piece class
void remove_trivial_seg(cell_piece *piece){	

	//cout << "remove_trivial_seg called" << endl;
	
	vector<Polyhedron*> *temp = piece->intersect_segments;
	
	//cout << "temp has size: " << temp->size() << endl;
	
	for(int i = 0, cap = temp->size(); i < cap; ++i, cap = temp->size()){
	
		//cout << "remove_trivial_seg processing edge " << i << endl;
	
		Segment *seg = static_cast<Segment*>(temp->at(i));
		
		//cout << "seg declared" << endl;
		
		vert_intersection *v1 = static_cast<vert_intersection*>
									(precise_verts.at(seg->list_of_verts[0]));
		vert_intersection *v2 = static_cast<vert_intersection*>
									(precise_verts.at(seg->list_of_verts[1]));
	
		//cout << "verts declared" << endl;
		
		//edge_stat(seg);
		
		if( ((v1->status != IN) && (v2->status != IN)) &&
			((v1->status == ON_BOUND) || (v2->status == ON_BOUND)) ){
				
			//cout << "suspicious edge, checking" << endl;
			
			body_edge *check = trivial_check(piece, seg);
			
			if(check != NULL){
				delete check;
			}
			else{
				//trivial intersection
				//cout << "%%% trivial intersection %%%" << endl;
				temp->at(i) = NULL;
				temp->erase(temp->begin() + i);
				i--;
				continue;
			}
		}
	}
}

//generates list of points/edges fields for a cell
void generate_cell_lists(CartCell* cell){
	//indice numbers
	int **temp_array = new int*[1];
	temp_array[0] = &CARTCELL_NUM_INDICES;
	cell->num_indices = temp_array;
	
	int v_newID = precise_verts.size();
	//list of verts and edges
	vert_precise *v2 = new vert_precise(precise_verts.at(cell->pos)->x,
										precise_verts.at(cell->pos)->y + cell->length, v_newID);
	vert_precise *v3 = new vert_precise(precise_verts.at(cell->pos)->x + cell->length,
										precise_verts.at(cell->pos)->y + cell->length, v_newID+1);
	vert_precise *v4 = new vert_precise(precise_verts.at(cell->pos)->x + cell->length,
										precise_verts.at(cell->pos)->y, v_newID+2);
	precise_verts.push_back(v2);
	precise_verts.push_back(v3);
	precise_verts.push_back(v4);
	
	
	edge *e1 = new edge(cell->pos, v_newID);	//pos - v2
	edge *e2 = new edge(v_newID, v_newID+1);	//v2 - v3
	edge *e3 = new edge(v_newID+1, v_newID+2);	//v3 - v4
	edge *e4 = new edge(v_newID+2, cell->pos);	//v4 - pos
	
	
	int e_newID = precise_edges.size();
	precise_edges.push_back(e1);
	precise_edges.push_back(e2);
	precise_edges.push_back(e3);
	precise_edges.push_back(e4);
	
	int *vert = new int[4];
	int *edge = new int[4];
	vert[0] = cell->pos;
	vert[1] = v_newID;
	vert[2] = v_newID+1;
	vert[3] = v_newID+2;
	edge[0] = e_newID;
	edge[1] = e_newID+1;
	edge[2] = e_newID+2;
	edge[3] = e_newID+3;
	
	cell->list_of_verts = vert;
	cell->list_of_edges = edge;
	/*
	//cout << "verts: " << *(precise_verts.at(cell->pos)) << " "
					  << *(precise_verts.at(v_newID)) << " "
					  << *(precise_verts.at(v_newID+1)) << " "
					  << *(precise_verts.at(v_newID+2)) << " " << endl;
							  
	//cout << "edges: " << *(precise_edges.at(e_newID)) << " "
					  << *(precise_edges.at(e_newID+1)) << " "
					  << *(precise_edges.at(e_newID+2)) << " "
					  << *(precise_edges.at(e_newID+3)) << " " << endl;
	*/
}

//produces the ADT for segments (init geometry)
ADTNode* make_segment_ADT(vector<Segment*> vec_seg, int depth, pt4D *region){
	//inserts a vector of segments into an ADT
	//base case
	if(vec_seg.size() == 0){
		return NULL;
	}
	
	int dir = depth % DIM;	//find direction
	
	//sizes of the children vectors
	int med_index = (vec_seg.size() / 2);
	int left_size = med_index;
	int right_size = med_index + (vec_seg.size() % 2 - 1);
	

	switch(dir){	
	case 0:	//if x
		sort(vec_seg.begin(), vec_seg.end(), comp_in_x);	//sort in dir of x				
		break;
	case 1:	//if y
		sort(vec_seg.begin(), vec_seg.end(), comp_in_y);	//sort in dir of y
		break;
	}
	
	vector<Segment*>::iterator it = vec_seg.begin();	//iterator for dividing the vector
	it += med_index;	//move iterator to middle of the vector
	double new_bound;
	
	//push until new_bound does not intersect anything
	switch(dir){	
	case 0:	//if x
		new_bound = max(precise_verts.at((*it)->list_of_verts[0])->x, 
						precise_verts.at((*it)->list_of_verts[1])->x);
							
		for(vector<Segment*>::iterator itr = it + 1; itr != vec_seg.end(); itr++){
			if( min(precise_verts.at((*itr)->list_of_verts[0])->x, 
					precise_verts.at((*itr)->list_of_verts[1])->x) <= new_bound ){
				
				left_size++;
				right_size--;
				new_bound = max(precise_verts.at((*itr)->list_of_verts[0])->x, 
								precise_verts.at((*itr)->list_of_verts[1])->x);
				it++;
			}
		}
		break;
		
	case 1:	//if y
		new_bound = max(precise_verts.at((*it)->list_of_verts[0])->y, 
						precise_verts.at((*it)->list_of_verts[1])->y);
						
		for(vector<Segment*>::iterator itr = it + 1; itr != vec_seg.end(); itr++){
			if( min(precise_verts.at((*itr)->list_of_verts[0])->y, 
					precise_verts.at((*itr)->list_of_verts[1])->y) <= new_bound ){
				
				left_size++;
				right_size--;
				new_bound = max(precise_verts.at((*itr)->list_of_verts[0])->y, 
								precise_verts.at((*itr)->list_of_verts[1])->y);
				it++;
			}
		}
		break;
	}
	
	//cout << "left size: " << left_size << endl;
	//cout << "right size: " << right_size << endl;
	
	//preparing the children vectors
	vector<Segment*> left;
	vector<Segment*> right;	
	left.resize(left_size);
	right.resize(right_size);
	
	//copying values to the children vectors
	copy(vec_seg.begin(), it, left.begin());
	copy(it+1, vec_seg.end(), right.begin());
	
	pt4D *l_bound = new pt4D(region->x_min, region->x_max, region->y_min, region->y_max);
	pt4D *r_bound = new pt4D(region->x_min, region->x_max, region->y_min, region->y_max);
	
	switch(dir){	
	case 0:	//if x
		l_bound->x_max = new_bound;
		r_bound->x_min = new_bound;
		break;
	case 1:	//if y
		l_bound->y_max = new_bound;
		r_bound->y_min = new_bound;
		break;
	}	
	
	//recusive call 
	ADTNode *temp = new ADTNode(*it, depth, region);
	temp->L_branch = make_segment_ADT(left, depth+1, l_bound);
	temp->R_branch = make_segment_ADT(right, depth+1, r_bound);
	
	return temp;
}

//finds all intersection candidates for a region by search through an ADT of the init geometry
vector<Polyhedron*> *find_intersections(pt4D *box, ADTNode *surface, bool (*overlap_func)(pt4D*, pt4D*) ){
	//find possible intersections of a box
	//check if objects on surface is in the box
	ADTNode *temp = surface;
	stack<ADTNode*> right_branches;	//holds the right branches
	vector<Polyhedron*> *candidates = new vector<Polyhedron*>;
	
	SEARCH:
		
		//if the regions overlap, is a candidate
		if( overlap_func(box, temp->poly->bounding_box) ){
			candidates->push_back(temp->poly);
		}

		//if right branch exist and overlaps, stock it up for later
		if( (temp->R_branch != NULL) && (overlap_x_or_y(box, temp->R_branch->region) || 
										 overlap_x_or_y(box, temp->R_branch->poly->bounding_box)) ){
			right_branches.push(temp->R_branch);
		}
		
		//if left branch exist and overlaps, check it next
		if( (temp->L_branch != NULL) && (overlap_x_or_y(box, temp->L_branch->region) || 
										 overlap_x_or_y(box, temp->L_branch->poly->bounding_box)) ){
			temp = temp->L_branch;
			goto SEARCH;
		}
		
		//once no more left branch, start checking the stocked right branches
		if(!right_branches.empty()){
			temp = right_branches.top();
			right_branches.pop();
			goto SEARCH;
		}
		
	if(candidates->size() == 0){
		delete candidates;
		return NULL;
	}
	return candidates;
}

//find_intersection that shows all steps, for debugging purposes
vector<Polyhedron*> *find_intersections_show(pt4D *box, ADTNode *surface, bool (*overlap_func)(pt4D*, pt4D*) ){
	//find possible intersections of a box
	//check if objects on surface is in the box
	ADTNode *temp = surface;
	stack<ADTNode*> right_branches;	//holds the right branches
	vector<Polyhedron*> *candidates = new vector<Polyhedron*>;
	
	//cout << "finding intersections" << endl;
	
	SEARCH:
		
		//cout << *box << *temp->region << *temp->poly->bounding_box
		//	 << overlap_func(box, temp->poly->bounding_box) << endl << endl;
		
		if( overlap_func(box, temp->poly->bounding_box) ){
			candidates->push_back(temp->poly);
		}
		//}
		if( (temp->R_branch != NULL) && (overlap_x_or_y(box, temp->R_branch->region) || 
										 overlap_x_or_y(box, temp->R_branch->poly->bounding_box)) ){
		
			//cout << "pushed in right branch stack: " << endl;
			//cout << *temp->R_branch->poly->bounding_box
			//	 << *temp->R_branch->region << endl;
		
			right_branches.push(temp->R_branch);
		}
		
		/////////
		if( (temp->L_branch != NULL) ){
			//cout << "left is: " << *temp->L_branch->region;
			//cout << overlap_x_or_y(box, temp->L_branch->region) << endl;
			//cout << overlap_x_or_y(box, temp->L_branch->poly->bounding_box) << endl;
		}
		/////////
		
		if( (temp->L_branch != NULL) && (overlap_x_or_y(box, temp->L_branch->region) ||
										 overlap_x_or_y(box, temp->L_branch->poly->bounding_box)) ){
		
			//cout << "left branch next: " << endl;
			//cout << *temp->L_branch->region << endl;
		
			temp = temp->L_branch;
			goto SEARCH;
		}
		if(!right_branches.empty()){
			temp = right_branches.top();
			right_branches.pop();
			goto SEARCH;
		}
		
	if(candidates->size() == 0){
		delete candidates;
		return NULL;
	}
	return candidates;
}

//initial mesh generating functions for Manual and Min_size mode
ADTNode* generate_dummy_mesh(double size, double min_size, int pos_ID){
	//generates dummy square mesh with all depth 0
	//improvement: use goto to avoid recursion
	
	pt4D *sBox = new pt4D(precise_verts.at(pos_ID)->x, precise_verts.at(pos_ID)->x + size,
						  precise_verts.at(pos_ID)->y, precise_verts.at(pos_ID)->y + size);
	CartCell *current_cell = new CartCell(sBox, size, pos_ID, *INTERNAL);
	pt4D *current_bound = sBox;
	ADTNode* mesh = new ADTNode(current_cell, 0, sBox);
	
	if(size <= min_size){
		return mesh;	//do not divide
	}
	
	//<--- divide cell in 4
	double med = size / 2.0;
		//<--- divide cell in 2
		pt4D *l_bound = new pt4D(current_bound->x_min, current_bound->x_max - med,
								 current_bound->y_min, current_bound->y_max);
		pt4D *r_bound = new pt4D(current_bound->x_min + med, current_bound->x_max,
								 current_bound->y_min, current_bound->y_max);
		Polyhedron *rect1 = new Polyhedron();
		Polyhedron *rect2 = new Polyhedron();
		rect1->bounding_box = l_bound;
		rect2->bounding_box = r_bound;
		mesh->L_branch = new ADTNode(rect1, 0, l_bound);
		mesh->R_branch = new ADTNode(rect2, 0, r_bound);
		// end of divide cell in 2 --->
	pt4D *b_c1 = new pt4D(l_bound->x_min, l_bound->x_max,
						  l_bound->y_min, l_bound->y_max - med);
	pt4D *b_c2 = new pt4D(l_bound->x_min, l_bound->x_max,
						  l_bound->y_min + med, l_bound->y_max);
	pt4D *b_c3 = new pt4D(r_bound->x_min, r_bound->x_max,
						  r_bound->y_min, r_bound->y_max - med);
	pt4D *b_c4 = new pt4D(r_bound->x_min, r_bound->x_max,
						  r_bound->y_min + med, r_bound->y_max);
			
	int newID = precise_verts.size();	//ID of the next added entry
	vert_precise *v2 = new vert_precise(precise_verts.at(current_cell->pos)->x,
										precise_verts.at(current_cell->pos)->y + med, newID);
	vert_precise *v3 = new vert_precise(precise_verts.at(current_cell->pos)->x + med,
										precise_verts.at(current_cell->pos)->y, newID+1);
	vert_precise *v4 = new vert_precise(precise_verts.at(current_cell->pos)->x + med,
										precise_verts.at(current_cell->pos)->y + med, newID+2);										
	precise_verts.push_back(v2);
	precise_verts.push_back(v3);
	precise_verts.push_back(v4);
	//recursive call to generate the children
	mesh->L_branch->L_branch = generate_dummy_mesh(med, min_size, current_cell->pos);
	mesh->L_branch->R_branch = generate_dummy_mesh(med, min_size, newID);
	mesh->R_branch->L_branch = generate_dummy_mesh(med, min_size, newID+1);
	mesh->R_branch->R_branch = generate_dummy_mesh(med, min_size, newID+2);	
	
	return mesh;
}
ADTNode* crop_mesh(ADTNode* mesh, pt4D *keep_box){
	//discards nodes from a dummy mesh that's not in keep_box
	//improvement: use goto to avoid recursion
	if(mesh == NULL){
		return NULL;
	}
	if(!overlap(keep_box, mesh->region)){
		free_ADTNode(mesh);
		return NULL;
	}
	if(contained(keep_box, mesh->region)){
		return mesh;
	}
	else{
		mesh->L_branch = crop_mesh(mesh->L_branch, keep_box);
		mesh->R_branch = crop_mesh(mesh->R_branch, keep_box);
	}
	return mesh;
}

//sets the pointer ptr to null
void set_null_ptr(ADTNode **ptr){
    *ptr = NULL;
}

//crop mesh function for Adaptive mode, sets new border edges' outside pointers to NULL
void new_crop_mesh_helper(ADTNode* mesh, pt4D *keep_box){
	if(mesh == NULL){
		return;
	}
	
	CartCell *cell = static_cast<CartCell*>(mesh->poly);
	
	ADTNode *top = find_other(cell->top, mesh);
	ADTNode *bot = find_other(cell->bot, mesh);
	ADTNode *left = find_other(cell->left, mesh);
	ADTNode *right = find_other(cell->right, mesh);
	
	if( top != NULL && !overlap(keep_box, top->region)){
		set_null_ptr(find_other_ref(cell->top, mesh));
	}
	if( bot != NULL && !overlap(keep_box, bot->region)){
		set_null_ptr(find_other_ref(cell->bot, mesh));
	}
	if( left != NULL && !overlap(keep_box, left->region)){
		set_null_ptr(find_other_ref(cell->left, mesh));
	}
	if( right != NULL && !overlap(keep_box, right->region)){
		set_null_ptr(find_other_ref(cell->right, mesh));
	}
	
	if(mesh->L_branch != NULL){
		new_crop_mesh_helper(mesh->L_branch->L_branch, keep_box);
		new_crop_mesh_helper(mesh->L_branch->R_branch, keep_box);
	}
	if(mesh->R_branch != NULL){
		new_crop_mesh_helper(mesh->R_branch->L_branch, keep_box);
		new_crop_mesh_helper(mesh->R_branch->R_branch, keep_box);
	}
	return;
}
ADTNode* new_crop_mesh(ADTNode* mesh, pt4D *keep_box){
	//discards nodes from a dummy mesh that's not in keep_box
	//improvement: use goto to avoid recursion
	if(mesh == NULL){
		return NULL;
	}
	if(!overlap(keep_box, mesh->region)){
		return NULL;
	}
	if(contained(keep_box, mesh->region)){
		return mesh;
	}
	else{
		mesh->L_branch = new_crop_mesh(mesh->L_branch, keep_box);
		mesh->R_branch = new_crop_mesh(mesh->R_branch, keep_box);
	}
	return mesh;
}

//check if two lines form an acute angle
bool is_acute(vert_precise *a, vert_precise *b, vert_precise *c, vert_precise *d){
				
	double ty1 = angle_between(0, 0, 0, 1, a, b);
	double ty2 = angle_between(0, 0, 0, 1, c, d);
	double tx1 = angle_between(0, 0, 1, 0, a, b);
	double tx2 = angle_between(0, 0, 1, 0, c, d);
	
	double ty3 = angle_between(0, 0, 0, 1, a, b);
	double ty4 = angle_between(0, 0, 0, 1, c, d);
	double tx3 = angle_between(0, 0, 1, 0, a, b);
	double tx4 = angle_between(0, 0, 1, 0, c, d);
	
	if( (fabs(tx1 - tx2) < PI/2) ||
		(fabs(ty1 - ty2) < PI/2) ||
		(fabs(tx3 - tx4) < PI/2) ||
		(fabs(ty3 - ty4) < PI/2) ){
		return true;
	}
	
	return false;
}
bool is_acute(Segment *seg1, Segment *seg2){
	return is_acute(precise_verts.at(seg1->list_of_verts[0]),	
					precise_verts.at(seg1->list_of_verts[1]),
					precise_verts.at(seg2->list_of_verts[0]),
					precise_verts.at(seg2->list_of_verts[1]));
}

//check if a vertice or CartCell is inside or outside of the geometry
bool inside_geometry(CartCell* cell, ADTNode *surface, int shift = 0){
	
	//account for the shift (zero on first call)
	double x1 = *GBOX_X_MIN + shift * (*INIT_MIN / 20);
	double x2 = precise_verts.at(cell->pos)->x;
	
	//position of cell to position on global box
	pt4D test_box = pt4D(min(x1, x2), max(x1, x2),
						 *GBOX_Y_MIN, precise_verts.at(cell->pos)->y);
						 
	vector<Polyhedron*> *intersect_cands;
	vector<Segment*> intersect_check;
	int num_intersections = 0;
	//int num_double_counts = 0;
	
	intersect_cands = find_intersections(&test_box, surface, overlap_nonstrict);
	if(intersect_cands == NULL){
		if(*INTERNAL){
			return true;
		}
		else{
			return false;
		}
	}
	for(int i = 0; i < intersect_cands->size(); i++){
		Segment *seg = static_cast<Segment*>(intersect_cands->at(i));
		if( do_lines_intersect(x1, *GBOX_Y_MIN,	//line one
									x2, 
									precise_verts.at(cell->pos)->y,
									precise_verts.at(seg->list_of_verts[0])->x,	//line two
									precise_verts.at(seg->list_of_verts[0])->y,
									precise_verts.at(seg->list_of_verts[1])->x,
									precise_verts.at(seg->list_of_verts[1])->y) ){
			num_intersections++;
			intersect_check.push_back(seg);
		}
	}
	
	//if intersects as boundary point, shifts the ray and try again
	while(intersect_check.size() > 0){
		Segment *check_seg = intersect_check.back();
		
		for(int i = 0; i < (intersect_check.size() - 1); i++){
			if(do_lines_intersect(intersect_check.at(i), check_seg)){
					// && !is_acute(intersect_check.at(i), check_seg)
				//num_double_counts++;
				return inside_geometry(cell, surface, shift + 1);
			}
		}
		intersect_check.back() = NULL;
		intersect_check.pop_back();
	}
	
	delete intersect_cands;
	
	if(!*INTERNAL)
		return static_cast<bool>(num_intersections % 2);
	else
		return !(static_cast<bool>(num_intersections % 2));
}
bool inside_geometry(vert_precise *pt, ADTNode *surface, int shift = 0){
	
	//account for the shift (zero on first call)
	double x1 = *GBOX_X_MIN + shift * (*INIT_MIN / 20);
	double x2 = pt->x;
	
	//position of cell to position on global box
	pt4D test_box = pt4D(min(x1, x2), max(x1, x2),
						 *GBOX_Y_MIN, pt->y);
						 
	vector<Polyhedron*> *intersect_cands;
	vector<Segment*> intersect_check;
	int num_intersections = 0;
	//int num_double_counts = 0;
	
	intersect_cands = find_intersections(&test_box, surface, overlap);
	if(intersect_cands == NULL){
		if(*INTERNAL){
			return true;
		}
		else{
			return false;
		}
	}
	
	for(int i = 0; i < intersect_cands->size(); i++){
		Segment *seg = static_cast<Segment*>(intersect_cands->at(i));
		
		if( do_lines_intersect(x1, *GBOX_Y_MIN,	//line one
									x2, 
									pt->y,
									precise_verts.at(seg->list_of_verts[0])->x,	//line two
									precise_verts.at(seg->list_of_verts[0])->y,
									precise_verts.at(seg->list_of_verts[1])->x,
									precise_verts.at(seg->list_of_verts[1])->y) ){
			num_intersections++;	
			intersect_check.push_back(seg);
		}
	}
	
	//if intersects as boundary point, shifts the ray and try again
	while(intersect_check.size() != 0){
		Segment *check_seg = intersect_check.back();
		for(int i = 0; i < (intersect_check.size() - 1); i++){
			if(do_lines_intersect(intersect_check.at(i), check_seg)){
					// && !is_acute(intersect_check.at(i), check_seg)
				//num_double_counts++;
				return inside_geometry(pt, surface, shift + 1);
			}
		}
		intersect_check.back() = NULL;
		intersect_check.pop_back();
	}
	
	delete intersect_cands;
	
	//cout << "num_intersections: " << num_intersections << endl;
	
	if(!*INTERNAL)
		return static_cast<bool>(num_intersections % 2);
	else
		return !(static_cast<bool>(num_intersections % 2));
}

//initial refinement for Manual and Minsize Mode
ADTNode* mesh_refine(ADTNode *surface, Refinement *ref, ADTNode *mesh, int depth){
	//helper for manual version of generate_mesh
	ADTNode* current_node = mesh;
	if(current_node == NULL){
		return mesh;	//base case
	}
	if((current_node->L_branch == NULL) && (current_node->R_branch == NULL)){	//if is leaf
		CartCell *current_cell = static_cast<CartCell*>(mesh->poly);
		pt4D* current_bound = current_cell->bounding_box;	//convenience pointers
		//no more divide, find intersection and generate verts and edges for this cell
		
		//cout << "the cell_bound in mesh_refine" << endl;
		//cout << *current_bound << endl;
		
		current_cell->intersect_candidates = find_intersections(current_bound, surface, overlap_nonstrict);
		
		if(current_cell->intersect_candidates != NULL){
			//cout << "intersect_candidates's size: " << current_cell->intersect_candidates->size() << endl;
		}
		
		//ref->table->at(depth) is the depth+1 list of pt4D refinement boxes
		if(depth >= ref->table->size()){
			//if the depth is more than the deepest level of refinement table
			
			if(current_cell->intersect_candidates != NULL){
				cut_cell_nodes.push_back(current_node);	//if has intersections, add to cut_cells
			}
			else{
				non_cut_cell_nodes.push_back(current_node);
				//check if inside geometry
				/*
				if( inside_geometry(current_cell, surface) ){
						//if the cell is inside
						//delete the precise vert and sub in zero, need dummy verts
						return NULL;
				}
				*/
			}
			
			return mesh; //make no modiification
		}		
		for(int j = 0; j < ref->table->at(depth)->size(); j++){	
			//table->at(depth)->at(j) iterates though all the pt4D of that level
			pt4D *ref_box = ref->table->at(depth)->at(j);			
			if(overlap(ref_box, current_bound)){	//if refine is needed
				//<---divide cell in 4
				double med = current_cell->length / 2.0;
					//<--- divide cell in 2
					pt4D *l_bound = new pt4D(current_bound->x_min, current_bound->x_max - med,
											 current_bound->y_min, current_bound->y_max);
					pt4D *r_bound = new pt4D(current_bound->x_min + med, current_bound->x_max,
											 current_bound->y_min, current_bound->y_max);
					Polyhedron *rect1 = new Polyhedron();
					Polyhedron *rect2 = new Polyhedron();
					rect1->bounding_box = l_bound;
					rect2->bounding_box = r_bound;
					mesh->L_branch = new ADTNode(rect1, depth, l_bound);
					mesh->R_branch = new ADTNode(rect2, depth, r_bound);
					// end of divide cell in 2 --->
				// 2 4
				// 1 3
				pt4D *b_c1 = new pt4D(l_bound->x_min, l_bound->x_max,
											 l_bound->y_min, l_bound->y_max - med);
				pt4D *b_c2 = new pt4D(l_bound->x_min, l_bound->x_max,
											 l_bound->y_min + med, l_bound->y_max);
				pt4D *b_c3 = new pt4D(r_bound->x_min, r_bound->x_max,
											 r_bound->y_min, r_bound->y_max - med);
				pt4D *b_c4 = new pt4D(r_bound->x_min, r_bound->x_max,
											 r_bound->y_min + med, r_bound->y_max);
				//current_cell->pos;
				int newID = precise_verts.size();	//ID of the next added entry
				vert_precise *v2 = new vert_precise(precise_verts.at(current_cell->pos)->x,
													precise_verts.at(current_cell->pos)->y + med, newID);
				vert_precise *v3 = new vert_precise(precise_verts.at(current_cell->pos)->x + med,
													precise_verts.at(current_cell->pos)->y, newID+1);
				vert_precise *v4 = new vert_precise(precise_verts.at(current_cell->pos)->x + med,
													precise_verts.at(current_cell->pos)->y + med, newID+2);
				precise_verts.push_back(v2);
				precise_verts.push_back(v3);
				precise_verts.push_back(v4);
				CartCell* c1 = new CartCell(b_c1, med, current_cell->pos, *INTERNAL);
				CartCell* c2 = new CartCell(b_c2, med, newID, *INTERNAL);
				CartCell* c3 = new CartCell(b_c3, med, newID+1, *INTERNAL);
				CartCell* c4 = new CartCell(b_c4, med, newID+2, *INTERNAL);
				mesh->L_branch->L_branch = new ADTNode(c1, depth+1, b_c1);
				mesh->L_branch->R_branch = new ADTNode(c2, depth+1, b_c2);
				mesh->R_branch->L_branch = new ADTNode(c3, depth+1, b_c3);
				mesh->R_branch->R_branch = new ADTNode(c4, depth+1, b_c4);
				//recursive call to refine the children
				mesh->L_branch->L_branch = mesh_refine(surface, ref, mesh->L_branch->L_branch, depth+1);
				mesh->L_branch->R_branch = mesh_refine(surface, ref, mesh->L_branch->R_branch, depth+1);
				mesh->R_branch->L_branch = mesh_refine(surface, ref, mesh->R_branch->L_branch, depth+1);
				mesh->R_branch->R_branch = mesh_refine(surface, ref, mesh->R_branch->R_branch, depth+1);
				break; //only needs to get refined once
			}
		}
	}
	else{
		//directly modified, return value ignored
		mesh_refine(surface, ref, mesh->L_branch, depth);
		mesh_refine(surface, ref, mesh->R_branch, depth);
	}
	return mesh;
}

//#initial refinement for Adaptive Mode
ADTNode* mesh_refine(ADTNode *surface, bool (*ref_cond)(ADTNode*, void*), ADTNode *mesh, int depth, void *ref_rule){
	//helper for adaptive (automatic) version of generate_mesh
	ADTNode* current_node = mesh;
	if(current_node == NULL){
		return mesh;	//base case
	}
	if((current_node->L_branch == NULL) && (current_node->R_branch == NULL)){	//if is leaf
		CartCell *current_cell = static_cast<CartCell*>(mesh->poly);
		pt4D* current_bound = current_cell->bounding_box;	//convenience pointers
		//no more divide, find intersection and generate verts and edges for this cell
		
		//cout << "the cell_bound in mesh_refine" << endl;
		//cout << *current_bound << endl;
		
		current_cell->intersect_candidates = find_intersections(current_bound, surface, overlap_nonstrict);
		
		if(!ref_cond(current_node, ref_rule)){
			if(current_cell->intersect_candidates != NULL){
				cut_cell_nodes.push_back(current_node);	//if has intersections, add to cut_cells
			}
			else{
				non_cut_cell_nodes.push_back(current_node);	//if not, add to non_cut_cells
			}			
			return mesh; //make no modiification
		}
		
		refine(current_node, false);
		//recursive call to refine the children
		mesh->L_branch->L_branch = mesh_refine(surface, ref_cond, mesh->L_branch->L_branch, depth+1, ref_rule);
		mesh->L_branch->R_branch = mesh_refine(surface, ref_cond, mesh->L_branch->R_branch, depth+1, ref_rule);
		mesh->R_branch->L_branch = mesh_refine(surface, ref_cond, mesh->R_branch->L_branch, depth+1, ref_rule);
		mesh->R_branch->R_branch = mesh_refine(surface, ref_cond, mesh->R_branch->R_branch, depth+1, ref_rule);
	}
	else{
		//directly modified, return value ignored
		mesh_refine(surface, ref_cond, mesh->L_branch, depth, ref_rule);
		mesh_refine(surface, ref_cond, mesh->R_branch, depth, ref_rule);
	}
	return mesh;
}

//resets the is_walked status of all edges in a cell_piece, preparing for a walk starting at prestart
void reset_edges(cell_piece *piece, vert_intersection *prestart){
	for(int i = 0; i < *(piece->num_indices[0]); i++){
		if(prestart != NULL){
			if( is_equal(precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt1), 
										  precise_verts.at(prestart->host->pt1)) && 
				is_equal(precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt2),
										  precise_verts.at(prestart->host->pt2)) ){					  
				precise_edges.at(piece->list_of_edges[i])->walked = true;
			}
			continue;
		}
		precise_edges.at(piece->list_of_edges[i])->walked = false;
	}
	return;
}

//walk algorithm to find path from a point intersection to another
void walk(vert_intersection *start, cell_piece *piece, vector<vert_intersection*> *pts,
		  vector<int> *edge_bound, vector<int> *vert_bound,
		  bool (*cond_func)(vert_precise*, void*), 
		  vert_precise *terminate, edge* ref_edge){
	//terminate is only used with cond_closed	
	if(terminate == NULL){
		
		//cout << "terminate is NULL" << endl;
		
		if(cond_func(start, terminate)){
			//walk completed
			//cout << "terminating cond place 1 satisified" << endl;
			
			//cout << "--@@@@@@@@@@@@@@--" << endl;
			//edge_stat(start->host);			
			//cout << *start << endl;
			//cout << "--@@@@@@@@@@@@@@--" << endl;
			//cout << "walk complete" << endl;
			return;
		}
	
		//finding cutting boundary, no need to divide intersection_pts
		vector<Polyhedron*> *temp = piece->intersect_segments;
		vector<vert_intersection*> *temp_vert_vec = piece->list_of_intersecton_points;
		for(int i = 0, cap = temp->size(); i < cap; ++i, cap = temp->size()){
		
			if(temp->at(i) == NULL){
				continue;
			}
			
			//find the next point in the chain
			if(is_pt_on_line(start, static_cast<Segment*>(temp->at(i)))){
				//cout << "not a lone point" << endl;
				
				vert_intersection *next;
				int lst_index;
				
				if(is_equal(start, precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0]))){
					next = static_cast<vert_intersection*>(precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[1]));		
					lst_index = 1;
				}
				else{
					next = static_cast<vert_intersection*>(precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0]));
					lst_index = 0;
				}
				
				//cout << "next vert before cutting" << endl;
				//cout << *next << endl;
				
				
				if(next->status == OUT){
					//cout << "next is out" << endl;
					//if OUT, find the onbound vert_intersection
					for(int j = 0, cap_v = temp_vert_vec->size(); j < cap_v; ++j, cap_v = temp_vert_vec->size()){
						if(is_pt_on_line(temp_vert_vec->at(j), static_cast<Segment*>(temp->at(i)))){
							next = temp_vert_vec->at(j);
							
							remove_pt(next, piece);
							break;
						}
					}
				}
				
				//walk until reaches boundary
				if(cond_func(next, terminate)){
					//walk completed
					//cout << "terminating cond place 2 satisified" << endl;
					//cout << *(next) << endl;
					//edge_stat(next->host);
					
					//cout << "--@@@@@@@@@@@@@@--" << endl;
					//cout << *(precise_verts.at(start->ID)) << endl;
					//cout << *(precise_verts.at(next->ID)) << endl;
					//cout << "--@@@@@@@@@@@@@@--" << endl;
									
					int pos_e = precise_edges.size();
					body_edge *e_last = new body_edge(start->ID, next->ID, NULL, false);
					precise_edges.push_back(e_last);
					
					edge_bound->push_back(pos_e);
					vert_bound->push_back(next->ID);
					temp->at(i) = NULL;
					temp->erase(temp->begin() + i);
										
					//cout << "walk complete" << endl;
					return;
				}
				
				//add found points to list
				vert_bound->push_back(next->ID);
				edge_bound->push_back(static_cast<Segment*>(temp->at(i))->list_of_edges[0]);
				
				//remove this segment from pool
				temp->at(i) = NULL;
				temp->erase(temp->begin() + i);
				i--;
				
				//cout << "before recursive call on walk" << endl;
				walk(next, piece, pts, edge_bound, vert_bound, cond_func, terminate, ref_edge);
				//cout << "after recursive call on walk" << endl;
				
				return;
			}
		}
	}
	else{
		//closing the polygon, iterate thru the edges
		
		//cout << "terminate is not NULL" << endl;
		//cout << "$$$$$$" << endl;
		//cout << "current node in walk is: " << endl;
		//cout << *(static_cast<vert_intersection*>(start)) << endl;
		//cout << "end node for walk is: " << endl;
		//cout << *(static_cast<vert_intersection*>(terminate)) << endl;
		
		if(cond_func(terminate, ref_edge)){
			//walk completed
			//cout << "terminating cond satisified place 1" << endl;
			//cout << "~~~~~~~~~~~~~~~~~walk complete~~~~~~~~~~~~~~~~~" << endl;
			return;
		}
	
		//walk on the outer walls of the cell/cell_piece
		for(int i = 0; i < *(piece->num_indices[0]); i++){
			body_edge *temp_eg = static_cast<body_edge*>(precise_edges.at(piece->list_of_edges[i]));
			
			//cout << "edge during walk" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt1)->x << ", "
			//	 << precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt1)->y << ")" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt2)->x << ", "
			//	 << precise_verts.at(precise_edges.at(piece->list_of_edges[i])->pt2)->y << ")" << endl;
				 
			//cout << "temp_eg walked?: " << temp_eg->walked << " " << is_pt_on_line(start, temp_eg) << endl;
						 
			//if temp_eg is not walked yet
			if((!temp_eg->walked) && is_pt_on_line(start, temp_eg)){
				//cout << "not a lone point" << endl;
				
				vert_intersection *next;
				if(is_equal(start, precise_verts.at(temp_eg->pt1))){
					next = static_cast<vert_intersection*>(precise_verts.at(temp_eg->pt2));
				}
				else{
					next = static_cast<vert_intersection*>(precise_verts.at(temp_eg->pt1));
				}
				
				body_edge *ref_edge = temp_eg;
				//cout << "&&&&&&&&&&&&ref_edge will be: " << endl;
				//cout << "(" << precise_verts.at(start->ID)->x << ", "
				// << precise_verts.at(start->ID)->y << ")" << endl;
				//cout << "(" << precise_verts.at(next->ID)->x << ", "
				// << precise_verts.at(next->ID)->y << ")" << endl;
				//cout << "&&&&&&&&&&&&" << endl;
				
				//cout << "end node for walk is: " << endl;
				//cout << *(static_cast<vert_intersection*>(terminate)) << endl;
				
				//once back to the starting point, close off and get one polygon
				if(cond_func(terminate, ref_edge)){
					//walk completed
					//cout << "terminating cond satisified place 2" << endl;
					
					//close the polygon
					int pos_e = precise_edges.size();
					body_edge *e_last = new body_edge(start->ID, terminate->ID, NULL, ref_edge->is_cart);
					precise_edges.push_back(e_last);
					next->host = e_last;
					edge_bound->push_back(pos_e);
					vert_bound->push_back(start->ID);
										
					//find intersection pts for last edge
					
					for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
						itr != piece->list_of_intersecton_points->end(); ++itr){
					//for(vector<vert_intersection*>::iterator itr = temp_eg->list_of_intersecton_points->begin(); 
					//	itr != temp_eg->list_of_intersecton_points->end(); ++itr){
						if(is_pt_on_line(*itr, e_last)){
							pts->push_back(*itr);
						}
					}			
					//cout << "~~~~~~~~~~~~~~~~~walk complete~~~~~~~~~~~~~~~~~" << endl;
					//delete ref_edge;
					return;
				}
				
				if( !is_equal(precise_verts.at(vert_bound->back()), start) ){
					edge_bound->push_back(piece->list_of_edges[i]);
					vert_bound->push_back(start->ID);
				}
								
				//takes intersection pts
				for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
						itr != piece->list_of_intersecton_points->end(); ++itr){
				//for(vector<vert_intersection*>::iterator itr = temp_eg->list_of_intersecton_points->begin(); 
				//	itr != temp_eg->list_of_intersecton_points->end(); ++itr){
					if(is_pt_on_line(*itr, temp_eg)){
						pts->push_back(*itr);
					}
				}
				
				//mark the edge as walked
				temp_eg->walked = true;

				walk(next, piece, pts, edge_bound, vert_bound, cond_func, terminate, ref_edge);
				//delete ref_edge;
				return;
			}
		}
	}
	//not found
	//cout << "err: lone point" << endl;
}

//divides the cell into different piece using the paths that cut the cell
void divide_cell_piece(cell_piece *piece, stack<cell_piece*> *parts){

	/*//cout << "%%%%%%%%%%%%%%%%  PIECE's list_of_intersecton_points list  %%%%%%%%%%%%%%%%%" << endl;
		for(int j = 0; j < piece->list_of_intersecton_points->size(); j++){
			//cout << *(piece->list_of_intersecton_points->at(j)) << endl;
		}
	//cout << "%%%%%%%%%%%%%%%%  END OF PIECE's list_of_intersecton_points list  %%%%%%%%%%%%%%%%%" << endl;
	//cout << "%%%%%%%%%%%%%%%%  PIECE's intersect_segments list  %%%%%%%%%%%%%%%%%" << endl;
		for(int i = 0; i < piece->intersect_segments->size(); i++){	
			edge_stat(static_cast<Segment*>(piece->intersect_segments->at(i)));
		}
	//cout << "%%%%%%%%%%%%%%%%  PIECE's intersect_segments list  %%%%%%%%%%%%%%%%%" << endl;*/
	
	removeDuplicates(piece->list_of_intersecton_points);
	
	if(piece->list_of_intersecton_points->size() < 2){
		//add the shape to parts
		//cout << "part was added" << endl;
		parts->push(piece);
		return;
	}
	
	//else divide the cell_piece
	
	//make backup of intersect segment
	vector<Polyhedron*> *intersect_segments_backup1 = new vector<Polyhedron*>;
	vector<Polyhedron*> *intersect_segments_backup2 = new vector<Polyhedron*>;
	for(int i = 0; i < piece->intersect_segments->size(); i++){	
		intersect_segments_backup1->push_back(piece->intersect_segments->at(i));
		intersect_segments_backup2->push_back(piece->intersect_segments->at(i));
	}
	
	
	//make a boundary container (vector of edgeIDs)
	vector<int> *edge_bound1 = new vector<int>;
	vector<int> *vert_bound1 = new vector<int>;
	
	vert_intersection *start = piece->list_of_intersecton_points->front();

	//use first point, toss out the OUT point
	
	//cout << "@@@@@--starting walk for dividing bound--@@@@@" << endl;
	
	vector<Polyhedron*> *temp = piece->intersect_segments;
	for(int i = 0, cap = temp->size(); i < cap; ++i, cap = temp->size()){	
		
		if(temp->at(i) == NULL){
			continue;
		}
		
		//cout << "being compared:" << endl;
		//cout << *(start) << endl;
		//cout << "vert1 " << precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0])->x << " "
		//	 << precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0])->y << endl;
		//cout << "vert2 " << precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[1])->x << " "
		//	 << precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[1])->y << endl;
	
		if(is_pt_on_line(start, static_cast<Segment*>(temp->at(i)))){
		
			//cout << "point is found on segment: " << i << endl;
			
			vert_intersection *next;
			vert_intersection *vert1 = static_cast<vert_intersection*>
										 (precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0]));
			vert_intersection *vert2 = static_cast<vert_intersection*>
										 (precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[1]));
										 
			
			//cout << ">>>>>>>>vert1<<<<<<<<" << endl;
			//cout << *vert1;
			//cout << ">>>>>>>>vert2<<<<<<<<" << endl;
			//cout << *vert2;
										 
			//checking if the intersection is trivial							
			if( (vert1->status != IN) && (vert2->status != IN) ){
			
				//if the number of vert intersection on this segment is less than 2, trivial
				body_edge *check = trivial_check(piece, static_cast<Segment*>(temp->at(i)));
				if(check != NULL){
					vert_bound1->push_back(check->pt1);
					vert_bound1->push_back(check->pt2);
					edge_bound1->push_back(precise_edges.size());
					precise_edges.push_back(check);
					next = static_cast<vert_intersection*>(precise_verts.at(check->pt2));
					remove_pt(next, piece);
				}
				else{
					//trivial intersection
					//cout << "%%% trivial intersection %%%" << endl;
					temp->at(i) = NULL;
					temp->erase(temp->begin() + i);
					i--;
					continue;
				}
			}
			
			int pos_v = precise_verts.size();
			int pos_e = precise_edges.size();
			vert_intersection *v_first = new vert_intersection(start->x, start->y, NULL, ON_BOUND, false, precise_verts.size());
			precise_verts.push_back(v_first);
			
			//find the next vertice, whichever's inside the cell out of the two on the edge hosting the intersection point
			if(vert1->status == IN){
				//0 is in, add start<->0 edge, add to edge_bound
				body_edge *e_first = new body_edge(pos_v, static_cast<Segment*>(temp->at(i))->list_of_verts[0], NULL, false);
				precise_edges.push_back(e_first);
				v_first->host = e_first;
				edge_bound1->push_back(pos_e);
				vert_bound1->push_back(pos_v);
				vert_bound1->push_back(static_cast<Segment*>(temp->at(i))->list_of_verts[0]);
				next = vert1;
			}
			else if(vert2->status == IN){
				//1 is in, add start<->1 edge, add to boundary
				body_edge *e_first = new body_edge(pos_v, static_cast<Segment*>(temp->at(i))->list_of_verts[1], NULL, false);
				precise_edges.push_back(e_first);
				v_first->host = e_first;
				edge_bound1->push_back(pos_e);
				vert_bound1->push_back(pos_v);
				vert_bound1->push_back(static_cast<Segment*>(temp->at(i))->list_of_verts[1]);
				next = vert2;
			}
			
			temp->at(i) = NULL;
			temp->erase(temp->begin() + i);
			i--;
			
			//cout << ">>>>>>>>next<<<<<<<<" << endl;
			//cout << *next;
			
			if( (piece->list_of_intersecton_points->size() == 1) && (next->status != ON_BOUND) ){
				//cout << "part was added at place 2" << endl;
				parts->push(piece);
				return;
			}
			
			//call walk once to find boundary (cond_onbound)
			//boundary that divides the cell/cell_piece
			walk(next, piece, NULL, edge_bound1, vert_bound1, cond_onbound, NULL, NULL);
			break;
		}
	}
	
	//duplicate the boundary (container)
	if(vert_bound1->size() == 0){
		//cout << "vert_bound1 size is 0, original" << endl;
		parts->push(piece);
		return;
	}
	
	int pos_e;
	vert_intersection *end = static_cast<vert_intersection*>(precise_verts.at(vert_bound1->back()));
	vector<int> *old_vert_bound;
	vector<int> *old_edge_bound;
	
	
	vector<int> *new_verts;
	vector<int> *new_edges = new vector<int>;
	
	old_vert_bound = vert_bound1;
	old_edge_bound = edge_bound1;
		
	if(vert_bound1->size() > boundary_cap){
		new_verts = max_dist_vert(vert_bound1, start, end, boundary_cap);
		pos_e = precise_edges.size();
		
		//generate edges
		for(int i = 0; i < (new_verts->size() - 1); i++){
			body_edge *eg = new body_edge(new_verts->at(i), new_verts->at(i+1), NULL, false);
			precise_edges.push_back(eg);
			print_boundary.push_back(eg);
			new_edges->push_back(pos_e+i);
		}
		
		vert_bound1 = new_verts;
		edge_bound1 = new_edges;
	}
	else{
		for(int i = 0; i < edge_bound1->size(); i++){
			print_boundary.push_back(precise_edges.at(edge_bound1->at(i)));
		}
	}
	
	/*
	if(vert_max_dist != -1){
		//if more than 2 lines, replace with 2 lines
		vector<int> *new_verts = new vector<int>;
		vector<int> *new_edges = new vector<int>;
		pos_e = precise_edges.size();
		body_edge *e1 = new body_edge(start->ID, vert_max_dist, NULL, false);
		body_edge *e2 = new body_edge(vert_max_dist, end->ID, NULL, false);
		precise_edges.push_back(e1);
		precise_edges.push_back(e2);
		
		print_boundary.push_back(e1);
		print_boundary.push_back(e2);
		
		new_verts->push_back(start->ID);
		new_verts->push_back(vert_max_dist);
		new_verts->push_back(end->ID);
		new_edges->push_back(pos_e);
		new_edges->push_back(pos_e+1);
		
		old_vert_bound = vert_bound1;
		old_edge_bound = edge_bound1;
		
		vert_bound1 = new_verts;
		edge_bound1 = new_edges;
	}
	else{
		//cout << "ERROR" << endl;
	}
	*/
	
	
	//cout << "@@@@@--starting walk on the cell_piece boundary--@@@@@" << endl;
	//walk to complete the polygon
	
	vector<int> *vert_bound2 = new vector<int>;
	vector<int> *edge_bound2 = new vector<int>;
	for(vector<int>::iterator itr = vert_bound1->begin(); itr != vert_bound1->end(); ++itr){
		vert_bound2->push_back(*itr);
	}
	for(vector<int>::iterator itr = edge_bound1->begin(); itr != edge_bound1->end(); ++itr){
		edge_bound2->push_back(*itr);
	}
	
	//cout << "last vert added was: " << endl;
	//cout << *(precise_verts.at(vert_bound1->back())) << endl;
	
	vector<vert_intersection*> *pts1 = new vector<vert_intersection*>;
	vector<vert_intersection*> *pts2 = new vector<vert_intersection*>;
	
	//remove the vert_instersection from detection list
	remove_pt(start, piece);
	
	
	
	//call walk again twice to find two polygons (is_equal)
	start = static_cast<vert_intersection*>(precise_verts.at(vert_bound1->back()));	
	end = static_cast<vert_intersection*>(precise_verts.at(vert_bound1->front()));

	vert_intersection *path1;
	vert_intersection *path2;
	int part_elenum = *(piece->num_indices[0]);
	
	bool same_side = false;
	vert_intersection *ss_vert;
	body_edge *ss_edge;
	vector<vert_intersection*> *ss_pts;
	
	//check if there are two on a side
	for(int i = 0; i < part_elenum; i++){
		body_edge *temp_beg = static_cast<body_edge*>(precise_edges.at(piece->list_of_edges[i]));
		if(is_pt_on_line(start, temp_beg) && is_pt_on_line(end, temp_beg)){
			same_side = true;
			path1 = end;
			
			//make end<->pt1, if start is NOT on, use pt1, else use pt2
			edge test = edge(end->ID, temp_beg->pt1);
			if(!is_pt_on_line(start, &test)){
				path2 = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt2));
				ss_vert = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt1));
			}
			else{
				path2 = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt1));
				ss_vert = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt2));
			}
			
			ss_pts = new vector<vert_intersection*>;
			ss_edge = new body_edge(ss_vert->ID, end->ID, ss_pts, temp_beg->is_cart);
			
			for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
				itr != piece->list_of_intersecton_points->end(); ++itr){
				if(is_pt_on_line(*itr, ss_edge)){
					ss_pts->push_back(*itr);
				}
			}
			
			break;
		}
	}
		
	
	//if not, use next's edge
	if(!same_side){
		for(int i = 0; i < part_elenum; i++){
			body_edge *temp_beg = static_cast<body_edge*>(precise_edges.at(piece->list_of_edges[i]));
			if(is_pt_on_line(start, temp_beg)){
				path1 = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt1));
				path2 = static_cast<vert_intersection*>(precise_verts.at(temp_beg->pt2));
				temp_beg->walked = true;
				break;
			}
		}
	}
	
	//cout << "checking the paths" << endl;
	//cout << "size vert_bound1: " << vert_bound1->size() << endl;
	//cout << "size edge_bound1: " << edge_bound1->size() << endl;
	/*
	if( (vert_bound1->size() - edge_bound1->size()) != 1 ){
		//cout << "err: sizes are wrong" << endl;
		
		for(int i = 0; i < vert_bound1->size(); i++){
			//cout << *(precise_verts.at(vert_bound1->at(i))) << endl;
		}
	}
	
	
	for(int i = 0; i < vert_bound1->size(); i++){
			//cout << *(precise_verts.at(vert_bound1->at(i))) << endl;
	}
	
	for(int i = 0; i < edge_bound1->size(); i++){
			edge_stat(precise_edges.at(edge_bound1->at(i)));
	}
	*/
	//cout << *path1 << endl;
	//cout << *path2 << endl;
	
	
	vert_intersection *prestart = new vert_intersection(start->x, start->y, start->host, 
														start->status, start->on_cell_bound, start->ID);
	
	//make edge start<->path1, start<->path2
	pos_e = precise_edges.size();
	
	if(!is_equal( precise_verts.at(vert_bound1->back()), precise_verts.at(path1->ID) )){	
		body_edge *e_mid1 = new body_edge(vert_bound1->back(), path1->ID, pts1, static_cast<body_edge*>(prestart->host)->is_cart);
		precise_edges.push_back(e_mid1);
		static_cast<vert_intersection*>(precise_verts.at(vert_bound1->back()))->host = e_mid1;
		
		//edge_stat(e_mid1);
		//edge_stat(precise_edges.at(pos_e));
		edge_bound1->push_back(pos_e);
		
		if(same_side){
			//path1
			e_mid1->is_cart = true;
		}
		
		//find intersection pts for mid (first) edge
		for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
			itr != piece->list_of_intersecton_points->end(); ++itr){
			if(is_pt_on_line(*itr, e_mid1)){
				pts1->push_back(*itr);
			}
		}
	}
	else{
		//pop the last vert
		vert_bound1->pop_back();
	}
	
	if(!is_equal( precise_verts.at(vert_bound2->back()), precise_verts.at(path2->ID) )){	
		body_edge *e_mid2 = new body_edge(vert_bound2->back(), path2->ID, pts2, static_cast<body_edge*>(prestart->host)->is_cart);
		precise_edges.push_back(e_mid2);
		static_cast<vert_intersection*>(precise_verts.at(vert_bound2->back()))->host = e_mid2;
		
		//edge_stat(e_mid2);
		//edge_stat(precise_edges.at(pos_e+1));
		edge_bound2->push_back(pos_e+1);
		
		//find intersection pts for mid (first) edge
		for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
			itr != piece->list_of_intersecton_points->end(); ++itr){
			if(is_pt_on_line(*itr, e_mid2)){
				pts2->push_back(*itr);
			}
		}
	}
	else{
		//pop the last vert
		vert_bound2->pop_back();
	}	
		
	//cout << "before walking" << endl;
	//cout << "size vert_bound1: " << vert_bound1->size() << endl;
	//cout << "size edge_bound1: " << edge_bound1->size() << endl;
	
	//cout << *prestart << endl;
	//edge_stat(prestart->host);
	
	reset_edges(piece, prestart);
	walk(path1, piece, pts1, edge_bound1, vert_bound1, is_in, end, NULL);
	
	//cout << *prestart << endl;
	//edge_stat(prestart->host);
	
	reset_edges(piece, prestart);
	walk(path2, piece, pts2, edge_bound2, vert_bound2, is_in, end, NULL);
	
	if(same_side){
		//add the end edge and vert, path2
		pos_e = precise_edges.size();
		precise_verts.push_back(ss_vert);
		precise_edges.push_back(ss_edge);
		vert_bound2->push_back(ss_vert->ID);
		edge_bound2->push_back(pos_e);
	}
	
	delete prestart;
	
	//cout << "after walking" << endl;
	//cout << "size vert_bound1: " << vert_bound1->size() << endl;
	//cout << "size edge_bound1: " << edge_bound1->size() << endl;
	
	/*
	//close the polygon
	pos_e = precise_edges.size();
	body_edge *e_last1 = new body_edge(vert_bound1->back(), vert_bound1->front(), NULL, false);
	precise_edges.push_back(e_last1);
	static_cast<vert_intersection*>(precise_verts.at(vert_bound1->back()))->host = e_last1;
	edge_bound1->push_back(pos_e);
	
	body_edge *e_last2 = new body_edge(vert_bound2->back(), vert_bound2->front(), NULL, false);
	precise_edges.push_back(e_last2);
	static_cast<vert_intersection*>(precise_verts.at(vert_bound2->back()))->host = e_last2;
	edge_bound2->push_back(pos_e+1);
	
	//find intersection pts for last edge
	for(vector<vert_intersection*>::iterator itr = piece->list_of_intersecton_points->begin(); 
		itr != piece->list_of_intersecton_points->end(); ++itr){
		if(is_pt_on_line(*itr, e_last1)){
			pts1->push_back(*itr);
		}
		if(is_pt_on_line(*itr, e_last2)){
			pts2->push_back(*itr);
		}
	}*/
	
	//make cell pieces out of these
	assert(edge_bound1->size() == vert_bound1->size());
	assert(edge_bound2->size() == vert_bound2->size());
	
	/*
	//cout << "~~~~!!!@!@#!@#!@~~~~piece 1: ~~~~!!!@!@#!@#!@~~~~" << endl;
	for(int temp_index = 0; temp_index < vert_bound1->size(); temp_index++){
		//cout << *(precise_verts.at(vert_bound1->at(temp_index))) << endl;
		edge_stat(precise_edges.at(edge_bound1->at(temp_index)));

		//print_edges.push_back(precise_edges.at(edge_bound1->at(temp_index)));
		
		//cout << endl;
	}
	
	//cout << "~~~~!!!@!@#!@#!@~~~~piece 2: ~~~~!!!@!@#!@#!@~~~~" << endl;
	for(int temp_index = 0; temp_index < vert_bound2->size(); temp_index++){
		//cout << *(precise_verts.at(vert_bound2->at(temp_index))) << endl;
		edge_stat(precise_edges.at(edge_bound2->at(temp_index)));
		
		//print_edges.push_back(precise_edges.at(edge_bound2->at(temp_index)));
		
		//cout << endl;
	}*/
	
	//prepare new fields for the children
	int *i1 = new int;
	int *i2 = new int;
	*i1 = vert_bound1->size();
	*i2 = vert_bound2->size();
	int **indice1 = new int*[1];
	int **indice2 = new int*[1];
	indice1[0] = i1;
	indice2[0] = i2;
	
	int *vert_list1 = new int[*(indice1[0])];
	int *edge_list1 = new int[*(indice1[0])];
	int *vert_list2 = new int[*(indice2[0])];
	int *edge_list2 = new int[*(indice2[0])];
	for(int i = 0; i < **indice1; i++){
		vert_list1[i] = vert_bound1->at(i);
		edge_list1[i] = edge_bound1->at(i);
	}
	for(int i = 0; i < **indice2; i++){
		vert_list2[i] = vert_bound2->at(i);
		edge_list2[i] = edge_bound2->at(i);
	}
		
	/*
	//cout << "before making children" << endl;
	//cout << piece->list_of_intersecton_points->size() << endl;
	
	for(int j = 0; j < piece->list_of_intersecton_points->size(); j++){
		//cout << *(piece->list_of_intersecton_points->at(j)) << endl;
	}
	
	//cout << "---" << endl;
	//cout << pts1->size() << endl;
	
	for(int j = 0; j < pts1->size(); j++){
		//cout << *(pts1->at(j)) << endl;
	}
	
	//cout << "@@@" << endl;
	//cout << pts2->size() << endl;
	
	for(int j = 0; j < pts2->size(); j++){
		//cout << *(pts2->at(j)) << endl;
	}
	*/
	
	cell_piece *child1 = new cell_piece(vert_list1, edge_list1, indice1, 
										pts1, intersect_segments_backup1);
	cell_piece *child2 = new cell_piece(vert_list2, edge_list2, indice2, 
										pts2, intersect_segments_backup2);	
	
	child1->old_vert_bound = old_vert_bound;
	child2->old_vert_bound = old_vert_bound;
	child1->old_edge_bound = old_edge_bound;
	child2->old_edge_bound = old_edge_bound;
	
	
	reset_edges(child1, NULL);
	reset_edges(child2, NULL);
	
	/*
	//cout << "%%%" << endl;
	//cout << child1->intersect_segments->size() << endl;
	for(int i = 0; i < child1->intersect_segments->size(); i++){	
		edge_stat(static_cast<Segment*>(child1->intersect_segments->at(i)));
	}
	
	//cout << child2->intersect_segments->size() << endl;
	for(int i = 0; i < child2->intersect_segments->size(); i++){	
		edge_stat(static_cast<Segment*>(child2->intersect_segments->at(i)));
	}
	*/
	
	//recursive calls
	divide_cell_piece(child1, parts);
	divide_cell_piece(child2, parts);
	
}

//makes a CutCell class from a processed cell_piece class, discarding inside parts
CutCell *make_cut_cell(cell_piece *piece, ADTNode *surface){
	stack<cell_piece*> *parts = new stack<cell_piece*>;
	
	//cout << "---------before dividing cell---------" << endl;
	divide_cell_piece(piece, parts);
	//cout << "---------after dividing cell---------" << endl;
	
	int *entries = new int;
	*entries = parts->size();
	int **indices = new int*[*entries];
	
	vector<int> *vert_vec = new vector<int>;
	vector<int> *edge_vec = new vector<int>;
		
	//assemble parts
	//check if in geometry during make, add everything at this stage
	int count = 0;
	
	//cout << "parts size: " << parts->size() << endl;
	
	while(!parts->empty()){
		//add each cell_piece to lists that makes a cut cell
		
		cell_piece *temp = parts->top();
		
		int part_elenum = *(temp->num_indices[0]);
		bool inside = false;
		
		
		for(int i = 0; i < part_elenum; i++){
			//vert_intersection *temp_interv = static_cast<vert_intersection*>(precise_verts.at(temp->list_of_verts[i]));
			body_edge *temp_beg = static_cast<body_edge*>(precise_edges.at(temp->list_of_edges[i]));
					
			if(temp_beg->is_cart){
				
				//cout << "testing edge:" << endl; 
				//edge_stat(temp_beg);
				
				vert_precise *midpoint = find_midpoint( precise_verts.at(temp_beg->pt1),
														precise_verts.at(temp_beg->pt2) );
														
				//cout << *midpoint << endl;
				
				//if inside, discard
				if(inside_geometry(midpoint, surface)){
					inside = true;
					delete midpoint;
					//cout << "midpoint set inside" << endl;
					break;
				}
				delete midpoint;
				
				////cout << "inside is currently: " << inside << endl;
			}
		}
		
		//cout << "---cell piece info---" << endl;
		cell_piece_info(temp);
		//cout << inside << endl;
		
		if(!inside){
		
			//cout << "element added to print" << endl;
			cell_piece_info(temp);
			//cout << "element was added to print" << endl;
			
			//keep track of the sizes of individual pieces
			indices[count] = temp->num_indices[0];
			count++;
			
			//putting it together
			for(int i = 0; i < part_elenum; i++){
				vert_vec->push_back( temp->list_of_verts[i] );
				edge_vec->push_back( temp->list_of_edges[i] );
			}
		}
		
		parts->top() = NULL;
		parts->pop();
	}
	
	//make and return the cut cell
	*entries = count;
	int* vert_lst = new int[vert_vec->size()];
	int* edge_lst = new int[edge_vec->size()];
	
	for(int i = 0; i < vert_vec->size(); i++){
    	vert_lst[i] = vert_vec->at(i);
    }
    
    for(int i = 0; i < edge_vec->size(); i++){
    	edge_lst[i] = edge_vec->at(i);
	}
	
	CutCell *result = new CutCell(vert_lst, edge_lst, entries, indices);
	
	delete vert_vec;
	delete edge_vec;
	
	//cout << "---CutCell info---" << endl;
	cell_piece_info(result);
	
	return result;
}

//refinement condition for Min_size, cells satisifying this condition are divided
bool ref_cond_minsize(ADTNode *node, void* min_size_ref_in){
	CartCell *cell = static_cast<CartCell*>(node->poly);
	vector<Minsize_Refinement*> *min_size_ref = static_cast<vector<Minsize_Refinement*>*>(min_size_ref_in);
	
	pt4D test_box = pt4D(precise_verts.at(cell->pos)->x, precise_verts.at(cell->pos)->x + cell->length,
						 precise_verts.at(cell->pos)->y, precise_verts.at(cell->pos)->y + cell->length);
	
	for(int i = 0; i < min_size_ref->size(); i++){
		if(overlap_nonstrict(min_size_ref->at(i)->box, &test_box) && (cell->length > min_size_ref->at(i)->min_size)){
			return true;
		}
	}
	return false;
}

//checks if the curvature within a cell is greater than the input parameter ALPHA
bool check_cell_curvature(CartCell *cell, Adaptive_variable *adaptive_var){
	//cout << "check_cell_curvature called" << endl;	
	
	
	//bool only_holes = true;
	
	vector<Polyhedron*> *new_candidates = new vector<Polyhedron*>;
	for(int i = 0; i < cell->intersect_candidates->size(); i++){
		new_candidates->push_back(cell->intersect_candidates->at(i));
		
		/*
		//a line is contain if both points are contained
		if(!contained_strict(precise_verts.at(static_cast<Segment*>(cell->intersect_candidates->at(i))->list_of_verts[0]), 
							cell->bounding_box) ||
		   !contained_strict(precise_verts.at(static_cast<Segment*>(cell->intersect_candidates->at(i))->list_of_verts[1]), 
							cell->bounding_box) ){
			only_holes = false;
		}
		*/
	}
	
	/*
	if(only_holes){
		return true;
	}
	*/	
	
	//cell preparation, replace vertices and edges with vert_intersection and body_edge classes
	//removes trivial intersection candidates
	generate_cell_lists(cell);
	vector<vert_intersection*> *full_list_intersect_pts = new vector<vert_intersection*>;
	
	for(int j = 0; j < CARTCELL_NUM_INDICES; j++){
		vector<vert_intersection*> *list_intersect_pts = new vector<vert_intersection*>;
	
		for(int i = 0; i < new_candidates->size(); i++){
			
			if(new_candidates->at(i) == NULL){
				continue;
			}
			
			Segment *seg_temp = static_cast<Segment*>(new_candidates->at(i));
			
			pt4D *cell_bound = new pt4D(precise_verts.at(cell->pos)->x,
										precise_verts.at(cell->pos)->x + cell->length,
										precise_verts.at(cell->pos)->y,
										precise_verts.at(cell->pos)->y + cell->length);
			
			if(!contained(cell->bounding_box, seg_temp->bounding_box)
			 && do_lines_intersect(seg_temp, precise_edges.at(cell->list_of_edges[j]))){
			 
				//the lines intersect, create vert_intersect
				edge *edge_temp = precise_edges.at(cell->list_of_edges[j]);
				double x1 = precise_verts.at(seg_temp->list_of_verts[0])->x;
				double y1 = precise_verts.at(seg_temp->list_of_verts[0])->y;
				double x2 = precise_verts.at(seg_temp->list_of_verts[1])->x;
				double y2 = precise_verts.at(seg_temp->list_of_verts[1])->y;
				vert_intersection *inter_vert;
				double x, y;
				
				if(precise_verts.at(edge_temp->pt1)->x == precise_verts.at(edge_temp->pt2)->x){
					//x is same, bound is vert line
					x = precise_verts.at(edge_temp->pt1)->x;
					if( fabs(x2 - x1) > *EPSILON ){
						//segment not vert line
						const double M = (y2 - y1) / (x2 - x1);	//slope
						const double B = -(M * x1) + y1;	//y-intercept
						y = (M * x + B);
						inter_vert = new vert_intersection(x, y, edge_temp, ON_BOUND, false, precise_verts.size());
						precise_verts.push_back(inter_vert);
						
						list_intersect_pts->push_back(inter_vert);
						full_list_intersect_pts->push_back(inter_vert);
					}
				}
				else{
					//y is same, bound is horz line
					y = precise_verts.at(edge_temp->pt1)->y;
					if( fabs(y2 - y1) > *EPSILON ){
						//segment not horz line			
						if( fabs(x2 - x1) <= *EPSILON ){
							//vert line, x from here, y from boundary
							x = x1;
						}
						else{
							const double M = (y2 - y1) / (x2 - x1);	//slope
							const double B = -(M * x1) + y1;	//y-intercept
							x = (y - B) / M;
						}
						inter_vert = new vert_intersection(x, y, edge_temp, ON_BOUND, false, precise_verts.size());
						precise_verts.push_back(inter_vert);
						
						list_intersect_pts->push_back(inter_vert);
						full_list_intersect_pts->push_back(inter_vert);
					}						
				}
			}
			//replace all the verts in segments with the subclass version (check in/out)
			vert_cond cond1 = is_inside(cell, precise_verts.at(seg_temp->list_of_verts[0]));
			vert_cond cond2 = is_inside(cell, precise_verts.at(seg_temp->list_of_verts[1]));
			
			bool b1 = (cond1 == ON_BOUND);
			bool b2 = (cond2 == ON_BOUND);		
			body_edge *eg = new body_edge(precise_edges.at(seg_temp->list_of_edges[0])->pt1, 
										  precise_edges.at(seg_temp->list_of_edges[0])->pt2, 
										  NULL, false);	
			vert_intersection *v1 = new vert_intersection(
											precise_verts.at(seg_temp->list_of_verts[0])->x,
											precise_verts.at(seg_temp->list_of_verts[0])->y,
											eg, cond1, b1, seg_temp->list_of_verts[0]);
			vert_intersection *v2 = new vert_intersection(
											precise_verts.at(seg_temp->list_of_verts[1])->x,
											precise_verts.at(seg_temp->list_of_verts[1])->y,
											eg, cond2, b2, seg_temp->list_of_verts[1]);	
													
			delete precise_verts.at(seg_temp->list_of_verts[0]);
			delete precise_verts.at(seg_temp->list_of_verts[1]);
			precise_verts.at(seg_temp->list_of_verts[0]) = v1;
			precise_verts.at(seg_temp->list_of_verts[1]) = v2;
		}

		//replace all the edges and verts in cartcell to subclass version
		vert_intersection *v_temp = new vert_intersection(
										precise_verts.at(cell->list_of_verts[j])->x,
										precise_verts.at(cell->list_of_verts[j])->y,
										NULL, ON_BOUND, true, 
										precise_verts.at(cell->list_of_verts[j])->ID);
		
		body_edge *e_temp = new body_edge(precise_edges.at(cell->list_of_edges[j])->pt1, 
										  precise_edges.at(cell->list_of_edges[j])->pt2, 
										  list_intersect_pts, true);
													  
		precise_verts.at(cell->list_of_verts[j]) = v_temp;
		precise_edges.at(cell->list_of_edges[j]) = e_temp;
		
		removeDuplicates(list_intersect_pts);
	}
	removeDuplicates(full_list_intersect_pts);
	//done preparing the cell

	//if no intersecting points but have intersecting segments, there is a hole
	if((full_list_intersect_pts->size() == 0) && (full_list_intersect_pts->size() != 0)){
		return true;
	}
		
	//base case
	if(new_candidates->size() < 2){
		//cout << "candidate size less than 2" << endl;
		return false;
	}
	
	/*
	if(new_candidates->size() > adaptive_var->max_intersections){
		return true;
	}*/
	
	cell_piece *piece = new cell_piece(cell->list_of_verts, cell->list_of_edges, 
										   cell->num_indices, full_list_intersect_pts, 
										   new_candidates);	   
	remove_trivial_seg(piece);
	
	//cout << "ref_cond_adaptive before dividing piece" << endl;
	
	double max_angle_x = -2.0 * PI;
	double min_angle_x = 2 * PI;
	double max_angle_y = -2.0 * PI;
	double min_angle_y = 2 * PI;
	
	while(piece->list_of_intersecton_points->size() > 1){	//while more boundaries are there
		
		//cout << "piece->list_of_intersecton_points->size() is: " 
		//	 << piece->list_of_intersecton_points->size() << endl;
		
		//find boundary
		vector<int> *edge_bound1 = new vector<int>;
		vector<int> *vert_bound1 = new vector<int>;
		vert_intersection *start = piece->list_of_intersecton_points->front();
			//use first point, toss out the OUT point
		vector<Polyhedron*> *temp = piece->intersect_segments;
		for(int i = 0, cap = temp->size(); i < cap; ++i, cap = temp->size()){	
			
			if(temp->at(i) == NULL){
				continue;
			}
			
			if(is_pt_on_line(start, static_cast<Segment*>(temp->at(i)))){
				vert_intersection *next;
				vert_intersection *vert1 = static_cast<vert_intersection*>
											 (precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[0]));
				vert_intersection *vert2 = static_cast<vert_intersection*>
											 (precise_verts.at(static_cast<Segment*>(temp->at(i))->list_of_verts[1]));
			
				if( (vert1->status != IN) && (vert2->status != IN) ){
					//if the number of vert intersection on this segment is less than 2, trivial
					body_edge *check = trivial_check(piece, static_cast<Segment*>(temp->at(i)));
					if(check != NULL){
						vert_bound1->push_back(check->pt1);
						vert_bound1->push_back(check->pt2);
						edge_bound1->push_back(precise_edges.size());
						precise_edges.push_back(check);
						next = static_cast<vert_intersection*>(precise_verts.at(check->pt2));
						//remove_pt(next, piece);
					}
					else{
						//trivial intersection
						temp->at(i) = NULL;
						temp->erase(temp->begin() + i);
						i--;
						continue;
					}
				}
				
				int pos_v = precise_verts.size();
				int pos_e = precise_edges.size();
				vert_intersection *v_first = new vert_intersection(start->x, start->y, NULL, ON_BOUND, false, precise_verts.size());
				precise_verts.push_back(v_first);
				
				if(vert1->status == IN){
					//0 is in, add start<->0 edge, add to edge_bound
					body_edge *e_first = new body_edge(pos_v, static_cast<Segment*>(temp->at(i))->list_of_verts[0], NULL, false);
					precise_edges.push_back(e_first);
					v_first->host = e_first;
					edge_bound1->push_back(pos_e);
					vert_bound1->push_back(pos_v);
					vert_bound1->push_back(static_cast<Segment*>(temp->at(i))->list_of_verts[0]);
					next = vert1;
				}
				else if(vert2->status == IN){
					//1 is in, add start<->1 edge, add to boundary
					body_edge *e_first = new body_edge(pos_v, static_cast<Segment*>(temp->at(i))->list_of_verts[1], NULL, false);
					precise_edges.push_back(e_first);
					v_first->host = e_first;
					edge_bound1->push_back(pos_e);
					vert_bound1->push_back(pos_v);
					vert_bound1->push_back(static_cast<Segment*>(temp->at(i))->list_of_verts[1]);
					next = vert2;
				}
				
				temp->at(i) = NULL;
				temp->erase(temp->begin() + i);
				i--;
				
				//call walk once to find boundary (cond_onbound)
				//cout << "before walk" << endl;
				walk(next, piece, NULL, edge_bound1, vert_bound1, cond_onbound, NULL, NULL);
				//cout << "after walk" << endl;
				break;
			}
		}
		
		vert_precise *a;
		vert_precise *b;
		
		if(vert_bound1->size() != 0){
			vert_intersection *end = static_cast<vert_intersection*>(precise_verts.at(vert_bound1->back()));
			
			//cout << "removing start" << endl;
			//cout << *start << endl;
			if(piece->list_of_intersecton_points->size() != 0){
				remove_pt(start, piece);
			}
			
			//cout << "removing end" << endl;
			//cout << *end << endl;
			
			//cout << piece->list_of_intersecton_points->size() << endl;
			if(piece->list_of_intersecton_points->size() != 0){
				remove_pt(end, piece);
			}
				
			a = precise_verts.at(vert_bound1->back());
			vert_bound1->pop_back();
		}
		else{
			//cout << "removing start" << endl;
			//cout << *start << endl;
			remove_pt(start, piece);
		}
		
		while(vert_bound1->size() != 0){
		
			//cout << "vert_bound1->size() is: " << vert_bound1->size() << endl;
			
			b = precise_verts.at(vert_bound1->back());
			
			double angle_x = angle_between(0, 0, 1, 0, a, b);
			double angle_y = angle_between(0, 0, 0, 1, a, b);
			
			/*
			if(angle_x < 0){
				angle_x = angle_x + PI;
			}
			if(angle_y < 0){
				angle_y = angle_y + PI;
			}
			if(is_slope_neg(eg)){
				angle_x = -angle_x;
				angle_y = -angle_y;
			}
			*/
			
			//cout << "---" << endl;
			//cout << "angle_x is: " << angle_x << endl;
			//cout << "angle_y is: " << angle_y << endl;
			//cout << "---" << endl;
			
			if(angle_x > max_angle_x){
				max_angle_x = angle_x;
			}
			if(angle_x < min_angle_x){
				min_angle_x = angle_x;
			}
			
			if(angle_y > max_angle_y){
				max_angle_y = angle_y;
			}
			if(angle_y < min_angle_y){
				min_angle_y = angle_y;
			}	
			
			a = b;
			vert_bound1->pop_back();
		}
		
		//cout << *cell->bounding_box << endl;
		//cout << "max_angle_x is: " << max_angle_x << endl;
		//cout << "min_angle_x is: " << min_angle_x << endl;
		//cout << "max_angle_y is: " << max_angle_y << endl;
		//cout << "min_angle_y is: " << min_angle_y << endl;		
		
		/*
			fmod((max_angle_x - min_angle_x), PI)
			fmod((max_angle_y - min_angle_y), PI)
		*/
		
		//cout << "delta_x is: " << fmod((max_angle_x - min_angle_x), 2*PI) << endl;
		//cout << "delta_y is: " << fmod((max_angle_y - min_angle_y), 2*PI) << endl;
		//cout << "alpha is: " << adaptive_var->alpha << endl;
		
		//if the angle is greater than alpha
		if( (fmod((max_angle_x - min_angle_x), 2*PI) > adaptive_var->alpha) ||
			(fmod((max_angle_y - min_angle_y), 2*PI) > adaptive_var->alpha) ){
			
			delete piece;
			return true;
		}
	}
	
	delete piece;
	return false;
}

//refinement condition for Adaptive, cells satisifying this condition are divided
bool ref_cond_adaptive(ADTNode *node, void* adaptive_var_in){
	CartCell *cell = static_cast<CartCell*>(node->poly);
	Adaptive_variable *adaptive_var = static_cast<Adaptive_variable*>(adaptive_var_in);
	
	int depth = log_base2(*INIT_MIN/cell->length);
	if(depth >= adaptive_var->max_refs){
		//cout << "max_ref reached" << endl;
		//cout << depth << " >= " << adaptive_var->max_refs << endl;
		return false;
	}
	
	if(cell->intersect_candidates == NULL){
		//cout << "trivial" << endl;
		return false;
	}
	
	if(cell->intersect_candidates != NULL){
		return check_cell_curvature(cell, adaptive_var);
	}
	
	return false;
}

//<---removeDuplicates for node vectors
bool is_node_equal(ADTNode *n1, ADTNode *n2){
	return (n1 == n2);
}
bool comp_cell(ADTNode *n1, ADTNode *n2){
	return (n1->depth > n2->depth);
}
void removeDuplicates(vector<ADTNode*> *vec){
   sort(vec->begin(), vec->end(), comp_cell);
   vec->erase(unique(vec->begin(), vec->end(), is_node_equal), vec->end());
}
//end of removeDuplicates for node vectors--->

//check if the cell associated with the node is a cutcell
void check_cutcell(ADTNode *node){
	vector<Polyhedron*> *temp = static_cast<CartCell*>(node->poly)->intersect_candidates;
	if(temp != NULL){
		cut_cell_nodes_backup.push_back(node);
	}
	else{
		non_cut_cell_nodes_backup.push_back(node);
	}
}

//reloads the proper cells into the lists of cutcells and non cutcells
void reload(){

	//cout << "cut_cell_nodes_backup: " << cut_cell_nodes_backup.size() << endl;
	
	for(int i = 0; i < cut_cell_nodes_backup.size(); i++){	
		cut_cell_nodes.push_back(cut_cell_nodes_backup.at(i));
		cut_cell_nodes_backup.at(i) = NULL;
	}
	cut_cell_nodes_backup.clear();
	
	
	
	//cout << "non_cut_cell_nodes_backup: " << non_cut_cell_nodes_backup.size() << endl;
	
	for(int i = 0; i < non_cut_cell_nodes_backup.size(); i++){	
		non_cut_cell_nodes.push_back(non_cut_cell_nodes_backup.at(i));
		non_cut_cell_nodes_backup.at(i) = NULL;
	}
	non_cut_cell_nodes_backup.clear();
	return;
}

//divides all cells passed without questioning for a set number of times
void divide_cut_cells(vector<ADTNode*> *to_be_processd, int num_times){
	if(num_times == 0){
		//replace cutcells (backup > cut, clear backup)
		reload();
		return;
	}
	
	vector<ADTNode*> *new_to_be_processed = new vector<ADTNode*>;
	
	while(!to_be_processd->empty()){
	
		ADTNode *node = to_be_processd->back();

		if(node == NULL){
			//cout << "base_case" << endl;
		}
		
		CartCell *cell = static_cast<CartCell*>(node->poly);

		vector<Polyhedron*> *temp = static_cast<CartCell*>(node->poly)->intersect_candidates;
		
		refine(node, false);
		new_to_be_processed->push_back(node->L_branch->L_branch);
		new_to_be_processed->push_back(node->L_branch->R_branch);
		new_to_be_processed->push_back(node->R_branch->L_branch);
		new_to_be_processed->push_back(node->R_branch->R_branch);
			

		CartCell *tempc1 = static_cast<CartCell*>(node->L_branch->L_branch->poly);
		CartCell *tempc2 = static_cast<CartCell*>(node->L_branch->R_branch->poly);
		CartCell *tempc3 = static_cast<CartCell*>(node->R_branch->L_branch->poly);
		CartCell *tempc4 = static_cast<CartCell*>(node->R_branch->R_branch->poly);
		
		vector<Polyhedron*> *t1 = new vector<Polyhedron*>;
		vector<Polyhedron*> *t2 = new vector<Polyhedron*>;
		vector<Polyhedron*> *t3 = new vector<Polyhedron*>;
		vector<Polyhedron*> *t4 = new vector<Polyhedron*>;
		for(int i = 0; i < temp->size(); i++){
			t1->push_back(temp->at(i));
			t2->push_back(temp->at(i));
			t3->push_back(temp->at(i));
			t4->push_back(temp->at(i));
		}
		
		tempc1->intersect_candidates = t1;
		tempc2->intersect_candidates = t2;
		tempc3->intersect_candidates = t3;
		tempc4->intersect_candidates = t4;


		/*
		//find_intersections(current_bound, surface, overlap_nonstrict);
		tempc1->intersect_candidates = find_intersections(tempc1->bounding_box, segments_ADT, overlap_nonstrict);
		tempc2->intersect_candidates = find_intersections(tempc2->bounding_box, segments_ADT, overlap_nonstrict);
		tempc3->intersect_candidates = find_intersections(tempc3->bounding_box, segments_ADT, overlap_nonstrict);
		tempc4->intersect_candidates = find_intersections(tempc4->bounding_box, segments_ADT, overlap_nonstrict);
		*/
		
		check_cutcell(node->L_branch->L_branch);
		check_cutcell(node->L_branch->R_branch);
		check_cutcell(node->R_branch->L_branch);
		check_cutcell(node->R_branch->R_branch);
		
		to_be_processd->back() = NULL;
		to_be_processd->pop_back();
	}
	
	divide_cut_cells(new_to_be_processed, num_times-1);
}

//check if the neighbour of a cell (node) is of higher than two refinement levels apart
bool check_neighbour(ADTNode *node, int reach = 1){
	CartCell *cell = static_cast<CartCell*>(node->poly);
	
	return ((edge_depth(cell->top) > reach) ||
			(edge_depth(cell->bot) > reach) ||
			(edge_depth(cell->left) > reach) ||
			(edge_depth(cell->right) > reach) );
}

//smoothen the current mesh by making sure all cells passed satisify smoothness criteria
void mesh_smoothing(vector<ADTNode*> *to_be_processd){	
	vector<ADTNode*> *new_to_be_processed = new vector<ADTNode*>;
	
	//cout << "smoothing called, size is: " << to_be_processd->size() << endl;
	
	while(!to_be_processd->empty()){
	
		ADTNode *node = to_be_processd->back();

		if(node == NULL){
			//cout << "base_case" << endl;
		}
		
		CartCell *cell = static_cast<CartCell*>(node->poly);
				
		vector<Polyhedron*> *temp = static_cast<CartCell*>(node->poly)->intersect_candidates;
		
		//if one of the neighbours is two levels higher
		if( check_neighbour(node) ){
			
			//cout << "refining" << endl;
			
			//refine cell (node) and queue new subcells
			refine(node, false);
			new_to_be_processed->push_back(node->L_branch->L_branch);
			new_to_be_processed->push_back(node->L_branch->R_branch);
			new_to_be_processed->push_back(node->R_branch->L_branch);
			new_to_be_processed->push_back(node->R_branch->R_branch);			
			
			if(temp != NULL){
			
				CartCell *tempc1 = static_cast<CartCell*>(node->L_branch->L_branch->poly);
				CartCell *tempc2 = static_cast<CartCell*>(node->L_branch->R_branch->poly);
				CartCell *tempc3 = static_cast<CartCell*>(node->R_branch->L_branch->poly);
				CartCell *tempc4 = static_cast<CartCell*>(node->R_branch->R_branch->poly);
			
				vector<Polyhedron*> *t1 = new vector<Polyhedron*>;
				vector<Polyhedron*> *t2 = new vector<Polyhedron*>;
				vector<Polyhedron*> *t3 = new vector<Polyhedron*>;
				vector<Polyhedron*> *t4 = new vector<Polyhedron*>;
				for(int i = 0; i < temp->size(); i++){
					t1->push_back(temp->at(i));
					t2->push_back(temp->at(i));
					t3->push_back(temp->at(i));
					t4->push_back(temp->at(i));
				}
				
				tempc1->intersect_candidates = t1;
				tempc2->intersect_candidates = t2;
				tempc3->intersect_candidates = t3;
				tempc4->intersect_candidates = t4;

				/*		
				//find_intersections(current_bound, surface, overlap_nonstrict);
				tempc1->intersect_candidates = find_intersections(tempc1->bounding_box, segments_ADT, overlap_nonstrict);
				tempc2->intersect_candidates = find_intersections(tempc2->bounding_box, segments_ADT, overlap_nonstrict);
				tempc3->intersect_candidates = find_intersections(tempc3->bounding_box, segments_ADT, overlap_nonstrict);
				tempc4->intersect_candidates = find_intersections(tempc4->bounding_box, segments_ADT, overlap_nonstrict);
				*/
			}
			
		}
		else{
			if(temp != NULL){
				cut_cell_nodes_backup.push_back(node);
			}
			else{
				non_cut_cell_nodes_backup.push_back(node);
			}
		}
		
		
		to_be_processd->back() = NULL;
		to_be_processd->pop_back();
	}
	
	if(!new_to_be_processed->empty()){
		removeDuplicates(new_to_be_processed);
		mesh_smoothing(new_to_be_processed);
	}
	return;
}

//buffering algorithm to create layers of padding around, refines all cells of a certain depth
bool buffer_helper(int i, int lvl_to_be_refined){
	
	return (i == lvl_to_be_refined); 
	
}
void buffer(int lvl_to_be_refined){
	reload();
	while(!cut_cell_nodes.empty()){
	
		ADTNode *node = cut_cell_nodes.back();

		if(node == NULL){
			//cout << "base_case" << endl;
		}
		
		CartCell *cell = static_cast<CartCell*>(node->poly);
				
		vector<Polyhedron*> *temp = static_cast<CartCell*>(node->poly)->intersect_candidates;
		
		//if cell is at the depth we want refined
		if( buffer_helper(node->depth, lvl_to_be_refined) ){
			
			//refine cell (node) and queue new subcells
			refine(node, false);						
			CartCell *tempc1 = static_cast<CartCell*>(node->L_branch->L_branch->poly);
			CartCell *tempc2 = static_cast<CartCell*>(node->L_branch->R_branch->poly);
			CartCell *tempc3 = static_cast<CartCell*>(node->R_branch->L_branch->poly);
			CartCell *tempc4 = static_cast<CartCell*>(node->R_branch->R_branch->poly);
										
			//find_intersections(current_bound, surface, overlap_nonstrict);
			tempc1->intersect_candidates = find_intersections(tempc1->bounding_box, segments_ADT, overlap_nonstrict);
			tempc2->intersect_candidates = find_intersections(tempc2->bounding_box, segments_ADT, overlap_nonstrict);
			tempc3->intersect_candidates = find_intersections(tempc3->bounding_box, segments_ADT, overlap_nonstrict);
			tempc4->intersect_candidates = find_intersections(tempc4->bounding_box, segments_ADT, overlap_nonstrict);
			
			check_cutcell(node->L_branch->L_branch);
			check_cutcell(node->L_branch->R_branch);
			check_cutcell(node->R_branch->L_branch);
			check_cutcell(node->R_branch->R_branch);
		
		}
		else{
			cut_cell_nodes_backup.push_back(node);
		}
		
		cut_cell_nodes.back() = NULL;
		cut_cell_nodes.pop_back();
	}
	while(!non_cut_cell_nodes.empty()){
		
		ADTNode *node = non_cut_cell_nodes.back();

		if(node == NULL){
			//cout << "base_case" << endl;
		}
		
		CartCell *cell = static_cast<CartCell*>(node->poly);
				
		vector<Polyhedron*> *temp = static_cast<CartCell*>(node->poly)->intersect_candidates;
		
		if( buffer_helper(node->depth, lvl_to_be_refined) ){
			
			//cout << "refining" << endl;
			
			refine(node, false);
			non_cut_cell_nodes_backup.push_back(node->L_branch->L_branch);
			non_cut_cell_nodes_backup.push_back(node->L_branch->R_branch);
			non_cut_cell_nodes_backup.push_back(node->R_branch->L_branch);
			non_cut_cell_nodes_backup.push_back(node->R_branch->R_branch);					
		}
		else{
			non_cut_cell_nodes_backup.push_back(node);
		}
		
		non_cut_cell_nodes.back() = NULL;
		non_cut_cell_nodes.pop_back();
	}
}

//painting algorithm to quickly mark cells as in or out of the geomrtry
void paint_helper(ADTNode *node, bool brush){
	if(node == NULL){
		return;
	}
	
	CartCell *cell = static_cast<CartCell*>(node->poly);
	
	//if cutcell or painted, do nothing
	if( (cell->intersect_candidates != NULL) || (cell->painted) ){
		return;
	}
	
	//cout << "painted: in? "  << brush << " " << *cell->bounding_box;
	
	cell->painted = true;
	cell->in_geo = brush;
	
	ADTNode *top = find_other(cell->top, node);
	ADTNode *bot = find_other(cell->bot, node);
	ADTNode *left = find_other(cell->left, node);
	ADTNode *right = find_other(cell->right, node);
	
	//recursively paint neighbours that are not cutcells
	paint_helper(top, brush);
	paint_helper(bot, brush);
	paint_helper(left, brush);
	paint_helper(right, brush);

	return;
}
void paint(vector<ADTNode*> *to_be_painted){
	
	for(int i = 0; i < to_be_painted->size(); i++){
		ADTNode *node = to_be_painted->at(i);
		CartCell *cell = static_cast<CartCell*>(node->poly);
		if(cell->painted){
			continue;	//skip painted cells
		}
		//get the brush type and paint the neighbours
		paint_helper(node, inside_geometry(cell, segments_ADT));
	}
}

//prints a cell (node)'s info
void print_cell(ADTNode *node){
	if(node == NULL){
		//cout << "NULL" << endl;
		return;
	}
	//cout << *(static_cast<CartCell*>(node->poly)->bounding_box);	
}

//prints all edge structural info of all cells in the mesh
void check_structure(){
	for(int i = 0; i < non_cut_cell_nodes.size(); i++){
		ADTNode *node = non_cut_cell_nodes.at(i);
		CartCell *cell = static_cast<CartCell*>(node->poly);
		
		//cout << "Cell is: ";
		print_cell(node);
		//cout << "Cell's top is: ";
		print_cell(find_other(cell->top, cell));
		//cout << "Cell's bot is: ";
		print_cell(find_other(cell->bot, cell));
		//cout << "Cell's left is: ";
		print_cell(find_other(cell->left, cell));
		//cout << "Cell's right is: ";
		print_cell(find_other(cell->right, cell));
		//cout << endl;
	}
	return;
}

//generates the mesh based on the input parameters
ADTNode* generate_mesh(ADTNode *surface, Refinement *ref, ADTNode *init_mesh, bool (*ref_cond)(ADTNode*, void*), void *ref_rule){
	//generates mesh with manual parameters
	
	ADTNode* result;
	
	if(*PROGRAM_GENE_TYPE == MANUAL){
		result = mesh_refine(surface, ref, init_mesh, 0);
		paint(&non_cut_cell_nodes_backup);
		paint(&non_cut_cell_nodes);
	}
	else{
		result = mesh_refine(surface, ref_cond, init_mesh, 0, ref_rule);
		
		if(*PROGRAM_GENE_TYPE == ADAPTIVE){
			//check_structure();
			
			//cout << "calling divide_cut_cells" << endl;
			Adaptive_variable *adaptive_var = static_cast<Adaptive_variable*>(ref_rule);

			
			//divide_cut_cells(&cut_cell_nodes, adaptive_var->max_refs - 4);
			
			
			//cout << "calling mesh_smoothing" << endl;
			mesh_smoothing(&cut_cell_nodes);
			
			for(int i = 0; i < cut_cell_nodes_backup.size(); i++){	
				cut_cell_nodes.push_back(cut_cell_nodes_backup.at(i));
				cut_cell_nodes_backup.at(i) = NULL;
			}
			cut_cell_nodes_backup.clear();
			
			mesh_smoothing(&non_cut_cell_nodes);
						
			//create paddings, by building layers of certain depth cells from deepest up
			for(int j = adaptive_var->max_refs; j > 0; j--){
				for(int i = 0; i < 2; i++){
					//buffering
					buffer(j);
					reload();		
					mesh_smoothing(&cut_cell_nodes);
					mesh_smoothing(&non_cut_cell_nodes);
				}
			}
			
			
			
			reload();
			
			//cout << "calling paint" << endl;
			paint(&non_cut_cell_nodes_backup);
			paint(&non_cut_cell_nodes);
		}
	}

	//deal with stack of cutcells
	//cout << "cut_cell process starts" << endl;
	
	int count = 0;
	//int cut_celL_limit = 14;
	
	
	while(!cut_cell_nodes.empty()){
		CartCell *cell = static_cast<CartCell*>(cut_cell_nodes.back()->poly);
		
		
		/*
		vector<Polyhedron*> *test = find_intersections_show(cell->bounding_box, surface, overlap_nonstrict);
		//cout << test->size() << endl;
		
		
		//THIS SHOULD BE NOT NULL
		//cout << "%%%%%%%%%%%%%%%%  CUT CELL's intersect_candidates list  %%%%%%%%%%%%%%%%%" << endl;
		for(int j = 0; j < cell->intersect_candidates->size(); j++){
			edge_stat(static_cast<Segment*>(cell->intersect_candidates->at(j)));
		}
		//cout << "%%%%%%%%%%%%%%%%  END OF CUT CELL's intersect_candidates list  %%%%%%%%%%%%%%%%%" << endl;
		*/
		
		generate_cell_lists(cell);
		
		//cell preparation, swapping vertices and edges with vert_intersection and body_edge classes
		vector<vert_intersection*> *full_list_intersect_pts = new vector<vert_intersection*>;
		
		//cout << "before preparing cell" << endl;
		
		for(int j = 0; j < CARTCELL_NUM_INDICES; j++){
			vector<vert_intersection*> *list_intersect_pts = new vector<vert_intersection*>;
			
			//cout << "before processing edge " << j << endl;
			
			for(int i = 0; i < cell->intersect_candidates->size(); i++){
				
				if(cell->intersect_candidates->at(i) == NULL){
					continue;
				}
				
				//cout << "before processing segment " << j << endl;
				
				Segment *seg_temp = static_cast<Segment*>(cell->intersect_candidates->at(i));
				
				pt4D *cell_bound = new pt4D(precise_verts.at(cell->pos)->x,
											precise_verts.at(cell->pos)->x + cell->length,
											precise_verts.at(cell->pos)->y,
											precise_verts.at(cell->pos)->y + cell->length);
				
				//cout << "before checking if not contained nor intersect" << endl;
				
				/*
				edge_stat(seg_temp);
				edge_stat(precise_edges.at(cell->list_of_edges[j]));
				//cout << do_lines_intersect(seg_temp, precise_edges.at(cell->list_of_edges[j])) << endl;
				//cout << *(cell->bounding_box) << endl;
				//cout << *(seg_temp->bounding_box) << endl;
				//cout << contained(cell->bounding_box, seg_temp->bounding_box) << endl;
				//cout << "wtf" << endl;
				*/
				
				if(!contained(cell->bounding_box, seg_temp->bounding_box)
				 && do_lines_intersect(seg_temp, precise_edges.at(cell->list_of_edges[j]))){
				 
					//the lines intersect, create vert_intersect
					//cout << "the lines intersect, create vert_intersect" << endl;
					
					edge *edge_temp = precise_edges.at(cell->list_of_edges[j]);
					double x1 = precise_verts.at(seg_temp->list_of_verts[0])->x;
					double y1 = precise_verts.at(seg_temp->list_of_verts[0])->y;
					double x2 = precise_verts.at(seg_temp->list_of_verts[1])->x;
					double y2 = precise_verts.at(seg_temp->list_of_verts[1])->y;
					vert_intersection *inter_vert;
					double x, y;
					
					if(precise_verts.at(edge_temp->pt1)->x == precise_verts.at(edge_temp->pt2)->x){
						//x is same, bound is vert line
						//cout << "x is same, bound is vert line" << endl;
						
						x = precise_verts.at(edge_temp->pt1)->x;
						if( fabs(x2 - x1) > *EPSILON ){
							//segment not vert line
							const double M = (y2 - y1) / (x2 - x1);	//slope
							const double B = -(M * x1) + y1;	//y-intercept
							y = (M * x + B);
							inter_vert = new vert_intersection(x, y, edge_temp, ON_BOUND, false, precise_verts.size());
							precise_verts.push_back(inter_vert);
							
							//cout << endl;
							//cout << "Segment: (" << x1 << ", " << y1 << ") -- (" << x2 << ", " << y2 << ")" << endl;
							//cout << *inter_vert << endl;
							
							list_intersect_pts->push_back(inter_vert);
							full_list_intersect_pts->push_back(inter_vert);
						}
					}
					else{
						//y is same, bound is horz line
						//cout << "y is same, bound is horz line" << endl;
						
						y = precise_verts.at(edge_temp->pt1)->y;
						if( fabs(y2 - y1) > *EPSILON ){
							//segment not horz line			
							if( fabs(x2 - x1) <= *EPSILON ){
								//vert line, x from here, y from boundary
								x = x1;
							}
							else{
								const double M = (y2 - y1) / (x2 - x1);	//slope
								const double B = -(M * x1) + y1;	//y-intercept
								x = (y - B) / M;
							}
							inter_vert = new vert_intersection(x, y, edge_temp, ON_BOUND, false, precise_verts.size());
							precise_verts.push_back(inter_vert);
							
							//cout << endl;
							//cout << "Segment: (" << x1 << ", " << y1 << ") -- (" << x2 << ", " << y2 << ")" << endl;
							//cout << *inter_vert << endl;
							//cout << *(precise_verts.at(inter_vert->ID)) << endl;
							
							list_intersect_pts->push_back(inter_vert);
							full_list_intersect_pts->push_back(inter_vert);
						}						
					}
				}
				//replace all the verts in segments with the subclass version (check in/out)
				//cout << "replacing s all the verts in segments with the subclass version" << endl;
				vert_cond cond1 = is_inside(cell, precise_verts.at(seg_temp->list_of_verts[0]));
				vert_cond cond2 = is_inside(cell, precise_verts.at(seg_temp->list_of_verts[1]));
				
				
				//cout << ">>>>>>>>>>the verts are<<<<<<<<<" << endl;
				//cout << *(precise_verts.at(seg_temp->list_of_verts[0])) << endl;
				//cout << *(precise_verts.at(seg_temp->list_of_verts[1])) << endl;
				//cout << cond1 << endl;
				//cout << cond2 << endl;
				
				bool b1 = (cond1 == ON_BOUND);
				bool b2 = (cond2 == ON_BOUND);		
				body_edge *eg = new body_edge(precise_edges.at(seg_temp->list_of_edges[0])->pt1, 
											  precise_edges.at(seg_temp->list_of_edges[0])->pt2, 
											  NULL, false);	
				vert_intersection *v1 = new vert_intersection(
												precise_verts.at(seg_temp->list_of_verts[0])->x,
												precise_verts.at(seg_temp->list_of_verts[0])->y,
												eg, cond1, b1, seg_temp->list_of_verts[0]);
				vert_intersection *v2 = new vert_intersection(
												precise_verts.at(seg_temp->list_of_verts[1])->x,
												precise_verts.at(seg_temp->list_of_verts[1])->y,
												eg, cond2, b2, seg_temp->list_of_verts[1]);	
												
				//cout << ">>>>>>>>>>the verts after<<<<<<<<<" << endl;
				//cout << *v1 << endl;
				//cout << *v2 << endl;
									
				delete precise_verts.at(seg_temp->list_of_verts[0]);
				delete precise_verts.at(seg_temp->list_of_verts[1]);
				precise_verts.at(seg_temp->list_of_verts[0]) = v1;
				precise_verts.at(seg_temp->list_of_verts[1]) = v2;
				
				//cout << "done processing segment " << j << endl;
			}
			
			//cout << "done processing edge " << j << endl;
			
			//replace all the edges and verts in cartcell to subclass version
			vert_intersection *v_temp = new vert_intersection(
											precise_verts.at(cell->list_of_verts[j])->x,
											precise_verts.at(cell->list_of_verts[j])->y,
											NULL, ON_BOUND, true, 
											precise_verts.at(cell->list_of_verts[j])->ID);
					
			//cout << "orginal edge" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt1)->x << ", "
			//	 << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt1)->y << ")" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt2)->x << ", "
			//	 << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt2)->y << ")" << endl;
			
											
			body_edge *e_temp = new body_edge(precise_edges.at(cell->list_of_edges[j])->pt1, 
											  precise_edges.at(cell->list_of_edges[j])->pt2, 
											  list_intersect_pts, true);
														  
			precise_verts.at(cell->list_of_verts[j]) = v_temp;
			precise_edges.at(cell->list_of_edges[j]) = e_temp;
			
			//cout << "new edge" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt1)->x << ", "
			//	 << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt1)->y << ")" << endl;
			//cout << "(" << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt2)->x << ", "
			//	 << precise_verts.at(precise_edges.at(cell->list_of_edges[j])->pt2)->y << ")" << endl;
			
			removeDuplicates(list_intersect_pts);
			//cout << "@@@--size of list_intersections: " << list_intersect_pts->size() << endl;
		}
		
		removeDuplicates(full_list_intersect_pts);
		//cout << "@@@--size of full_list_intersections: " << full_list_intersect_pts->size() << endl;
		
		//done preparing the cell
		
		//cout << "cut_cell prepared" << endl;
		
		cell_piece *piece = new cell_piece(cell->list_of_verts, cell->list_of_edges, 
										   cell->num_indices, full_list_intersect_pts, 
										   cell->intersect_candidates);
										   
		remove_trivial_seg(piece);
		
		/*
		//cout << "num_indices: " << **(cell->num_indices) << endl;
		for(int i = 0; i < **(cell->num_indices); i++){
			//cout << *(precise_verts.at(cell->list_of_verts[i])) << endl;
			//cout << *(precise_edges.at(cell->list_of_edges[i])) << endl;
		}
		
		//cout << "full_list_intersect_pts:" << endl;
		for(int j = 0; j < full_list_intersect_pts->size(); j++){
			//cout << *(full_list_intersect_pts->at(j)) << endl;
		}
		
		for(int j = 0; j < piece->intersect_segments->size(); j++){
			edge_stat(static_cast<Segment*>(piece->intersect_segments->at(j)));
		}
		////cout << *(piece) << endl;
		*/
		
		//make the CutCell in place to replace the orginal cell		
		CutCell *new_cutcell = make_cut_cell(piece, surface);		
		CutCell *temp_cell = static_cast<CutCell*>(cut_cell_nodes.back()->poly);
		
		*temp_cell = *new_cutcell;	
		
		cut_cell_nodes.back() = NULL;
		cut_cell_nodes.pop_back();
		
	}
	
	return result;	
}

//reads the input msh file
void read_geo_file(const char* shape_txt){
	ifstream infile;
	infile.open(shape_txt);
	
	string junk;
	int int_junk;
	
	//removing version number
	getline(infile, junk);
	getline(infile, junk);
	getline(infile, junk);
	infile >> junk; //Nodes
	
	int num_nodes;
	infile >> num_nodes;
	
	for(int i = 0; i < num_nodes; i++){
		int node_ID;
		double x;
		double y;
		double z;
	 	infile >> node_ID;
	 	infile >> x;
	 	infile >> y;
	 	infile >> z;
		vert_precise *vert_temp = new vert_precise(x, y, precise_verts.size());
	 	//add the new vertex
	 	precise_verts.push_back(vert_temp);
	}
 	
 	infile >> junk; //EndNodes
 	infile >> junk; //Edges

	int num_edges;
	infile >> num_edges;
	geo_num = num_edges;
	for(int i = 0; i < num_edges; i++){
		int edge_ID;
		int n1;
		int n2;
	 	infile >> edge_ID;
	 	
	 	//edge attributes
	 	infile >> int_junk;
	 	infile >> int_junk;
	 	//useless ^
	 	
	 	infile >> n1;
	 	infile >> n2;	 	
	 	edge *edge_temp = new edge(n1-1, n2-1);
	 	//indexing
	 	int pos = precise_edges.size();
	 	//add the new edges
	 	precise_edges.push_back(edge_temp);
	 	
	 	//add the new pt4D
	 	double x1 = min(precise_verts.at(n1-1)->x, precise_verts.at(n2-1)->x);
	 	double x2 = max(precise_verts.at(n1-1)->x, precise_verts.at(n2-1)->x);
	 	double y1 = min(precise_verts.at(n1-1)->y, precise_verts.at(n2-1)->y);
	 	double y2 = max(precise_verts.at(n1-1)->y, precise_verts.at(n2-1)->y);
	 					
	 	pt4D *box_temp = new pt4D(x1, x2, y1, y2);
	 
	 	//add the new segment	
	 	int *seg_verts = new int[2];
	 	seg_verts[0] = n1-1;
	 	seg_verts[1] = n2-1;
	 	
	 	int *seg_edges = new int[1];
	 	seg_edges[0] = pos;	 	
	 	Segment *seg_temp = new Segment(box_temp, seg_verts, seg_edges);
	 	
	 	segments.push_back(seg_temp);	 	
	}
	 	
	infile.close();
	return;
}

//returns the refinement table in ref_txt for Manual mode
Refinement *read_ref_file(const char* ref_txt){
	vector<vector<pt4D*>*> *ref_table = new vector<vector<pt4D*>*>;
	ifstream infile;
	infile.open(ref_txt);
	int ref_lvls;	//number of ref levels
	int num_refs;	//number of refs in each level
	infile >> ref_lvls;
	////cout << "ref_lvls: " << ref_lvls << endl;
	for(int i = 0; i < ref_lvls; i++){
		vector<pt4D*> *temp_ref_vec = new vector<pt4D*>;
		infile >> num_refs;
		////cout << "num_refs: " << num_refs << endl;
		for(int j = 0; j < num_refs; j++){
			double x1, x2, y1, y2;
			infile >> x1 >> x2 >> y1 >> y2;
			////cout << x1 << " " << x2 << " " << y1 << " " << y2 << endl;
			pt4D *temp_ref = new pt4D(x1, x2, y1, y2);
			temp_ref_vec->push_back(temp_ref);
		}
		ref_table->push_back(temp_ref_vec);
	}
	infile.close();
	Refinement *ref = new Refinement(ref_table);
	return ref;
}

//imports the refinement rules from ref_txt for Min_size mode
void min_size_ref_read(vector<Minsize_Refinement*> *ref_minsize, const char* ref_txt){
	ifstream infile;
	infile.open(ref_txt);
	
	int list_size;
	infile >> list_size;
	
	for(int i = 0; i < list_size; i++){
		double x1, x2, y1, y2, min_size;
		infile >> x1 >> x2 >> y1 >> y2 >> min_size;
		
		pt4D *temp_box = new pt4D(x1, x2, y1, y2);
		Minsize_Refinement *temp = new Minsize_Refinement(temp_box, min_size);
		
		ref_minsize->push_back(temp);
	}
}

//returns the variables used for Adaptive mode
Adaptive_variable *adapt_ref_read(const char* ref_txt){
	ifstream infile;
	infile.open(ref_txt);
	
	string junk;
	int max_refs;
	double alpha;
	
	infile >> junk;
	infile >> max_refs;
	infile >> junk;
	infile >> alpha;
	
	Adaptive_variable *ref = new Adaptive_variable(max_refs, alpha);
	
	return ref;	
}

//prints the mesh's ADT structure to file
void print_grid_helper(queue<ADTNode*> *print_queue, ofstream &outfile){
	//NULL branches are not printed
	int depth = 0;
	ADTNode *ptr_null = NULL;	
	print_queue->push(ptr_null);
	
	////cout << "before entering GRIFD_PRINT" << endl;
	GRID_PRINT:
		if(print_queue->empty()){
			return;
		}
		if(print_queue->front() == NULL){
			print_queue->pop();
			goto GRID_PRINT;
		}
		outfile << "####################################################################" << endl;
		outfile << "Level: " << depth << endl;
		print_queue->push(ptr_null);	//end of level
		LEVEL_PRINT:
			ADTNode *temp = print_queue->front();
			CartCell *temp_cell = new CartCell(static_cast<CartCell*>(temp->poly)->bounding_box,
											   static_cast<CartCell*>(temp->poly)->length,
											   static_cast<CartCell*>(temp->poly)->pos, *INTERNAL);
			outfile << *temp_cell << endl;
			if(temp->L_branch != NULL){
				print_queue->push(temp->L_branch);
			}
			if(temp->R_branch != NULL){
				print_queue->push(temp->R_branch);
			}
			print_queue->pop();
			if(print_queue->front() == NULL){
				//this level is completed
				depth++;
				print_queue->pop();
				goto GRID_PRINT;
			}
		goto LEVEL_PRINT;
}
void print_grid_to_file(ADTNode *mesh, const char* filename){
	ofstream outfile;
	outfile.open(filename);
	outfile << "Alternating Digital Tree" << endl;
	outfile << "Even levels are CartCells" << endl;
	queue<ADTNode*> *print_queue = new queue<ADTNode*>;
	print_queue->push(mesh);
	//NULL branches are not printed
	print_grid_helper(print_queue, outfile);
	outfile.close();
	return;
}

//algorithm for finding a vertice using binary search
int bin_search(vector<vert_precise*> vec, vert_precise *value, int low, int hi){
	if(hi < low){
		return -1;
	}
	int mid = (low + hi) / 2;
	if( more_than(vec.at(mid), value) ){
		return bin_search(vec, value, low, mid-1);
	}
	else if( less_than(vec.at(mid), value) ){
		return bin_search(vec, value, mid+1, hi);
	}
	else{
		return mid;
	}
}
int find(int ID){
	return bin_search(precise_verts, precise_verts_backup.at(ID), 0, precise_verts.size()-1);
}
int find(vert_precise *vert){
	return bin_search(precise_verts, vert, 0, precise_verts.size()-1);
}

//algorithm for finding new boundary elements (20000 elements)
bool is_present(vector<edge*> *edges, edge *eg){
	for(int i = 0; i < edges->size(); i++){
		if(eg->pt1 == edges->at(i)->pt1 && eg->pt2 == edges->at(i)->pt2){
			return true;
		}
	}
	return false;
}
void keep_dups(){
	//print_boundary
	vector<edge*> removed;
	
	for(vector<edge*>::iterator itr = print_boundary.begin(); itr != print_boundary.end(); ++itr){
		if(!is_present(&removed, *itr)){
			removed.push_back(*itr);
			print_boundary.erase(itr);
			itr--;
		}
	}

	return;
}

//writes the n-sided polygon into the mesh, helper function for print_msh
void print_sol_ngon(n_gon *ngon, ofstream &outfile){
	
	/*
	for(int i = 0; i < ngon->n; i++){
		outfile << ID << " 1 0 " 
				<< find(ngon->edges->at(i)->pt1) + 1 << " "
				<< find(ngon->edges->at(i)->pt2) + 1 << endl;
	}*/
	
	int *verts;
	verts = new int[ngon->n];
	
	int next = ngon->edges->at(0)->pt1;
	
	if(next == ngon->edges->at(1)->pt1 || next == ngon->edges->at(1)->pt2){
		next = ngon->edges->at(0)->pt2;
	}
	
	for(int i = 0; i < ngon->n; i++){
		verts[i] = next;
		
		if(next == ngon->edges->at(i)->pt1){
			next = ngon->edges->at(i)->pt2;
		}
		else{
			next = ngon->edges->at(i)->pt1;
		}
	}
	
	
	/*
	outfile << ID << " " << (ngon->n + 1) << " 601 600 " << ngon->n;
	
	for(int i = 0; i < ngon->n; i++){
		
		outfile << " " << find(verts[i]) + 1;

	}
	outfile << endl;
	*/
	
	int i,j,k;
	int count = 0;
   	double z;
   	
   	for (i = 0 ; i < ngon->n ; i++){
      	j = (i + 1) % ngon->n;
      	k = (i + 2) % ngon->n;
      	z  = (precise_verts_backup.at(verts[j])->x - precise_verts_backup.at(verts[i])->x) * 
      		 (precise_verts_backup.at(verts[k])->y - precise_verts_backup.at(verts[j])->y);
      	z -= (precise_verts_backup.at(verts[j])->y - precise_verts_backup.at(verts[i])->y) * 
      		 (precise_verts_backup.at(verts[k])->x - precise_verts_backup.at(verts[j])->x);
      	if (z < 0)
         	count--;
      	else if (z > 0)
         	count++;
   	}
   	
   	
   	int code = (ngon->n + 1);
	
	if(ngon->n == 4){
		code = 3;
	}
	
	if(ngon->n == 3){
		code = 2;
	}
	
   	if (count > 0){
     	outfile << code << " 601 600 " << ngon->n;
		for(int i = 0; i < ngon->n; i++){
			outfile << " " << find(verts[i]) + 1;
		}
		outfile << endl;
	}
   	else if (count < 0){
    	outfile << code << " 601 600 " << ngon->n;
		for(int i = ngon->n - 1; i >= 0; i--){
			outfile << " " << find(verts[i]) + 1;
		}
		outfile << endl;
    }
    
    delete [] verts;
    
    /*
	int min_x = INT_MAX;
	int ID_pt;
	int ID_prev;
	int ID_next;
	
	for(int i = 0; i < ngon->n; i++){
		if(precise_verts_backup.at(verts[i])->x < min_x){
			min_x = precise_verts_backup.at(verts[i])->x;
			ID_pt = i;
			
			if(i == 0){
				ID_prev = ngon->n - 1;
				ID_next = i + 1;
			}
			else if(i == (ngon->n - 1)){
				ID_prev = i - 1;
				ID_next = 0;				
			}
			else{
				ID_prev = i - 1;
				ID_next = i + 1;
			}
			
		}
	}
	
	double v1x = precise_verts_backup.at(ID_pt)->x - precise_verts_backup.at(ID_prev)->x;
	double v1y = precise_verts_backup.at(ID_pt)->y - precise_verts_backup.at(ID_prev)->y;
	double v2x = precise_verts_backup.at(ID_next)->x - precise_verts_backup.at(ID_pt)->x;
	double v2y = precise_verts_backup.at(ID_next)->y - precise_verts_backup.at(ID_pt)->y;
	
	double z = (v1x * v2y) - (v1y * v2x);
	*/	
}

//prints the mesh to a msh file, queue used to avoid recursion and stack overflow
void print_msh_helper(queue<ADTNode*> *print_queue, ofstream &outfile){
	ADTNode *ptr_null = NULL;	
	print_queue->push(ptr_null);
	MESH_PRINT:
		//cout << "MESH_PRINT" << endl;
		if(print_queue->empty()){
			return;
		}
		//NULL is a sentinel value
		if(print_queue->front() == NULL){
			print_queue->pop();
			goto MESH_PRINT;
		}
		print_queue->push(ptr_null);	//end of level
		MESH_LEVEL_PRINT:	
			//cout << "MESH_LEVEL_PRINT" << endl;
			if(print_queue->front() == NULL){
				print_queue->pop();
				goto MESH_PRINT;
			}
			ADTNode *temp = print_queue->front();
			Graph* temp_poly = static_cast<Graph*>(temp->poly);
			
			if((temp->L_branch == NULL) && (temp->R_branch == NULL)){
				if(temp->poly->generic){
					delete print_queue->front();
					print_queue->pop();
					goto MESH_PRINT;
				}
				
				CartCell *temp_cell = static_cast<CartCell*>(temp->poly);
											   
				//cout << "arrived at a leaf" << endl;
				if(temp_poly->num_entries == NULL){
					//if not specified then its not a bodycell, cartcell used
					
					if(temp_poly->is_cut){
						//cout << "err: the cutcell has NULL for num entries" << endl;
					}
					
					temp_poly->num_entries = &CARTCELL_NUM_ENTRIES;
					
					int **temp_array = new int*[1];
					temp_array[0] = &CARTCELL_NUM_INDICES;
					temp_poly->num_indices = temp_array;
					generate_cell_lists(temp_cell);
					
					if(temp_cell->in_geo){
						delete print_queue->front();
						print_queue->pop();
						goto MESH_LEVEL_PRINT;
					}
				}
				
				int part_sum = 0;
				
				for(int i = 0; i < *(temp_poly->num_entries); i++){
				
					n_gon *temp_piece = new n_gon(*(temp_poly->num_indices[i]));
					temp_piece->edges = new vector<edge*>;
					
					for(int j = 0; j < *(temp_poly->num_indices[i]); j++){
						//add to print verts and edges
						//int id = print_edges.size();
						
						if(temp->poly->is_cut){
							
							//cout << "in print_msh_helper" << endl;
							CutCell *temp_cut = static_cast<CutCell*>(temp->poly);
							cell_piece_info(temp_cut);
							
							temp_piece->edges->push_back(precise_edges.at(temp_cut->list_of_edges[part_sum+j]));
							
							print_boundary.push_back(precise_edges.at(temp_cut->list_of_edges[part_sum+j]));
						}
						else{
							print_vertIDs.push_back(temp_cell->list_of_verts[j]);
						}
					}
					
					part_sum += *(temp_poly->num_indices[i]);
					
					if(temp->poly->is_cut){
						print_cutcells.push_back(temp_piece);
					}
					else{
						delete temp_piece;
					}
				}
			}
			else{
				//cout << "adding to queue" << endl;
				if(temp->L_branch != NULL){
					print_queue->push(temp->L_branch);
				}
				if(temp->R_branch != NULL){
					print_queue->push(temp->R_branch);
				}
			}
			delete print_queue->front();
			print_queue->pop();
			if(print_queue->front() == NULL){
				print_queue->pop();
				goto MESH_PRINT;
			}
		goto MESH_LEVEL_PRINT;
}
void print_msh(ADTNode *mesh, const char* filename){
	ofstream outfile;
	outfile.open(filename);

	queue<ADTNode*> *print_queue = new queue<ADTNode*>;
	print_queue->push(mesh);
	
	//cout << "before print_ms_helper" << endl;
	print_msh_helper(print_queue, outfile);
	
	//printing starts
	//copy, removeDup
	for(int i = 0; i < precise_verts.size(); i++){
		precise_verts_backup.push_back(precise_verts.at(i));
	}
	removeDuplicates(&precise_verts);


	int box_num_x = ((*GBOX_X_MAX - *GBOX_X_MIN) / *INIT_MIN);
	int box_num_y = ((*GBOX_Y_MAX - *GBOX_Y_MIN) / *INIT_MIN);
	int box_ele_num = 2 * box_num_x + 2 * box_num_y;	
	//always print precise_verts.at(bin_search(precise_verts, precise_verts_backup.at(val), 0, precise_verts.size()-1)
	//or find(val)
	
	/*
	outfile << "$MeshFormat" << endl;
	outfile << "2.2 0 8" << endl;
	outfile << "$EndMeshFormat" << endl;
	*/
	
	cout << "[Card_2d] printing nodes..." << endl;
	
	outfile << "$NOD" << endl;
	outfile << precise_verts.size() << endl;
	for(int i = 0; i < precise_verts.size(); i++){
		outfile << i+1 << " " << precise_verts.at(i)->x << " "
							  << precise_verts.at(i)->y << " 0" << endl;
	}
	outfile << "$ENDNOD" << endl;
	
	cout << "[Card_2d] printing elements..." << endl;
	
	outfile << "$ELM" << endl;
	
	//only keep dups for print_boundary
	keep_dups();
	
	outfile << print_cutcells.size() + (print_vertIDs.size() / 4) + box_ele_num + print_boundary.size() << endl;
	
	int current_ID = 1;
	
	for(int i = 0; i < print_boundary.size(); i++){
		//remember to change the element count
		outfile << current_ID + i << " 1 20000 10 2 " 
								  << find(print_boundary.at(i)->pt1) + 1 << " "
								  << find(print_boundary.at(i)->pt2) + 1 << endl;
	}
	current_ID += print_boundary.size();
	
	cout << "[Card_2d] print_boundary done" << endl;
	
	int CartID = 0;
	for(int i = 0; i < print_vertIDs.size(); i+=4){
		outfile << current_ID + CartID << " 3 601 600 4 " 
									   << find(print_vertIDs.at(i+3)) + 1 << " "
									   << find(print_vertIDs.at(i+2)) + 1 << " "
									   << find(print_vertIDs.at(i+1)) + 1 << " "
									   << find(print_vertIDs.at(i)) + 1 << endl;
		CartID++;
	}
	current_ID += print_vertIDs.size() / 4;
	
	cout << "[Card_2d] print_cartcells done" << endl;
	
	//CutCells
	for(int i = 0; i < print_cutcells.size(); i++){
		outfile << current_ID + i << " ";
		print_sol_ngon(print_cutcells.at(i), outfile);
	}
	current_ID += print_cutcells.size();
	
	cout << "[Card_2d] print_cutcells done" << endl;
	
	/*
	for(int i = 0; i < print_edges.size(); i++){
		outfile << i + current_ID << " 1 0 " 
											<< find(print_edges.at(i)->pt1) + 1 << " "
								      		<< find(print_edges.at(i)->pt2) + 1 << endl;
	}
	current_ID += print_edges.size();*/
	
	//Box	
	int j = 0;
	
	vert_precise top1 = vert_precise(*GBOX_X_MIN, *GBOX_Y_MIN, -1);
	vert_precise bot1 = vert_precise(*GBOX_X_MIN, *GBOX_Y_MAX, -1);
	int top_start = find(&top1);
	int bot_start = find(&bot1);
	for(int i = 1; i <= box_num_x; i++){
		vert_precise top2 = vert_precise(*GBOX_X_MIN + (i * *INIT_MIN), *GBOX_Y_MIN, -1);
		vert_precise bot2 = vert_precise(*GBOX_X_MIN + (i * *INIT_MIN), *GBOX_Y_MAX, -1);
		int top_end = find(&top2);
		int bot_end = find(&bot2);
		
		outfile << current_ID + j << " 1 30000 1 2 " 
								  << top_start + 1 << " "
								  << top_end + 1 << endl;
		j++;
		outfile << current_ID + j << " 1 30000 1 2 " 
								  << bot_start + 1 << " "
								  << bot_end + 1 << endl;
		j++;
		
		top_start = top_end;
		bot_start = bot_end;
	}
	
	vert_precise left1  = vert_precise(*GBOX_X_MIN, *GBOX_Y_MIN, -1);
	vert_precise right1 = vert_precise(*GBOX_X_MAX, *GBOX_Y_MIN, -1);
	int left_start = find(&left1);
	int right_start = find(&right1);
	for(int i = 1; i <= box_num_y; i++){
		vert_precise left2  = vert_precise(*GBOX_X_MIN, *GBOX_Y_MIN + (i * *INIT_MIN), -1);
		vert_precise right2 = vert_precise(*GBOX_X_MAX, *GBOX_Y_MIN + (i * *INIT_MIN), -1);
		int left_end = find(&left2);
		int right_end = find(&right2);
		
		outfile << current_ID + j << " 1 30000 1 2 " 
								  << left_start + 1 << " "
								  << left_end + 1 << endl;
		j++;						  
		outfile << current_ID + j << " 1 30000 1 2 " 
								  << right_start + 1 << " "
								  << right_end + 1 << endl;
		j++;
		
		left_start = left_end;
		right_start = right_end;		
	}
	
	/*
	for(int i = 0; i < print_box_edges.size(); i++){
		outfile << current_ID + i << " 1 0 " 
								  << find(print_box_edges.at(i)->pt1) + 1 << " "
								  << find(print_box_edges.at(i)->pt2) + 1 << endl;
	}
	*/
	
	//cout << "print_outer_bound done" << endl;
	
	
	outfile << "$ENDELM" << endl;

	outfile.close();
	return;
}

//printing the refinement table obtained (Manual mode)
void print_ref_table(Refinement *ref){
	//cout << "num levels: " << ref->table->size() << endl;
	for(int i = 0; i < ref->table->size(); i++){
		//cout << "Refinement level: " << i+1 << endl;
		for(int j = 0; j < ref->table->at(i)->size(); j++){
			//cout << "Table#: " << j+1 << endl;
			//cout << *(ref->table->at(i)->at(j)) << endl;
		}
	}
}

//deallocate memory used
void deallocate(){
	//deallocate
	/*for(int i = 0; i < precise_verts.size(); i++){
		////cout << *precise_verts.at(i) << endl;	//testing if input is correct
		if(precise_verts.at(i) != NULL)
			delete precise_verts.at(i);
	}*/
	
	precise_verts.clear();
	precise_verts_backup.clear();
	
	/*for(int i = 0; i < precise_edges.size(); i++){
		////cout << *precise_edges.at(i) << endl;	//testing if input is correct
		if(precise_edges.at(i) != NULL)
			delete precise_edges.at(i);
	}*/
	
	precise_edges.clear();
	
	/*for(int i = 0; i < segments.size(); i++){
		////cout << *segments.at(i) << endl;	//testing if input is correct
		if(segments.at(i) != NULL)
			delete segments.at(i);
	}*/
	segments.clear();
	
	delete GBOX_X_MIN;
	delete GBOX_X_MAX;
	delete GBOX_Y_MIN;
	delete GBOX_Y_MAX;
	delete INIT_SIZE;
	delete INIT_MIN;
	delete EPSILON;
	delete INTERNAL;
	delete PROGRAM_GENE_TYPE;

	return;
}

void read_parameters(const char* para_txt){
	ifstream infile;
	infile.open(para_txt);
	
	string junk;
	double temp;
	bool keep_mode;
	int gene_mode;
	
	//GBOX_X_MIN
	infile >> junk;
	infile >> temp;
	*GBOX_X_MIN = temp;
	
	//GBOX_X_MAX
	infile >> junk;
	infile >> temp;
	*GBOX_X_MAX = temp;
	
	//GBOX_Y_MIN
	infile >> junk;
	infile >> temp;
	*GBOX_Y_MIN  = temp;
	
	//GBOX_Y_MAX
	infile >> junk;
	infile >> temp;
	*GBOX_Y_MAX = temp;
	
	//INIT_MIN
	infile >> junk;
	infile >> temp;
	*INIT_MIN = temp;
	
	//INIT_SIZE
	int multipler = max( (*GBOX_X_MAX-*GBOX_X_MIN), (*GBOX_Y_MAX-*GBOX_Y_MIN) ) / *INIT_MIN + 1;
	*INIT_SIZE = *INIT_MIN;

	for(int i = 0; ((*INIT_MIN) * pow(4.0, i)) < multipler; i++){
		*INIT_SIZE = (*INIT_MIN) * pow(4.0, i+1);
	}
	
	//PROGRAM_GENE_TYPE
	infile >> junk;
	infile >> gene_mode;
	
	switch(gene_mode){
	case 0:
		*PROGRAM_GENE_TYPE = MANUAL;
		break;
	case 1:
		*PROGRAM_GENE_TYPE = MIN_SIZE;
		break;
	case 2:
		*PROGRAM_GENE_TYPE = ADAPTIVE;
		break;
	default:
		cout << "invalid program_gene_type" << endl;
		exit(1);
	}

	//INTERNAL
	infile >> junk;
	infile >> keep_mode;
	*INTERNAL = keep_mode;

	//EPSILON
	*EPSILON = 1e-10 * (*INIT_MIN);
	
	infile.close();
}

int main(int argc, char *argv[]){

	if(argc < 4){
		cout << "[Card_2d] Usage: InputGmshFile InputRefinementFile OutputMeshFile InputParameters" << endl;
		exit(1);
	}
		
	const char *INF_NAME = argv[1];
	const char *REF_NAME = argv[2];
	const char *OUTF_NAME = argv[3];
	const char *PARAMETERS = argv[4];
	
	/*
	cout << INF_NAME << endl;
	cout << REF_NAME << endl;
	cout << OUTF_NAME << endl;
	cout << PARAMETERS << endl;
	cout << "reading parameters" << endl;
	*/
	
	read_parameters(PARAMETERS);
	
	/*
	cout << *GBOX_X_MIN << endl;
	cout << *GBOX_X_MAX << endl;
	cout << *GBOX_Y_MIN << endl;
	cout << *GBOX_Y_MAX << endl;
	cout << *INIT_SIZE << endl;
	cout << *INIT_MIN << endl;
	cout << *INTERNAL << endl;
	cout << *EPSILON << endl;
	cout << *PROGRAM_GENE_TYPE << " TYPE: 0 is MANUAL, 1 is MIN_SIZE, 2 is ADATPIVE" << endl;
	*/
	
	//<--- initializing input parameters
	read_geo_file(INF_NAME);
	
	pt4D *GBOX = new pt4D(*GBOX_X_MIN, *GBOX_X_MAX, *GBOX_Y_MIN, *GBOX_Y_MAX);
	//end of initializing parameters --->
	
	//<--- gene 16*16 grid with min size 1 then crop it for init mesh
	mesh_coord = new vert_precise(*GBOX_X_MIN, *GBOX_Y_MIN, precise_verts.size());
	precise_verts.push_back(mesh_coord);
	
	//cout << "before segment ADT" << endl;
	segments_ADT = make_segment_ADT(segments, 0, GBOX);
	
	/*
	//print the input geo to cout
	//cout << "precise_verts size is: " << precise_verts.size() << endl;
	//cout << "segments size is: " << segments.size() << endl;
	for(int i = 0; i < segments.size(); i++){
		edge_stat(segments.at(i));
	}
	print_ADT(segments_ADT);
	*/
	
	final_mesh = new_generate_dummy_mesh(*INIT_SIZE, *INIT_MIN, mesh_coord->ID);
	
	if(*PROGRAM_GENE_TYPE == ADAPTIVE){
		new_crop_mesh_helper(final_mesh, GBOX);
	}
	
	final_mesh = new_crop_mesh(final_mesh, GBOX);
	
	//end of crop mesh --->
	
	
	//cout << "end of reading input files" << endl;
	
	//<--- generate mesh with ref_tables and surface
	//print_ref_table(ref_table);	//ref_table check
	
	//cout << "[Card_2d] before generating mesh" << endl;
	
	if(*PROGRAM_GENE_TYPE == MANUAL){
		Refinement *ref_table;
		ref_table = read_ref_file(REF_NAME);
		
		final_mesh = generate_mesh(segments_ADT, ref_table, final_mesh, NULL, NULL);
	}
	else if(*PROGRAM_GENE_TYPE == MIN_SIZE){
		vector<Minsize_Refinement*> *ref_minsize = new vector<Minsize_Refinement*>;
		min_size_ref_read(ref_minsize, REF_NAME);
		
		final_mesh = generate_mesh(segments_ADT, NULL, final_mesh, ref_cond_minsize, ref_minsize);
		
		while(ref_minsize->size() != 0){
			ref_minsize->pop_back();
		}
	}
	else if(*PROGRAM_GENE_TYPE == ADAPTIVE){
		Adaptive_variable *ref_adapt = adapt_ref_read(REF_NAME);
		
		final_mesh = generate_mesh(segments_ADT, NULL, final_mesh, ref_cond_adaptive, ref_adapt);
		
		delete ref_adapt;
	}
	//end of generate mesh --->

	//cout << "after gene mesh, before printing" << endl;
	
	print_msh(final_mesh, OUTF_NAME);
	
	deallocate();
	return 0;
}


