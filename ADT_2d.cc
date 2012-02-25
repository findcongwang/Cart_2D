/*
 * ADT.cc
 *
 *  Created on: Sep 7, 2011
 *      Author: c233wang
 */

#include <iostream>
#include "ADT_2d.h"

ADTNode::ADTNode(Polyhedron* poly_in, int depth_in, pt4D *region_in){
	poly = poly_in;
	depth = depth_in;
	region = region_in;
	L_branch = NULL;
	R_branch = NULL;
}
