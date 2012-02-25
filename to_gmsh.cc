#include <iostream>
#include <fstream>
#include <vector>
#include <stack>
#include <string>
#include <queue>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <algorithm>

using namespace std;

void read_geo_file(const char* shape_txt, const char* filename, const char* out_txt){
	ifstream infile;
	infile.open(shape_txt);
	
	ofstream outfile;
	ofstream tempfile;
	ofstream tempfile2;
	outfile.open(filename);
	tempfile.open("to_gmsh.temp");
	tempfile2.open("to_gmsh.temp2");
	
	string junk;
	int int_junk;
	int type;
	
	outfile << "$MeshFormat" << endl;
	outfile << "2.2 0 8" << endl;
	outfile << "$EndMeshFormat" << endl;
	
	infile >> junk; //Nodes
	outfile << "$Nodes" << endl;
	
	int num_nodes;
	infile >> num_nodes;
	outfile << num_nodes;
	
	for(int i = 0; i <= num_nodes; i++){
		getline(infile, junk);
		outfile << junk << endl;
	}
 	
 	infile >> junk; //EndNodes
 	infile >> junk; //Edges
 	outfile << "$EndNodes" << endl;
 	outfile << "$Elements" << endl;
	outfile << "ELEMENTSIZE" << endl;
	
	int num_edges;
	int new_ID = num_nodes + 1;
	infile >> num_edges;
	
	for(int i = 0; i < num_edges; i++){
		int edge_ID;
		int n1;
		int n2;
	 	infile >> edge_ID;
	 	
	 	//edge attributes
	 	infile >> type;
	 	
	 	if(type == 2){
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3;
	 		infile >> v1 >> v2 >> v3;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	
	 	if(type == 3){
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		outfile << new_ID << " 3 0";
	 		for(int i = 0; i < 4; i++){
	 			infile >> int_junk;
	 			outfile <<  " " << int_junk;
	 		}
	 		outfile << endl;
	 		new_ID++;
	 	}
	 	
	 	else if(type == 5){ //4
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4;
	 		infile >> v1 >> v2 >> v3 >> v4;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	else if(type == 6){ //5
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v1 << endl;
	 		new_ID++;
	 	}
 		else if(type == 7){ //6
 			infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5, v6;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v6 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v6 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	else if(type == 8){ //7
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5, v6, v7;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v6 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v6 << " " << v7 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v7 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	
	 	else if(type == 1){ //edge
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2;
	 		infile >> v1 >> v2;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 	}
	 	else{
	 		getline(infile, junk);
	 	}
	 
	}
	
	outfile << "$EndElements" << endl;
	
	int num_ele = new_ID - num_nodes - 1;
	
	tempfile << num_ele;
	tempfile2 << out_txt;
	
	infile.close();
	outfile.close();
	tempfile.close();
	return;
}
/*
void read_geo_file2(const char* shape_txt, const char* filename){
	ifstream infile;
	infile.open(shape_txt);
	
	ofstream outfile;
	outfile.open(filename);
	
	string junk;
	int int_junk;
	int type;
	
	outfile << "$MeshFormat" << endl;
	outfile << "2.2 0 8" << endl;
	outfile << "$EndMeshFormat" << endl;
	
	infile >> junk; //Nodes
	outfile << "$Nodes" << endl;
	
	int num_nodes;
	infile >> num_nodes;
	outfile << num_nodes;
	
	for(int i = 0; i <= num_nodes; i++){
		getline(infile, junk);
		outfile << junk << endl;
	}
 	
 	infile >> junk; //EndNodes
 	infile >> junk; //Edges
 	outfile << "$EndNodes" << endl;
 	outfile << "$Elements" << endl;

	int num_edges;
	int new_ID = num_nodes + 1;
	infile >> num_edges;
	
	for(int i = 0; i < num_edges; i++){
		int edge_ID;
		int n1;
		int n2;
	 	infile >> edge_ID;
	 	
	 	//edge attributes
	 	infile >> type;
	 	
	 	if(type == 3){
	 		infile >> int_junk;
	 		//infile >> int_junk;
	 		//infile >> int_junk;
	 		outfile << new_ID << " 3 0";
	 		for(int i = 0; i < 4; i++){
	 			infile >> int_junk;
	 			outfile <<  " " << int_junk;
	 		}
	 		outfile << endl;
	 		new_ID++;
	 	}
	 	
	 	else if(type == 5){ //4
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4;
	 		infile >> v1 >> v2 >> v3 >> v4;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	else if(type == 6){ //5
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v1 << endl;
	 		new_ID++;
	 	}
 		else if(type == 7){ //6
 			infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5, v6;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v6 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v6 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	else if(type == 8){ //7
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		infile >> int_junk;
	 		int v1, v2, v3, v4 ,v5, v6, v7;
	 		infile >> v1 >> v2 >> v3 >> v4 >> v5 >> v6 >> v7;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v2 << " " << v3 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v3 << " " << v4 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v4 << " " << v5 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v5 << " " << v6 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v6 << " " << v7 << endl;
	 		new_ID++;
	 		outfile << new_ID << " 1 0 " << v7 << " " << v1 << endl;
	 		new_ID++;
	 	}
	 	else if(type == 1){ //edge
	 		infile >> int_junk;
	 		//infile >> int_junk;
	 		//infile >> int_junk;
	 		int v1, v2;
	 		infile >> v1 >> v2;
	 		outfile << new_ID << " 1 0 " << v1 << " " << v2 << endl;
	 		new_ID++;
	 	}
	 	else{
	 		getline(infile, junk);
	 	}
	 
	}
	
	outfile << "$EndElements" << endl;
	 	
	infile.close();
	return;
}
void read_geo_file3(const char* shape_txt, const char* filename){
	ifstream infile;
	infile.open(shape_txt);
	
	ofstream outfile;
	outfile.open(filename);
	
	string junk;
	int int_junk;
	int type;

	
	infile >> junk; //Nodes
	outfile << "$NOD" << endl;
	
	int num_nodes;
	infile >> num_nodes;
	outfile << num_nodes;
	
	for(int i = 0; i <= num_nodes; i++){
		getline(infile, junk);
		outfile << junk << endl;
	}
 	
 	infile >> junk; //EndNodes
 	infile >> junk; //Edges
 	outfile << "$ENDNOD" << endl;
 	outfile << "$ELM" << endl;

	int num_edges;
	int new_ID = num_nodes + 1;
	infile >> num_edges;
	
	int currentID = 1;
	
	for(int i = 0; i <= num_edges; i++){
		int edge_ID;
		int n1;
		int n2;
	 	infile >> edge_ID;
	
		getline(infile, junk);
		
		outfile << currentID << " " << junk << endl;
		
		currentID++;
	 
	}
	
	outfile << "$ENDELM" << endl;
	 	
	infile.close();
	return;
}
void shift_file(const char* shape_txt, const char* filename){
	ifstream infile;
	infile.open(shape_txt);
	
	ofstream outfile;
	outfile.open(filename);
	
	string junk;
	int int_junk;
	
	//removing version number
	getline(infile, junk);
	getline(infile, junk);
	getline(infile, junk);
	infile >> junk; //Nodes
	
	double shift = 200;
	
	int num_nodes;
	infile >> num_nodes;
	outfile << num_nodes << endl;
	
	for(int i = 0; i < num_nodes; i++){
		int node_ID;
		double x;
		double y;
		double z;
	 	infile >> node_ID;
	 	infile >> x;
	 	infile >> y;
	 	infile >> z;
		
		y += 0.35;
	 	outfile << (node_ID+shift) << " " << x << " " << y << " " << z << endl;
	}
	
	infile >> junk; //EndNodes
 	infile >> junk; //Edges

	int num_edges;
	infile >> num_edges;
	
	outfile << num_edges << endl;
	
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
	 	
	 	outfile << (edge_ID+shift) << " 1 0 " << (n1+shift) << " " << (n2+shift) << endl;
	}
	 	 	
	 	
	infile.close();
	outfile.close();
	return;
}
*/

int main(int argc, char *argv[]){
	
	if(argc < 2){
		cout << "Usage: InputFileName OutputFileName" << endl;
		exit(1);
	}
	
	const char *OUTF_NAME = argv[1];
	const char *OUTF_NAME_GMSH = argv[2];
	
	read_geo_file(OUTF_NAME, "to_gmsh_unprocess", OUTF_NAME_GMSH);
	
	//read_geo_file2("output_wing_adaptive.msh", "output.msh");
	//read_geo_file3("output_wing_mesh.msh", "output.msh");
	//shift_file("input.msh", "out.msh");
}