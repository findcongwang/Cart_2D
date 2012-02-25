//dimensions of the global box:
const double GBOX_X_MIN = -15.8;
const double GBOX_X_MAX = 16.2;
const double GBOX_Y_MIN = -15.8;
const double GBOX_Y_MAX = 16.2;

const double INIT_SIZE = 32;
//size of initial square mesh, 
//INIT_SIZE = INIT_MIN * 4^j > max(X_MAX-X_MIN, Y_MAX-T_MIN), for the smallest j

const double INIT_MIN = 2;	//size of initial mesh cells
const char *INF_NAME = "two_wing.msh";	//input msh file for geometry (gmsh format)
const char *REF_NAME = "adaptive_conds.txt";	//parameters for the generation mode (Adaptive)
const char *OUTF_NAME = "output_wings.msh";	//name for the output msh file
const char *OUTF_NAME_GMSH = "output_wings_gmsh.msh";	//name of the gmsh format output file
const bool INTERNAL = false;	//false: external cells are kept; true: internal cells are kept

enum gene_type {MANUAL, MIN_SIZE, ADAPTIVE, TYPE_COUNT};
const gene_type PROGRAM_GENE_TYPE = ADAPTIVE;	//generation mode
