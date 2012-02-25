const double GBOX_X_MIN = -15.8;
const double GBOX_X_MAX = 16.2;
const double GBOX_Y_MIN = -15.8;
const double GBOX_Y_MAX = 16.2;
const double INIT_SIZE = 32;
const double INIT_MIN = 2;
const char *INF_NAME = "two_wing.msh";
const char *REF_NAME = "adaptive_conds.txt";
const char *OUTF_NAME = "output_wings.msh";
const char *OUTF_NAME_GMSH = "output_wings_gmsh.msh";
const bool INTERNAL = false;

enum gene_type {MANUAL, MIN_SIZE, ADAPTIVE, TYPE_COUNT};
const gene_type PROGRAM_GENE_TYPE = ADAPTIVE;
