const double GBOX_X_MIN = -4.2;
const double GBOX_X_MAX = 3.8;
const double GBOX_Y_MIN = -4.2;
const double GBOX_Y_MAX = 3.8;
const double INIT_SIZE = 8;
const double INIT_MIN = 0.5;
const char *INF_NAME = "toy.msh";
const char *REF_NAME = "multi.txt";
const char *OUTF_NAME = "output_wing_adaptive.msh";
const char *OUTSOLF_NAME = "output_wing_mesh.msh";
const bool INTERNAL = false;

enum gene_type {MANUAL, MIN_SIZE, ADAPTIVE, TYPE_COUNT};
const gene_type PROGRAM_GENE_TYPE = ADAPTIVE;
