const double GBOX_X_MIN = -8000.2;
const double GBOX_X_MAX = 23999.8;
const double GBOX_Y_MIN = -8000.2;
const double GBOX_Y_MAX = 23999.8;
const double INIT_SIZE = 32000;
const double INIT_MIN = 8000;
const char *INF_NAME = "lake_withholes.msh";
const char *REF_NAME = "adaptive_conds.txt";
const char *OUTF_NAME = "output_lake.msh";
const char *OUTF_NAME_GMSH = "output_lake_gmsh.msh";
const bool INTERNAL = true;

enum gene_type {MANUAL, MIN_SIZE, ADAPTIVE, TYPE_COUNT};
const gene_type PROGRAM_GENE_TYPE = ADAPTIVE;
