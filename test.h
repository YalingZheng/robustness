#ifndef neutralheader
#define neutralheader

// local includes
#include "stocc.h"
#include "./Mersenne-1.1/MersenneTwister.h"
#include "red_black_tree.h"

// namespace
using namespace std;

// useful constants
#define TRUE 1
#define FALSE 0

#define LINESIZE 128
#define PROTEINSIZE 512

#define GEN 1024
#define POPN 1024

enum HASHTAGS { HASHUPDATE=100, HASHQUERY, ENERGYTRANSFER, HASHREQUEST, SERVEREXIT  };
enum ROBTAGS { DNATRANSFER=10, IDLEWORKER, NEXTDNA, PARTIALENERGY, LOOPEXIT };
enum DISPATCHTAGS { TRANSFER=50, STABLETRANSFER, IDLE, NEXT, EXIT };

// random number generator
MTRand mtrand1;

// constant for amino acid order in delta(G) calculation.
const char AA[20] = {'a','e','q','d','n','l','g','k','s','v','r','t','p','i','m','f','y','c','w','h'};

// constant contact energy for delta(G) calculation.
float contactenergy[20][20] = {
    {-0.0479, 0.1346, -0.0457, 0.1018, 0.1049, -0.1711, 0.1844, 0.0691, 0.0464, -0.1431, 0.1049, 0.0310, 0.1462, -0.0737, -0.0847, -0.1119, -0.1408, -0.1085, -0.0880, 0.0266},
    {0.1346, 0.1259, -0.0413, 0.1581, 0.1146, 0.0802, 0.2311, -0.2403, 0.0823, 0.1010, -0.3511, 0.0675, 0.2241, 0.1103, 0.0637, 0.0885, -0.0522, 0.1550, -0.0967, -0.0827},
    {-0.0457, -0.0413, -0.0550, -0.0728, -0.0050, -0.0172, 0.1710, -0.0735, 0.1169, 0.1061, 0.0059, -0.0243, 0.1127, -0.0480, -0.1038, -0.0171, -0.1431, 0.0715, -0.0540, -0.0125},
    {0.1018, 0.1581, -0.0728, 0.0840, 0.0192, 0.2673, 0.1115, -0.1154, 0.0424, 0.2728, -0.1859, 0.1043, 0.2386, 0.1892, -0.0197, 0.0827, -0.1165, 0.1169, -0.0124, -0.0749},
    {0.1049, 0.1146, -0.0050, 0.0192, -0.0917, 0.0890, 0.1196, -0.0381, 0.1452, 0.1180, -0.0150, 0.0155, 0.1560, 0.1485, 0.0124, 0.0018, -0.1149, -0.0844, -0.0250, 0.0386},
    {-0.1711, 0.0802, -0.0172, 0.2673, 0.0890, -0.5067, 0.0782, 0.0543, 0.0959, -0.4593, -0.0651, -0.0316, 0.0745, -0.5112, -0.1822, -0.5450, -0.2614, -0.1305, -0.2639, -0.0169},
    {0.1844, 0.2311, 0.1710, 0.1115, 0.1196, 0.0782, 0.2219, 0.1963, 0.1075, 0.1859, -0.0251, 0.1763, 0.2131, 0.1174, -0.0573, 0.0789, -0.0176, -0.0982, -0.1567, 0.0979},
    {0.0691, -0.2403, -0.0735, -0.1154, -0.0381, 0.0543, 0.1963, 0.1216, 0.1690, 0.0609, 0.0839, 0.0467, 0.1099, 0.0682, 0.0866, -0.0416, -0.1120, -0.0330, -0.1152, 0.0390},
    {0.0464, 0.0823, 0.1169, 0.0424, 0.1452, 0.0959, 0.1075, 0.1690, 0.0941, 0.1766, 0.0442, 0.0228, 0.1626, 0.0332, 0.0185, 0.0398, 0.0214, -0.0132, -0.0802, -0.0005},
    {-0.1431, 0.1010, 0.1061, 0.2728, 0.1180, -0.4593, 0.1859, 0.0609, 0.1766, -0.5193, 0.0475, 0.0119, 0.0868, -0.4223, -0.2127, -0.4001, -0.2792, -0.2349, -0.2898, -0.0039},
    {0.1049, -0.3511, 0.0059, -0.1859, -0.0150, -0.0651, -0.0251, 0.0839, 0.0442, 0.0475, 0.0306, -0.0210, -0.0614, -0.0266, -0.0163, -0.0904, -0.1369, 0.0544, -0.2070, -0.0184},
    {0.0310, 0.0675, -0.0243, 0.1043, 0.0155, -0.0316, 0.1763, 0.0467, 0.0228, 0.0119, -0.0210, 0.0150, 0.1908, -0.0700, 0.0018, -0.1120, -0.1445, -0.0013, 0.0052, 0.0681},
    {0.1462, 0.2241, 0.1127, 0.2386, 0.1560, 0.0745, 0.2131, 0.1099, 0.1626, 0.0868, -0.0614, 0.1908, 0.1077, 0.0882, -0.0069, -0.0604, -0.1326, 0.0545, -0.0910, 0.0295},
    {-0.0737, 0.1103, -0.0480, 0.1892, 0.1485, -0.5112, 0.1174, 0.0682, 0.0332, -0.4223, -0.0266, -0.0700, 0.0882, -0.5852, -0.2137, -0.3791, -0.3164, -0.2235, -0.1961, -0.0326},
    {-0.0847, 0.0637, -0.1038, -0.0197, 0.0124, -0.1822, -0.0573, 0.0866, 0.0185, -0.2127, -0.0163, 0.0018, -0.0069, -0.2137, -0.1059, -0.1785, -0.1621, -0.0557, -0.0775, -0.0345},
    {-0.1119, 0.0885, -0.0171, 0.0827, 0.0018, -0.5450, 0.0789, -0.0416, 0.0398, -0.4001, -0.0904, -0.1120, -0.0604, -0.3791, -0.1785, -0.3088, -0.4212, -0.3262, -0.3405, -0.1250},
    {-0.1408, -0.0522, -0.1431, -0.1165, -0.1149, -0.2614, -0.0176, -0.1120, 0.0214, -0.2792, -0.1369, -0.1445, -0.1326, -0.3164, -0.1621, -0.4212, -0.2793, -0.2444, -0.3209, -0.1976},
    {-0.1085, 0.1550, 0.0715, 0.1169, -0.0844, -0.1305, -0.0982, -0.0330, -0.013, -0.2349, 0.0544, -0.0013, 0.0545, -0.2235, -0.0557, -0.3262, -0.2444, -1.0442, -0.1176, -0.0701},
    {-0.0880, -0.0967, -0.0540, -0.0124, -0.0250, -0.2639, -0.1567, -0.1152, -0.0802, -0.2898, -0.2070, 0.0052, -0.0910, -0.1961, -0.0775, -0.3405, -0.3209, -0.1176, -0.1066, -0.0200},
    {0.0266, -0.0827, -0.0125, -0.0749, 0.0386, -0.0169, 0.0979, 0.0390, -0.0005, -0.0039, -0.0184, 0.0681, 0.0295, -0.0326, -0.0345, -0.1250, -0.1976, -0.0701, -0.0200, 0.0005}
};

// atom structure
struct Atom {
    float x;
    float y;
    float z;
    char seq;
    int index;
};

// hash structure
struct HashTable {
    char dna[PROTEINSIZE*3];
    unsigned long hash;
    int index;
    int duplicates;
    float energy;
    struct HashTable* next;
};

// options structure
struct Opt {
    int gen;
    int pop;
    char aafile[256];
    char pdbfile[256];
    bool robustness;
    char fdna[512];
};

// type definitions
typedef struct HashTable hashtable_t;
typedef struct Atom atom_t;
typedef struct Opt opt_t;

// global hashtable
hashtable_t *head;
hashtable_t *tail;

// global stable hashtable
hashtable_t *stablehead;
hashtable_t* stabletail;

//global options
opt_t* options;

// function prototypes
float energy( atom_t*, int, int& );
void print_atom( atom_t* ) ;
void print_contactenergy( float [][20] );
void read_pdbfile( const char*, atom_t*, int& );
void read_dnafile( const char* , char* );
void write_dnafile( const char*, char* );
void Mutate( char[], int );
char translate( char* );
void translate_protein( char[], char[] , int );
float robustness( hashtable_t* );
void submit_robustness( char* );
float dna_energy( char* );

hashtable_t* inithash( char[], unsigned long int, bool );
void add( hashtable_t*, hashtable_t**, hashtable_t** );
hashtable_t* search( unsigned long int, hashtable_t**, bool* );
hashtable_t* update( char[], hashtable_t**, hashtable_t**, bool );
void print_hash( hashtable_t* );
void print_hashtable( hashtable_t **, bool );
unsigned long sdbm( unsigned char[] );
void clearhash( hashtable_t**, hashtable_t** );

hashtable_t* mostcommondna( hashtable_t** );
int get_robustnesstable( hashtable_t*, char[][PROTEINSIZE*3] );
void HashServer();

opt_t* cmd_options(int, char**);
void usage(char*);

#endif
