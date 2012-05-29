#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <ctime>
#include <cctype>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

//#include "mpi.h"

#include "neutralrobustness.h"

using namespace std;

unsigned int energy_calls=0;

int main ( int argc, char* argv[] ) {

    atom_t atoms[PROTEINSIZE*10];
    int atomsize, aasize, diff;

    float eng, energyoutput;

    int g, i, j, k, m, z, mut[POPN];

    char dnainput[PROTEINSIZE*3] = {'\0'};
    char dna[POPN][PROTEINSIZE*3] = {'\0'};
    char stable[POPN][PROTEINSIZE*3] = {'\0'};
    char cmd[512] = {'\0'};

    hashtable_t* sh;
    hashtable_t* h;

    int seed = time(0);
    TRandomMersenne rg(seed);

    options = cmd_options(argc, argv);

    if( options->robustness ) {
        FILE* fstable = fopen( options->fdna, "r" );
        while( fgets( dnainput, PROTEINSIZE*3, fstable ) ){
            for( i=0; i < PROTEINSIZE*3; i++ ) {
                if( dnainput[i] == '\n' ) {
                    dnainput[i] = '\0';
                    break;
                }
            }
            h = update( dnainput, &head, &tail, TRUE );
            print_hash( h );
            printf( "robustness:: %5.3f\n\n", (float) robustness( h ) );
        }
        return 0;
    }

    // calculate original protein contact energy
    sprintf( cmd, "./Scwrl4 -p ./Scwrl4.ini -i %s -o refined.pdb -s %s -h -v  > logfile1", options->pdbfile, options->aafile );
    system( cmd );
    read_pdbfile( (const char*) "refined.pdb", atoms, atomsize );
    eng = energy( atoms, atomsize, aasize );
    printf( "%d Generations, %d Individuals\n\n", options->gen, options->pop );
    printf( "initial energy = %f\n\n", eng );

    /**** start evolution ****/


    // read input dna and replenish first population.
    read_dnafile( (const char*) "dna", dnainput );
    for( i=0; i < options->pop; i++ ) {
        strcpy( dna[i], dnainput );
    }

    for( g=0; g < options->gen; g++ ) {
        printf( "\nGeneration #%d:: \n", g+1 );

        k=0;
        seed = time(0);
        rg.RandomInit(seed);
        StochasticLib1 stochastic(seed + rg.IRandom(0,options->pop));
        for( j=0; j < options->pop; j++ ) {
            // mutate organism
            mut[j] = stochastic.Poisson(1);
            Mutate( dna[j], mut[j] );

            // a hash table is used to avoid recalculating energies
            h = update( dna[j], &head, &tail, TRUE );

            // cull
            if ( (h->energy <= eng + 0.055) && (h->energy >= eng - 0.055) ) {
                // individuals that reach here are stable
                strcpy( stable[k++], h->dna );
                sh = update( h->dna, &stablehead, &stabletail, FALSE );
                sh->energy = h->energy;
                print_hash( h );

                /** robustness could be submitted from here for this individual **/
            }
        }

        /** robustness could be submitted from here for this stable population **/

        // replenish
        diff = options->pop-k;
        printf("\nstablesize = %d, populationsize = %d, diff = %d\n", k, options->pop, diff);
        for( m=0; m < diff; m++ ) {
            z = rg.IRandom( 0, k-1 );
            strcpy( stable[k++], stable[z] );
        }

        for( i=0; i < options->pop; i++ ) {
            strcpy( dna[i], stable[i] );
        }
    }

    /** robustness could be submitted from here for all stable individuals **/
    printf("\nRobustness ----\n");
    int total=0;
    for( h=stablehead; h != NULL; h=h->next ) {
        total++;
    }
    printf("total number of individuals = %d\n\n", total);

    FILE* fstable = NULL;
    char fname[512];
    char hostname[512];

    i=0; j=0;
    gethostname( hostname, 512 );

   // struct stat st;
   // if(stat("./cmd",&st) != 0) {
   //     mkdir( "./cmd",  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
   //     mkdir( "./cmd/cmd",  S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH );
   // }

    for( h=stablehead; h != NULL; h=h->next ) {
        if( i%100 == 0 ) {
            if( fstable != NULL ) {
                fclose( fstable );
                submit_robustness( fname );
            }
            sprintf( fname, "./%s-%d.dna", hostname, j++ );
            fstable = fopen( fname, "w+" );
        }
        print_hash( h );
        fprintf( fstable, "%s\n", h->dna );
        i++;
    }
    submit_robustness( fname );
    fclose( fstable );

    printf("\n\nthere were %d calls to dna_energy\n", energy_calls);

    return 0;
}

void submit_robustness( char* fdna ) {
    FILE* condor;
    char fname[64];
    char cmd[96];

    sprintf( fname, "%s.cmd", fdna );
    condor = fopen( fname, "w" );
    if ( condor != NULL ) {
        fprintf( condor, "####################\n" );
        fprintf( condor, "##\n" );
        fprintf( condor, "## Robustness command file\n" );
        fprintf( condor, "##\n" );
        fprintf( condor, "####################\n" );
        fprintf( condor, "\n" );
        fprintf( condor, "universe        = vanilla\n" );
        fprintf( condor, "+ProjectName = \"TG-STA120004S\"\n");
	fprintf( condor, "executable      = neutralrobustness\n" );
        fprintf( condor, "output          = %s.out\n", fname);
        fprintf( condor, "error           = %s.err\n", fname);
        fprintf( condor, "log             = %s.log\n", fname);
        fprintf( condor, "arguments       = -r %s\n", fdna);
	fprintf( condor, "requirements = (Memory>=1)&&(Arch==\"X86_64\")\n");

        fprintf( condor, "transfer_input_files = %s,Scwrl4,bbDepRotLib.bin,Scwrl4.ini,refined.pdb, temp.pdb\n", fdna);
	fprintf( condor, "should_transfer_files = IF_NEEDED\n" );
        fprintf( condor, "when_to_transfer_output = ON_EXIT\n" );
        fprintf( condor, "queue\n" );
        fclose( condor );
    }

//sprintf( cmd, "condor_submit %s", fname );
  //  system( cmd );
}

float robustness( hashtable_t* h ) {

    char permutations[PROTEINSIZE*3][PROTEINSIZE*3] = {'\0'};
    float energyoutput=0.0;
    int count = get_robustnesstable( h, permutations );

    for( int r=0; r < count; r++ ) {
        energyoutput += dna_energy( permutations[r] );
    }

    return energyoutput/count;
}

float energy( atom_t* atoms, int size, int& seqsize ) {

    int aacontacts[PROTEINSIZE][PROTEINSIZE] = {0};

    int i, j, l, s;

    float dist;
    float eng = 0.0;
    int f=0, m=0;

    float xdiffsquared;
    float ydiffsquared;
    float zdiffsquared;

    char aminoacidseq[PROTEINSIZE] = {0};

    seqsize=0;
    for(int q=0;q<size;q++)
    {
        if(atoms[q].index!=atoms[q+1].index || size < q + 2 )
            aminoacidseq[seqsize++] = atoms[q].seq;
    } /* end for */

    for( i = 0; i < size; i++ ) {
        for( j = i; j < size; j++ ) {
            if( abs(atoms[i].index - atoms[j].index) >= 2 ) {
                xdiffsquared = (atoms[i].x - atoms[j].x) * (atoms[i].x - atoms[j].x);
                ydiffsquared = (atoms[i].y - atoms[j].y) * (atoms[i].y - atoms[j].y);
                zdiffsquared = (atoms[i].z - atoms[j].z) * (atoms[i].z - atoms[j].z);

                dist = xdiffsquared + ydiffsquared + zdiffsquared;

                if ( dist <= 20.25 ) {
                    if ( aacontacts[atoms[i].index][atoms[j].index] == 0 ) {
                        aacontacts[atoms[i].index][atoms[j].index] = 1;
                        for( l=0; l < 20; l++ ) {
                            if(AA[l] == aminoacidseq[atoms[i].index-1]) {
                                m=l;
                                break;
                            }
                        }
                        for( s=0; s < 20; s++ ) {
                            if(AA[s] == aminoacidseq[atoms[j].index-1]) {
                                f=s;
                                break;
                            }
                        }

                        eng += contactenergy[m][f];
                        m=0;
                        f=0;
                    }
                }
            }
        } /* end inner for */
    } /* end outer for */

    return eng;
}

void print_atom( atom_t* atom ) {
    printf( "%c, %d, %f, %f, %f\n", atom->seq, atom->index, atom->x, atom->y, atom->z );
}

void print_contactenergy( float contactenergy[][20] ) {

    int i, j;
    for ( i = 0; i < 20; i++ ) {
        for ( j = 0; j < 20; j++ ) {
            printf( "%f\n", contactenergy[i][j] );
        }
    }
}

void read_dnafile( const char* filename, char* dna ) {

    int i;
    char line[PROTEINSIZE*3];
    FILE* dnafile;

    dnafile = fopen( filename, "r" );
    if ( dnafile != NULL ) {
        if ( fgets( line, PROTEINSIZE*3, dnafile ) != NULL ) {
            if( line != NULL ) {
                strncpy( dna, line, PROTEINSIZE*3 );
            }
        }
    } else {
        printf( "No dna file found.. aborting.\n" );
        exit(-1);
    }

    for( i=0; i < strlen(dna); i++) {
        if( dna[i] == '\n' ) {
            dna[i] = '\0';
        }
    }

    return;
}

void write_dnafile( const char* filename, char* dna ) {

    FILE* dnafile;

    dnafile = fopen( filename, "w" );
    if ( dnafile != NULL ) {
        fprintf( dnafile, "%s%c", dna, EOF );
    } else {
        printf( "No dna file found.. aborting.\n" );
        exit(-1);
    }
    fclose(dnafile);

    return;
}

void read_pdbfile( const char* filename, atom_t* atoms, int& aasize ) {

    char aa;
    char line[LINESIZE];
    char* substr;
    FILE* pdbfile;

    aasize = 0;
    pdbfile = fopen( filename, "r" );
    if ( pdbfile != NULL ) {
        while ( fgets( line, LINESIZE, pdbfile ) != NULL ) {
            if ( line != NULL && strstr( line, (const char*)"ATOM") && line[21] == 'A' ) {
                // find aminoacid one-letter.
                substr = &line[17];
                if ( strncmp( substr, (const char*)"ALA", 3 ) == 0 ) {
                    aa = 'a';
                } else if ( strncmp( substr, (const char*)"ARG", 3 ) == 0 ) {
                    aa = 'r';
                } else if ( strncmp( substr, (const char*)"ASN", 3 ) == 0 ) {
                    aa = 'n';
                } else if ( strncmp( substr, (const char*)"ASP", 3 ) == 0 ) {
                    aa = 'd';
                } else if ( strncmp( substr, (const char*)"CYS", 3 ) == 0 ) {
                    aa = 'c';
                } else if ( strncmp( substr, (const char*)"GLN", 3 ) == 0 ) {
                    aa = 'q';
                } else if ( strncmp( substr, (const char*)"GLU", 3 ) == 0 ) {
                    aa = 'e';
                } else if ( strncmp( substr, (const char*)"GLY", 3 ) == 0 ) {
                    aa = 'g';
                } else if ( strncmp( substr, (const char*)"HIS", 3 ) == 0 ) {
                    aa = 'h';
                } else if ( strncmp( substr, (const char*)"ILE", 3 ) == 0 ) {
                    aa = 'i';
                } else if ( strncmp( substr, (const char*)"LEU", 3 ) == 0 ) {
                    aa = 'l';
                } else if ( strncmp( substr, (const char*)"LYS", 3 ) == 0 ) {
                    aa = 'k';
                } else if ( strncmp( substr, (const char*)"MET", 3 ) == 0 ) {
                    aa = 'm';
                } else if ( strncmp( substr, (const char*)"PHE", 3 ) == 0 ) {
                    aa = 'f';
                } else if ( strncmp( substr, (const char*)"PRO", 3 ) == 0 ) {
                    aa = 'p';
                } else if ( strncmp( substr, (const char*)"SER", 3 ) == 0 ) {
                    aa = 's';
                } else if ( strncmp( substr, (const char*)"THR", 3 ) == 0 ) {
                    aa = 't';
                } else if ( strncmp( substr, (const char*)"TRP", 3 ) == 0 ) {
                    aa = 'w';
                } else if ( strncmp( substr, (const char*)"TYR", 3 ) == 0 ) {
                    aa = 'y';
                } else if ( strncmp( substr, (const char*)"VAL", 3 ) == 0 ) {
                    aa = 'v';
                } else {
                    printf( "Unknown aminoacid name.. aborting.\n" );
                    exit(0);
                }
                atoms[aasize].seq = aa;

                // find aminoacid number.
                substr = &line[23];
                sscanf( substr, "%d", &atoms[aasize].index );

                // find atom x-coordinate.
                substr = &line[31];
                sscanf( substr, "%f", &atoms[aasize].x );

                // find atom y-coordinate.
                substr = &line[39];
                sscanf( substr, "%f", &atoms[aasize].y );

                // find atom z-coordinate.
                substr = &line[47];
                sscanf( substr, "%f", &atoms[aasize].z );

                //print_atom(&atoms[aasize]);
                aasize++;
            }
        }

        //print_contactenergy(contactenergy);

    } else {
        printf( "No refined pdb (%s) file found.. aborting.\n", filename );
        exit(-1);
    }

    fclose(pdbfile);
    return;
}

void Mutate(char d[], int m) {

    int i, nuc, nu, trans;

    if( m != 0 ) {
        for( i=0; i<m; i++ ) {
            nuc = mtrand1.randInt( (int)strlen(d) - 1 );
            nu = mtrand1.randInt( 2 );
            trans = mtrand1.randInt( 2 );

            if(d[nuc] == 'T') {
                if(nu == 0 || nu == 1)
                    d[nuc] = 'C';
                else if(trans == 0)
                    d[nuc] = 'A';
                else
                    d[nuc] = 'G';
            } else if(d[nuc] == 'A') {
                if(nu == 0 || nu == 1)
                    d[nuc] = 'G';
                else if(trans == 0)
                    d[nuc] = 'T';
                else
                    d[nuc] = 'C';
            } else if(d[nuc] == 'C') {
                if(nu == 0 || nu == 1)
                    d[nuc] = 'T';
                else if(trans == 0)
                    d[nuc] = 'A';
                else
                    d[nuc] = 'G';
            } else if(d[nuc] == 'G') {
                if(nu == 0 || nu == 1)
                    d[nuc] = 'A';
                else if(trans == 0)
                    d[nuc] = 'T';
                else
                    d[nuc] = 'C';
            }
        }
    }

    return;
}

char translate(char* codon) {

    if      (strncmp( codon, (const char*) "TCA", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "TCC", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "TCG", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "TCT", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "TTC", 3 ) == 0)      return 'f';
    else if (strncmp( codon, (const char*) "TTT", 3 ) == 0)      return 'f';
    else if (strncmp( codon, (const char*) "TTA", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "TTG", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "TAC", 3 ) == 0)      return 'y';
    else if (strncmp( codon, (const char*) "TAT", 3 ) == 0)      return 'y';
    else if (strncmp( codon, (const char*) "TAA", 3 ) == 0)      return '_';
    else if (strncmp( codon, (const char*) "TAG", 3 ) == 0)      return '_';
    else if (strncmp( codon, (const char*) "TGC", 3 ) == 0)      return 'c';
    else if (strncmp( codon, (const char*) "TGT", 3 ) == 0)      return 'c';
    else if (strncmp( codon, (const char*) "TGA", 3 ) == 0)      return '_';
    else if (strncmp( codon, (const char*) "TGG", 3 ) == 0)      return 'w';
    else if (strncmp( codon, (const char*) "CTA", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "CTC", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "CTG", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "CTT", 3 ) == 0)      return 'l';
    else if (strncmp( codon, (const char*) "CCA", 3 ) == 0)      return 'p';
    else if (strncmp( codon, (const char*) "CCC", 3 ) == 0)      return 'p';
    else if (strncmp( codon, (const char*) "CCG", 3 ) == 0)      return 'p';
    else if (strncmp( codon, (const char*) "CCT", 3 ) == 0)      return 'p';
    else if (strncmp( codon, (const char*) "CAC", 3 ) == 0)      return 'h';
    else if (strncmp( codon, (const char*) "CAT", 3 ) == 0)      return 'h';
    else if (strncmp( codon, (const char*) "CAA", 3 ) == 0)      return 'q';
    else if (strncmp( codon, (const char*) "CAG", 3 ) == 0)      return 'q';
    else if (strncmp( codon, (const char*) "CGA", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "CGC", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "CGG", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "CGT", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "ATA", 3 ) == 0)      return 'i';
    else if (strncmp( codon, (const char*) "ATC", 3 ) == 0)      return 'i';
    else if (strncmp( codon, (const char*) "ATT", 3 ) == 0)      return 'i';
    else if (strncmp( codon, (const char*) "ATG", 3 ) == 0)      return 'm';
    else if (strncmp( codon, (const char*) "ACA", 3 ) == 0)      return 't';
    else if (strncmp( codon, (const char*) "ACC", 3 ) == 0)      return 't';
    else if (strncmp( codon, (const char*) "ACG", 3 ) == 0)      return 't';
    else if (strncmp( codon, (const char*) "ACT", 3 ) == 0)      return 't';
    else if (strncmp( codon, (const char*) "AAC", 3 ) == 0)      return 'n';
    else if (strncmp( codon, (const char*) "AAT", 3 ) == 0)      return 'n';
    else if (strncmp( codon, (const char*) "AAA", 3 ) == 0)      return 'k';
    else if (strncmp( codon, (const char*) "AAG", 3 ) == 0)      return 'k';
    else if (strncmp( codon, (const char*) "AGC", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "AGT", 3 ) == 0)      return 's';
    else if (strncmp( codon, (const char*) "AGA", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "AGG", 3 ) == 0)      return 'r';
    else if (strncmp( codon, (const char*) "GTA", 3 ) == 0)      return 'v';
    else if (strncmp( codon, (const char*) "GTC", 3 ) == 0)      return 'v';
    else if (strncmp( codon, (const char*) "GTG", 3 ) == 0)      return 'v';
    else if (strncmp( codon, (const char*) "GTT", 3 ) == 0)      return 'v';
    else if (strncmp( codon, (const char*) "GCA", 3 ) == 0)      return 'a';
    else if (strncmp( codon, (const char*) "GCC", 3 ) == 0)      return 'a';
    else if (strncmp( codon, (const char*) "GCG", 3 ) == 0)      return 'a';
    else if (strncmp( codon, (const char*) "GCT", 3 ) == 0)      return 'a';
    else if (strncmp( codon, (const char*) "GAC", 3 ) == 0)      return 'd';
    else if (strncmp( codon, (const char*) "GAT", 3 ) == 0)      return 'd';
    else if (strncmp( codon, (const char*) "GAA", 3 ) == 0)      return 'e';
    else if (strncmp( codon, (const char*) "GAG", 3 ) == 0)      return 'e';
    else if (strncmp( codon, (const char*) "GGA", 3 ) == 0)      return 'g';
    else if (strncmp( codon, (const char*) "GGC", 3 ) == 0)      return 'g';
    else if (strncmp( codon, (const char*) "GGG", 3 ) == 0)      return 'g';
    else if (strncmp( codon, (const char*) "GGT", 3 ) == 0)      return 'g';
    else {
        printf( "Unrecognized codon\n" );
        exit(-1);
    }
}

void translate_protein( char dna[], char prot[], int size ) {

    int k;
    char codon[4] = {'\0'};
    for( k=0; k < size; k++ ) {
        prot[k] = '\0';
    }

    for( k=0; k < strlen(dna)-2; k+=3 ) {
        codon[0] = dna[k];
        codon[1] = dna[k+1];
        codon[2] = dna[k+2];
        prot[k/3] = translate(codon);
    }
}

float dna_energy( char* dna ) {

    energy_calls++;

    char mutatedprot[PROTEINSIZE*3] = {'\0'};
    atom_t atoms[PROTEINSIZE*10];
    int atomsize, aasize, diff;
    float energyoutput=0.0;
    int count;

    int rank;
    char filename[50] = {'\0'};
    char hostname[25] = {'\0'}; 
    char logfilename[50] = {'\0'};
    
    char cmd[512] = {'\0'};

    //printf("going to calculate energy.. ok?\n");
    translate_protein( dna, mutatedprot, sizeof(mutatedprot) );
    if( strrchr( mutatedprot, '_' ) == NULL ) {
        //printf("mutatedprot = %s\n", mutatedprot);

 	gethostname(hostname, 512);
	pid_t p = getpid();
	
        sprintf( filename, "%s_%d.out\0", hostname, p );
        
	char* pdbfilename = "temp.pdb"; 
        sprintf( logfilename, "%s.scwrl4out\0", filename );
        sprintf( cmd, "./Scwrl4 -p ./Scwrl4.ini -i %s -o %s -s %s -h -v  > %s\0", options->pdbfile, pdbfilename, filename, logfilename );

        write_dnafile( (const char*) filename, mutatedprot );
        system( cmd );
        read_pdbfile( (const char*) pdbfilename, atoms, atomsize );
        
        energyoutput = energy( atoms, atomsize, aasize );
        //printf( "energy = %f, prot = %s, hash = %u\n", energyoutput, mutatedprot, sdbm((unsigned char*) dna) );
    } else {
        energyoutput = 0.0;
    }

    return energyoutput;
}

hashtable_t* mostcommondna( hashtable_t** nhead ) {

    // find hash with more duplicates and return it.

    hashtable_t* h;
    hashtable_t* max;

    if( *nhead != NULL ) {
        max = *nhead;
        for ( h=*nhead; h != NULL; h=h->next ) {
            if( h->duplicates >= max->duplicates ) {
                max = h;
            }
        }
    }

    return max;
}

/**
 *  sdbm algorithm copied and modified from www.cse.yorku.ca/~oz/hash.pdf on April 14, 2011
 *
 *  this algorithm was created for sdbm (a public-domain reimplementation of ndbm)
 *  database library. it was found to do well in scrambling bits, causing better
 *  distribution of the keys and fewer splits. it also happens to be a good general
 *  hashing function with good distribution. the actual function is:
 *
 *      hash(i) = hash(i - 1) * 65599 + str[i];
 *
 *  what is included below is the faster version used in gawk. [there is even a faster, duff-device version]
 *  the magic constant 65599 was picked out of thin air while experimenting with different constants, and
 *  turns out to be a prime. this is one of the algorithms used in berkeley db (see sleepycat) and elsewhere.
 * */
unsigned long sdbm( unsigned char str[] ) {

    unsigned long hash = 0;
    int c=0;

    while( str[c] != '\0' ) {
        //printf( "c = %d, hash = %d, char = %c\n", c, hash, str[c] );
        hash = str[c++] + (hash << 6) + (hash << 16) - hash;
    }
    //printf( "hash = %u ", hash );

    return hash;
}

hashtable_t* inithash( char dna[], unsigned long int hash, bool calcenergy ) {

    hashtable_t *n = (hashtable_t*) malloc(sizeof(hashtable_t));

    if( n != NULL ) {
        strcpy( n->dna, dna );
        n->hash = hash;
        n->duplicates = 1;
        n->energy = 0.0;
        n->next = NULL;
    } else {
        printf( "error allocating memory..\n");
    }
    //print_hash( n );
    return n;
}

void add( hashtable_t* h, hashtable_t** nhead, hashtable_t** ntail ) {

    if( *nhead == NULL ) {
        *nhead = h;
        *ntail = *nhead;
    } else {
        (*ntail)->next = h;
        h->next = NULL;
        *ntail = h;
    }
}

hashtable_t* search( unsigned long int hash, hashtable_t** nhead, bool* found ) {

    hashtable_t *h;
    hashtable_t *prev=*nhead;

    *found = 0;

    //printf( "hash = %lu\n", hash );
    if( *nhead != NULL ) {
        prev=*nhead;
        for( h=*nhead; h != NULL; h=h->next ) {
            if( h->hash == hash ) {
                *found = 1;
                return h;
            } else if( prev->hash < hash && h->hash > hash ) {
                return prev;
            }
            prev=h;
        }
    }

    return NULL;
}

hashtable_t* update( char dna[], hashtable_t** nhead, hashtable_t** ntail, bool calcenergy ) {

    hashtable_t* h;
    hashtable_t* n;
    unsigned long int hash;
    bool found = 0;

    char mutatedprot[PROTEINSIZE*3] = {'\0'};
    translate_protein( dna, mutatedprot, sizeof(mutatedprot) );

    // compute current dna hash.
    hash = sdbm( (unsigned char*)mutatedprot );
    //hash = sdbm( (unsigned char*)dna );
    //printf( "hash = %u\n", hash );

    h = search( hash, nhead, &found );
    if( h != NULL && found ) {
        h->duplicates++;
    } else if ( h != NULL && !found) {
        n = inithash( dna, hash, calcenergy );
        n->next = h->next;
        h->next = n;
        if (calcenergy) {
            n->energy = dna_energy( n->dna );
        }
        h=n;
    } else {
        h = inithash( dna, hash, calcenergy );
        add( h, nhead, ntail );
        if (calcenergy) {
            h->energy = dna_energy( h->dna );
        } else {
            h->energy = 0.0;
        }
    }

    return h;
}

hashtable_t* updateByKey( unsigned long int hash, hashtable_t** nhead, hashtable_t** ntail, bool calcenergy ) {

    hashtable_t* h;
    hashtable_t* n;
    bool found = 0;
    char dna[1] = {'\0'};

    h = search( hash, nhead, &found );
    if( h != NULL && found ) {
        h->duplicates++;
    } else if ( h != NULL && !found) {
        n = inithash( dna, hash, calcenergy );
        n->next = h->next;
        h->next = n;
        if (calcenergy) {
            n->energy = dna_energy( n->dna );
        }
        h=n;
    } else {
        h = inithash( dna, hash, calcenergy );
        add( h, nhead, ntail );
        if (calcenergy) {
            h->energy = dna_energy( h->dna );
        }
    }

    return h;
}

void print_hash( hashtable_t* h ) {
    if( h != NULL ) {
        if( strlen(h->dna) != 0 )
            printf( "dna = %s, energy = %f, hash = %lu, duplicates = %d\n", h->dna, h->energy, h->hash, h->duplicates );
        else
            printf( "energy = %f, hash = %lu, duplicates = %d\n", h->energy, h->hash, h->duplicates );
    }
}

void print_hashtable( hashtable_t** nhead, bool printh ) {

    hashtable_t* h;

    int hashtablesize = 0;
    int total = 0;

    printf( "\nHash Table:\n" );
    for( h=*nhead; h != NULL; h=h->next ) {
        if( printh ) {
            print_hash( h );
        }
        total += h->duplicates + 1;
        hashtablesize++;
    }
    printf( "hashsize = %d, total = %d\n", hashtablesize, total );

    return;
}

int get_hashtablesize( hashtable_t** nhead ) {

    hashtable_t* h;

    int hashtablesize = 0;

    //printf( "Hash Table: " );
    for( h=*nhead; h != NULL; h=h->next ) {
        //print_hash( h );
        hashtablesize++;
    }
    //printf( "hashsize = %d\n", hashtablesize );

    return hashtablesize;
}

void clearhash( hashtable_t** nhead, hashtable_t** ntail ) {

    hashtable_t* curr;
    while( *nhead != NULL ) {
        curr = *nhead;
        *nhead = (*nhead)->next;
        free(curr);
    }
    *nhead = NULL;
    *ntail = NULL;
    return;
}

int get_robustnesstable( hashtable_t* h, char permutations[][PROTEINSIZE*3] ) {

    int k;
    int length;
    char aa;
    char* dna;
    int psize=0;
    unsigned long int hash;

    if( h != NULL && h->dna != NULL ) {
        dna = h->dna;

        length = strlen(dna);
        for( k=0; k < length; k++ ) {
            aa = dna[k];
            switch( dna[k] ) {
                case 'A':
                    dna[k] = 'C';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'T';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'G';
                    strcpy( permutations[psize++], dna );
                    break;
                case 'T':
                    dna[k] = 'C';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'A';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'G';
                    strcpy( permutations[psize++], dna );
                    break;
                case 'G':
                    dna[k] = 'C';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'T';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'A';
                    strcpy( permutations[psize++], dna );
                    break;
                case 'C':
                    dna[k] = 'A';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'T';
                    strcpy( permutations[psize++], dna );
                    dna[k] = 'G';
                    strcpy( permutations[psize++], dna );
                    break;
            }
            dna[k] = aa;
        }
    }

    return psize;
}

struct RB_Node {
    unsigned long hash;
    float energy;
    int duplicates;
};

typedef struct RB_Node rb_node_t;

void HashDest(void* a) {
    free((rb_node_t*)a);
}

int HashComp(const void* a,const void* b) {
    rb_node_t* h1 = (rb_node_t*)a;
    rb_node_t* h2 = (rb_node_t*)b;
    if( h1->hash > h2->hash) return(1);
    if( h1->hash < h2->hash) return(-1);
    return(0);
}

void HashPrint( const void* a ) {
    rb_node_t* h = (rb_node_t*)a;
    if( h != NULL ) {
            printf( "[[energy = %f, hash = %lu, duplicates = %d]]\n", h->energy, h->hash, h->duplicates );
    }
}

void InfoPrint(void* a) {
    ;
}

void InfoDest(void *a){
    ;
}

rb_node_t* initnode( unsigned long int hash ) {

    rb_node_t *n = (rb_node_t*) malloc(sizeof(rb_node_t));

    if( n != NULL ) {
        n->hash = hash;
        n->duplicates = 0;
        n->energy = 0.0;
    } else {
        printf( "error allocating memory..\n");
    }
    return n;
}

void usage(char* program) {
    printf("Usage: %s [options]\n\n", program);
    printf("Options:\n");
    printf("  -g\tnumber of generations (default: 500)\n");
    printf("  -p\tnumber of populations (default: 1000)\n");
    printf("  -aa\tinital aminoacid file (default: aa)\n");
    printf("  -pdb\tpdb file (default: 1EGL.pdb)\n");
    printf("  -h\tprint help dialog\n");
    printf("\n\n");
}

opt_t* cmd_options(int argc, char** argv) {
    int i = 1;
    opt_t *options = (opt_t*) malloc(sizeof(opt_t));

    options->gen = 500;
    options->pop = 1000;
    strcpy(options->aafile, "aa\0");
    strcpy(options->pdbfile, "1EGL.pdb\0");
    options->robustness = FALSE;

    while( i < argc ) {
        if( strcmp( argv[i], (const char*)"-r" ) == 0 ) {
            options->robustness = TRUE;
            strcpy( options->fdna, argv[++i] );
        } else if( strcmp( argv[i], (const char*)"-g" ) == 0 ) {
            options->gen = atoi(argv[++i]);
        } else if( strcmp( argv[i], (const char*)"-p" ) == 0 ) {
            options->pop = atoi(argv[++i]);
        } else if( strcmp( argv[i], (const char*)"-aa" ) == 0 ) {
            strcpy( options->aafile, argv[++i] );
        } else if( strcmp( argv[i], (const char*)"-pdb" ) == 0 ) {
            strcpy( options->pdbfile, argv[++i] );
        } else {
            usage(argv[0]);
            exit(1);
        }
        i++;
    }

    return options;
}

