/**
 * Lipid Lattice @version 0.0.3
 * 
 * Simulates a simple Monte-Carlo simulations of a lattice system 
 * of lipids and CHOL.
 * 
 * Based on interactions derived by MD simulations performed by Dr. Fabian Keller.
 * The first version of a squared lattice system for pure lipid mixtures was implementend by Davit Hokabyan & Andreas Heuer.
 *
 * Example Input and Output Files can be found in "/Check" for
 * comparison after changing the code.
 * 
 * Changelog:
 * 
 * 0.0.4:   * Changed the algorithm in "mc_move_order" (see there)
 * 
 * 0.0.3:   * Changed "Phi_P" (see there)
 * 
 * 0.0.2:   * Changed output from ".dat" to ".csv"
 *          * Added example files in "/Check" 
 * 
 * 0.0.1:   * Initial release.
 * /

/* Header Files */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <sys/stat.h>

/* Defintions */

#define N 10001                                         /** max. number of lattice sites */
#define M 21                                            /** max. number of lipid types */
#define SWAP(a,b) {double temp=(a);(a)=(b);(b)=temp;}   /** swaps two values a and b (duh)*/

/* Functions */

/* Main functions */

void main();                                            

void read_input();                                      

void make_neighbors();                                  
void make_particle_numbers();                           
void make_initialization();                             
void make_interaction();                                

/* Monte-Carlo functions */

void calc_all();                                        
void calc_evolution();                                  
void lattice_init();                                    

void mc_move_chol(int);                                 
void mc_move_order(int);                                
void mc_move_particle(int);                             

/* Calculating output variables */

void count_neighbors();                                 
int count_cc(int);                                      
int count_cc_moore(int);                                
int count_dc(int);                                      
int count2_dd(int,int);                                 
int count2_cd(int,int);                                 
int count2_cc(int,int);                                 
void count_order();                                     

/* Printing output files */

void plot_order(int,char[]);                            
void plot_order_corr(int,char[]);                       
void plot_order_neigh(int,char[]);                      
void plot_lattice(int,char[]);                          
void plot_energies(int,char[]);                         
void plot_info(int,char[]);                                 

/* Calculating energies within the lattice */

double ener_total();                                    
double ener_dd(int,int);                                
double ener_cd(int,int);                                
double ener_cc(int,int);                                
double ener_self(int);                                  

/* Interaction functions as determined by MD simulations */

double NN_CD_L(int);                                    
double NN_CD_P(int);                                    
double NN_DD_L(int);                                    
double NN_DD_P(int);                                    

double H_self_L(int);                                   
double H_self_P(int);                                   

double E_L(int);                                        
double E_P(int);                                        
double E_init_L(int);                                   
double E_init_P(int);                                   

double Phi_L(int,int);                                  
double Phi_P(int,int);                                  
double Phi_init_L(int,int);                             
double Phi_init_P(int,int);                             

double H_DD_LL_avg(int,int);                            
double H_DD_PP_avg(int,int);                            
double H_DD_LP_avg(int,int);                            

double H_DD_LL_min(int,int);                            
double H_DD_PP_min(int,int);                            
double H_DD_LP_min(int,int);                            

double H_CC_L(int);                                     
double H_CC_P(int);                                     

double H_CD_L(int,int);                                 
double H_CD_P(int,int);                                 

void MD_order();                                        

/* Technical details */

void transfer1(int,int*);                               
int transfer2(int*);                                    

int per(int,int);                                       
void sort(int,double *,int *);                           
void fisher_yate_shuffle(int*);                         

int min_order(int,int);                                 
int avg_order(int,int);                                 
void norm_order();                                      

double norm_dist(double,double,double);                  

void print_progress(double);                            


/* Variables */

/* Setting up the simulations */

char project_name[1024];                                /** name used for output directories (s. input.txt) */ 

int lat_type;                                           /** defines the simulation's mode (0: MC-simulation, 1: IBI, 2: CHOL's entropy) */
int int_type;                                           /** defines the interaction type used (0: avg, 1: minimal order parameter) */ 
int mv_type;                                            /** defines the jump type used (0: random site, 1: to the nearest site) */

int dim;                                                /** dimension of the lattice (should be 2 in most cases) */
int lattice_size;                                       /** defines the number of sites in every dimension */

double read_temperature;                                /** temperature in K as described in input */
double temperature;                                     /** temperature in kJ as used in the simulation */
int nr_types;                                           /** number of different lipids */
int lipid_type;                                         /** identity of the lipids  */
double c_lipid[M];                                      // concentration of lipid of type M
double c_chol;                                          // concentration of CHOL
int orderwidth;                                         // defines the change of the order parameter

int nr_runs;                                            // number of runs
int nr_steps;                                           // number of MC steps
int nr_count;                                           // number of steps counted for the calculation of output parameters
int nr_data;                                            // TO-DO!
double pequil;                                          // percentage after which the output parameters are calculated

int reverse[4];                                         // array to reverse the directions in the lattice
int length[4];                                          // array for the lattice length in every dimension
double quotf;                                           // quotient for calculations

int nr_lattice;                                         // number of lattice sites
int lattice_sites[N+1];                                 // array of all lattice sites (for Fisher-Yates-Shuffle)
int nr_chol;                                            // number of CHOL particles
int nrt[M];                                             // number of lipid particles

/* Describing the lattice */

int particlet[N];                                       // denotes the particletyp of lattice site N
int order[N];                                           // denotes the order of lipid on site N
int cholpos[N];                                         // denotes if a CHOL is on lattice site N

/* Output parameters */

int neightype[M][M];                                    // contacts between PL and PL
int neightype_dc[M];                                    // contacts between PL and CHOL
int neightype_cc;                                       // contacts between CHOL and CHOL
int neightype_cc_moore;                                 // contacts between CHOL and CHOL (Moore neighborhood)

/* Neighbor lists */

int neigh[N][7], neigh_cd[N][7], neigh_dc[N][7];                    // neighboring sites of lattice site N
int neigh2_dd[N][4][10],neigh2_cd[N][4][10],neigh2_cc[N][4][10];    // neighboring sites in all directions for pairs with lattice site N
int neigh_cd_2[N][12], neigh_cd_3[N][20];                           // second and third "coordination sphere" (Von Neumann)
int neigh_moore[N][3];                                              // neighboring CHOL sites of CHOL lattice site N (Moore)
int cmax_dd,cmax_cd,cmax_cc;                                        // max number of CHOL neighbors for specific pairs

/* Output parameters */

double p_MC[M][151];                                    // order parameter distribution
double p_MC_neigh[M][151][5];                           // order parameter distribution in dependence on the CHOL neighbors
double p_MC_1[M][151];                                  // order parameter distribution in the second coordination sphere around CHOL
double p_MC_2[M][151];                                  // order parameter distribution in the third coordination sphere around CHOL
double p_MC_3[M][151];                                  // order parameter distribution in the fourth coordination sphere around CHOL

int cc_contacts;                                        // CHOL-CHOL contacts (von Neumann)
int cc_moore;                                           // CHOL-CHOL contacts (Moore)
double neighpair_cd[10];                                // number of CHOL around cd pairs

double ener_mem;                                        // stores the total energy
double contacts_ener[5];                                // energy of CHOL-CHOL contacts depending on CHOL neighbors
double contacts_num[5];                                 // number of CHOL-CHOL contacts depending on CHOL neighbors

double H_dd_avg_track[2][301][10];                      // tracks the energy between two lipids based on their average order parameter
double H_dd_min_track[2][151][10];                      // tracks the energy between two lipids based on their minimal order parameter

/* Interaction functions */

double nneigh_cd[M][5], nneigh_dd[M][5];                // nearest neighbor functions of lipids around CHOL/lipids

double H_self[M][151];                                  // self energy of lipids
double H_dd_min[M][M][151][10];                         // interaction energy between PL and PL based on the minimal order parameter
double H_dd_avg[M][M][151][10];                         // interaction energy between PL and PL based on the average order parameter
double H_cd[M][151][10];                                // interaction energy between PL and CHOL
double H_cc[10];                                        // interaction energy between CHOL and CHOL

double Entro[M][151];                                   // entropic function determined in the lattice model (Iterative Boltzmann Inversion)
double Phi[M][151][5];                                  // entropic function for the effect of CHOL on the lipid's chain entropy

double p_MD[151];                                       // order parameter distribution determined in the MD simulations


/**
 * Main function performing all the simulations.
 * 
 * Start reading here. 
 */
void main()
{
    char ufname[1024];                                  // placeholder for several filenames
    int seed = 128;                                     // random seed

    /* Printing Developer Information*/

    printf("\n------------------------------------------------------\n");
    printf("         _   _\n");
    printf("        (.)_(.)\n");
    printf("     _ (   _   ) _\n");
    printf("    / \\/`-----'\\/ \\   Welcome to Lipid Lattice 0.0.4\n"); 
    printf("  __\\ ( (     ) ) /__\n");
    printf("  )   /\\ \\._./ /\\   (\n");
    printf("   )_/ /|\\   /|\\ \\_(\n");
    printf("\n------------------------------------------------------\n\n");

    /* Reading the Input-File */
    read_input();                                       // reads the input file named "input.txt"
    printf("Input read successfully...\n");             // confirms the correct reading of the input



    /* Setting up Variables */
    MD_order();
    reverse[0] = 1, reverse[1] = 0, reverse[2] = 3, reverse[3] = 2;
    if (nr_types == 1) c_lipid[1] = 1.0;                  // if there is only lipid its concentration is 1

    if (c_chol > 0.5) printf("ERROR: CHOL-Conc. is too high (0.5 is the possible Maximum).\n");

    c_chol = (c_chol)/(1-c_chol);                       // calculating the concentration of CHOL in its sub-lattice
    temperature = read_temperature/120.0;               // converting temperature from K to kJ
    length[1] = length[2] = lattice_size;               // defining the length of the lattice in every dimension
    nr_data = 1;                                        // TO-DO!

    if (lipid_type != 1 && nr_types == 2) printf("ERROR: Too many types of lipids without DLiPC.\n");
    if (nr_steps%(nr_data*nr_count) != 0) printf("ERROR: Wrong number of steps given.\n");
    if (nr_count > nr_steps/nr_data) printf("ERROR: Data Points too high!\n");
    if (lat_type != 0 && nr_types == 2) printf("ERROR: Choosen simulation mode only for single lipid type.\n");

    if (lat_type == 1 || lat_type == 2) printf("ERROR: The choosen simulation is not implementend! Please refer to older scripts.\n");

    if (nr_types == 2) printf("WARNING: Simulations with two lipids is not yet properly implementend!\n");

    srandom(seed);                                      // random with the same seed
    //srandom(clock());                                   // random with the clock to have different starting structures
    quotf = pow(2.0,-31.0);                             // quotient used for several calculations

    printf("Variables defined...\n");

    /* Setting up the Lattice */

    make_neighbors();
    make_particle_numbers();
    make_initialization();
    make_interaction();

    char run_folder[1024];
    sprintf(run_folder, "%s", project_name);
    struct stat st = {0};
    if (stat(run_folder,&st) == -1){
        mkdir(run_folder,0700);
    }

    printf("System initialized...\n \n");

    calc_all();

}

/**
 * Function gathering all functions performing the Monte-Carlo-Simulation. 
 */
void calc_all()
{
    calc_evolution();
}

/**
 * Function performing the Monte-Carlo-Simulation.
 */
void calc_evolution()
{
    FILE *fconfig;
    int l,n,t;
    int m;
    char fname[1024], str[12];
    char ufname[1024];
    double count = 0.1;

    for (n=0;n<=nr_runs;n++){

        lattice_init();
        count_neighbors();

        if (lat_type == 0) {
            char run_folder[1024];
            sprintf(run_folder, "%s/%d", project_name,n);
            struct stat st = {0};
            if (stat(run_folder,&st) == -1){
                mkdir(run_folder,0700);
            }
        }
        sprintf(fname,"%s/%d/contacts.dat",project_name,n);
        fconfig = fopen(fname,"w");

        if (nr_types == 1) fprintf(fconfig, "t dd cc cd\n");
        if (nr_types == 2) fprintf(fconfig, "t 11 22 12 cc c1 c2\n");

        for (t=1;t<=nr_steps;t++){
            
            if (t%(nr_steps/(nr_data*nr_count)) == 0){

                count_neighbors();

                if (nr_types == 1){
                    fprintf(fconfig, "%d %d %d %d\n",
                    t,neightype[1][1],neightype_cc,neightype_dc[1]);
                }

                if (nr_types == 2){
                    fprintf(fconfig, "%d %d %d %d %d %d %d\n",
                    t, neightype[1][1], neightype[2][2], neightype[2][1],
                    neightype_cc, neightype_dc[1], neightype_dc[2]);
                }
            }

            if ((double)t > pequil*(double)nr_steps){

                count_order();
                count_neighbors();
                cc_contacts += neightype_cc;
                cc_moore += neightype_cc_moore;
                ener_mem += ener_total();

                for (l=0;l<nr_lattice;l++){
                    if (cholpos[l] == 1){
                        int num = count_cc(l);
                        contacts_num[num] ++;
                        for (m=0;m<4;m++){
                            int cc = count2_cd(l,m);
                            neighpair_cd[cc] ++;
                            contacts_ener[num] += ener_cd(l,m);
                        }
                    }
                }
            }

            fisher_yate_shuffle(lattice_sites);

            for (l=0;l<nr_lattice;l++){
                int lm = lattice_sites[l];
                if (cholpos[lm] == 1) mc_move_chol(lm);
                mc_move_order(lm);
                mc_move_particle(lm);
            }

            for (l=0;l<nr_lattice;l++){
                for (m=0;m<2*dim;m++){

                    int x = neigh[l][m];
                    int num_neighbor = count2_dd(l,m);

                    int S_min = min_order(l,x);
                    int S_avg = order[l] + order[x];

                    double ener_min = H_dd_min[1][1][S_min][num_neighbor];
                    double ener_avg = H_dd_avg[1][1][S_avg][num_neighbor];

                    H_dd_avg_track[0][S_avg][num_neighbor] ++;
                    H_dd_avg_track[1][S_avg][num_neighbor] += ener_min;

                    H_dd_min_track[0][S_min][num_neighbor] ++;
                    H_dd_min_track[1][S_min][num_neighbor] += ener_avg;
                }
            }

            if ((double)t == round(count*(double)nr_steps)){
                print_progress(count);
                count += 0.1;
            }

        }

        if (lat_type == 0){

            plot_order(n,ufname);
            plot_order_neigh(n,ufname);
            plot_order_corr(n,ufname);
            plot_lattice(n,ufname);
            plot_energies(n,ufname);

            plot_info(n,ufname);
        }

        fclose(fconfig);

    }

    printf("------------------------------------------------------\n \n");
}

/**
 * Moves the CHOL on a given...
 * @param l lattice site. 
 */
void mc_move_chol(int l)
{
    int o,k;
    int m;
    int x,y,g,h;
    double delta,jc_old,jc_new;
    int l_new;

    int inter_cc[nr_lattice][4],inter_cd[nr_lattice][4],inter_dd[nr_lattice][4];
    int inter_self[nr_lattice];

    if (mv_type == 0){
        l_new = random()%(nr_lattice);
        if (l_new = l) l_new = (l_new + 1)%nr_lattice;
    }
    if (mv_type == 1){
        m = random()%(2*dim);
        l_new = neigh[l][m];
    }

    if (cholpos[l_new] == 0){

        for (o=0;o<nr_lattice;o++){
            inter_self[o] = 0;
            for (m=0;m<4;m++){
                inter_cc[o][m] = 0;
                inter_cd[o][m] = 0;
                inter_dd[o][m] = 0;
            }
        }
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            inter_self[g] ++;   // self_interaction
            inter_self[h] ++;   // self_interaction

            inter_cd[l][m] ++;

            if (cholpos[neigh_dc[g][m]] == 1){
                inter_cd[neigh_dc[g][m]][reverse[m]] ++;
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                inter_cd[neigh_dc[h][m]][reverse[m]] ++;
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){ 
                    inter_cd[x][k] ++;
                    if (cholpos[neigh[x][k]] == 1){
                        inter_cc[x][k] ++;
                        inter_cc[neigh[x][k]][reverse[k]] ++;
                    }
                }

                if (cholpos[y] == 1){
                    inter_cd[y][k] ++;
                    if (cholpos[neigh[y][k]] == 1){
                        inter_cc[y][k] ++;
                        inter_cc[neigh[y][k]][reverse[k]] ++;
                    }
                }

                if (m<2){
                    inter_dd[g][k] ++;
                    inter_dd[neigh[g][k]][reverse[k]] ++;

                    inter_dd[h][k] ++;
                    inter_dd[neigh[h][k]][reverse[k]] ++;
                }
            }

            if (m == 2){
                inter_dd[g][2] ++;
                inter_dd[neigh[g][2]][reverse[2]] ++;
                inter_dd[g][1] ++;
                inter_dd[neigh[g][1]][reverse[1]] ++;

                inter_dd[h][2] ++;
                inter_dd[neigh[h][2]][reverse[2]] ++;
                inter_dd[h][1] ++;
                inter_dd[neigh[h][1]][reverse[1]] ++;
            }

            if (m == 3){
                inter_dd[g][3] ++;
                inter_dd[neigh[g][3]][reverse[3]] ++;
                inter_dd[g][0] ++;
                inter_dd[neigh[g][0]][reverse[0]] ++;

                inter_dd[h][3] ++;
                inter_dd[neigh[h][3]][reverse[3]] ++;
                inter_dd[h][0] ++;
                inter_dd[neigh[h][0]][reverse[0]] ++;
            }
        }

        jc_old = 0;
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            jc_old += ener_self(g);
            jc_old += ener_self(h);

            jc_old += ener_cd(l,m)/(double)(inter_cd[l][m]);

            if (cholpos[neigh_dc[g][m]] == 1){
                jc_old += ener_cd(neigh_dc[g][m],reverse[m])/(double)(inter_cd[neigh_dc[g][m]][reverse[m]]);
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                jc_old += ener_cd(neigh_dc[h][m],reverse[m])/(double)(inter_cd[neigh_dc[h][m]][reverse[m]]);
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){
                    jc_old += ener_cd(x,k)/(double)(inter_cd[x][k]);
                    if (cholpos[neigh[x][k]] == 1) jc_old += ener_cc(x,k)/(double)(inter_cc[x][k]);
                }

                if (cholpos[y] == 1){
                    jc_old += ener_cd(y,k)/(double)(inter_cd[y][k]);
                    if (cholpos[neigh[y][k]] == 1) jc_old += ener_cc(y,k)/(double)(inter_cc[y][k]);
                }

                if (m<2){
                    jc_old += ener_dd(g,k)/(double)(inter_dd[g][k]);
                    jc_old += ener_dd(h,k)/(double)(inter_dd[h][k]);
                }
            }

            if (m == 2){
                jc_old += ener_dd(g,2)/(double)(inter_dd[g][2]);
                jc_old += ener_dd(g,1)/(double)(inter_dd[g][1]);

                jc_old += ener_dd(h,2)/(double)(inter_dd[h][2]);
                jc_old += ener_dd(h,1)/(double)(inter_dd[h][1]);
            }

            if (m == 3){
                jc_old += ener_dd(g,3)/(double)(inter_dd[g][3]);
                jc_old += ener_dd(g,0)/(double)(inter_dd[g][0]);

                jc_old += ener_dd(h,3)/(double)(inter_dd[h][3]);
                jc_old += ener_dd(h,0)/(double)(inter_dd[h][0]);
            }
        }
        
        SWAP(cholpos[l],cholpos[l_new]);

        for (o=0;o<nr_lattice;o++){
            inter_self[o] = 0;
            for (m=0;m<4;m++){
                inter_cc[o][m] = 0;
                inter_cd[o][m] = 0;
                inter_dd[o][m] = 0;
            }
        }
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            inter_self[g] ++;   // ener_self
            inter_self[h] ++;   // ener_self

            inter_cd[l_new][m] ++;

            if (cholpos[neigh_dc[g][m]] == 1){
                inter_cd[neigh_dc[g][m]][reverse[m]] ++;
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                inter_cd[neigh_dc[h][m]][reverse[m]] ++;
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){ 
                    inter_cd[x][k] ++;
                    if (cholpos[neigh[x][k]] == 1){
                        inter_cc[x][k] ++;
                        inter_cc[neigh[x][k]][reverse[k]] ++;
                    }
                }

                if (cholpos[y] == 1){
                    inter_cd[y][k] ++;
                    if (cholpos[neigh[y][k]] == 1){
                        inter_cc[y][k] ++;
                        inter_cc[neigh[y][k]][reverse[k]] ++;
                    }
                }

                if (m<2){
                    inter_dd[g][k] ++;
                    inter_dd[neigh[g][k]][reverse[k]] ++;

                    inter_dd[h][k] ++;
                    inter_dd[neigh[h][k]][reverse[k]] ++;
                }
            }

            if (m == 2){
                inter_dd[g][2] ++;
                inter_dd[neigh[g][2]][reverse[2]] ++;
                inter_dd[g][1] ++;
                inter_dd[neigh[g][1]][reverse[1]] ++;

                inter_dd[h][2] ++;
                inter_dd[neigh[h][2]][reverse[2]] ++;
                inter_dd[h][1] ++;
                inter_dd[neigh[h][1]][reverse[1]] ++;
            }

            if (m == 3){
                inter_dd[g][3] ++;
                inter_dd[neigh[g][3]][reverse[3]] ++;
                inter_dd[g][0] ++;
                inter_dd[neigh[g][0]][reverse[0]] ++;

                inter_dd[h][3] ++;
                inter_dd[neigh[h][3]][reverse[3]] ++;
                inter_dd[h][0] ++;
                inter_dd[neigh[h][0]][reverse[0]] ++;
            }
        }

        jc_new = 0;
        for (m=0;m<2*dim;m++){

            x = neigh[l][m];
            y = neigh[l_new][m];
            g = neigh_cd[l][m];
            h = neigh_cd[l_new][m];

            jc_new += ener_self(h);
            jc_new += ener_self(g);

            jc_new += ener_cd(l_new,m)/(double)(inter_cd[l_new][m]);

            if (cholpos[neigh_dc[g][m]] == 1){
                jc_new += ener_cd(neigh_dc[g][m],reverse[m])/(double)(inter_cd[neigh_dc[g][m]][reverse[m]]);
            }
            if (cholpos[neigh_dc[h][m]] == 1){
                jc_new += ener_cd(neigh_dc[h][m],reverse[m])/(double)(inter_cd[neigh_dc[h][m]][reverse[m]]);
            }

            for (k=0;k<2*dim;k++){

                if (cholpos[x] == 1){
                    jc_new += ener_cd(x,k)/(double)(inter_cd[x][k]);
                    if (cholpos[neigh[x][k]] == 1) jc_new += ener_cc(x,k)/(double)(inter_cc[x][k]);
                }

                if (cholpos[y] == 1){
                    jc_new += ener_cd(y,k)/(double)(inter_cd[y][k]);
                    if (cholpos[neigh[y][k]] == 1) jc_new += ener_cc(y,k)/(double)(inter_cc[y][k]);
                }

                if (m<2){
                    jc_new += ener_dd(g,k)/(double)(inter_dd[g][k]);
                    jc_new += ener_dd(h,k)/(double)(inter_dd[h][k]);
                }
            }

            if (m == 2){
                jc_new += ener_dd(g,2)/(double)(inter_dd[g][2]);
                jc_new += ener_dd(g,1)/(double)(inter_dd[g][1]);

                jc_new += ener_dd(h,2)/(double)(inter_dd[h][2]);
                jc_new += ener_dd(h,1)/(double)(inter_dd[h][1]);
            }

            if (m == 3){
                jc_new += ener_dd(g,3)/(double)(inter_dd[g][3]);
                jc_new += ener_dd(g,0)/(double)(inter_dd[g][0]);

                jc_new += ener_dd(h,3)/(double)(inter_dd[h][3]);
                jc_new += ener_dd(h,0)/(double)(inter_dd[h][0]);
            }
        }

        delta = jc_new - jc_old;

        if (exp(-delta/temperature) < random()*quotf){
            SWAP(cholpos[l],cholpos[l_new]);
        }
    }
}

/**
 * Moves the order parameter of a lipid on a given...
 *
 * @param l lattice site.
 * 
 * @version 0.0.4 Changed algorithm so that illegal moves are not reattempted.
 */
void mc_move_order(int l)
{
    double delta,jc_old,jc_new;
    int delta_order;
    int m;

    jc_old = 0;
    for (m=0;m<2*dim;m++){
        jc_old += ener_dd(l,m);
        if (cholpos[neigh_dc[l][m]] == 1){
            jc_old += ener_cd(neigh_dc[l][m],reverse[m]);
        }
    }
    jc_old += ener_self(l);

    /** In version 0.0.3 moves changing the order parameter to an illegal
     * value were not automatically rejected but just "rerolled".
     */

    /*while (0 == 0){
        delta_order = random()%orderwidth - (orderwidth-1)/2;
        if (lipid_type == 1 && particlet[l] == 1){ 
            if (order[l] + delta_order >= 30 && order[l] + delta_order <= 120){
                break;
            }
        }
        else{
            if (order[l] + delta_order >= 30 && order[l] + delta_order <= 150){
                break;
            }
        }
    }*/

    delta_order = random()%orderwidth - (orderwidth-1)/2;
    order[l] += delta_order;
    if (lipid_type == 1 && particlet[l] == 1){
        if (order[l] < 30 || order[l] > 120){
            order[l] -= delta_order;
            return;
        }
    }
    else{
        if (order[l] < 30 || order[l] > 150){
            order[l] -= delta_order;
            return;
        }
    }

    jc_new = 0;
    for (m=0;m<2*dim;m++){
        jc_new += ener_dd(l,m);
        if (cholpos[neigh_dc[l][m]] == 1){
            jc_new += ener_cd(neigh_dc[l][m],reverse[m]);
        }
    }
    jc_new += ener_self(l);

    delta = jc_new - jc_old;
    if (exp(-delta/temperature) < random()*quotf){
        order[l] -= delta_order;
    }
}

/**
 * Moves the particle (lipid) on a given...
 *
 * @param l lattice site. 
 */
void mc_move_particle(int l)
{
    double delta,jc_old,jc_new;
    int m,l_new;

    if (mv_type == 0){
        l_new = random()%(nr_lattice);
        if (l_new == l) l_new = (l_new + 1)%nr_lattice;
    }
    if (mv_type == 1){
        m = random()%(2*dim);
        l_new = neigh[l][m];
    }

    if (particlet[l] != particlet[l_new]){

        jc_old = 0;
        for (m=0;m<2*dim;m++){
            if (neigh[l][m] != l_new) jc_old += ener_dd(l,m);
            else jc_old += 0.5*ener_dd(l,m);

            if (neigh[l_new][m] != l) jc_old += ener_dd(l_new,m);
            else jc_old += 0.5*ener_dd(l_new,m);

            if (cholpos[neigh_dc[l_new][m]] == 1){
                jc_old += ener_cd(neigh_dc[l_new][m],reverse[m]);
            }
            if (cholpos[neigh_dc[l][m]] == 1){
                jc_old += ener_cd(neigh_dc[l][m],reverse[m]);
            }
        }

        SWAP(particlet[l],particlet[l_new]);
        SWAP(order[l],order[l_new]);

        jc_new = 0;
        for (m=0;m<2*dim;m++){
            if (neigh[l][m] != l_new) jc_new += ener_dd(l,m);
            else jc_new += 0.5*ener_dd(l,m);

            if (neigh[l_new][m] != l) jc_new += ener_dd(l_new,m);
            else jc_new += 0.5*ener_dd(l_new,m);

            if (cholpos[neigh_dc[l_new][m]] == 1){
                jc_new += ener_cd(neigh_dc[l_new][m],reverse[m]);
            }
            if (cholpos[neigh_dc[l][m]] == 1){
                jc_new += ener_cd(neigh_dc[l][m],reverse[m]);
            }
        }

        delta = jc_new - jc_old;

        if (exp(-delta/temperature) < random()*quotf){
            SWAP(particlet[l],particlet[l_new]);
            SWAP(order[l],order[l_new]);
        }
    }
}

/**
 * Counts the number of different contacts throughout
 * the lattice and saves them in the corresponding variables.
 */
void count_neighbors()
{
    int i,j;
    int l,x;
    int m;

    for (i=1;i<=nr_types;i++){
        neightype_dc[i] = 0;
        for (j=1;j<=i;j++){
            neightype[i][j] = 0;
        }
    }
    neightype_cc = 0;
    neightype_cc_moore = 0;

    for (l=0;l<nr_lattice;l++){
        for (m=0;m<2*dim;m+=2){
            x = neigh[l][m];
            if (particlet[l] >= particlet[x]){ 
                neightype[particlet[l]][particlet[x]] ++;
            }
            else{
                neightype[particlet[x]][particlet[l]] ++;
            }
        }
    }

    for (l=0;l<nr_lattice;l++){
        neightype_dc[particlet[l]] += count_dc(l);
    }

    for (l=0;l<nr_lattice;l++){
        if (cholpos[l] == 1) neightype_cc += count_cc(l);
        if (cholpos[l] == 1) neightype_cc_moore += count_cc_moore(l);
    }
    neightype_cc /= 2;
    neightype_cc_moore /= 2;
}

/**
 * Counts and returns the number of CHOL around a lipid on a
 * given lattice site, according to a Von Neumann neighborhood.
 * 
 * @param l Lattice Site
 */
int count_dc(int l)
{
    int m,n=0;
    for (m=0;m<=3;m++) n += cholpos[neigh_dc[l][m]];
    return(n);
}

/**
 * Counts and returns the number of CHOL-CHOL contacts on a
 * given lattice site, according to a Von Neumann neighborhood.
 * 
 * @param l Lattice Site
 */
int count_cc(int l)
{
    int m,n=0;
    for (m=0;m<=3;m++) n += cholpos[neigh[l][m]];
    return (n);
}

/**
 * Counts and returns the number of CHOL-CHOL contacts on a
 * given lattice site, according to a Moore neighborhood.
 * 
 * @param l Lattice Site
 */
int count_cc_moore(int l)
{
    int m,n=0;
    for (m=0;m<=3;m++) n+= cholpos[neigh_moore[l][m]];
    return (n);
}

/** 
 * Counts and returns the number of CHOL around a Lipid-Lipid pair,
 * starting on a given lattice site in a given direction. 
 * 
 * @param l Lattice Site
 * @param m Direction
 */
int count2_dd(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_dd;c++) n += cholpos[neigh2_dd[l][m][c]];
    return (n);
}

/** 
 * Counts and returns the number of CHOL around a CHOL-Lipid pair,
 * starting on a given lattice site in a given direction. 
 * 
 * @param l Lattice Site
 * @param m Direction
 */
int count2_cd(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_cd;c++) n += cholpos[neigh2_cd[l][m][c]];
    return (n);
}

/** 
 * Counts and returns the number of CHOL around a CHOL-CHOL pair,
 * starting on a given lattice site in a given direction. 
 * 
 * @param l Lattice Site
 * @param m Direction
 */
int count2_cc(int l,int m)
{
    int c,n=0;
    for (c=1;c<=cmax_cc;c++) n += cholpos[neigh2_cc[l][m][c]];
    return (n);
}

/** 
 * Calculates the order parameter distributions and saves them in
 * the corresponding arrays. 
 */
void count_order()
{
    int c,l;
    int d,e,b;
    int x;

    for (l=0;l<nr_lattice;l++){
        c = count_dc(l);
        p_MC[particlet[l]][order[l]] ++;
        p_MC_neigh[particlet[l]][order[l]][c] ++;
        if (cholpos[l] == 1){
            for (b = 0; b <= 3; b++){
                x = neigh_cd[l][c];
                p_MC_1[particlet[x]][order[x]] ++;
            }
            for (d = 0; d <= 11; d++){
                x = neigh_cd_2[l][d];
                p_MC_2[particlet[x]][order[x]] ++;
            }
            for (e=0; e <= 19; e++){
                x = neigh_cd_3[l][e];
                p_MC_3[particlet[x]][order[x]] ++;
            }
        }
    }
}

/** 
 * Calculates and returns all energies within the lattice. 
 */
double ener_total()
{
    int l,m;
    double enertot;

    enertot = 0;
    for (l=0;l<nr_lattice;l++){
        for (m=0;m<2*dim;m++){
            enertot += ener_dd(l,m);
            if (cholpos[neigh_dc[l][m]] == 1){
                enertot += ener_cd(neigh_dc[l][m],reverse[m]);
            }
            if (cholpos[l] == 1){
                enertot += ener_cd(l,m);
                if (cholpos[neigh[l][m]] == 1){
                    enertot += ener_cc(l,m);
                }
            }
        }
    }
    enertot /= 2.0;
    for (l=0;l<nr_lattice;l++){
        enertot += ener_self(l);
    }
    return(enertot);
}

/**
 * Calculates and returns of a Lipid-Lipid pair, 
 * starting on a given lattice site in a specific direction.
 * 
 * @param l Lattice Site
 * @param m Direction
 */
double ener_dd(int l, int m)
{
    int i = particlet[l];
    int x = neigh[l][m];
    int j = particlet[x];

    int c = count2_dd(l,m);
    int ci = count_dc(l);
    int cj = count_dc(x);

    if (int_type == 0){
        return ((1.0/8.0)*(nneigh_dd[i][ci] + nneigh_dd[j][cj])*
            H_dd_avg[i][j][avg_order(l,x)][c]);
    }

    if (int_type == 1){
        return ((1.0/8.0)*(nneigh_dd[i][ci] + nneigh_dd[j][cj])*
            H_dd_min[i][j][min_order(l,x)][c]);
    }
}

/**
 * Calculates and returns of a CHOL-Lipid pair, 
 * starting on a given lattice site in a specific direction.
 * 
 * @param l Lattice Site
 * @param m Direction
 */
double ener_cd(int l, int m)
{
    int x = neigh_cd[l][m];
    int i = particlet[x];

    int c = count2_cd(l,m);
    int cm = count_cc(l);

    return ((1/4.0)*(nneigh_cd[i][cm])*H_cd[i][order[x]][c]);
}

/**
 * Calculates and returns of a CHOL-CHOL pair, 
 * starting on a given lattice site in a specific direction.
 * 
 * @param l Lattice Site
 * @param m Direction
 */
double ener_cc(int l, int m)
{
    int c = count2_cc(l,m);
    int lm = neigh[l][m];
    
    int cl = count_cc(l);
    int cm = count_cc(lm);

    return(H_cc[c]); // Note: In some cases it's useful to use H_cc[c]*(c*0.25 + 1.00) instead
}

/**
 * Calculates and returns the self energy of the lipid on a given...
 * 
 * @param l lattice site.
 */
double ener_self(int l)
{
    int i = particlet[l];
    int c = count_dc(l);

    return (H_self[i][order[l]] - temperature*(Entro[i][order[l]]+Phi[i][order[l]][c]));
}

/**
 * Initializes the lattice structure.
 */
void lattice_init()
{
    int sum,i,j,k;
    double x[N];
    for (i=1 ; i <= nr_lattice ; i++)
    x[i] = random();

    i = 1;
    for (j=1 ; j <= nr_types ; j++)
    for (k=1 ; k <= nrt[j] ; k++)
    {
        particlet[i] = j;
        i++;
    }

    sort(nr_lattice,x,particlet);
    for (i=0 ; i < nr_lattice ; i++) particlet[i] = particlet[i+1];

    for (i=0 ; i < nr_lattice; i++){
        if (lipid_type == 1 && particlet[i] == 1){
            while (order[i] > 120 || order[i] < 30){ 
                order[i] = 50 + (int)(30*random()*quotf);
            }  
        }
        else{
            while (order[i] > 150 || order[i] < 30){
                order[i] = 30 + (int)(30*(random()*quotf)); 
            }
        } 
    }

    for (i=1 ; i <= nr_chol; i++){
        cholpos[i] = 1; x[i] = random(); 
    }
    for (i=nr_chol+1 ; i <= nr_lattice; i++){
        cholpos[i] = 0; x[i] = random(); 
    }
    sort(nr_lattice,x,cholpos);
    for (i=0 ; i < nr_lattice ; i++) cholpos[i] = cholpos[i+1];
}

/**
 * Normalizes all order parameter distributions.
 */
void norm_order()
{
    int s,c,i;
    double p_MC_norm;
    double p_MC_1_norm, p_MC_2_norm, p_MC_3_norm;

    for (i=1;i<=nr_types;i++){
        p_MC_norm = 0;
        for (s=0;s<=150;s++){
            p_MC_norm += p_MC[i][s];
            p_MC_1_norm += p_MC_1[i][s];
            p_MC_2_norm += p_MC_2[i][s];
            p_MC_3_norm += p_MC_3[i][s];
        }
        for (s=0;s<=150;s++){
            p_MC[i][s] /= p_MC_norm;
            p_MC_1[i][s] /= p_MC_1_norm;
            p_MC_2[i][s] /= p_MC_2_norm;
            p_MC_3[i][s] /= p_MC_3_norm;
            for (c=0;c<=4;c++){
                p_MC_neigh[i][s][c] /= p_MC_norm;
            }
        }
    }
}

/**
 * Returns the average of two order parameters.
 * (Simply the sum, to still be able to use int).
 * 
 * @param l First order parameter.
 * @param o Second order parameter.
 */
int avg_order(int l, int o)
{
    return (order[l]+order[o]);
}

/**
 * Returns the minimal of two order parameters.
 * 
 * @param l First order parameter.
 * @param o Second order parameter.
 */
int min_order(int l, int o)
{
    if (order[l] < order[o]){
        return order[l];
    }
    else{
        return order[o];
    }
}

/**
 * Returns the value of nearest neighbor
 * function of DLiPC around CHOL
 * for a given... 
 * 
 * @param ord number of CHOL neighbors.
 */
double NN_CD_L(int neighbor)
{
    if (neighbor == 0)
    {
      return (6.380169972613336 - 0.006504015670675491*read_temperature);
    }
    if (neighbor == 1)
    {
      return (6.191307287280158 -0.006244334905148628*read_temperature);
    }
    if (neighbor == 2)
    {
      return (6.023217470971802 -0.006163443775914513*read_temperature);
    }
    if (neighbor == 3)
    {
      return (5.688301398563751 -0.005732889536537225*read_temperature);
    }
    if (neighbor == 4)
    {
      return (4.953838034411138 -0.0041862507153207894*read_temperature);
    }
}

/**
 * Returns the value of nearest neighbor
 * function of DPPC around CHOL
 * for a given... 
 * 
 * @param ord number of CHOL neighbors.
 */
double NN_CD_P(int neighbor)
{
    if (neighbor == 0)
    {
      return (6.149307668194549 + (-1.2123337139100219/ 
      (1 + exp(-0.279420884957051 * (read_temperature - 323.5622172239552)))));
    }
    if (neighbor == 1)
    {
      return (6.132031203205369 + (-1.3391102930661072 / 
      (1 + exp(-0.2848108023841632 * (read_temperature - 323.4614496375261)))));
    }
    if (neighbor == 2)
    {
      return (6.09118413190892 + (-1.431194943529606 / 
      (1 + exp(-0.2632957680847098 * (read_temperature - 323.22226054097797)))));
    }
    if (neighbor == 3)
    {
      return (6.0063999161659885 + (-1.5233708107172301 / 
      (1 + exp(-0.26322941684422435 * (read_temperature - 323.3131242555826)))));
    }
    if (neighbor == 4)
    {
      return (5.836996032158986 + (-1.4812685885334167	 / 
      (1 + exp(-0.23159639245144456 * (read_temperature - 322.96436681509624)))));
    }
}

/**
 * Returns the value of the nearest neighbor
 * function of DLiPC around DLiPC
 * for a given... 
 * 
 * @param ord number of CHOL neighbors.
 */
double NN_DD_L(int neighbor)
{
    if (neighbor == 0)
    {
      return (5.103242510197753 - 0.005085225361896779*read_temperature);
    }
    if (neighbor == 1)
    {
      return (5.113452890355549 - 0.005357476571283622*read_temperature);
    }
    if (neighbor == 2)
    {
      return (5.041564491528647 - 0.005448578551134591*read_temperature);
    }
    if (neighbor == 3)
    {
      return (4.907827350161882 - 0.005455112431505342*read_temperature);
    }
    if (neighbor == 4)
    {
      return (4.437067046342299 - 0.004487121970234406*read_temperature);
    }
}  

/**
 * Returns the value of the nearest neighbor
 * function of DPPC around DPPC
 * for a given... 
 * 
 * @param ord number of CHOL neighbors.
 */
double NN_DD_P(int neighbor)
{
    if (neighbor == 0)
    {
      return (5.29864598362317 + (-1.046036898442689/ 
      (1 + exp(-0.25679510115746684 * (read_temperature - 324.0478762880865)))));
    }
    if (neighbor == 1)
    {
      return (5.291045158728617 + (-1.155530326754457 / 
      (1 + exp(-0.28882041884126536 * (read_temperature - 324.2842194308451)))));
    }
    if (neighbor == 2)
    {
      return (5.2916826656954115 + (-1.2845977425568815 / 
      (1 + exp(-0.28664779185612405 * (read_temperature - 323.9102671158312)))));
    }
    if (neighbor == 3)
    {
      return (5.241010005432148 + (-1.3888048522693892 / 
      (1 + exp(-0.2828101614730355 * (read_temperature - 323.715223040918)))));
    }
    if (neighbor == 4)
    {
      return (5.1776041414015666 + (-1.4834873657631906 / 
      (1 + exp(-0.27016244921668536 * (read_temperature - 323.4582319860251)))));
    }
}

/**
 * Returns the value of DLiPC's self interaction
 * for a given... 
 * 
 * @param ord order parameter.
 */
double H_self_L(int ord)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    return (-198.2144241389774
    - 5.767114535779552*o 
    - 21.392325238025208*o*o
    + 126.763192813133*o*o*o
    - 168.20406177882367*o*o*o*o);
}

/**
 * Returns the value of DPPC's self interaction
 * for a given... 
 * 
 * @param ord order parameter.
 */
double H_self_P(int ord)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    return (-223.3148767446365
    + 11.182283360518415*o 
    - 24.886674151265726*o*o
    - 23.449813607231356*o*o*o
    + 9.052427099706513*o*o*o*o);
}

/**
 * Returns the value of the DLiPC-DLiPC interaction 
 * function for a given...
 * 
 * @param ord MINIMAL order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_DD_LL_min(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.2385777507702)
      + (-85.0/(1+exp(-5.5152388615934305*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-61.907988669194644)
      + (-85.0/(1+exp(-5.434468636915715*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.51808320647846)
      + (-85.0/(1+exp(-5.7222222085142835*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.80010184207919)
      + (-85.0/(1+exp(-6.046112173136968*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-54.10623762430083)
      + (-85.0/(1+exp(-5.73206006922409*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-52.25671633964124)
      + (-85.0/(1+exp(-5.7767043757695635*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-49.57073396130829)
      + (-65.0/(1+exp(-5.917923875494118*(o-1.0)))));
    }
}

/**
 * Returns the value of the DLiPC-DLiPC interaction 
 * function for a given...
 * 
 * @param ord AVERAGE order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_DD_LL_avg(int ord, int neighbor)
{
    double p = (double)ord/2.0;
    double o = (p-50.0)/100.0;

    double m = o*o;

    if (o > (0.0)) return (-51.34502*m*o + 11.80470*m - 0.20345*o - 64.72701 + 2.5*neighbor);

    else return (-64.35 + 2.5*neighbor);
}

/**
 * Returns the value of the DPPC-DPPC interaction 
 * function for a given...
 * 
 * @param ord MINIMAL order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_DD_PP_min(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.01505834779738)
      + (-85.0/(1+exp(-5.790501529508871*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-62.140061908466556)
      + (-71.26511634433881/(1+exp(-5.738073754577456*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.13452426429441)
      + (-67.70009843013932/(1+exp(-5.5794520397459575*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.25342719515398)
      + (-67.07890880720473/(1+exp(-5.6396362558218724*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-52.77151128708946)
      + (-70.34178084725163/(1+exp(-5.87616410257819*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-49.71869099631439)
      + (-65.0/(1+exp(-6.040477808247997*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-43.935239780557815)
      + (-76.95590795060127/(1+exp(-5.0*(o-1.0)))));
    }
}

/**
 * Returns the value of the DPPC-DPPC interaction 
 * function for a given...
 * 
 * @param ord AVERAGE order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_DD_PP_avg(int ord, int neighbor)
{
    double p = (double)ord/2.0;
    double o = (p-50.0)/100.0;

    double m = o*o;

    if (neighbor == 0){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286);
        else return (-64.0);
    }

    if (neighbor == 1){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90);
        else return (-64.0 + 1.90);
    }

    if (neighbor == 2){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90 + 3.00);
        else return (-64.0 + 1.90 + 3.00);
    }

    if (neighbor == 3){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90 + 3.00*2);
        else return (-64.0 + 1.90 + 3.00*2);
    }

    if (neighbor == 4){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90 + 3.00*2 + 3.51);
        else return (-64.0 + 1.90 + 3.00*2 + 3.51);
    }

    if (neighbor == 5){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90 + 3.00*3 + 3.51);
        else return (-64.0 + 1.90 + 3.00*3 + 3.51);
    }

    if (neighbor == 6){
        if (o > (-0.1)) return (-34.14372*m*o + 4.44322*m - 0.29837*o - 64.32286 + 1.90 + 3.00*4 + 3.51);
        else return (-64.0 + 1.90 + 3.00*4 + 3.51);
    }
}

/**
 * Returns the value of the DPPC-DLiPC interaction 
 * function for a given...
 * 
 * @param ord MINIMAL order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_DD_LP_min(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.00634092896922)
      + (-85.0/(1+exp(-5.3300771721508475*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-61.50693262332394)
      + (-85.0/(1+exp(-5.8495526115388845*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.08429503488477)
      + (-85.0/(1+exp(-6.234155020540408*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.3906908602111)
      + (-68.39937312113554/(1+exp(-6.0663592559087*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-52.790928365243936)
      + (-85.0/(1+exp(-5.627236498118878*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-49.083301346156844)
      + (-85.0/(1+exp(-5.589864043019059*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-47.150888458064365)
      + (-85.0/(1+exp(-5.621171915591369*(o-1.0)))));
    }
}

/**
 * Returns the value of the DPPC-DLiPC interaction 
 * function for a given...
 * 
 * @param ord AVERAGE order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 * 
 * @warning WIP!
 */
double H_DD_LP_avg(int ord,int neighbor) 
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
      return ((-64.00634092896922)
      + (-85.0/(1+exp(-5.3300771721508475*(o-1.0)))));
    }
    if (neighbor == 1)
    {
      return ((-61.50693262332394)
      + (-85.0/(1+exp(-5.8495526115388845*(o-1.0)))));
    }
    if (neighbor == 2)
    {
      return ((-59.08429503488477)
      + (-85.0/(1+exp(-6.234155020540408*(o-1.0)))));
    }
    if (neighbor == 3)
    {
      return ((-56.3906908602111)
      + (-68.39937312113554/(1+exp(-6.0663592559087*(o-1.0)))));
    }
    if (neighbor == 4)
    {
      return ((-52.790928365243936)
      + (-85.0/(1+exp(-5.627236498118878*(o-1.0)))));
    }
    if (neighbor == 5)
    {
      return ((-49.083301346156844)
      + (-85.0/(1+exp(-5.589864043019059*(o-1.0)))));
    }
    if (neighbor == 6)
    {
      return ((-47.150888458064365)
      + (-85.0/(1+exp(-5.621171915591369*(o-1.0)))));
    }
}

/**
 * Returns the value of the CHOL-CHOL interaction 
 * function derived from MD simulations with DLiPC 
 * for a given...
 * 
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_CC_L(int neighbor)
{
    return (-23.73548301296039 + 1.7338424224038076*neighbor);
}

/**
 * Returns the value of the CHOL-CHOL interaction 
 * function derived from MD simulations with DPPC 
 * for a given...
 * 
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_CC_P(int neighbor)
{
    return (-21.406267414625454 + 1.4939079023403155*neighbor);
}

/**
 * Returns the value of the CHOL-DLiPC interaction 
 * function for a given...
 * 
 * @param ord order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_CD_L(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    double n = (double)neighbor;
    //return (-20*o-28+2*n);
    if (neighbor == 0)
    {
      double q = pow(((0.7539372160118417-1.5)/(o-1.5)),1.6592956112169561);
      return (36.483321544385326*q*(q-2.0));
    }
    if (neighbor == 1)
    {
      double q = pow(((0.7641800265641312-1.5)/(o-1.5)),1.6255609292682194);
      return (34.13691226888765*q*(q-2.0));
    }
    if (neighbor == 2)
    {
      double q = pow(((0.7699158189596492-1.5)/(o-1.5)),1.8291013645164085);
      return (33.68025196550303*q*(q-2.0));
    }
    if (neighbor == 3)
    {
      double q = pow(((0.7228854333900939-1.5)/(o-1.5)),2.4111347746661385);
      return (30.4277779459905*q*(q-2.0));
    }
    if (neighbor == 4)
    {
      double q = pow(((0.6740667747083389-1.5)/(o-1.5)),2.5);
      return (26.601605192448396*q*(q-2.0));
    }
    if (neighbor == 5)
    {
      double q = pow(((0.7419107199339535-1.5)/(o-1.5)),2.2050855814147825);
      return (26.949507831362077*q*(q-2.0));
    }
}

/**
 * Returns the value of the CHOL-DPPC interaction 
 * function for a given...
 * 
 * @param ord order parameter.
 * @param neighbor number of CHOL neighbors around the pair.
 */
double H_CD_P(int ord,int neighbor)
{
    double p = (double)ord;
    double o = (p-50.0)/100;
    if (neighbor == 0)
    {
        double q = pow(((0.8015034812730596-1.5)/(o-1.5)),1.2104586776529458);
        return (34.162896996372645*q*(q-2.0));
    }
    if (neighbor == 1)
    {
        double q = pow(((0.7469858015529757-1.5)/(o-1.5)),1.4201657855084255);
        return (31.950030406520444*q*(q-2.0));
    }
    if (neighbor == 2)
    {
        double q = pow(((0.7529462715453211-1.5)/(o-1.5)),1.4532308249852621);
        return (29.894418246485863*q*(q-2.0));
    }
    if (neighbor == 3)
    {
        double q = pow(((0.7569774323440122-1.5)/(o-1.5)),1.3923986491593072);
        return (27.70235738706128*q*(q-2.0));
    }
    if (neighbor == 4)
    {
        double q = pow(((0.7613793423645616-1.5)/(o-1.5)),1.6463683453577225);
        return (25.6221899051367*q*(q-2.0));
    }
    if (neighbor == 5)
    {
        double q = pow(((0.8164971196090005-1.5)/(o-1.5)),1.3387806908751303);
        return (23.2981162468298*q*(q-2.0));
    }
}

/**
 * Returns the value of the entropic function of DLiPC for a...
 * 
 * @param ord given order parameter.
 */
double E_L(int ord)
{
    double entro_avg[] = {-32.585663	,
-30.866203	,
-29.232518	,
-27.681936	,
-26.211767	,
-24.819306	,
-23.501839	,
-22.256655	,
-21.081049	,
-19.972329	,
-18.927823	,
-17.944885	,
-17.020896	,
-16.153274	,
-15.339477	,
-14.577006	,
-13.86341	,
-13.196286	,
-12.57329	,
-11.992131	,
-11.450581	,
-10.946472	,
-10.4777	,
-10.042228	,
-9.638087	,
-9.263374	,
-8.916259	,
-8.59498	,
-8.297849	,
-8.023245	,
-5.879377	,
-5.774855	,
-5.680293	,
-5.587116	,
-5.495373	,
-5.40958	,
-5.256765	,
-5.124981	,
-4.841459	,
-4.879505	,
-4.76921	,
-4.672118	,
-4.584211	,
-4.367889	,
-4.280518	,
-4.185753	,
-4.147191	,
-4.068844	,
-4.038433	,
-3.994268	,
-3.952651	,
-3.897551	,
-3.86577	,
-3.823503	,
-3.795968	,
-3.759546	,
-3.742374	,
-3.714611	,
-3.702783	,
-3.690262	,
-3.667644	,
-3.66403	,
-3.649278	,
-3.645	,
-3.644549	,
-3.651427	,
-3.655617	,
-3.665521	,
-3.682499	,
-3.701393	,
-3.724392	,
-3.75123	,
-3.783737	,
-3.820218	,
-3.864088	,
-3.909537	,
-3.960334	,
-4.018678	,
-4.082179	,
-4.151228	,
-4.2264	,
-4.30541	,
-4.394138	,
-4.486408	,
-4.58696	,
-4.693939	,
-4.806593	,
-4.92647	,
-5.054721	,
-5.188704	,
-5.332174	,
-5.484733	,
-5.641997	,
-5.809017	,
-5.987276	,
-6.174251	,
-6.366742	,
-6.573127	,
-6.789823	,
-7.018619	,
-7.259292	,
-7.51496	,
-7.782613	,
-8.063208	,
-8.362974	,
-8.677154	,
-9.008832	,
-9.362476	,
-9.73028	,
-10.127916	,
-10.547617	,
-10.9955	,
-11.468903	,
-11.974876	,
-12.512817	,
-13.090008	,
-13.629436	,
-14.195804	,
-14.805544	,
-15.447538	,
-16.133729	,
-13.942953	,
-14.676721	,
-15.473254	,
-16.337913	,
-17.276442	,
-18.294987	,
-19.400117	,
-20.598843	,
-21.898641	,
-23.307475	,
-24.833817	,
-26.486672	,
-28.275603	,
-30.210754	,
-32.302876	,
-34.563352	,
-37.004226	,
-39.638228	,
-42.478806	,
-45.54015	,
-48.837225	,
-52.3858	,
-56.20248	,
-60.304737	,
-64.710944	,
-69.440406	,
-74.513399	,
-79.951198	,
-85.776118	,
-92.011551};
double entro_min[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, -2.12944, -2.037605, -1.937061, -1.846972, -1.753329, -1.67774, -1.530248, -1.384189, -1.267692, -1.154616, -1.052285, -0.949279, -0.857928, -0.781794, -0.697712000000001, -0.630328, -0.565137, -0.502534, -0.448644, -0.391088, -0.342475, -0.297835, -0.25408, -0.218283, -0.182523, -0.14943, -0.119477, -0.0923769999999999, -0.0690029999999999, -0.0513080000000001, -0.0345270000000002, -0.0205630000000001, -0.012143, -0.00445800000000007, 0, -0.00146200000000007, -0.00513300000000028, -0.0130000000000003, -0.0229170000000001, -0.0382600000000002, -0.0579700000000001, -0.0829219999999999, -0.110859, -0.143277, -0.181157, -0.224968, -0.272803, -0.326512, -0.384787, -0.449207, -0.519846, -0.595263000000001, -0.676419000000001, -0.7646, -0.858887, -0.958931, -1.066007, -1.179457, -1.299697, -1.426658, -1.561501, -1.703715, -1.853455, -2.012468, -2.177832, -2.352923, -2.537326, -2.733284, -2.935633, -3.150797, -3.376641, -3.615728, -3.866462, -4.13049, -4.41024, -4.704347, -5.017327, -5.347598, -5.696178, -6.068123, -6.461887, -6.880999, -7.324958, -7.80177, -8.310836, -8.852782, -9.360283, -9.893097, -10.469662, -11.07636, -11.71812, 121, 122, 123, 124, 125, 126, 127, 128, 129, 130, 131, 132, 133, 134, 135, 136, 137, 138, 139, 140, 141, 142, 143, 144, 145, 146, 147, 148, 149, 150};
    
    if (int_type == 0) return entro_avg[ord];
    if (int_type == 1) return entro_min[ord];
}


/**
 * Returns the value of the entropic function of DPPC for a...
 * 
 * @param ord given order parameter.
 */
double E_P(int ord)
{
    double entro_avg[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,-4.118937,-3.852361,-3.772671,-3.697232,-3.546854,-3.496151,-3.340507,-3.16266,-2.961434,-2.842948,-2.713178,-2.532163,-2.389425,-2.249913,-2.112968,-1.978765,-1.802305,-1.673868,-1.367041,-1.26247,-1.153913,-1.039294,-1.019863,-0.871165,-0.708574,-0.632807,-0.549904,-0.469855,-0.402968,-0.350183,-0.28764,-0.243385,-0.194656,-0.177376,-0.125448,-0.122756,-0.078265,-0.0648940000000002,-0.057086,-0.077941,0,-0.0793690000000002,-0.0952250000000001,-0.126744,-0.0963850000000002,-0.144892,-0.250412,-0.246676,-0.350902,-0.404134,-0.4623,-0.572485,-0.649316,-0.753291,-0.866856,-0.987921,-1.12086,-1.259962,-1.409906,-1.56693,-1.691124,-1.888587,-2.077237,-2.256544,-2.464971,-2.680907,-2.909033,-3.147398,-3.399947,-3.659804,-3.929231,-4.191412,-4.489291,-4.793285,-5.095529,-5.440075,-5.76721,-6.116733,-6.484325,-6.852125,-7.257355,-7.63348,-8.052163,-8.486288,-8.930542,-9.379204,-9.853483,-10.340207,-10.842307,-11.354913,-11.886182,-12.42949,-12.990063,-13.562656,-14.152452,-14.755464,-15.375649,-16.008366,-16.660715,-17.328951,-18.008317,-18.708756,-19.420442,-20.15205,-20.900657,-21.669096,-22.447564,-23.244881,-24.051944,-24.88605,-25.729152,-26.597938,-27.497828,-28.38938,-29.287956,-30.240541,-31.08282,-31.97013,-32.841195,-33.703116,-34.599741};
    double entro_min[] = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,-4.118937,-3.852361,-3.772671,-3.697232,-3.546854,-3.496151,-3.340507,-3.16266,-2.961434,-2.842948,-2.713178,-2.532163,-2.389425,-2.249913,-2.112968,-1.978765,-1.802305,-1.673868,-1.367041,-1.26247,-1.153913,-1.039294,-1.019863,-0.871165,-0.708574,-0.632807,-0.549904,-0.469855,-0.402968,-0.350183,-0.28764,-0.243385,-0.194656,-0.177376,-0.125448,-0.122756,-0.078265,-0.0648940000000002,-0.057086,-0.077941,0,-0.0793690000000002,-0.0952250000000001,-0.126744,-0.0963850000000002,-0.144892,-0.250412,-0.246676,-0.350902,-0.404134,-0.4623,-0.572485,-0.649316,-0.753291,-0.866856,-0.987921,-1.12086,-1.259962,-1.409906,-1.56693,-1.691124,-1.888587,-2.077237,-2.256544,-2.464971,-2.680907,-2.909033,-3.147398,-3.399947,-3.659804,-3.929231,-4.191412,-4.489291,-4.793285,-5.095529,-5.440075,-5.76721,-6.116733,-6.484325,-6.852125,-7.257355,-7.63348,-8.052163,-8.486288,-8.930542,-9.379204,-9.853483,-10.340207,-10.842307,-11.354913,-11.886182,-12.42949,-12.990063,-13.562656,-14.152452,-14.755464,-15.375649,-16.008366,-16.660715,-17.328951,-18.008317,-18.708756,-19.420442,-20.15205,-20.900657,-21.669096,-22.447564,-23.244881,-24.051944,-24.88605,-25.729152,-26.597938,-27.497828,-28.38938,-29.287956,-30.240541,-31.08282,-31.97013,-32.841195,-33.703116,-34.599741};
    
    if (int_type == 0) return entro_avg[ord];
    if (int_type == 1) return entro_min[ord];
}

/**
 * Returns an initial guess for the entropic function for DLiPC.
 * (According to the Boltzmann distribution.)
 * 
 * @param ord Given order Parameter.
 */
double E_init_L(int ord)
{
    double temp = temperature*120;
    double init = (log(p_MD[ord])) + (1/(temp*0.00831448664))*(0.5*H_dd_avg[1][1][ord][0]*nneigh_dd[1][0] + H_self[1][ord]);
    return init;
}

/**
 * Returns an initial guess for the entropic function for DPPC.
 * (According to the Boltzmann distribution.)
 * 
 * @param ord Given order Parameter.
 */
double E_init_P(int ord)
{
    double temp = temperature*120;
    double init = (log(p_MD[ord])) + (1/(temp*0.00831448664))*(0.5*H_dd_avg[1][1][ord][0]*nneigh_dd[1][0] + H_self[1][ord]);
    return init;
}

/**
 * Returns the value of Phi_DLiPC for a...
 * 
 * @param ord given order parameter.
 * @param neighbor given number of CHOL around the lipid.
 */
double Phi_L(int ord,int neighbor) 
{
    double s = ((double)ord-50.0)/100.0;
    double n = (double)neighbor;

    return (-27.324707*s*s+10.250992*s)*(n-0.00432*n*n);
}

/**
 * Returns the value of Phi_DPPC for a...
 * 
 * @param ord given order parameter.
 * @param neighbor given number of CHOL around the lipid.
 * 
 * @version 0.0.3 Added restriction for four neighbors below a certain S value.
 */
double Phi_P(int ord,int neighbor)
{
    double s = ((double)ord-50.0)/100.0;
    double n = (double)neighbor;

    if(n == 4 && ord <= 92) return (-30.344);
    return (-57.274922*s*s+48.593694*s)*(n-0.434*n*n);
}

/**
 * Initalizes a first guess for Phi = (a*S+b*S*S)*(N_Chol + c*N_Chol*N_Chol)
 * @warning WIP!
 */
double Phi_init_L(int ord,int neighbor) // WIP!
{
    return 0.0;
}

/**
 * Initalizes a first guess for Phi = (a*S+b*S*S)*(N_Chol + c*N_Chol*N_Chol)
 * @warning WIP!
 */
double Phi_init_P(int ord,int neighbor) // WIP!
{
    return 0.0;
}

/**
 * Calls the order parameter distribution from
 * the MD simulations.
 */
void MD_order()
{
    double norm_sum = 0.0;

    for (int s=0; s<=150; s++){
        double p = (double)s;
        double o = (p-50.0)/100.0;
        if (lipid_type != 0){
            if (c_chol == 0.0){
                p_MD[s] = (exp(-0.14122 + 7.51277*o -9.36903*o*o -4.43679*o*o*o -97.86418*o*o*o*o +192.92704*o*o*o*o*o+19.37517*o*o*o*o*o*o -168.20577*o*o*o*o*o*o*o));
            }
            if (c_chol == 0.1){
                p_MD[s] = norm_dist(o,0.275,0.115);
            }
            if (c_chol == 0.2){
                p_MD[s] = norm_dist(o,0.31,0.1);
            }
            if (c_chol == 0.3){
                p_MD[s] = norm_dist(o,0.36,0.13);
            }
        }
        if (lipid_type == 0){
            if (c_chol == 0.0){
                p_MD[s] = (exp(-0.9767356 + 8.69286553*o -12.7808724*o*o +12.12000201*o*o*o -21.41776641*o*o*o*o + 7.14478559*o*o*o*o*o));
            }
            if (c_chol == 0.1){
                p_MD[s] = norm_dist(o,0.45,0.1);
            }
            if (c_chol == 0.2){
                p_MD[s] = norm_dist(o,0.55,0.15);
            }
            if (c_chol == 0.3){
                p_MD[s] = norm_dist(o,0.7,0.05);
            }
        }
    }

    for (int s=0; s<=150; s++){
        norm_sum += p_MD[s];
    }

    for (int s=0; s<=150; s++){
      p_MD[s] /= (norm_sum);
    }
}

/**
 * Initializes the nearest neighbor lists.
 * (A little bit complicated, tinker at your own risk.)
 */
void make_neighbors()
{
    int i,j,k,l,vec[4],vecn[4],m,keep; 
    int c_1[4][10],c_2[4][10] ;

    nr_lattice = 1;
    for (i=1 ; i <= dim ; i++) 
    {
      nr_lattice *= length[i]; //number of lattice squares is product of lengths in all directions
    }
    nr_chol = (int)(c_chol * nr_lattice);
    if (nr_lattice > N) printf("ERROR: More Particles than Lattice Sites. \n");
    c_chol = (double)(nr_chol)/(double)(nr_lattice);

    for (i=0 ; i < nr_lattice ; i++) 
    {
        transfer1(i,vec);
        l = 0;
        for (j=1 ; j<= dim ; j++)
        for (k=-1 ; k <= 1 ; k+=2)
        {
            keep = vec[j];
            vec[j] = per(vec[j] +  k,j);
            m = transfer2(vec);
            neigh[i][l] = m;   // counts the neighbors of site i, i.e. neigh[i][0],...,neigh[i][2*dim-1]
            vec[j] = keep;
            l++;
        }

        transfer1(i,vec);
        vecn[1] = per(vec[1] -1,1); vecn[2] = per(vec[2] -1,2); neigh_dc[i][0] = transfer2(vecn);
        vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] ,2); neigh_dc[i][1] = transfer2(vecn);
        vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] -1,2); neigh_dc[i][2] = transfer2(vecn);
        vecn[1] = per(vec[1] - 1 ,1); vecn[2] = per(vec[2] ,2); neigh_dc[i][3] = transfer2(vecn);

        transfer1(i,vec);
        vecn[1] = per(vec[1] ,1); vecn[2] = per(vec[2] ,2); neigh_cd[i][0] = transfer2(vecn);
        vecn[1] = per(vec[1]+1 ,1); vecn[2] = per(vec[2]+1 ,2); neigh_cd[i][1] = transfer2(vecn);
        vecn[1] = per(vec[1]+1 ,1); vecn[2] = per(vec[2] ,2); neigh_cd[i][2] = transfer2(vecn);
        vecn[1] = per(vec[1]  ,1); vecn[2] = per(vec[2]+1 ,2); neigh_cd[i][3] = transfer2(vecn);

        transfer1(i,vec);
        vecn[1] = per(vec[1],1); vecn[2] = per(vec[2]-1,2); neigh_cd_2[i][0] = transfer2(vecn);
        vecn[1] = per(vec[1],1); vecn[2] = per(vec[2]+2,2); neigh_cd_2[i][1] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]-1,2); neigh_cd_2[i][2] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2],2); neigh_cd_2[i][3] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]+1,2); neigh_cd_2[i][4] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]+2,2); neigh_cd_2[i][5] = transfer2(vecn);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]-1,2); neigh_cd_2[i][6] = transfer2(vecn);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]+2,2); neigh_cd_2[i][7] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2]-1,2); neigh_cd_2[i][8] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2],2); neigh_cd_2[i][9] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2]+1,2); neigh_cd_2[i][10] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2]+2,2); neigh_cd_2[i][11] = transfer2(vecn);

        transfer1(i,vec);
        vecn[1] = per(vec[1],1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][0] = transfer2(vecn);
        vecn[1] = per(vec[1],1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][1] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][2] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][3] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][4] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2]-1,2); neigh_cd_3[i][5] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2],2); neigh_cd_3[i][6] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2]+1,2); neigh_cd_3[i][7] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2]+2,2); neigh_cd_3[i][8] = transfer2(vecn);
        vecn[1] = per(vec[1]-2,1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][9] = transfer2(vecn);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][10] = transfer2(vecn);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][11] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][12] = transfer2(vecn);
        vecn[1] = per(vec[1]+2,1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][13] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2]-2,2); neigh_cd_3[i][14] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2]-1,2); neigh_cd_3[i][15] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2],2); neigh_cd_3[i][16] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2]+1,2); neigh_cd_3[i][17] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2]+2,2); neigh_cd_3[i][18] = transfer2(vecn);
        vecn[1] = per(vec[1]+3,1); vecn[2] = per(vec[2]+3,2); neigh_cd_3[i][19] = transfer2(vecn);

        transfer1(i,vec);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]+1,2); neigh_moore[i][0] = transfer2(vecn);
        vecn[1] = per(vec[1]+1,1); vecn[2] = per(vec[2]-1,2); neigh_moore[i][1] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]+1,2); neigh_moore[i][2] = transfer2(vecn);
        vecn[1] = per(vec[1]-1,1); vecn[2] = per(vec[2]-1,2); neigh_moore[i][3] = transfer2(vecn);

    }
  // Starting point is always  (2,2); 

    cmax_dd = 6;
    c_1[0][1] = 0;  c_2[0][1] = 1;    
    c_1[0][2] = 0;  c_2[0][2] = 2;    
    c_1[0][3] = 1;  c_2[0][3] = 1;    
    c_1[0][4] = 1;  c_2[0][4] = 2;    
    c_1[0][5] = 2;  c_2[0][5] = 1;    
    c_1[0][6] = 2;  c_2[0][6] = 2;    

    c_1[1][1] = 1;  c_2[1][1] = 1;    
    c_1[1][2] = 1;  c_2[1][2] = 2;    
    c_1[1][3] = 2;  c_2[1][3] = 1;    
    c_1[1][4] = 2;  c_2[1][4] = 2;    
    c_1[1][5] = 3;  c_2[1][5] = 1;    
    c_1[1][6] = 3;  c_2[1][6] = 2;    

    c_1[2][1] = 1;  c_2[2][1] = 0;    
    c_1[2][2] = 2;  c_2[2][2] = 0;    
    c_1[2][3] = 1;  c_2[2][3] = 1;    
    c_1[2][4] = 2;  c_2[2][4] = 1;    
    c_1[2][5] = 1;  c_2[2][5] = 2;    
    c_1[2][6] = 2;  c_2[2][6] = 2;    

    c_1[3][1] = 1;  c_2[3][1] = 1;    
    c_1[3][2] = 2;  c_2[3][2] = 1;    
    c_1[3][3] = 1;  c_2[3][3] = 2;    
    c_1[3][4] = 2;  c_2[3][4] = 2;    
    c_1[3][5] = 1;  c_2[3][5] = 3;    
    c_1[3][6] = 2;  c_2[3][6] = 3;    


    for (i=0 ; i < nr_lattice ; i++) 
    { 
      for (j=0 ; j < 2*dim ; j++)
      {
        for (l = 1 ; l <= cmax_dd; l++)
        {
          transfer1(i,vec); 
          vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
          vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
          neigh2_dd[i][j][l] = transfer2(vecn);
        }
      }
    }

    cmax_cd = 5;
    c_1[0][1] = 1;  c_2[0][1] = 2;    
    c_1[0][2] = 2;  c_2[0][2] = 1;    
    c_1[0][3] = 3;  c_2[0][3] = 2;    
    c_1[0][4] = 2;  c_2[0][4] = 3;    
    c_1[0][5] = 1;  c_2[0][5] = 1;    

    c_1[1][1] = 1;  c_2[1][1] = 2;    
    c_1[1][2] = 2;  c_2[1][2] = 1;    
    c_1[1][3] = 3;  c_2[1][3] = 2;    
    c_1[1][4] = 2;  c_2[1][4] = 3;    
    c_1[1][5] = 3;  c_2[1][5] = 3;   

    c_1[2][1] = 1;  c_2[2][1] = 2;    
    c_1[2][2] = 2;  c_2[2][2] = 1;    
    c_1[2][3] = 3;  c_2[2][3] = 2;    
    c_1[2][4] = 2;  c_2[2][4] = 3;    
    c_1[2][5] = 3;  c_2[2][5] = 1;    

    c_1[3][1] = 1;  c_2[3][1] = 2;    
    c_1[3][2] = 2;  c_2[3][2] = 1;    
    c_1[3][3] = 3;  c_2[3][3] = 2;    
    c_1[3][4] = 2;  c_2[3][4] = 3;    
    c_1[3][5] = 1;  c_2[3][5] = 3; 

    for (i=0 ; i < nr_lattice ; i++) 
    { 
      for (j=0 ; j < 2*dim ; j++)
      {
        for (l = 1 ; l <= cmax_cd; l++)
        {
          transfer1(i,vec); 
          vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
          vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
          neigh2_cd[i][j][l] = transfer2(vecn);
        }
      }
    }


    cmax_cc = 6;
    c_1[0][1] = 0;  c_2[0][1] = 2;    
    c_1[0][2] = 1;  c_2[0][2] = 1;    
    c_1[0][3] = 1;  c_2[0][3] = 3;    
    c_1[0][4] = 2;  c_2[0][4] = 1;    
    c_1[0][5] = 2;  c_2[0][5] = 3;    
    c_1[0][6] = 3;  c_2[0][6] = 2;    

    c_1[1][1] = 1;  c_2[1][1] = 2;    
    c_1[1][2] = 2;  c_2[1][2] = 1;    
    c_1[1][3] = 2;  c_2[1][3] = 3;    
    c_1[1][4] = 3;  c_2[1][4] = 1;    
    c_1[1][5] = 3;  c_2[1][5] = 3;    
    c_1[1][6] = 4;  c_2[1][6] = 2;    

    c_1[2][1] = 2;  c_2[2][1] = 0;    
    c_1[2][2] = 1;  c_2[2][2] = 1;    
    c_1[2][3] = 3;  c_2[2][3] = 1;    
    c_1[2][4] = 1;  c_2[2][4] = 2;    
    c_1[2][5] = 3;  c_2[2][5] = 2;    
    c_1[2][6] = 2;  c_2[2][6] = 3;    

    c_1[3][1] = 2;  c_2[3][1] = 1;    
    c_1[3][2] = 1;  c_2[3][2] = 2;    
    c_1[3][3] = 3;  c_2[3][3] = 2;    
    c_1[3][4] = 1;  c_2[3][4] = 3;    
    c_1[3][5] = 3;  c_2[3][5] = 3;    
    c_1[3][6] = 2;  c_2[3][6] = 4;

    for (i=0 ; i < nr_lattice ; i++) 
    { 
      for (j=0 ; j < 2*dim ; j++)
      {
        for (l = 1 ; l <= cmax_dd; l++)
        {
          transfer1(i,vec); 
          vecn[1] = per(c_1[j][l]+vec[1]-2,1);  
          vecn[2] = per(c_2[j][l]+vec[2]-2,2);  
          neigh2_cc[i][j][l] = transfer2(vecn);
        }
      }
    }
}

/**
 * Calculates the number of specific 
 * particles for the simulation.
 */
void make_particle_numbers()
{
    int i;
    int sum = 0;
    double num = 0.0;

    for (i=1;i<nr_types;i++){
        num = (double)(nr_lattice)*c_lipid[i];
        nrt[i] = (int)round(num);
        sum += nrt[i];
    }
    nrt[nr_types] = nr_lattice - sum; 
}

/**
 * Initializes the necessary output parameters.
 */
void make_initialization()
{
    int i,s,c;
    for (i=1;i<=nr_types;i++){
        for (s=0;s<=150;s++){
            p_MC[i][s] = 0;
            p_MC_1[i][s] = 0;
            p_MC_2[i][s] = 0;
            p_MC_3[i][s] = 0;
            for (c=0;c<=4;c++){
                p_MC_neigh[i][s][c] = 0;
            }
        }
    }
    cc_contacts = 0;
    cc_moore = 0;
    ener_mem = 0;

    for (int j=0; j<=4; j++){
        contacts_ener[j] = 0;
        contacts_num[j] = 0;
    }

    for (int j=0; j<=10; j++){
        neighpair_cd[j] = 0;
    }

    for (int l=0; l<nr_lattice;l++){
        lattice_sites[l] = l;
    }

    for (i=0;i<=9;i++){
        for (s=0;s<=150;s++){
            H_dd_min_track[0][s][i] = 0.0;
            H_dd_min_track[1][s][i] = 0.0;
        }
        for (s=0;s<=300;s++){
            H_dd_avg_track[0][s][i] = 0.0;
            H_dd_avg_track[1][s][i] = 0.0;
        }
    }
}

/**
 * Initializes the interaction and nearest neighbor
 * functions.
 */
void make_interaction()
{
    int a,b,c;
    int s;

    if (lipid_type == 1 && nr_types == 1)
    {
        for (a=0;a<=4;a++){
            nneigh_dd[1][a] = NN_DD_L(a);
            nneigh_cd[1][a] = NN_CD_L(a);
        }
        for (s=0;s<=150;s++){
            H_self[1][s] = H_self_L(s);

            if (lat_type == 1) Entro[1][s] = E_init_L(s);
            else Entro[1][s] = E_L(s);

            for (c=0;c<=4;c++){
                if (lat_type == 2) Phi[1][s][c] = Phi_init_L(s,c);
                else Phi[1][s][c] = Phi_L(s,c);
            }

            for (b=0;b<=6;b++){
                H_dd_min[1][1][s][b] = H_DD_LL_min(s,b);
                H_cc[b] = H_CC_L(b);
            }

            for (c=0;c<=5;c++) H_cd[1][s][c] = H_CD_L(s,c);
        }
        for (s=0;s<=300;s++){
            for (b=0;b<=6;b++){
                H_dd_avg[1][1][s][b] = H_DD_LL_avg(s,b);
            }
        }
    }

    if (lipid_type == 0 && nr_types == 1)
    {
        for (a=0;a<=4;a++){
            nneigh_dd[1][a] = NN_DD_P(a);
            nneigh_cd[1][a] = NN_CD_P(a);
        }
        for (s=0;s<=150;s++){
            H_self[1][s] = H_self_P(s);

            if (lat_type == 1) Entro[1][s] = E_init_P(s);
            else Entro[1][s] = E_P(s);

            for (c=0;c<=4;c++){
                if (lat_type == 2) Phi[1][s][c] = Phi_init_P(s,c);
                else Phi[1][s][c] = Phi_P(s,c);
            }

            for (b=0;b<=6;b++){
                H_dd_min[1][1][s][b] = H_DD_PP_min(s,b);
                H_cc[b] = H_CC_P(b);
            }

            for (c=0;c<=5;c++) H_cd[1][s][c] = H_CD_P(s,c);
        }
        for (s=0;s<=300;s++){
            for (b=0;b<=6;b++){
                H_dd_avg[1][1][s][b] = H_DD_PP_avg(s,b);
            }
        }
    }

    if (lipid_type == 1 && nr_types == 2)
    {
        for (a=0;a<=4;a++){
            nneigh_dd[1][a] = NN_DD_L(a);
            nneigh_cd[1][a] = NN_CD_L(a);

            nneigh_dd[2][a] = NN_DD_P(a);
            nneigh_cd[2][a] = NN_CD_P(a);
        }
        for (s=0;s<=150;s++){
            H_self[1][s] = H_self_L(s);
            H_self[2][s] = H_self_P(s);

            Entro[1][s] = E_L(s);
            Entro[2][s] = E_P(s);

            
            for (c=0;c<=4;c++){
                Phi[1][s][c] = Phi_L(s,c);
                Phi[2][s][c] = Phi_P(s,c);
            }

            for (b=0;b<=6;b++){
                H_dd_min[1][1][s][b] = H_DD_LL_min(s,b);
                H_dd_min[2][2][s][b] = H_DD_PP_min(s,b);
                H_dd_min[1][2][s][b] = H_dd_min[2][1][s][b] = H_DD_LP_min(s,b);
                H_cc[b] = H_CC_L(b);
            }

            for (c=0;c<=5;c++){
                H_cd[1][s][c] = H_CD_L(s,c);
                H_cd[2][s][c] = H_CD_P(s,c);
            }
        }

        for (s=0;s<=300;s++){
            for (b=0;b<=6;b++){
                H_dd_avg[1][1][s][b] = H_DD_LL_avg(s,b);
                H_dd_avg[2][2][s][b] = H_DD_PP_avg(s,b);
                H_dd_avg[1][2][s][b] = H_dd_min[2][1][s][b] = H_DD_LP_avg(s,b);
            }
        }
    }
}

/**
 * Returns the vector corresponding
 * to the given lattice site.
 * 
 * @param i Lattice site to transfer.
 * @param vec The given vector.
 */
void transfer1(int l, int* vec)
{
    int ll,m;

    ll = nr_lattice/length[1];
    for (m=1;m<=dim;m++){
        vec[m] = l/ll;
        l -= (ll*vec[m]);
        if (m<dim) ll /= length[m+1];
    }
}

/**
 * Returns the lattice site corresponding
 * to the given vector.
 * 
 * @param vec Vector to transfer.
 */
int transfer2(int *vec)
{
    int m;
    int l = 0;
    int ll = 1;

    for (m=1;m<=dim;m++){
        l += vec[dim+1-m]*ll;
        ll *= length[dim+1-m];
    }
    return(l);
}

/**
 * Periodic boundary conditions, when shifting
 * lattice sites. Returns the new site.
 * 
 * @param i A given lattice site.
 * @param d Direction of the lattice.
 */
int per(int i, int d) 
{
   if (i >= length[d] ) i -= length[d]; 
   if (i < 0) i += length[d]; 
   return(i);
}

/**
 * Returns the value of a normal distribution
 * for a given x (ord).
 * 
 * @param ord The value for the calculation.
 * @param mean Mean of the normal distribution.
 * @param sigma Standard deviation.
 */
double norm_dist(double ord, double mean, double sigma)
{
    double vorfaktor = sqrt(2*3.14159265359)*sigma;
    double exponent = (-0.5)*(pow((ord-mean),2))*(1/pow(sigma,2));
    double e_funktion = exp(exponent);

    return e_funktion/vorfaktor;
}

/**
 * Sorts the element of two arrays by their sizes.
 * 
 * @param n The range of the arrays that is sort. 
 * @param arr One of the arrays to sort.
 * @param brr The other array.
 */
void sort(int n, double *arr, int *brr)
{
    int l,i;
    double a,b;

    for (l=2;l<=n;l++){
        a = arr[l];
        b = brr[l];
        i = l-1;
        while (i > 0 && arr[i] > a){
            arr[i+1] = arr[i];
            brr[i+1] = brr[i];
            i--;
        }
        arr[i+1] = a;
        brr[i+1] = b;
    }
}

/**
 * Randomizes the elements of a given 
 * array via the Fisher-Yates-Shuffle.
 * 
 * @param array The array to shuffle.
 */
void fisher_yate_shuffle(int *array)
{
    int size = nr_lattice;
    for (int i = size - 1; i>0; i--){
        int j = rand() % (i+1);
        int temp = array[i];
        array[i] = array[j];
        array[j] = temp;
    }
}

/**
 * Reads the input file (conveniently named input.dat).
 */
void read_input()
{
    FILE *input;                                        // defines the input file
    char line[100];                                     // necessary for looping through the lines of the input file

    input = fopen("input.txt","r");                     // you can change the name of the input file here
    if (input == NULL) printf("ERROR: Input not found.\n");

    while (fgets(line,sizeof(line),input))              // now the input parameters will be read
    {
        if (strstr(line, "Project Name: ") != NULL) sscanf(line, "Project Name: %s",&project_name);

        if (strstr(line, "Lattice Type: ") != NULL) sscanf(line, "Lattice Type: %d",&lat_type);
        if (strstr(line, "Interaction Type: ") != NULL) sscanf(line, "Interaction Type: %d",&int_type);
        if (strstr(line, "Movement Type: ") != NULL) sscanf(line, "Movement Type: %d",&mv_type);

        if (strstr(line, "Dimension: ") != NULL) sscanf(line, "Dimension: %d",&dim);
        if (strstr(line, "Lattice Size: ") != NULL) sscanf(line, "Lattice Size: %d",&lattice_size);

        if (strstr(line, "Temperature: ") != NULL) sscanf(line, "Temperature: %lf",&read_temperature);
        if (strstr(line, "Number of Lipids: ") != NULL) sscanf(line, "Number of Lipids: %d",&nr_types);
        if (strstr(line, "Lipid Type: ") != NULL) sscanf(line, "Lipid Type: %d",&lipid_type);
        if (strstr(line, "Concentration (DLiPC): ") != NULL) sscanf(line, "Concentration (DLiPC): %lf", &c_lipid[1]);
        if (strstr(line, "Concentration (CHOL): ") != NULL) sscanf(line, "Concentration (CHOL): %lf", &c_chol);
        if (strstr(line, "Orderwidth: ") != NULL) sscanf(line, "Orderwidth: %d", &orderwidth);

        if (strstr(line, "Number of Runs: ") != NULL) sscanf(line, "Number of Runs: %d", &nr_runs);
        if (strstr(line, "Number of Steps: ") != NULL) sscanf(line, "Number of Steps: %d", &nr_steps);
        if (strstr(line, "Data Points: ") != NULL) sscanf(line, "Data Points: %d", &nr_count);
        if (strstr(line, "Equilibrium: ") != NULL) sscanf(line, "Equilibrium: %lf", &pequil);
    }
    fclose(input);                                      // closes the input file
}

/**
 * Prints a progress bar for the simulation to the console.
 * 
 * @param progress Progress of the simulation time. 
 */
void print_progress(double progress) {
    const int barwidth = 45;
    const int frog_width = 1; 
    int pos = (int)((barwidth - frog_width) * progress);
    printf("[");
    for (int i = 0; i < barwidth; i++) {
        if (i == pos) {
            printf("X"); 
            i += frog_width - 1; 
        } else {
            printf(i < pos ? "~" : "-");
        }
    }
    printf("] %.1f%%\r", progress * 100);
    fflush(stdout);
}

/**
 * Writes the order parameter distribution.
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_order(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024];
    int s,i;

    norm_order();

    sprintf(fname, "%s/%d/order.csv", project_name,n);

    fconfig = fopen(fname,"w");

    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d",s);
        for (i=1;i<=nr_types;i++){
            fprintf(fconfig,",%lf",p_MC[i][s]);
        }
        fprintf(fconfig,",%lf\n",p_MD[s]);
    }
    fclose(fconfig);
}

/**
 * Writes the order parameter distribution 
 * in dependence on the number of neighboring
 * CHOL to the respective lipids to order_neigh.dat. 
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_order_neigh(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024];
    int s,c,i;

    norm_order();

    sprintf(fname, "%s/%d/order_neigh.csv", project_name,n);

    fconfig = fopen(fname,"w");
    
    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d",s);
        for (i=1;i<=nr_types;i++){
            for (c=0;c<=4;c++){
                fprintf(fconfig,",%lf",p_MC_neigh[i][s][c]);
            }
        }
        fprintf(fconfig,"\n");
    }
    fclose(fconfig);
}

/**
 * Writes the order parameter distribution 
 * in dependence on the lipids distance to
 * CHOL to order_corr.dat. 
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_order_corr(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024];
    int s,c,i;

    norm_order();

    sprintf(fname, "%s/%d/order_corr.csv", project_name,n);

    fconfig = fopen(fname,"w");

    fprintf(fconfig, "S");
    for (i=1; i<=nr_types;i++){
        for (c=1;c<=3;c++){
            fprintf(fconfig, ",%d_%d",i,c);
        }
    }
    fprintf(fconfig, "\n");

    for (s=0;s<=150;s++){
        fprintf(fconfig,"%d",s);
        for (i=1;i<=nr_types;i++){
            fprintf(fconfig,",%lf,%lf,%lf",p_MC_1[i][s],p_MC_2[i][s],p_MC_3[i][s]);
        }
        fprintf(fconfig,"\n");
    }
    fclose(fconfig);
}

/**
 * Writes the energies of the simulation to energies.dat. 
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_energies(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024];
    int s,i;

    sprintf(fname, "%s/%d/energies.csv", project_name,n);

    fconfig = fopen(fname,"w");

    fprintf(fconfig, "S_min");
    for (i=0; i<=6;i++){
        fprintf(fconfig, ",#C_min_%d", i);
    }
    fprintf(fconfig, "\n");

    for (s=0;s<=150;s++){
        double p = ((double)s-50)/100;
        fprintf(fconfig,"%lf",p);
        for (i=0;i<=6;i++){
            fprintf(fconfig,",%lf",(H_dd_min_track[1][s][i]/H_dd_min_track[0][s][i]));
        }
        fprintf(fconfig,"\n");
    }

    fprintf(fconfig,"\n\n");

    fprintf(fconfig, "S_avg");
    for (i=0; i<=6;i++){
        fprintf(fconfig, ",#C_avg_%d", i);
    }
    fprintf(fconfig, "\n");
    
    for (s=0;s<=300;s++){
        double p = (((double)s/2.0)-50)/100;
        fprintf(fconfig,"%lf",p);
        for (i=0;i<=6;i++){
            fprintf(fconfig,",%lf",(H_dd_avg_track[1][s][i]/H_dd_avg_track[0][s][i]));
        }
        fprintf(fconfig,"\n");
    }

    fprintf(fconfig,"\n\n");

    fprintf(fconfig, "S_avg");
    for (i=0; i<=6;i++){
        fprintf(fconfig, ",#C_avg_%d", i);
    }
    fprintf(fconfig, "\n");

    for (s=0;s<=150;s++){
        double p = ((double)s-50)/100;
        fprintf(fconfig,"%lf",p);
        for (i=0;i<=6;i++){
            fprintf(fconfig,",%lf",(H_dd_avg_track[1][s*2][i]/H_dd_avg_track[0][s*2][i]));
        }
        fprintf(fconfig,"\n");
    }
    fclose(fconfig);
}

/**
 * Writes the final lattice structure to lattice.dat. 
 * Can be visualized with Lattice_Plot.py.
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_lattice(int n, char ufname[])
{
    FILE *fconfig;
    char fname[1024], str[100];
    int l,x,vec[4];
    double vec_chol[4];

    sprintf(fname, "%s/%d/lattice.csv", project_name,n);

    fconfig = fopen(fname, "w");
    for (vec[1]=0;vec[1]<length[1];vec[1]++){
        for (vec[2]=0;vec[2]<length[2];vec[2]++){
            x = transfer2(vec);
            fprintf(fconfig, "%d,%d,%d,%d\n", vec[1],vec[2],particlet[x],order[x]);
        }
    }

    for (l=0;l<nr_lattice;l++){
        if (cholpos[l] == 1){
            transfer1(l,vec);
            vec_chol[1] = (double)(vec[1]) + 0.5;
            vec_chol[2] = (double)(vec[2]) + 0.5;
            if (vec_chol[1] > (double)(length[1])) vec_chol[1] -= (double)(length[1]);
            if (vec_chol[2] > (double)(length[2])) vec_chol[2] -= (double)(length[2]);
            fprintf(fconfig, "%lf,%lf\n",vec_chol[1],vec_chol[2]);
        }
    }
    fclose(fconfig);
}

/**
 * Writes an overview of the most important results of the simulation to data.dat.
 * 
 * @param n Number of the respective run. 
 * @param ufname Placeholder for the name of the output file.
 */
void plot_info(int n,char ufname[])
{
    FILE *fconfig;
    char fname[1024];

    int l,i,s;
    double enertot,contacts;
    double contacts_moore,contacts_sum;

    double c_system = (c_chol)/(1+c_chol);

    sprintf(fname,"%s/%d/info.dat",project_name,n);

    fconfig = fopen(fname,"w");

    if (nr_types == 1 && c_chol == 0.0) fprintf(fconfig,"Unitary System - ");
    if (nr_types == 1 && c_chol != 0.0) fprintf(fconfig,"Binary System - ");
    if (nr_types == 2 && c_chol == 0.0) fprintf(fconfig,"Binary System - ");
    if (nr_types == 2 && c_chol != 0.0) fprintf(fconfig,"Tertiary System - ");
    if (lipid_type == 1) fprintf(fconfig,"DLiPC");
    if (lipid_type == 0) fprintf(fconfig,"DPPC");
    if (nr_types == 2) fprintf(fconfig,"/DPPC");
    if (c_chol != 0.0) fprintf(fconfig,"/CHOL");
    fprintf(fconfig,"\n \n");

    if (lat_type == 0) fprintf(fconfig,"Simulation Type: Simple MC Run\n");    
    if (lat_type == 1) fprintf(fconfig,"Simulation Type: Iterative Boltzmann Inversion\n");
    if (lat_type == 2) fprintf(fconfig,"Simulation Type: Determining CHOL's Effect on Entropy\n");
    if (int_type == 0) fprintf(fconfig,"Interaction Type: Average Order Parameter\n");
    if (int_type == 1) fprintf(fconfig,"Interaction Type: Minimal Order Parameter\n");
    if (mv_type == 0) fprintf(fconfig,"Movement Type: Random Jumps\n");
    if (mv_type == 1) fprintf(fconfig,"Movement Type: Adjacent Jumps\n");
    fprintf(fconfig,"\n");
    
    if (c_chol != 0.0) fprintf(fconfig,"CHOL-Conc.: %lf (%d Particles)\n",c_system,nr_chol);
    if (nr_types == 2){
        fprintf(fconfig,"DLiPC-Conc.: %lf (%d Particles)\n",c_lipid[1],nrt[1]);
        fprintf(fconfig,"DPPC-Conc.: %lf (%d Particles)\n",c_lipid[2],nrt[2]);
    }
    fprintf(fconfig,"Temperature: %d K\n",(int)read_temperature);

    fprintf(fconfig,"Average Order Parameter: ");
    for (i=1;i<=nr_types;i++){
        double mean = 0;
        double norm = 0;
        for (s=0;s<=150;s++){
            mean += (s-50)/100.0*p_MC[i][s];
            norm += p_MC[i][s];
        }
        for (s=0;s<=150;s++){
            p_MC[i][s] /= norm;
        }
        fprintf(fconfig," %lf",mean/norm);
    }
    fprintf(fconfig,"\n");
    enertot = ener_mem/((double)nr_steps-pequil*(double)nr_steps);
    fprintf(fconfig,"Total Energy: %lf\n",enertot);
    contacts = 0.0;
    contacts = (double)cc_contacts/((double)nr_steps-pequil*(double)nr_steps);
    fprintf(fconfig,"#CC-Contacts (Equilibrium): %lf\n",contacts);
    contacts /= (double)nr_chol*2*(((double)nr_chol-1)/((double)(nr_lattice)-1));
    fprintf(fconfig,"#CC-Contacts (Normalized): %lf\n",contacts);
}