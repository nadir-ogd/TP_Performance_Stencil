#include <iostream>
#include <assert.h>
#include <math.h>
#include <sys/time.h>


/********************************************************
 * Kernels à optimizer dans le cadre du TP optimisation
 * du module "Architecture et Modèle de calcul"
 * 
 * NOM_1: Nom Prénom Binôme1
 * EMAIL_1 : Email Binôme1
 * NOM_2: Nom Prénom Binôme2
 * EMAIL_2 : Email Binôme2
 * 
 * Compilation :
 * > g++ -O3 stencil.cxx -o stencil
 * 
 ********************************************************/

double
dml_micros()
{
        static struct timezone tz;
        static struct timeval  tv;
        gettimeofday(&tv,&tz);
        return((tv.tv_sec*1000000.0)+tv.tv_usec);
}
using namespace std;

typedef unsigned long long ui64;

#define ITER 5
#define NBTEST 3
//#define DEBUG

const ui64 order=8;
ui64 DIMX,DIMY,DIMZ;
ui64 MAXX,MAXY,MAXZ;
ui64 xyplane,MATsize;

// retourne un offset dans le centre de la matrice les dimensions sont [0..DIM-1]
inline
ui64 DIMXYZ(ui64 x,ui64 y,ui64 z){
        return((z+order)*xyplane+(y+order)*MAXX+x+order);
}

// retourne un offset dans la matrice les dimensions sont [-order..DIM+order-1] mais en indices de [0..DIM+2*order-1]
inline
ui64 MATXYZ(ui64 x,ui64 y,ui64 z){
        return(x+ y*MAXX+z*xyplane);
}

double *matA;
double *matB;
double *matC;


void init(ui64 dim){
        // l initialisation ne fait pas partie de l exercise , elle peut etre optimisee mais n est pas mesuree car elle remplie de facon artificielle les matrices
        // les donnees n influent pas sur la performance

        DIMX  = dim;
        DIMY  = dim;
        DIMZ  = dim;
        MAXX=DIMX+2*order;
        MAXY=DIMY+2*order;
        MAXZ=DIMZ+2*order;
        xyplane=MAXX*MAXY;
        MATsize=MAXX*MAXY*MAXZ;

        // dynamically allocate memory of size DIMX*DIMY*DIMZ+ghost region on 6 faces
        matA = new double[MATsize];
        assert( matA!=NULL);
        matB = new double[MATsize];
        assert( matB!=NULL);
        matC = new double[MATsize];
        assert( matC!=NULL);

        // Initialisation centre et bords
        // Les matrices A et C sont mises a zero
        // A en la matrice d emtree et C la matrice de sortie
        // La matrice B est un stencil constant pour le run
        for (ui64 z = 0; z < MAXZ; z++) {
                for (ui64 y = 0; y < MAXY; y++){
                        for (ui64 x = 0; x < MAXX; x++){
                                matA[MATXYZ(x,y,z)] = 0.0;
                                matC[MATXYZ(x,y,z)] = 0.0;
                                matB[MATXYZ(x,y,z)] = sin(z*cos(x+0.311)*cos(y+.817)+.613);
                        }
                }
        }
        // Initialisation centre de A qui est la matrice de data
        for (ui64 z = 0; z < DIMZ; z++) {
                for (ui64 y = 0; y < DIMY; y++){
                        for (ui64 x = 0; x < DIMX; x++){
                                matA[DIMXYZ(x,y,z)] = 1.0;
                        }
                }
        }

}

void one_iteration()
{
                for (ui64 z = 0; z < DIMZ; z++) {
                        for (ui64 y = 0; y < DIMY; y++){
                                for (ui64 x = 0; x < DIMX; x++){
                                        matC[DIMXYZ(x,y,z)] = matA[DIMXYZ(x,y,z)]*matB[DIMXYZ(x,y,z)] ;
                                        for (ui64 o = 1; o <= order; o++){
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x+o,y,z)]*matB[DIMXYZ(x+o,y,z)] / pow(17.0,o);
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x-o,y,z)]*matB[DIMXYZ(x-o,y,z)] / pow(17.0,o);
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x,y+o,z)]*matB[DIMXYZ(x,y+o,z)] / pow(17.0,o);
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x,y-o,z)]*matB[DIMXYZ(x,y-o,z)] / pow(17.0,o);
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x,y,z+o)]*matB[DIMXYZ(x,y,z+o)] / pow(17.0,o);
                                               matC[DIMXYZ(x,y,z)]+= matA[DIMXYZ(x,y,z-o)]*matB[DIMXYZ(x,y,z-o)] / pow(17.0,o);
                                        }
                                }
                        }
                }
                //  A=C
                for (ui64 z = 0; z < DIMZ; z++) {
                        for (ui64 y = 0; y < DIMY; y++){
                                for (ui64 x = 0; x < DIMX; x++){
                                        matA[DIMXYZ(x,y,z)] = matC[DIMXYZ(x,y,z)];
                                }
                        }
                }
}

int test_val(double *val_ref, double *mat){
        for(ui64 i=0;i<5;i++)
                if (fabs((matA[DIMXYZ(DIMX/2+i,DIMY/2+i,DIMZ/2+i)] - val_ref[i])/val_ref[i])>ldexp(1,-10))
                        return 0;
                        
        return 1; 
}

int main(const int argc,char **argv){
        /* 
         * Pensez à ne faire les tests que sur un petit cas en mode développment
         */
        ui64 val_dim[NBTEST]={50,100,200};
        ui64 perf_ref[NBTEST]={164590,1249490,10277669};
        double val_ref[3][5]={
                { 0.633645047811189 , 0.540528349585631 , 0.588436820082514 , -2.147536753829282, -0.004316384996407},
                { 1.463043183363641 , -2.359859487644307, -0.017670676323675,  0.524241548985922, -1.223974699224165},
                { -0.575619469035569, -1.061195501985062, -1.064424865572597, -0.000530250045409,  1.016904895796148}
        };

        double speedup, avg_speedup=0.;

        for(uint32_t idim=0; idim<NBTEST; idim++){
                printf("Test de performance pour la dimension %lld \n", val_dim[idim]);

                init(val_dim[idim]);

                //phase1
                double t1=dml_micros();
                for (ui64 i = 0; i < ITER; i++) {
                        // calcule 1 iteration Jacobi   C=B@A
                        one_iteration();
                }
                double t2=dml_micros();
                printf("Temps            %10.3lf\n",(t2-t1) );
                printf("Nb Cycle / point %10.3lf\n",(t2-t1)*1000.0/DIMX/DIMY/DIMZ/ITER );
                speedup = (t2-t1)/(double)perf_ref[idim];
                printf("Speedup          %10.3lf\n\n", speedup);
                avg_speedup += speedup;
                if (test_val(val_ref[idim], matA) == 0){
                        printf("ERREUR: résultat au dessus de l'erreur limite \n");
                        printf("__REF__");  
                        for(ui64 i=0;i<5;i++)printf(" %18.15lf",val_ref[idim][i]);
                        printf("\n__YOU__");
                        for(ui64 i=0;i<5;i++)printf(" %18.15lf",matA[DIMXYZ(DIMX/2+i,DIMY/2+i,DIMZ/2+i)]);
                        printf("\n");                       
                }
#ifdef DEBUG
                        printf("_0_ ");
                        for(ui64 i=0;i<5;i++)printf(" %18.15lf",matA[DIMXYZ(DIMX/2+i,DIMY/2+i,DIMZ/2+i)]);
                        double ns_point=(t2-t1)*1000.0/DIMX/DIMY/DIMZ;
                        printf("  %10.0lf  %10.3lf %lld %lld %lld\n",t2-t1,ns_point,DIMX,DIMY,DIMZ);
#endif

                delete[] matA;
                delete[] matB;
                delete[] matC;
        }
        printf("AVG Speedup      %10.3lf\n", avg_speedup/NBTEST);
 
        return 0;
}
