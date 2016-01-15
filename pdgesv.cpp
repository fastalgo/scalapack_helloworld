#include <mpi.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <cstdlib>
//#include "Scalapack.h"
#include <mkl.h>
#include <mkl_scalapack.h>
#include "mkl_lapacke.h"
#include <mkl_cblas.h>

#define mat(matriz,coluna,i,j) (matriz[i*coluna+j])

#define p_of_i(i,bs,p) ( MKL_INT((i-1)/bs)%p)
#define l_of_i(i,bs,p) ( MKL_INT((i-1)/(p*bs)))
#define x_of_i(i,bs,p) (((i-1)%bs)+1)

//#define   numroc_      NUMROC

using namespace std;

extern "C" 
{
    /* BLACS C interface */
    void Cblacs_pinfo(int* mypnum, int* nprocs);
    void Cblacs_get( MKL_INT context, MKL_INT request, MKL_INT* value);
    int  Cblacs_gridinit( MKL_INT* context, char * order, MKL_INT np_row, MKL_INT np_col);
    void Cblacs_gridinfo( MKL_INT context, MKL_INT*  np_row, MKL_INT* np_col, MKL_INT*  my_row,
    MKL_INT*  my_col);
    int  numroc_( MKL_INT *n, MKL_INT *nb, MKL_INT *iproc, MKL_INT *isrcproc, MKL_INT *nprocs);
    void Cblacs_gridexit(MKL_INT ictxt);
    void Cblacs_barrier(MKL_INT ictxt, char * order);
}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol);
int getIndex(MKL_INT row, MKL_INT col,MKL_INT NCOLS) {return row*NCOLS+col;}

/*
CTEST_Scalapack::CTEST_Scalapack(void)
{
}

CTEST_Scalapack::~CTEST_Scalapack(void)
{
}
*/

int main(int argc, char ** argv) 
{

    int nprocs = 0;//MPI::COMM_WORLD.Get_size();
    int rank = 0;//MPI::COMM_WORLD.Get_rank();

    MPI_Init(&argc,&argv);    
    MPI_Comm_size(MPI_COMM_WORLD,&nprocs);
    MPI_Comm_rank(MPI_COMM_WORLD,&rank);

    std::cout<<"Returned: "<<" ";
    std::cout << "Hello World! I am " << rank << " of " << nprocs <<
    std::endl;

    srand(1);
    MKL_INT myrow=0;
    MKL_INT mycol=0;
    MKL_INT ictxt=0;
    MKL_INT nprow=0,npcol=0;

    MKL_INT BLOCK_SIZE =2; //this gonna be tricky - should be 64, but cannot be larger than the original matrix

    MKL_INT locR=0, locC=0;
    MKL_INT block = BLOCK_SIZE;
    MKL_INT izero = 0;
    MKL_INT matrix_size = 9;
   
    MKL_INT myone = 1;
    
    MKL_INT nrhs = 1;
   
    MKL_INT info=0;
  
    int i=0,j=0;
    double mone=(-1.e0),pone=(1.e0);
    double AnormF=0.e0, XnormF=0.e0, RnormF=0.e0, BnormF=0.e0, residF=0.e0,eps=0.e0;

    find_nps(nprocs,nprow,npcol);

    Cblacs_pinfo( &rank, &nprocs ) ;
    Cblacs_get(-1, 0, &ictxt);
    Cblacs_gridinit(&ictxt, "Row", nprow, npcol);
    Cblacs_gridinfo(ictxt, &nprow, &npcol, &myrow, &mycol);
    
    locR = numroc_(&matrix_size, &block, &myrow, &izero, &nprow);
    locC = numroc_(&matrix_size, &block, &mycol, &izero, &npcol);

   
    ////GLOBAL
    double * A = new double[matrix_size*matrix_size]();
    double * B = new double[matrix_size]();
    double * Acpy = new double[matrix_size*matrix_size]();
    double * Bcpy = new double[matrix_size]();
    
    //LOCAL
    double * local_know_vector = new double[locR]();
    double * local_matrix = new double[locR*locC]();
    
    MKL_INT* ipiv = new MKL_INT [locC*locR*block+1000000]();

    
    B[2] = 1;
    B[3] = 0;
    B[4] = 0;
    B[5] = 0;
    
    
    
    A[0] = 19;
    A[1] = 3;
    A[2] = 1;
    A[3] = 12;
    A[4] = 1;
    A[5] = 16;
    A[6] = 1;
    A[7] = 3;
    A[8] = 11;
    
    A[9] = -19;
    A[10] = 3;
    A[11] = 1;
    A[12] = 12;
    A[13] = 1;
    A[14] = 16;
    A[15] = 1;
    A[16] = 3;
    A[17] = 11;
    
    A[18] = -19;
    A[19] = -3;
    A[20] = 1;
    A[21] = 12;
    A[22] = 1;
    A[23] = 16;
    A[24] = 1;
    A[25] = 3;
    A[26] = 11;
    
    A[27] = -19;
    A[28] = -3;
    A[29] = -1;
    A[30] = 12;
    A[31] = 1;
    A[32] = 16;
    A[33] = 1;
    A[34] = 3;
    A[35] = 11;
    
    A[36] = -19;
    A[37] = -3;
    A[38] = -1;
    A[39] = -12;
    A[40] = 1;
    A[41] = 16;
    A[42] = 1;
    A[43] = 3;
    A[44] = 11;
    
    A[45] = -19;
    A[46] = -3;
    A[47] = -1;
    A[48] = -12;
    A[49] = -1;
    A[50] = 16;
    A[51] = 1;
    A[52] = 3;
    A[53] = 11;
    
    A[54] = -19;
    A[55] = -3;
    A[56] = -1;
    A[57] = -12;
    A[58] = -1;
    A[59] = -16;
    A[60] = 1;
    A[61] = 3;
    A[62] = 11;
    
    A[63] = -19;
    A[64] = -3;
    A[65] = -1;
    A[66] = -12;
    A[67] = -1;
    A[68] = -16;
    A[69] = -1;
    A[70] = 3;
    A[71] = 11;
    
    A[72] = -19;
    A[73] = -3;
    A[74] = -1;
    A[75] = -12;
    A[76] = -1;
    A[77] = -16;
    A[78] = -1;
    A[79] = -3;
    A[80] = 11;

    MKL_INT* descA  = new MKL_INT[9]();
    MKL_INT* descB  = new MKL_INT[9]();
   
    descA[0] = 1; // descriptor type
    descA[1] = ictxt; // blacs context
    descA[2] = matrix_size; // global number of rows
    descA[3] = matrix_size; // global number of columns
    descA[4] = block; // row block size
    descA[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descA[6] = 0; // initial process row(DEFINED 0)
    descA[7] = 0; // initial process column (DEFINED 0)
    descA[8] = locR; // leading dimension of local array

    descB[0] = 1; // descriptor type
    descB[1] = ictxt; // blacs context
    descB[2] = matrix_size; // global number of rows
    descB[3] = 1; // global number of columns
    descB[4] = block; // row block size
    descB[5] = block; // column block size (DEFINED EQUAL THAN ROW BLOCK SIZE)
    descB[6] = 0; // initial process row(DEFINED 0)
    descB[7] = 0; // initial process column (DEFINED 0)
    descB[8] = locR; // leading dimension of local array

    int il=0, jl=0;
    for(i=1; i< matrix_size+1; i++) 
    {
       for(j=1; j< matrix_size+1; j++) 
       {
    
        int pi = p_of_i(i,block,nprow);
        
        int li = l_of_i(i,block,nprow);

        int xi = x_of_i(i,block,nprow);
        //printf("i = %d, j = %d, pi = %d, li = %d\n",i,j,pi,li);;fflush(stdout);
        int pj = p_of_i(j,block,npcol);
        
        int lj = l_of_i(j,block,npcol);
        
        int xj = x_of_i(j,block,npcol);
        //printf("i = %d, j = %d, pj = %d, lj = %d, xj = %d\n",i,j,pj,lj,xj);;fflush(stdout);

        if( (pi == myrow) && (pj == mycol)) 
        {
            il = li*block+xi;
            jl = lj*block+xj;
            local_matrix[getIndex(il-1, jl-1, locC)] = A[getIndex(i-1,j-1,matrix_size)];
        }
    
        if(  (pi == myrow) &&(mycol==0)  )
        {
            local_know_vector[il-1] = B[i-1];
        }

       }
    
    }
      
    ////STARTING PDGESV
    pdgesv_(&matrix_size, &nrhs, local_matrix, &myone, &myone, descA, ipiv, local_know_vector, &myone, &myone, descB, &info);
    
    if(rank==0)
      {
        if(info != 0) cout <<"PDGESV problem! Info "<<info<<endl;
      }
    
    
    for(i=0; i< locR; i++)
    {
      cout<<"**\n"<<"rank "<<rank<<"  answer: "<<local_know_vector[i]<<endl;
    }

    if(NULL!=descA)                        {delete [] descA; descA=NULL;} 
    if(NULL!=descB)                        {delete [] descB; descB=NULL;} 
    if(NULL!=local_know_vector)            {delete [] local_know_vector; local_know_vector=NULL;} 
    if(NULL!=local_matrix)                {delete [] local_matrix; local_matrix=NULL;} 
    if(NULL!=Acpy)                        {delete [] Acpy; Acpy=NULL;} 
    if(NULL!=Bcpy)                        {delete [] Bcpy; Bcpy=NULL;} 
    if(NULL!=A)                            {delete [] A; A=NULL;} 
    if(NULL!=B)                            {delete [] B; B=NULL;} 
    

    Cblacs_gridexit(ictxt);

    return 0;

}

void find_nps(MKL_INT np, MKL_INT &nprow, MKL_INT & npcol) 
{
#if 1
    nprow = (int)sqrt( np );
    npcol = np / nprow;
    return;

#else

MKL_INT min_nprow=100000;
MKL_INT min_npcol=100000;

nprow = np;
npcol = np;

while(1) {

   npcol--;
  if(np%2==0   ) {
  if(npcol ==1){
   nprow --;
   npcol = nprow;
  }
  }else {
  if(npcol ==0){
   nprow --;
   npcol = nprow;
  }

  }

  if(nprow*npcol == np) {
    min_npcol = npcol;
    if(nprow < min_nprow)    min_nprow = nprow;
  }

    if(nprow ==1 ) break;

}

nprow = min_nprow;
npcol = min_npcol;

#endif
}
