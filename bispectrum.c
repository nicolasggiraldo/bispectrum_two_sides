#include           <stdio.h>
#include         <stdlib.h>
#include           <math.h>
#include          <fftw3.h>
#include         <string.h>
#include    <gsl/gsl_rng.h>
#include <gsl/gsl_spline.h>
#include           <time.h>
#include         <unistd.h>
#include         <limits.h>
#include            <mpi.h>

#include "constansts_and_structures.h"
#include                 "functions.h"

/*
  FFTW_MEASURE: find the optimal plan by actually computing several FFTs  
  FFTW_ESTIMATE: do not run any FFT and provide a "reasonable" plan
  FFTW_OUT_OF_PLACE: a plan assumes that the in and out are distinct 
  FFTW_IN_PLACE: a plan assumes that the in and out are same
*/

int main(int argc, char *argv[]){
  int i, j, k, l;                  // Counter in the X, Y, and Z axis.
  unsigned long int rand_i;        /* Index for the subsample of density 
				      contrast in the range between                
				      k1-DELTA_K/2 to k1+DELTA_K/2 */
  unsigned long int rand_j;
  //unsigned long int rand_k;
  long int id_cell;                // Counter for the id of the cell
  double kMag;                     /* Magnitude of the wavenumber vector */
  FILE *fout=NULL;                 // File handler to output
  double *kpos;                    /* Array of positions values according 
				      to FFTW k-position convention */ 
  int *indexpos=NULL;              /* Array of index values according 
				      to FFTW k-position convention */
  int indexcord[2];
  int m3[3];
  int rank, size;
  int Nbins;
  int *taskBin=NULL;
  double (*W_k_Re)(double) = NULL; // Addres to the window function
  double (*W_k_Im)(double) = NULL; // Addres to the window function
  double (*W2_k)(double) = NULL; // Addres to the window function
  struct densityContrast *q1=NULL; /* Array of structure with subset of 
				      density contrast with range between 
				      k-DELTA_K/2 to k+DELTA_K/2 */
  struct densityContrast *q2=NULL; 
  struct binStruct *bindata=NULL;  /* Data structure with information of
				      the measures in the taken bins */
  
  double denConCor_i[2],denConCor_j[2],denConCor_k[2];
  double WRe, WIm, WMag2;
  double D20_complexki0,D20_complexkj0,D20_complexkk0;
  double D20_complexki1,D20_complexkj1,D20_complexkk1;
    
  MPI_Status status;
  int chunk;

  /*
  gsl_rng *r=NULL;
  long seed;
  const gsl_rng_type *T;
  seed = time(NULL) * getpid();
  T = gsl_rng_gfsr4;
  r = gsl_rng_alloc(T);
  gsl_rng_set(r, seed);
  //*/
  
  
  
  /////////////////////
  //* MPI BEGGINING *//
  /////////////////////
  MPI_Init(&argc, &argv);
  
  // Number of processors
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  // Task numbering
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  
  
  //////////////////////////////
  //* READING PARAMETER FILE *//
  //////////////////////////////
  if(rank == 0){
    if(argc<2){
      printf("\n***********************************");
      printf("***********************************\n");
      printf("%s: You must specify the name of the parameter file\n",argv[0]);
      printf("For example: %s pathOfFile/parameterFile.txt\n",argv[0]);
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  

  // Reading parameter file and verifying there is no error.
  switch( read_parameters(argv[1], rank) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad path (or name) to the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: Bad settings in the parameter file.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////
  //* READING CELL BINARY FILE *//
  ////////////////////////////////
  switch( readBinaryFile(rank) ){
  case -1 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: The parameter file could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  case -2 :
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: FFTW arrays could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////
  //* EXECUTING FFTW3 ROUTINES *//
  ////////////////////////////////
  if(rank == 0){
    /* Do forward FFT */
    fftw_execute(forwardPlan);
    printf("\n-----------------------------------------------\n");
    printf("Fourier Transform succes.\n");
    fftw_free(denConX);
    fftw_destroy_plan(forwardPlan);
  }
  
  /* Waiting the processor 0 to do the FFT */
  MPI_Barrier(MPI_COMM_WORLD);
  
  
  
  ///////////////////////////////////////////////////
  //* SENDING ALL DELTAK GRID TO OTHER PROCESSORS *//
  ///////////////////////////////////////////////////
  chunk = (GV.NGRID3)/sizeof(fftw_complex);
  for(l=0; l<sizeof(fftw_complex); l++){
    MPI_Bcast( denConK+(l*chunk),
	       chunk*sizeof(fftw_complex),
	       MPI_BYTE,
	       0,
	       MPI_COMM_WORLD );
  }
  
  
  
  ///////////////////////////////////////
  //* SETTING FFTW3 COORDINATE SYSTEM *//
  ///////////////////////////////////////
  
  /* Position array for storing in the densityContrast */
  kpos = (double *) calloc(GV.NGRID, sizeof(double));
  if(kpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: kpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  indexpos = (int *) calloc(GV.NGRID, sizeof(int));
  if(indexpos == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: indexpos array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  /* Setting index and positions according to FFTW convention  */
  for( i=0; i<GV.NGRID; i++ ){
    /* REMEMBER kF is (2.0*PI)/L */
    kpos[i]     = (i<GV.NGRID/2) ? GV.KF*i : GV.KF*(i-GV.NGRID);
    indexpos[i] = (i<GV.NGRID/2) ?       i :       (i-GV.NGRID);
  }// for i

  /* These are the minimum and maximum integer vectors 
     possibles according to the grid taken  */
  indexcord[MIN] = -(GV.NGRID/2);
  indexcord[MAX] =  (GV.NGRID/2);
  
  
  
  ///////////////////////////////////////////
  //* DEFINING THE WINDOW FUNCTION TO USE *//
  ///////////////////////////////////////////
  if( strcmp(GV.SCHEME, "NGP") == 0 ){
    W_k_Re = W_NGP; //NGP
    W_k_Im = zero;
    W2_k = Sum_W2_NGP;
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){ 
    W_k_Re = W_CIC; //CIC
    W_k_Im = zero;
    W2_k = Sum_W2_CIC;
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    W_k_Re = W_TSC; // TSC
    W_k_Im = zero;
    W2_k = Sum_W2_TSC;
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){
    
    len_array_D20 = 399;
    int err;
    
    k_D20     = (double *) calloc( len_array_D20, sizeof(double));
    Re_D20    = (double *) calloc( len_array_D20, sizeof(double));
    Im_D20    = (double *) calloc( len_array_D20, sizeof(double));
    Kmag2_D20 = (double *) calloc( len_array_D20, sizeof(double));

    fin_D20 = fopen("./D20.txt", "r");

    for(i=0; i<len_array_D20; i++){
      err=fscanf( fin_D20,
		  "%lf %lf %lf %lf", 
		  &k_D20[i],
		  &Re_D20[i],
		  &Im_D20[i],
		  &Kmag2_D20[i] );
    }

    fclose(fin_D20);
    
    // GSL interpolation allocation
    acc = gsl_interp_accel_alloc(); // accelerator

    // spline interpolation
    splineRe   = gsl_spline_alloc( gsl_interp_cspline,
				   len_array_D20); 
    splineIm   = gsl_spline_alloc( gsl_interp_cspline,
				   len_array_D20);
    splineMag2 = gsl_spline_alloc( gsl_interp_cspline,
				   len_array_D20);
    
    // GSL init
    gsl_spline_init(splineRe,   k_D20, Re_D20,    len_array_D20);
    gsl_spline_init(splineIm,   k_D20, Im_D20,    len_array_D20);
    gsl_spline_init(splineMag2, k_D20, Kmag2_D20, len_array_D20);
    
    if(err){}
    
    W_k_Re = W_D20_Re; // D20
    W_k_Im = W_D20_Im;
    W2_k = Sum_W2_D20;
  }
  
  
  
  ////////////////////////////////////
  //* DECONVOLVING WINDOW FUNCTION *//
  ////////////////////////////////////

  /*
    if( strcmp(GV.SCHEME, "D20") == 0 ){
    double D20_complexi0, D20_complexj0,D20_complexk0;
    double D20_complexi1, D20_complexj1,D20_complexk1;
    double WRe, WIm, WMag2;
    double dRe, dIm;
    for(i=0; i<GV.NGRID; i++){
    
    //D20_complexi0=gsl_spline_eval(splineRe, kpos[i]/GV.KN, acc);
    //D20_complexi1=gsl_spline_eval(splineIm, kpos[i]/GV.KN, acc);
    D20_complexi0=gsl_spline_eval(splineRe, (1.0*indexpos[i])/(GV.NGRID/2), acc);
    D20_complexi1=gsl_spline_eval(splineIm, (1.0*indexpos[i])/(GV.NGRID/2), acc);
    //indexpos[i]
    
    for(j=0; j<GV.NGRID; j++){
    
    //D20_complexj0=gsl_spline_eval(splineRe, kpos[j]/GV.KN, acc);
    //D20_complexj1=gsl_spline_eval(splineIm, kpos[j]/GV.KN, acc);
    D20_complexj0=gsl_spline_eval(splineRe, (1.0*indexpos[j])/(GV.NGRID/2), acc);
    D20_complexj1=gsl_spline_eval(splineIm, (1.0*indexpos[j])/(GV.NGRID/2), acc);
    
    for(k=0; k<GV.NGRID; k++){
    
    D20_complexk0=gsl_spline_eval(splineRe, (1.0*indexpos[k])/(GV.NGRID/2), acc);
    D20_complexk1=gsl_spline_eval(splineIm, (1.0*indexpos[k])/(GV.NGRID/2), acc);
    
    WRe = ( +D20_complexi0 * D20_complexj0 * D20_complexk0
    -D20_complexi0 * D20_complexj1 * D20_complexk1
    -D20_complexi1 * D20_complexj0 * D20_complexk1
    -D20_complexi1 * D20_complexj1 * D20_complexk0);
    
    WIm = ( +D20_complexi0 * D20_complexj0 * D20_complexk1
    +D20_complexi0 * D20_complexj1 * D20_complexk0
    +D20_complexi1 * D20_complexj0 * D20_complexk0
    -D20_complexi1 * D20_complexj1 * D20_complexk1);
    
    //WRe=1.0;
    //WIm=0.0;
    
    WMag2 = (WRe*WRe) + (WIm*WIm);
    
    id_cell = INDEX(i,j,k);
    
    dRe = ( denConK[id_cell][0]*WRe + denConK[id_cell][1]*WIm)/WMag2;
    dIm = (-denConK[id_cell][0]*WIm + denConK[id_cell][1]*WRe)/WMag2;
    
    denConK[id_cell][0] = dRe;
    denConK[id_cell][1] = dIm;
    
    }// for k
    }// for j
    }// for i
    }
    else{ // NGP, CIC, TSC
    
    for(i=0; i<GV.NGRID; i++){
    for(j=0; j<GV.NGRID; j++){
    for(k=0; k<GV.NGRID; k++){
    
    id_cell = INDEX(i,j,k);
    
    denConK[id_cell][0] /= ( W_k(kpos[i])*W_k(kpos[j])*W_k(kpos[k]) );
    denConK[id_cell][1] /= ( W_k(kpos[i])*W_k(kpos[j])*W_k(kpos[k]) );
    
    }// for k
    }// for j
    }// for i
    
    }
  */
    
    
  ////////////////////////////////////
  //* GETTING SET K FOR THE K1 LEG *//
  ////////////////////////////////////
  q1 = (struct densityContrast *) calloc( floor(3.0*M_PI*GV.NGRID2*GV.S_KF*0.5), 
					  sizeof(struct densityContrast) );
  if(q1 == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: q1 array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ////////////////////////////////////
  //* GETTING SET K FOR THE K2 LEG *//
  ////////////////////////////////////
  q2 = (struct densityContrast *) calloc( floor(3.0*M_PI*GV.NGRID2*GV.S_KF*0.5), 
					  sizeof(struct densityContrast) );
  if(q2 == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: q2 array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  
  
  
  ///////////////////////////
  // BISPECTRUM ESTIMATION //
  ///////////////////////////
  
  Nbins = (int) ceil( GV.KN / GV.DELTA_K );
  
  if(rank == 0){
    printf("\n-----------------------------------------------\n");
    printf(" Domain decomposition... \n");
    printf(" %d Total bins\n", Nbins); fflush(stdout);
  }
  
  /* domain decomposition */
  // building bins array
  // declara arreglo que contenga la info de los dk
  // otro inicializado a cero que va a contener el valor de Bk
  bindata = (struct binStruct *) calloc( Nbins, sizeof(struct binStruct) );
  if(bindata == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: bindata array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }

  //Numero de bines por task = Numero total de bine/size;
  //binspertask = Nbins/size;
  // dar la lista de bines que va a calcular 
  
  taskBin = (int *) calloc( Nbins, sizeof(int) );
  if(taskBin == NULL){
    printf("\n***********************************");
    printf("***********************************\n");
    printf("Error: taskBin array could not be allocated.\n" );
    printf("***********************************");
    printf("***********************************\n\n");
    exit(0);
  }
  

  i=0;
  for(l=Nbins-1; l>=0; l--){
    taskBin[l] = i%size;
    
    bindata[l].k3 = (l+0.5)*GV.DELTA_K;
    bindata[l].Nk3 = 0;
    bindata[l].sumDelta2_k3 = 0.0;
    bindata[l].Pk3 = 0.0;         
    bindata[l].Pk3_Error = 0.0;
    
    bindata[l].Ntri = 0;       
    bindata[l].I_delta3 = 0.0; 
    bindata[l].Bk = 0.0;
    bindata[l].Bk_Error = 0.0; 
    bindata[l].Qk = 0.0;       
    bindata[l].Qk_Error = 0.0;
    bindata[l].Bk_shotnoise = 0.0;
    
    if(rank == 0){
      printf("rank %d has bin %d (k=%f)\n", 
	     taskBin[l], l, bindata[l].k3); fflush(stdout);
    }
    
    i++;
  }
  
  // Getting set q1 associated to k1
  for(i=0; i<GV.NGRID; i++){
    for(j=0; j<GV.NGRID; j++){
      for(k=1; k<(GV.NGRID/2+1); k++){
	
	kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	
	if( ( GV.K1-GV.DELTA_K*0.5 < kMag ) && ( kMag < GV.K1+GV.DELTA_K*0.5 ) ){
	  
	  id_cell = INDEX(i,j,k);
	  
	  q1[GV.Nk1].id   = id_cell;
	  q1[GV.Nk1].kMag = kMag;
	  
	  q1[GV.Nk1].triplex[X] = indexpos[i];
	  q1[GV.Nk1].triplex[Y] = indexpos[j];
	  q1[GV.Nk1].triplex[Z] = indexpos[k];
	  
	  GV.Nk1++;
	  
	  GV.sumDelta2_k1 += COMPLEXMAG(denConK, id_cell)/( W2_k(kpos[i]) *
							    W2_k(kpos[j]) *
							    W2_k(kpos[k]) );
	  
	}// if bindata[l].k1-GV.DELTA_K*0.5 < kMag < bindata[l].k1+GV.DELTA_K*0.5 )
	
      }// for k
    }// for j
  }// for i
    
  GV.Pk1  = GV.sumDelta2_k1 / GV.Nk1;
  GV.Pk1 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
  GV.Pk1 -= GV.SHOT_NOISE;
  GV.Pk1_Error = (GV.Pk1*(GV.DELTA_K/GV.K1))/sqrt(2.0*M_PI);


  // Getting sets associated to k2
  for(i=0; i<GV.NGRID; i++){
    for(j=0; j<GV.NGRID; j++){
      for(k=1; k<(GV.NGRID/2+1); k++){
	
	kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	
	if( ( GV.K2-GV.DELTA_K*0.5 < kMag ) && ( kMag < GV.K2+GV.DELTA_K*0.5 ) ){
	  
	  id_cell = INDEX(i,j,k);
	  
	  q2[GV.Nk2].id   = id_cell;
	  q2[GV.Nk2].kMag = kMag;
	  
	  q2[GV.Nk2].triplex[X] = indexpos[i];
	  q2[GV.Nk2].triplex[Y] = indexpos[j];
	  q2[GV.Nk2].triplex[Z] = indexpos[k];
	  
	  GV.Nk2++;
	  
	  GV.sumDelta2_k2 += COMPLEXMAG(denConK, id_cell)/( W2_k(kpos[i]) *
							    W2_k(kpos[j]) *
							    W2_k(kpos[k]) );
	  
	}// if bindata[l].k2-GV.DELTA_K*0.5 < kMag < bindata[l].k2+GV.DELTA_K*0.5
	
      }// for k
    }// for j
  }// for i
    
  GV.Pk2  = GV.sumDelta2_k2 / GV.Nk2;
  GV.Pk2 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
  GV.Pk2 -= GV.SHOT_NOISE;
  GV.Pk2_Error = (GV.Pk2*(GV.DELTA_K/GV.K2))/sqrt(2.0*M_PI);


  
  if(rank==0){
    printf("\n-----------------------------------------------\n");
    printf("Bispectrum calculation");
    printf("\n-----------------------------------------------\n");
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(l=0; l<Nbins; l++){
    if(rank != taskBin[l])
      continue;

    if( fabs(GV.K1-GV.K2)-(1.5*GV.DELTA_K)<=bindata[l].k3 &&  bindata[l].k3<=(GV.K1+GV.K2)+(1.5*GV.DELTA_K) )
      {
	printf("rank:%3d, k3 = %lf\t skipped!\n", rank, bindata[l].k3);
	fflush(stdout);
	continue;
      }
    
    printf("rank:%3d, k3 = %lf\n", rank, bindata[l].k3);
    fflush(stdout);
    
    // Getting sets associated to k3
    for(i=0; i<GV.NGRID; i++){
      for(j=0; j<GV.NGRID; j++){
	for(k=1; k<(GV.NGRID/2+1); k++){
	
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	
	  if( ( bindata[l].k3-GV.DELTA_K*0.5 < kMag ) && 
	      ( kMag < bindata[l].k3+GV.DELTA_K*0.5 ) ){
	    
	    id_cell = INDEX(i,j,k);
	    
	    bindata[l].Nk3++;
	    
	    bindata[l].sumDelta2_k3 += COMPLEXMAG(denConK, id_cell)/( W2_k(kpos[i]) *
								      W2_k(kpos[j]) *
								      W2_k(kpos[k]) );
	  
	  }// if bindata[l].k2-GV.DELTA_K*0.5 < kMag < bindata[l].k2+GV.DELTA_K*0.5 )
	  
	}// for k
      }// for j
    }// for i
    
    bindata[l].Pk3  = bindata[l].sumDelta2_k3 / bindata[l].Nk3;
    bindata[l].Pk3 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
    bindata[l].Pk3 -= GV.SHOT_NOISE;
    bindata[l].Pk3_Error = (bindata[l].Pk3*(GV.DELTA_K/bindata[l].k3))/sqrt(2.0*M_PI);
    
    
    /* Taking all possible triangle configurations */
    for( rand_i=0; rand_i<GV.Nk1; rand_i++){
      for( rand_j=0; rand_j<GV.Nk2; rand_j++){
	
	/*
	  do{
	  rand_i = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_j = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_k = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  }while( (rand_i==rand_j) || (rand_i==rand_k) || (rand_j==rand_k) );
	*/
	
	/* If both sides have the same id then continue */
	if(q1[rand_i].id == q2[rand_j].id)
	  continue;
	
	m3[X] = - q1[rand_i].triplex[X] - q2[rand_j].triplex[X];
	m3[Y] = - q1[rand_i].triplex[Y] - q2[rand_j].triplex[Y];
	m3[Z] = - q1[rand_i].triplex[Z] - q2[rand_j].triplex[Z];
	
	//if( indexcord[MIN]<=m3[X] && m3[X]<=indexcord[MAX] && 
	//indexcord[MIN]<=m3[Y] && m3[Y]<=indexcord[MAX] && 
	//indexcord[MIN]<=m3[Z] && m3[Z]<=indexcord[MAX] ){
	if( fabs(m3[X])<=indexcord[MAX] && 
	    fabs(m3[Y])<=indexcord[MAX] && 
	    fabs(m3[Z])<=indexcord[MAX] ){

	  kMag = VECTORMAG(m3[X]*GV.KF,m3[Y]*GV.KF,m3[Z]*GV.KF);
	  
	  //i = (m3[X]>=0) ? m3[X] : GV.NGRID+m3[X];
	  //j = (m3[Y]>=0) ? m3[Y] : GV.NGRID+m3[Y];
	  //k = (m3[Z]>=0) ? m3[Z] : GV.NGRID+m3[Z];
	  
	  //kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	  
	  if( ( bindata[l].k3-GV.DELTA_K*0.5 < kMag ) && 
	      ( kMag < bindata[l].k3+GV.DELTA_K*0.5 ) ){

	    if(m3[Z]==indexcord[MAX]){
	      m3[X] *= -1;
	      m3[Y] *= -1;
	      m3[Z] *= -1;
	    }
	    i = (m3[X]>=0) ? m3[X] : GV.NGRID+m3[X];
	    j = (m3[Y]>=0) ? m3[Y] : GV.NGRID+m3[Y];
	    k = (m3[Z]>=0) ? m3[Z] : GV.NGRID+m3[Z];
	    
	    id_cell = INDEX(i,j,k);

	    if( q1[rand_i].id==id_cell || q2[rand_j].id==id_cell)
	      continue;

	    /*
	      bindata[l].I_delta3 += 
	      (+denConK[q1[rand_i].id][0]*denConK[q2[rand_j].id][0]*denConK[id_cell][0]
	      -denConK[q1[rand_i].id][0]*denConK[q2[rand_j].id][1]*denConK[id_cell][1]
	      -denConK[q1[rand_i].id][1]*denConK[q2[rand_j].id][0]*denConK[id_cell][1]
	      -denConK[q1[rand_i].id][1]*denConK[q2[rand_j].id][1]*denConK[id_cell][0] );
	    */

	    if( strcmp(GV.SCHEME, "D20") == 0 ){ // D20
	      
	      D20_complexki0 = W_k_Re(q1[rand_i].triplex[X]*GV.KF);
	      D20_complexkj0 = W_k_Re(q1[rand_i].triplex[Y]*GV.KF);
	      D20_complexkk0 = W_k_Re(q1[rand_i].triplex[Z]*GV.KF);

	      D20_complexki1 = W_k_Im(q1[rand_i].triplex[X]*GV.KF);
	      D20_complexkj1 = W_k_Im(q1[rand_i].triplex[Y]*GV.KF);
	      D20_complexkk1 = W_k_Im(q1[rand_i].triplex[Z]*GV.KF);

	      WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki0 * D20_complexkj1 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj0 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
	      
	      WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
		      +D20_complexki0 * D20_complexkj1 * D20_complexkk0
		      +D20_complexki1 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk1);

	      WMag2 = (WRe*WRe) + (WIm*WIm);
	      
	      denConCor_i[0] = (denConK[q1[rand_i].id][0]*WRe + denConK[q1[rand_i].id][1]*WIm)/WMag2;
	      denConCor_i[1] = (denConK[q1[rand_i].id][1]*WRe - denConK[q1[rand_i].id][0]*WIm)/WMag2;

	      //////////////////////////////////////////////////////////////////////
	      
	      D20_complexki0 = W_k_Re(q2[rand_j].triplex[X]*GV.KF);
	      D20_complexkj0 = W_k_Re(q2[rand_j].triplex[Y]*GV.KF);
	      D20_complexkk0 = W_k_Re(q2[rand_j].triplex[Z]*GV.KF);

	      D20_complexki1 = W_k_Im(q2[rand_j].triplex[X]*GV.KF);
	      D20_complexkj1 = W_k_Im(q2[rand_j].triplex[Y]*GV.KF);
	      D20_complexkk1 = W_k_Im(q2[rand_j].triplex[Z]*GV.KF);

	      WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki0 * D20_complexkj1 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj0 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
	      
	      WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
		      +D20_complexki0 * D20_complexkj1 * D20_complexkk0
		      +D20_complexki1 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk1);

	      WMag2 = (WRe*WRe) + (WIm*WIm);
	      
	      denConCor_j[0] = (denConK[q2[rand_j].id][0]*WRe + denConK[q2[rand_j].id][1]*WIm)/WMag2;
	      denConCor_j[1] = (denConK[q2[rand_j].id][1]*WRe - denConK[q2[rand_j].id][0]*WIm)/WMag2;

	      //////////////////////////////////////////////////////////////////////

	      D20_complexki0 = W_k_Re(m3[X]*GV.KF);
	      D20_complexkj0 = W_k_Re(m3[Y]*GV.KF);
	      D20_complexkk0 = W_k_Re(m3[Z]*GV.KF);

	      D20_complexki1 = W_k_Im(m3[X]*GV.KF);
	      D20_complexkj1 = W_k_Im(m3[Y]*GV.KF);
	      D20_complexkk1 = W_k_Im(m3[Z]*GV.KF);

	      WRe = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki0 * D20_complexkj1 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj0 * D20_complexkk1
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk0);
	      
	      WIm = ( +D20_complexki0 * D20_complexkj0 * D20_complexkk1
		      +D20_complexki0 * D20_complexkj1 * D20_complexkk0
		      +D20_complexki1 * D20_complexkj0 * D20_complexkk0
		      -D20_complexki1 * D20_complexkj1 * D20_complexkk1);

	      WMag2 = (WRe*WRe) + (WIm*WIm);
	      	      
	      denConCor_k[0] = (denConK[id_cell][0]*WRe + denConK[id_cell][1]*WIm)/WMag2;
	      denConCor_k[1] = (denConK[id_cell][1]*WRe - denConK[id_cell][0]*WIm)/WMag2;
	      
	    }else{ // NGP, CIC, TSC
	      
	      denConCor_i[0] = denConK[q1[rand_i].id][0]/( W_k_Re(q1[rand_i].triplex[X]*GV.KF) *
							   W_k_Re(q1[rand_i].triplex[Y]*GV.KF) *
							   W_k_Re(q1[rand_i].triplex[Z]*GV.KF) );
	      denConCor_i[1] = denConK[q1[rand_i].id][1]/( W_k_Re(q1[rand_i].triplex[X]*GV.KF) *
							   W_k_Re(q1[rand_i].triplex[Y]*GV.KF) *
							   W_k_Re(q1[rand_i].triplex[Z]*GV.KF) );
	      
	      denConCor_j[0] = denConK[q2[rand_j].id][0]/( W_k_Re(q2[rand_j].triplex[X]*GV.KF) *
							   W_k_Re(q2[rand_j].triplex[Y]*GV.KF) *
							   W_k_Re(q2[rand_j].triplex[Z]*GV.KF) );
	      denConCor_j[1] = denConK[q2[rand_j].id][1]/( W_k_Re(q2[rand_j].triplex[X]*GV.KF) *
							   W_k_Re(q2[rand_j].triplex[Y]*GV.KF) *
							   W_k_Re(q2[rand_j].triplex[Z]*GV.KF) );
	      
	      denConCor_k[0] = denConK[id_cell][0]/( W_k_Re(m3[X]*GV.KF) *
						     W_k_Re(m3[Y]*GV.KF) *
						     W_k_Re(m3[Z]*GV.KF) );
	      denConCor_k[1] = denConK[id_cell][1]/( W_k_Re(m3[X]*GV.KF) *
						     W_k_Re(m3[Y]*GV.KF) *
						     W_k_Re(m3[Z]*GV.KF) );
	      
	    }

	    bindata[l].I_delta3 += (+denConCor_i[0]*denConCor_j[0]*denConCor_k[0]
				    -denConCor_i[0]*denConCor_j[1]*denConCor_k[1]
				    -denConCor_i[1]*denConCor_j[0]*denConCor_k[1]
				    -denConCor_i[1]*denConCor_j[1]*denConCor_k[0] );
	    	        
	    bindata[l].Ntri++;
	      
	  }// if bindata[l].k1-GV.DELTA_K*0.5 < kMag < bindata[l].k1+GV.DELTA_K*0.5

	}// if

      }// for rand_j
    }// for rand_i

    if(bindata[l].Ntri==0)
      continue;

    // Discrete bispectrum value
    bindata[l].Bk  = bindata[l].I_delta3 / (1.0*bindata[l].Ntri);
    bindata[l].Bk *= (GV.SIM_VOL/(1.0*GV.NGRID3)); 
    bindata[l].Bk *= (GV.SIM_VOL/(1.0*GV.NGRID3));
    bindata[l].Bk *= (1.0/(1.0*GV.NGRID3));
    
    // Stimating shot noise for bispectrum
    bindata[l].Bk_shotnoise = GV.SHOT_NOISE * (GV.Pk1+GV.Pk2+bindata[l].Pk3) + POW2(GV.SHOT_NOISE);
    
    /* Substracting shotnoise term  */
    bindata[l].Bk -= bindata[l].Bk_shotnoise;

    /* Bispectrum error */
    bindata[l].Bk_Error  = sqrt(M_PI/(GV.K1*GV.K2*bindata[l].k3*POW3(GV.S_KF)));
    bindata[l].Bk_Error *= sqrt(GV.Pk1*GV.Pk2*bindata[l].Pk3);

    /* Estimating dimensionless bispectrum */
    bindata[l].Qk  = bindata[l].Bk;
    bindata[l].Qk /= ( GV.Pk1*GV.Pk2 + 
		       GV.Pk2*bindata[l].Pk3 + 
		       bindata[l].Pk3*GV.Pk1 );

    bindata[l].Qk_Error  = bindata[l].Bk_Error;
    bindata[l].Qk_Error -= bindata[l].Qk*( (GV.Pk2+bindata[l].Pk3)*GV.Pk1_Error+
					   (bindata[l].Pk3+GV.Pk1)*GV.Pk2_Error+
					   (GV.Pk1+GV.Pk2)*bindata[l].Pk3_Error);
    bindata[l].Qk_Error /= ( GV.Pk1*GV.Pk2 + 
			     GV.Pk2*bindata[l].Pk3 + 
			     bindata[l].Pk3*GV.Pk1 );
    
  }// for l
  
  if(rank==0){
    printf("\n-----------------------------------------------\n");
    printf("Sending results to root process");
    printf("\n-----------------------------------------------\n");
    fflush(stdout);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  if(rank == 0){

    for(l=1; l<size; l++){
	
      for(j=0; j<Nbins; j++){

	if(l == taskBin[j]){
	  MPI_Recv(&bindata[j], 
		   sizeof(struct binStruct), 
		   MPI_BYTE, 
		   taskBin[j], 
		   j, 
		   MPI_COMM_WORLD, 
		   &status);
	  printf("Receiving bin %d from %d to 0\n",j, rank); fflush(stdout);
	}// if
      }// for j
    }// for l
  }// if rank==0
  else{
    
    for(l=1; l<size; l++){
      
      if(rank == l){
	for(j=0; j<Nbins; j++){
	  if(l == taskBin[j]){
	    printf("Sending bin %d from rank %d to 0\n",j, rank); fflush(stdout);
	    MPI_Send(&bindata[j], 
		     (int) sizeof(struct binStruct), 
		     MPI_BYTE, 
		     0, 
		     j, 
		     MPI_COMM_WORLD);
	  }// l==taskBin[j]
	}// for j
      }// rank==l
    }// for l
  }// else
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  if(rank == 0){
    
    /*
      sprintf(buff,"TodoElResultadoPaNicolas.%d",rank);
      fout = fopen(buff,"w");
      
      for(l=0; l<Nbins; l++)
      fprintf(fout,"%lf %e %e %e %e %ld %ld %e %e\n",
      bindata[l].k1, bindata[l].Pk1, bindata[l].Bk[RE], bindata[l].Bk[IM], bindata[l].Qk, 
      bindata[l].Ntri, bindata[l].Nk1, bindata[l].I[RE], bindata[l].I[IM]);
      
      fclose(fout);
      
      }
    */
  
    /* Saving data in the outfile */
    printf("\n-----------------------------------------------\n");
    printf("Saving data in %s\n", GV.OUTPUT);
  
    
    fout = fopen(GV.OUTPUT, "w");
    if(fout == NULL){
      printf("\n***********************************");
      printf("***********************************\n");
      printf("Error: Outfile could not be allocated.\n" );
      printf("***********************************");
      printf("***********************************\n\n");
      exit(0);
    }
    
    /* Writing header */
    fprintf(fout,"# NGRID          = %d\n",  GV.NGRID);
    fprintf(fout,"# GADGET VERSION = %d\n",  GV.GADGET_VERSION);
    fprintf(fout,"# L              = %lf\n", GV.L);
    fprintf(fout,"# SIM VOL        = %lf\n", GV.SIM_VOL);
    fprintf(fout,"# NP TOT         = %ld\n", GV.NP_TOT);
    fprintf(fout,"# TOTAL MASS     = %lf\n", GV.TOTAL_MASS);
    fprintf(fout,"# RHO MEAN       = %lf\n", GV.RHO_MEAN);
    fprintf(fout,"# VOL_CELL       = %lf\n", GV.VOL_CELL);
    fprintf(fout,"# H              = %lf\n", GV.H);
    fprintf(fout,"# DELTA k        = %lf\n", GV.DELTA_K);
    fprintf(fout,"# kF             = %lf\n", GV.KF);
    fprintf(fout,"# kN             = %lf\n", GV.KN);
    fprintf(fout,"# Shot Noise     = %lf\n", GV.SHOT_NOISE);
    fprintf(fout,"# SCHEME         = %s\n",  GV.SCHEME);
    fprintf(fout,"# OMEGA_M0       = %lf\n", GV.OMEGA_M0);
    fprintf(fout,"# OMEGA_L0       = %lf\n", GV.OMEGA_L0);
    fprintf(fout,"# ZRS            = %lf\n", GV.ZRS);
    fprintf(fout,"# HUBBLEPARAM    = %lf\n", GV.HUBBLEPARAM);
    fprintf(fout,"\n");
    
    fprintf(fout,"#%19s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n",
	    "k1", "k2", "k3", "P(k1)", "P(k2)", "P(k3)", "B(k1,k2,k3)", "Q(k1,k2,k3)", 
	    "Error_P(k1)", "Error_P(k2)", "Error_P(k3)", "Error_B123", "Error_Q123");
    
    /* Printing bispectrum data  */
    for(l=0; l<Nbins; l++){
      if(bindata[l].Ntri==0)
	continue;
      fprintf(fout,"%20lf %20lf %20lf %20e %20e %20e %20e %20e %20e %20e %20e %20e %20e\n",
	      GV.K1, GV.K2, bindata[l].k3, GV.Pk1, GV.Pk2, bindata[l].Pk3, 
	      bindata[l].Bk, bindata[l].Qk, 
	      GV.Pk1_Error, GV.Pk2_Error, bindata[l].Pk3_Error, 
	      bindata[l].Bk_Error, bindata[l].Qk_Error);
    }// for l
    fclose(fout);
    
  }// if rank==0
    
  
  
  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
  if( strcmp(GV.SCHEME, "D20") == 0 ){
    gsl_spline_free( splineRe );
    gsl_spline_free( splineIm );
    gsl_spline_free( splineMag2 );
    gsl_interp_accel_free( acc );

    free(k_D20);
    free(Re_D20);
    free(Im_D20);
    free(Kmag2_D20);
  }
  //fclose(fout);
  free(kpos);
  free(indexpos);
  fftw_free(denConK);
  free(q1);
  free(q2);
  free(bindata);
  free(taskBin);
  
  
  
  /////////////////////
  //* FINISHING MPI *//
  /////////////////////
  MPI_Finalize();
  
  
  
  return 0;

}
