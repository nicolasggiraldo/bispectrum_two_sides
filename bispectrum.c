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
  unsigned long int rand_i;        /* Index for the subsample of density contrast */
  unsigned long int rand_j;        /* in the range between                        */
  //unsigned long int rand_k; /* k1-DELTA_K/2 to k1+DELTA_K/2                */
  long int id_cell;                // Counter for the id of the cell
  double kMag;                     /* Magnitude of the wavenumber vector */
  FILE *fout=NULL;                 // File handler to output
  double *kpos;                    /* Array of positions values according to FFTW 
				      k-position convention */ 
  int *indexpos=NULL;              /* Array of index values according to FFTW 
				      k-position convention */
  int indexcord[2];
  int m3[3];
  int rank, size;
  int Nbins;
  int *taskBin=NULL;
  double (*W_k)(double) = NULL; // Memory addres to the window function
  struct densityContrast *q1=NULL; /* Array of structure with the subset of 
				      density contrast with range between 
				      k1-DELTA_K/2 to k1+DELTA_K/2 */
  struct densityContrast *q2=NULL; /* Array of structure with the subset of 
				      density contrast with range between 
				      k2-DELTA_K/2 to k2+DELTA_K/2 */
  struct binStruct *bindata=NULL;  /* Data structure with information of the 
				      measures in the taken bins */
  
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
  for(l=0; l<sizeof(fftw_complex); l++)
    MPI_Bcast(denConK+(l*chunk), chunk*sizeof(fftw_complex), MPI_BYTE, 0, MPI_COMM_WORLD); 
  
  
  
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
  indexcord[MAX] =  (GV.NGRID/2)-1;
  
  
  
  ///////////////////////////////////////////
  //* DEFINING THE WINDOW FUNCTION TO USE *//
  ///////////////////////////////////////////
  if( strcmp(GV.SCHEME, "NGP") == 0 ){
    W_k = W_NGP; //NGP
  }
  else if( strcmp(GV.SCHEME, "CIC") == 0 ){ 
    W_k = W_CIC; //CIC
  }
  else if( strcmp(GV.SCHEME, "TSC") == 0 ){
    W_k = W_TSC; // TSC
  }
  else if( strcmp(GV.SCHEME, "D20") == 0 ){

    len_array_D20 = 39;
    x_D20 = (double *) calloc( len_array_D20, sizeof(double));
    y_D20 = (double *) calloc( len_array_D20, sizeof(double));

    x_D20[0]  = 0.000000000000;
    x_D20[1]  = 0.105263157895;
    x_D20[2]  = 0.210526315789;
    x_D20[3]  = 0.315789473684;
    x_D20[4]  = 0.421052631579;
    x_D20[5]  = 0.526315789474;
    x_D20[6]  = 0.631578947368;
    x_D20[7]  = 0.736842105263;
    x_D20[8]  = 0.842105263158;
    x_D20[9]  = 0.947368421053;
    x_D20[10] = 1.052631578947;
    x_D20[11] = 1.157894736842;
    x_D20[12] = 1.263157894737;
    x_D20[13] = 1.368421052632;
    x_D20[14] = 1.473684210526;
    x_D20[15] = 1.578947368421;
    x_D20[16] = 1.684210526316;
    x_D20[17] = 1.789473684211;
    x_D20[18] = 1.894736842105;
    x_D20[19] = 2.000000000000;
    x_D20[20] = 2.105263157895;
    x_D20[21] = 2.210526315789;
    x_D20[22] = 2.315789473684;
    x_D20[23] = 2.421052631579;
    x_D20[24] = 2.526315789474;
    x_D20[25] = 2.631578947368;
    x_D20[26] = 2.736842105263;
    x_D20[27] = 2.842105263158;
    x_D20[28] = 2.947368421053;
    x_D20[29] = 3.052631578947;
    x_D20[30] = 3.157894736842;
    x_D20[31] = 3.263157894737;
    x_D20[32] = 3.368421052632;
    x_D20[33] = 3.473684210526;
    x_D20[34] = 3.578947368421;
    x_D20[35] = 3.684210526316;
    x_D20[36] = 3.789473684211;
    x_D20[37] = 3.894736842105;
    x_D20[38] = 4.000000000000;
    
    y_D20[0]  = 1.000000000000;
    y_D20[1]  = 0.999999986964;
    y_D20[2]  = 0.999999947844;
    y_D20[3]  = 0.999999865073;
    y_D20[4]  = 0.999996615502;
    y_D20[5]  = 0.999865850668;
    y_D20[6]  = 0.997888980263;
    y_D20[7]  = 0.983901583682;
    y_D20[8]  = 0.929829547217;
    y_D20[9]  = 0.801587086807;
    y_D20[10] = 0.597770493103;
    y_D20[11] = 0.367763547417;
    y_D20[12] = 0.178328576874;
    y_D20[13] = 0.064525672159;
    y_D20[14] = 0.016095961628;
    y_D20[15] = 0.002430216092;
    y_D20[16] = 0.000174403891;
    y_D20[17] = 0.000003474525;
    y_D20[18] = 0.000000003493;
    y_D20[19] = 0.000000000009;
    y_D20[20] = 0.000000002651;
    y_D20[21] = 0.000001909750;
    y_D20[22] = 0.000068979711;
    y_D20[23] = 0.000666947923;
    y_D20[24] = 0.002917323470;
    y_D20[25] = 0.007255209141;
    y_D20[26] = 0.011531053543;
    y_D20[27] = 0.012504254578;
    y_D20[28] = 0.009622917401;
    y_D20[29] = 0.005399117417;
    y_D20[30] = 0.002259675865;
    y_D20[31] = 0.000719121329;
    y_D20[32] = 0.000174034058;
    y_D20[33] = 0.000030590879;
    y_D20[34] = 0.000003474474;
    y_D20[35] = 0.000000201352;
    y_D20[36] = 0.000000003493;
    y_D20[37] = 0.000000000002;
    y_D20[38] = 0.000000000007;
    
    // GSL interpolation allocation
    acc    = gsl_interp_accel_alloc(); // accelerator
    spline = gsl_spline_alloc(gsl_interp_cspline, len_array_D20); // spline

    // GSL init
    gsl_spline_init(spline, x_D20, y_D20, len_array_D20);
    
    W_k = W_D20; // D20
  }
  
  
  
  ////////////////////////////////////
  //* DECONVOLVING WINDOW FUNCTION *//
  ////////////////////////////////////
  
  for(i=0; i<GV.NGRID; i++){
    for(j=0; j<GV.NGRID; j++){
      for(k=0; k<GV.NGRID; k++){
	
	id_cell = INDEX(i,j,k);
	
	denConK[id_cell][0] /= ( W_k(kpos[i])*W_k(kpos[j])*W_k(kpos[k]) );
	denConK[id_cell][1] /= ( W_k(kpos[i])*W_k(kpos[j])*W_k(kpos[k]) );
	
      }// for k
    }// for j
  }// for i
    
    
    
  ////////////////////////////////////
  //* GETTING SET K FOR THE K1 LEG *//
  ////////////////////////////////////
  q1 = (struct densityContrast *) calloc( floor(3.0*M_PI*GV.NGRID*GV.NGRID*GV.S_KF*0.5), 
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
  q2 = (struct densityContrast *) calloc( floor(3.0*M_PI*GV.NGRID*GV.NGRID*GV.S_KF*0.5), 
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
	  
	  GV.sumDelta2_k1 += COMPLEXMAG(denConK, id_cell) ;
	  
	}// if bindata[l].k1-GV.DELTA_K*0.5 < kMag < bindata[l].k1+GV.DELTA_K*0.5 )
	
      }// for k
    }// for j
  }// for i
    
  GV.Pk1  = GV.sumDelta2_k1 / GV.Nk1;
  GV.Pk1 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
  GV.Pk1 -= GV.SHOT_NOISE;


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
	  
	  GV.sumDelta2_k2 += COMPLEXMAG(denConK, id_cell) ;
	  
	}// if bindata[l].k2-GV.DELTA_K*0.5 < kMag < bindata[l].k2+GV.DELTA_K*0.5
	
      }// for k
    }// for j
  }// for i
    
  GV.Pk2  = GV.sumDelta2_k2 / GV.Nk2;
  GV.Pk2 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
  GV.Pk2 -= GV.SHOT_NOISE;


  
  if(rank==0){
    printf("\n-----------------------------------------------\n");
    printf("Bispectrum calculation");
    printf("\n-----------------------------------------------\n");
  }
  
  MPI_Barrier(MPI_COMM_WORLD);
  
  for(l=0; l<Nbins; l++){
    if(rank != taskBin[l])
      continue;
    
    printf("rank:%3d, k3 = %lf\n", rank, bindata[l].k3);
    fflush(stdout);
    
    // Getting sets associated to k3
    for(i=0; i<GV.NGRID; i++){
      for(j=0; j<GV.NGRID; j++){
	for(k=1; k<(GV.NGRID/2+1); k++){
	
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	
	  if( ( bindata[l].k3-GV.DELTA_K*0.5 < kMag ) && ( kMag < bindata[l].k3+GV.DELTA_K*0.5 ) ){
	    
	    id_cell = INDEX(i,j,k);
	    
	    bindata[l].Nk3++;
	    
	    bindata[l].sumDelta2_k3 += COMPLEXMAG(denConK, id_cell) ;
	  
	  }// if bindata[l].k2-GV.DELTA_K*0.5 < kMag < bindata[l].k2+GV.DELTA_K*0.5 )
	  
	}// for k
      }// for j
    }// for i
    
    bindata[l].Pk3  = bindata[l].sumDelta2_k3 / bindata[l].Nk3;
    bindata[l].Pk3 *= GV.SIM_VOL / (1.0 * GV.NGRID3 * GV.NGRID3);
    bindata[l].Pk3 -= GV.SHOT_NOISE;
    
    
    for( rand_i=0; rand_i<GV.Nk1; rand_i++){
      for( rand_j=0; rand_j<GV.Nk2; rand_j++){
	
	/*
	  do{
	  rand_i = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_j = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  rand_k = gsl_rng_uniform_int( r, (unsigned long int) Nk1 );
	  }while( (rand_i==rand_j) || (rand_i==rand_k) || (rand_j==rand_k) );
	*/
	
	/* Looking for closed triangles */
	if(q1[rand_i].id == q2[rand_j].id)
	  continue;
	
	m3[X] = - q1[rand_i].triplex[X] - q2[rand_j].triplex[X];
	m3[Y] = - q1[rand_i].triplex[Y] - q2[rand_j].triplex[Y];
	m3[Z] = - q1[rand_i].triplex[Z] - q2[rand_j].triplex[Z];
	
	if( indexcord[MIN]<=m3[X] && m3[X]<=indexcord[MAX] && 
	    indexcord[MIN]<=m3[Y] && m3[Y]<=indexcord[MAX] && 
	    indexcord[MIN]<=m3[Z] && m3[Z]<=indexcord[MAX] ){

	  i = (m3[X]>=0) ? m3[X] : GV.NGRID+m3[X];
	  j = (m3[Y]>=0) ? m3[Y] : GV.NGRID+m3[Y];
	  k = (m3[Z]>=0) ? m3[Z] : GV.NGRID+m3[Z];
	  
	  kMag = VECTORMAG(kpos[i],kpos[j],kpos[k]);
	  
	  if( ( bindata[l].k3-GV.DELTA_K*0.5 < kMag ) && 
	      ( kMag < bindata[l].k3+GV.DELTA_K*0.5 ) ){
	    
	    id_cell = INDEX(i,j,k);

	    if( q1[rand_i].id==id_cell || q2[rand_j].id==id_cell)
	      continue;
	    
	    bindata[l].I_delta3 += 
	      (+denConK[q1[rand_i].id][0]*denConK[q2[rand_j].id][0]*denConK[id_cell][0]
	       -denConK[q1[rand_i].id][0]*denConK[q2[rand_j].id][1]*denConK[id_cell][1]
	       -denConK[q1[rand_i].id][1]*denConK[q2[rand_j].id][0]*denConK[id_cell][1]
	       -denConK[q1[rand_i].id][1]*denConK[q2[rand_j].id][1]*denConK[id_cell][0] );
	    	        
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
    bindata[l].Bk_shotnoise = GV.SHOT_NOISE * (GV.Pk1+GV.Pk2+bindata[l].Pk3) + (GV.SHOT_NOISE*GV.SHOT_NOISE);
    
    /* Substracting shotnoise term  */
    bindata[l].Bk -= bindata[l].Bk_shotnoise;

    /* Estimating dimensionless bispectrum */
    bindata[l].Qk  = bindata[l].Bk;
    bindata[l].Qk /= ( GV.Pk1*GV.Pk2 + GV.Pk2*bindata[l].Pk3 + bindata[l].Pk3*GV.Pk1 );

    //Pk1_Error = sqrt( (GV.KF*GV.KF*GV.KF)/(4.0*M_PI*k1*k1*GV.DELTA_K) ) * Pk1;
    //Bk_Error  = sqrt( (GV.KF*GV.KF*GV.KF)/(4.0*M_PI*k1*k1*GV.DELTA_K) ) * Pk1;
    
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
	  MPI_Recv(&bindata[j], sizeof(struct binStruct), MPI_BYTE, taskBin[j], j, MPI_COMM_WORLD, &status);
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
	    MPI_Send(&bindata[j], (int) sizeof(struct binStruct), MPI_BYTE, 0, j, MPI_COMM_WORLD);
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
    
    fprintf(fout,"#%19s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s %20s\n",
	    "k1", "k2", "k3", "P(k1)", "P(k2)", "P(k3)", "B(k1,k2,k3)", "Q(k1,k2,k3)", 
	    "sumDelta2_k1", "sumDelta2_k2", "sumDelta2_k3", "Nk1", "Nk2", "Nk3", 
	    "I_delta3", "Ntri");
    
    /* Printing bispectrum data  */
    for(l=0; l<Nbins; l++){
      if(bindata[l].Ntri==0)
	continue;
      fprintf(fout,"%20lf %20lf %20lf %20e %20e %20e %20e %20e %20e %20e %20e %20ld %20ld %20ld %20e %20ld\n",
	      GV.K1, GV.K2, bindata[l].k3, GV.Pk1, GV.Pk2, bindata[l].Pk3, bindata[l].Bk, bindata[l].Qk, 
	      GV.sumDelta2_k1, GV.sumDelta2_k2, bindata[l].sumDelta2_k3, GV.Nk1, GV.Nk2, bindata[l].Nk3, 
	      bindata[l].I_delta3, bindata[l].Ntri);
    }// for l
    fclose(fout);
    
  }// if rank==0
    
  
  
  
  ///////////////////
  //* FREE MEMORY *//
  ///////////////////
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
