/*
 * Function:  read_parameters
 * --------------------
 * Reads the parameter file in which are the main parameters 
 * necessary to run the code.
 *
 * The information loaded in the same order are: 
 * FILE_NAME:      File name path of the GADGET binary file.
 * OUTPUT:         Path of the output file.
 * DELTA_K:        Width of space sampled for the calculation of 
 *                 the power spectrum. The value is given in terms 
 *                 of the fundamental frequency kF. The value
 *                 should be bigger or equal than 1.
 * 
 *  param_file_name: String with the name of the parameter file.
 *
 *  returns: Integer value.
 *            0 --> There is no error. 
 *           -1 --> There is an error loading the parameter file.
 *           -2 --> There is an error whith the settings of the 
 *                  parameter file.
 */
int read_parameters(char param_file_name[], int rank){

  FILE *cfg=NULL;  // Stream to the parameter (config) file
  int len = 200;   // Len of the read parameter
  char *buf =NULL; // buf variables to be used to read strings variables
  char *buf1=NULL;
  char *buf2=NULL; 
  char *dumb;
  
  
  if( (cfg=fopen(param_file_name,"r"))==NULL ){
    printf("%s not found.\n", param_file_name);
    // Value -1 means there is an error loading the param file
    return -1;
  }

  buf  = (char *) malloc( len*sizeof(char) );
  buf1 = (char *) malloc( len*sizeof(char) );
  buf2 = (char *) malloc( len*sizeof(char) );
  
  /* Reading FILE_NAME parameter */  
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    printf("No 'FILE_NAME' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.FILE_NAME = strdup(buf2);
    if(rank==0)
      printf("Reading from File: %s\n", GV.FILE_NAME);
  }
  
  /* Reading OUTPUT parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
    printf("No 'OUTPUT' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.OUTPUT = strdup(buf2);
    if(rank==0)
      printf("Output File: %s\n", GV.OUTPUT);
  }

  /* Reading K1 parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
     printf("No 'K1' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.K1=atof(buf2);
    if(GV.K1 > 0.0){
      if(rank==0)
	printf("K1 side: %lf\n", GV.K1);
      
      GV.sumDelta2_k1 = 0.0;
      GV.Nk1 = 0;
      GV.Pk1 = 0.0;
      GV.Pk1_Error = 0.0;
    }
    else{
      printf("Invalid 'K1' setting in configuration file.\n");
      return -2;
    }
  }

  /* Reading K2 parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
     printf("No 'K2' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.K2=atof(buf2);
    if(GV.K2 > 0.0){
      if(rank==0)
	printf("K2 side: %lf\n", GV.K2);

      GV.sumDelta2_k2 = 0.0;
      GV.Nk2 = 0;
      GV.Pk2 = 0.0;
      GV.Pk2_Error = 0.0;
    }
    else{
      printf("Invalid 'K2' setting in configuration file.\n");
      return -2;
    }
  }
      
  /* Reading S_KF parameter */
  //dumb=fgets(buf,len,cfg);
  do{dumb=fgets(buf, len, cfg);}while(dumb[0]=='#');
  if( sscanf(buf,"%s%s",buf1,buf2) < 2 ){
     printf("No 'S_KF' setting in configuration file.\n");
    return -2;
  }
  else{
    GV.S_KF=atof(buf2);
    if(GV.S_KF >= 1.0){
      if(rank==0)
	printf("Binning width in terms of the fundamental frequency kF: %lf\n", GV.S_KF);
    }
    else{
      printf("Invalid 'S_KF' setting in configuration file.\n");
      return -2;
    }
  }

  if(dumb==NULL){}
      
  fclose(cfg);
  free(buf);
  free(buf1);
  free(buf2);

  return 0;
}



/*
 * Function:  readBinaryFile
 * --------------------                  
 * Reads a binary file with the information of the cell and store 
 * the information in the data structure variable *part* it also 
 * returns the total number of particles.
 * 
 *  There are no arguments in the routiene.              
 *
 *  returns: Integer value.               
 *            0 --> There is no error. 
 *           -1 --> There is an error loading file
 *           -2 --> Structure cell could not be allocated. 
 */
int readBinaryFile(int rank){

  FILE *fdata;
  long int id_cell;
  double density_Contrast;
  int n[3]; // Number of grids in each axis for FFT estimation
  size_t err;

  if(rank == 0){
    printf("\n-----------------------------------------------\n");
    printf("Reading file:   %s\n", GV.FILE_NAME);
  }
  
  fdata = fopen(GV.FILE_NAME,"rb");
  if(fdata == NULL){
    printf("File %s cannot be open\n", GV.FILE_NAME);
    return -1;
  }
  

  /* Getting cosmological parameters of the simulation */
  err = fread(&GV.OMEGA_M0,    sizeof(double), 1, fdata);
  err = fread(&GV.OMEGA_L0,    sizeof(double), 1, fdata);
  err = fread(&GV.ZRS,         sizeof(double), 1, fdata);
  err = fread(&GV.HUBBLEPARAM, sizeof(double), 1, fdata);
  
  /* Getting simulation parameters */
  err = fread(&GV.NGRID,          sizeof(int),      1, fdata);
  err = fread(&GV.GADGET_VERSION, sizeof(int),      1, fdata);
  err = fread(&GV.L,              sizeof(double),   1, fdata);
  err = fread(&GV.NP_TOT,         sizeof(long int), 1, fdata);
  err = fread(&GV.TOTAL_MASS,     sizeof(double),   1, fdata);
  err = fread(&GV.RHO_MEAN,       sizeof(double),   1, fdata);
  err = fread(&GV.VOL_CELL,       sizeof(double),   1, fdata);
  err = fread(&GV.H,              sizeof(double),   1, fdata);
  err = fread(&(GV.SCHEME[0]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[1]),    sizeof(char),     1, fdata);
  err = fread(&(GV.SCHEME[2]),    sizeof(char),     1, fdata);
  GV.SCHEME[3] = '\0';

  GV.NGRID2     = (1L*GV.NGRID) * (1L*GV.NGRID);
  GV.NGRID3     = (1L*GV.NGRID) * (1L*GV.NGRID) * (1L*GV.NGRID);
  GV.SIM_VOL    = GV.L * GV.L * GV.L;
  GV.KF         = (2.0*M_PI) / GV.L;
  GV.DELTA_K    = GV.S_KF * GV.KF;
  GV.SHOT_NOISE = GV.VOL_CELL / GV.NP_TOT;
  GV.KN         = M_PI / GV.H;
  
  if( rank==0 ){
    printf("\n-----------------------------------------------\n");
    printf("The original snapshot has a total of %ld particles\n", GV.NP_TOT);
    printf("----------------------------------------\n");
    printf(" * Redshift...     %16.8lf\n", GV.ZRS);
    printf(" * Omega0...       %16.8lf\n", GV.OMEGA_M0);
    printf(" * OmageLa...      %16.8lf\n", GV.OMEGA_L0);
    printf(" * Hubbleparam...  %16.8lf\n", GV.HUBBLEPARAM);
    printf("----------------------------------------\n");
    printf(" * Boxsize...      %16.8lf\n", GV.L);
    printf(" * Ngrid...        %16d\n",    GV.NGRID);
    printf(" * SimMass...      %16.8e\n", GV.TOTAL_MASS);
    printf(" * Scheme...       %16s\n",    GV.SCHEME);
    printf("----------------------------------------\n");
    printf(" * kF...           %16.8lf\n", GV.KF);
    printf(" * kN...           %16.8lf\n", GV.KN);
    printf(" * DELTA_k...      %16.8lf\n", GV.DELTA_K);
    printf(" * PSshotNoise...  %16.8e\n", GV.SHOT_NOISE);
  }
  
  denConK = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
  if(denConK == NULL){
    printf("FFTW arrays could not be allocated\n");
    return -2;
  }//if
  
  if(rank == 0){
    /* Memory allocation for the input and output FFTW arrays,
       for the space, initialize input array to (1.,0.) */
    denConX = fftw_malloc( sizeof(fftw_complex) * GV.NGRID3 );
    if(denConX == NULL || denConK == NULL){
      printf("FFTW arrays could not be allocated\n");
      return -2;
    }//if
    
    /************************/
    /* FFTW plan allocation */
    /************************/
    
    /* Number of grids in each axis */
    n[X] = n[Y] = n[Z] = GV.NGRID;
    
    /* FFTW plan for a 3D Fourier transform */
    forwardPlan = fftw_plan_dft(3, n, denConX, denConK,	FFTW_FORWARD, FFTW_ESTIMATE);
    
    /****************************/
    /* Getting density contrast */
    /****************************/
    for(id_cell=0L; id_cell<GV.NGRID3; id_cell++){
      err = fread(&density_Contrast, sizeof(double), 1, fdata);
      denConX[id_cell][0] = density_Contrast;
      denConX[id_cell][1] = 0.0;
    }//for id_cell
    
  }// if rank == 0
  
  if(err){};

  fclose(fdata);
  
  return 0;
}

double sinc(double x){
  if(fabs(x)<1.0e-10)
    return 1.0;
  else
    return sin(M_PI*x)/(M_PI*x);
}

double pow2(double x){
  return x*x;
}

double pow3(double x){
  return x*x*x;
}

double pow4(double x){
  return x*x*x*x;
}



// Window function in Fourier space
double W_NGP(double k){
  return sinc( (0.5*k)/GV.KN );
}

double W_CIC(double k){
  return pow2( sinc( (0.5*k)/GV.KN ) );
}

double W_TSC(double k){
  return pow3( sinc( (0.5*k)/GV.KN ) );
}

double W_D20_Re(double k){
  return gsl_spline_eval(splineRe, (M_PI*k)/GV.KN, acc);
}

double W_D20_Im(double k){
  return gsl_spline_eval(splineIm, (M_PI*k)/GV.KN, acc);
}

double zero(double k){
  return 0.0;
}


//Sum square window function
double Sum_W2_NGP(double k){
  return 1.0;
}

double Sum_W2_CIC(double k){
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - 0.666666666666667*pow2( sine );
}

double Sum_W2_TSC(double k){
  double sine = sin( (M_PI*0.5*k)/GV.KN );
  return 1.0 - pow2( sine ) + 0.133333333333333*pow4( sine );
}

double Sum_W2_D20(double k){
  return 1.0;
}
