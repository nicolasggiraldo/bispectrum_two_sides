/////////////////////////////////////////////////////
// HEADER WITH ALL DATA SRUCTTURES FOR THE PROGRAM //
/////////////////////////////////////////////////////



/////////////////////////////
// PREPROCESSOR DIRECTIVES //
/////////////////////////////

/* Index preprocessor for the C-Order */
#define INDEX(i,j,k) (1L*k) + GV.NGRID * ( (1L*j) + GV.NGRID * (1L*i) )
#define COMPLEXMAG(A,i) ( (A[i][0] * A[i][0]) + (A[i][1] * A[i][1]) )
#define VECTORMAG(x,y,z) sqrt( (x * x) + (y * y) + (z * z) )

#define RE 0  // Macro to indicate the real part of a complex value
#define IM 1  // Macro to indicate the imaginary part of a complex value

#define MIN 0 // Macro to indicate the minimum value of a set
#define MAX 1 // Macro to indicate the maximum value of a set

#define X  0  // Macro to indicate the x component of a vector
#define Y  1  // Macro to indicate the y component of a vector
#define Z  2  // Macro to indicate the z component of a vector



////////////////////////////////////////////
// GLOBAL VARIABLES FOR THE FFTW ROUTINES //
////////////////////////////////////////////

// FFTW variables
fftw_complex *denConX = NULL; // Density contrast in X-space
fftw_complex *denConK = NULL; // Density contrast in K-space
fftw_plan forwardPlan;

// FFTW variables
//fftw_complex *D20_complex = NULL; // D20 array to do Fourier transform
//fftw_plan D20Plan;



////////////////////////////////////////////////////////////////////
// GLOBAL VARIABLES FOR INTERPOLATION IN DAUBECHIES-TYPES SCHEMES //
///////////////////////////////////////////////////////////////////
gsl_interp_accel *acc = NULL;
gsl_spline *splineRe = NULL;
gsl_spline *splineIm = NULL;
gsl_spline *splineMag2 = NULL;
int len_array_D20 = 0;
double *k_D20 = NULL;
double *Re_D20 = NULL;
double *Im_D20 = NULL;
double *Kmag2_D20 = NULL;
FILE *fin_D20 = NULL;



///////////////////////////////
// STRUCTURES OF THE PROGRAM //
///////////////////////////////

/* Density contrast */
struct densityContrast{
  long int id;
  double kMag;
  int triplex[3];
};

/* Global variales */
struct globalVariables{
  int      NGRID;          // Number of cell in each axis
  long int NGRID2;         // Axis cell number to the square NGRID^2
  long int NGRID3;         // Total number of cells (NGRID3 = NGRID^3)
  int      GADGET_VERSION; // GADGET version of the snapshot
  double   L;              // Lenght of the simulation in Mpc
  double   SIM_VOL;        // Volume of the simulation
  long int NP_TOT;         // Total number of particles in the simulation
  double   TOTAL_MASS;     // Total mass of all particles in the simulation
  double   RHO_MEAN;       // Mean density of ALL the simulation
  double   VOL_CELL;       // Volume of each cell
  double   H;              // Size of the cell
  double   S_KF;           // Constant to define the value of DELTA_K = S_KF * KF
  double   DELTA_K;        // Delta k for the binning
  double   KF;             // Fundamental frequency kF
  double   SHOT_NOISE;     // shotNoise of the Power Spectrum
  double   KN;             // Nyquist frequency
  char     *FILE_NAME;     // Path of the GADGET binary
  char     *OUTPUT;        // Name of the outputfile
  char     SCHEME[4];      // Scheme used for grid assignation

  /* COSMOLOGICAL PARAMETERS OF THE SIMULATION */
  double OMEGA_M0;         //Omega matter at present time
  double OMEGA_L0;         //Omega Lambda at present time
  double ZRS;              //Redshift of the simulation
  double HUBBLEPARAM;      //Hubble parameter of the simulation

  /* Parameters asociated to power spectrum estimation in k1,k2 */
  // k1 leg
  double K1;           /* Side of the triangle to evaluate the
			  Bispectrum equilateral configuration */
  double sumDelta2_k1; // Sum delta2 in the k1 range
  long int Nk1;        // Number of elements in array q1 
  double Pk1;          // Power spectrum at k1 
  double Pk1_Error;    // Power spectrum error at k1 
  
  // k2 leg
  double K2;           /* Side of the triangle to evaluate the                                                     
			  Bispectrum equilateral configuration */
  double sumDelta2_k2; // Sum delta2 in the k1 range
  long int Nk2;        // Number of elements in array q1 
  double Pk2;          // Power spectrum at k1 
  double Pk2_Error;    // Power spectrum error at k1 
  
}GV;

/* data for output */
struct binStruct{
  double k3;           /* Side of the triangle to evaluate the                                                     
			  Bispectrum equilateral configuration */
  double sumDelta2_k3; // Sum delta2 in the k1 range
  long int Nk3;        // Number of elements in array q1 
  double Pk3;          // Power spectrum at k1 
  double Pk3_Error;    // Power spectrum error at k1 

  double I_delta3;     /* Sum of elements denConk1*denConk2*denConk3 which 
			  form a triangle. Because denConk is a complex 
			  quantity the sum is also complex */
  long int Ntri;       // Number of counted triangles for bispectrum estimation  
  double Bk;           /* Bispectrum value to measure from the simulation 
			  snapshot in general the bispectrum is a complex                                   
			  quantity but as the avarage of density constrast                           
			  is taken the imaginary part goes to zero (see                       
			  Gil-Marin et al. 2012, J. Cosm. Astrop. Phys) */
  double Bk_Error;     // Bispectrum error at k1, k2, k3
  double Qk;           /* Reduced bispectrum, this value will be evaluated
			  with the real part of the bispectrum only */
  double Qk_Error;     // Reduced bispectrum error
  double Bk_shotnoise; // Bispectrum shotnoise
};
