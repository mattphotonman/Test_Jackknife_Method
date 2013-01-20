#include <iostream>
#include <sstream>
#include <fstream>
#include <math.h>
#include <cstdlib>
#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

using namespace std;

double ran3(long *idum);
void transform_randoms(double** rans,double** trans_rans,int n,int m,int nrv);
void calc_aves_jks(double** trans_rans,double* aves_trans_rans,double** jks_trans_rans,double*** db_jks_trans_rans,int n,int nrv);
void calc_aves_errs(double* aves_trans_rans,double** jks_trans_rans,double*** db_jks_trans_rans,double* aves_single,double* errs_single,double* err_errs_single,int n,int nrv,int nq);
void calc_quantities(double* ran_vars,double* quantities,int nrv,int nq);
void output_results(double** aves,double** errs,double** err_errs,int N,int nq);
void histogram(double list[],int N,string filename);
double round_nearest(double val,double interval);
double inv_erf(double x);

const double Pi=3.1415926535897932385;

int main(int argc, char *argv[])
{
  long idum=-1;
  //Use random numbers by calling ran3(&idum)
  
  const int m=2;  //Number of different uniform random numbers
                  //generated per measurement
  const int nrv=60; //Number of random variables obtained by transforming
                    //the uniform random variables
  const int nq=159; //Number of quantities calculated from random variables

  //Read command line arguments
  /*
    N=Number of experiments
    n=Number of repeated measurements in an experiment
  */
  if (argc!=3) {
    cout << "Usage: jackknife N n\n";
    return 1;
  }
  
  int ntmp=0, N, n;
  ntmp++;
  sscanf(argv[ntmp],"%d",&N);
  ntmp++;
  sscanf(argv[ntmp],"%d",&n);

  //Declare arrays needed
  double** rans=new double* [n];
  double** trans_rans=new double* [n];
  double*** db_jks_trans_rans=new double** [n];
  double** jks_trans_rans=new double* [n];
  double* aves_trans_rans=new double [nrv];
  for (int i=0; i<n; i++) {
    rans[i]=new double [m];
    trans_rans[i]=new double [nrv];
    jks_trans_rans[i]=new double [nrv];
    db_jks_trans_rans[i]=new double* [n];
    for (int j=0; j<n; j++)
      db_jks_trans_rans[i][j]=new double [nrv];
  }
  double** aves=new double* [N];
  double** errs=new double* [N];
  double** err_errs=new double* [N];
  for (int i=0; i<N; i++) {
    aves[i]=new double [nq];
    errs[i]=new double [nq];
    err_errs[i]=new double [nq];
  }

  //Loop over experiments
  for (int ii=0; ii<N; ii++) {
    //Generate n x m random numbers in [0,1]
    for (int i=0; i<n; i++)
      for (int j=0; j<m; j++)
	rans[i][j]=ran3(&idum);
    
    //Transform these random numbers to get the distribution you want.
    transform_randoms(rans,trans_rans,n,m,nrv);

    //Calculate average, jackknife values and double
    //jackknife values of random variables
    calc_aves_jks(trans_rans,aves_trans_rans,jks_trans_rans,db_jks_trans_rans,n,nrv);

    //Calculate the average, jackknife error and error
    //on the error of quantities calculated from the 
    //random variables.
    calc_aves_errs(aves_trans_rans,jks_trans_rans,db_jks_trans_rans,aves[ii],errs[ii],err_errs[ii],n,nrv,nq);

  }

  //Ouput information about the distribution
  //of aves and errs
  output_results(aves,errs,err_errs,N,nq);

  //Delete memory allocated
  for (int i=0; i<n; i++) {
    delete [] rans[i];
    delete [] trans_rans[i];
    delete [] jks_trans_rans[i];
    for (int j=0; j<n; j++)
      delete [] db_jks_trans_rans[i][j];
    delete [] db_jks_trans_rans[i];
  }
  delete [] rans;
  delete [] trans_rans;
  delete [] db_jks_trans_rans;
  delete [] jks_trans_rans;
  delete [] aves_trans_rans;
  for (int i=0; i<N; i++) {
    delete [] aves[i];
    delete [] errs[i];
    delete [] err_errs[i];
  }
  delete [] aves;
  delete [] errs;
  delete [] err_errs;
  
  return 0;
}

double ran3(long *idum)
{
  //Returns a uniform random deviate between 0.0 and 1.0.  Set idum to any
  //negative value to initialize or reinitialize the sequence.
  
  static int inext, inextp;
  static long ma[56];         //The value 56 (range ma[1..55]) is special and
                              //should not be modified; see Knuth.
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  
  if (*idum < 0 || iff == 0) {   //Initialization
    iff=1;
    mj=labs(MSEED-labs(*idum));  //Initialize ma[55] using the see idum and the
                                 //large number MSEED.
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1; i<=54; i++) {      //Now initialize the rest of the table, in a
                                 //slightly random order, with numbers that
                                 //are not especially random.
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1; k<=4; k++)         //We randomize them by "warming up the
                                 //generator".
      for (i=1; i<=55; i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;      //Prepare indices for our first generated number.
    inextp=31;    //The constant 31 is special; see Knuth.
    *idum=1;
  }
  //Here is where we start, except on initialization.
  if (++inext == 56) inext=1;    //Increment inext and inextp, wrapping around
                                 //56 to 1.
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];       //Generate a new random number subtractively.
  if (mj < MZ) mj += MBIG;       //Be sure that it is in range.
  ma[inext]=mj;                  //Store it,
  return mj*FAC;                 //and output the derived uniform deviate.
}

void transform_randoms(double** rans,double** trans_rans,int n,int m,int nrv)
{
  //Modify this function to get random variables with the distributions
  //that you want by transforming the uniform random variables.

  for (int i=0; i<n; i++) {
    //Uniform random variable
    trans_rans[i][0]=rans[i][0];
    //Gaussian random variable
    trans_rans[i][1]=sqrt(2.0)*inv_erf(2.0*rans[i][0]-1.0);
    //Random variable obtained by squaring uniform (distribution 1/(2*sqrt(x)))
    trans_rans[i][2]=rans[i][0]*rans[i][0];
    //Same as the first three, but statistically independent of them
    trans_rans[i][3]=rans[i][1];
    trans_rans[i][4]=sqrt(2.0)*inv_erf(2.0*rans[i][1]-1.0);
    trans_rans[i][5]=rans[i][1]*rans[i][1];
    //Variables correlated to the first two sets by taking linear combinations
    //of the two.
    int ntmp=5;
    double theta=Pi/4.0;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=cos(theta)*trans_rans[i][j]+sin(theta)*trans_rans[i][k];
      }
    theta=Pi/6.0;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=cos(theta)*trans_rans[i][j]+sin(theta)*trans_rans[i][k];
      }
    theta=2.0*Pi/3.0;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=cos(theta)*trans_rans[i][j]+sin(theta)*trans_rans[i][k];
      }
    theta=5.0*Pi/6.0;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=cos(theta)*trans_rans[i][j]+sin(theta)*trans_rans[i][k];
      }
    //Variables correlated to the first two sets by taking functions of the two.
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=0.75*trans_rans[i][j]*trans_rans[i][k]*trans_rans[i][k];
      }
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	ntmp++;
	trans_rans[i][ntmp]=0.4*sin(trans_rans[i][j])*cos(trans_rans[i][k]);
      }
    //Check that I allocated the correct number of trans_rans.
    ntmp++;
    if (ntmp!=nrv) {
      cout << "Error: nrv=" << nrv << " but tried to define " << ntmp << " transformed random variables.\n";
      exit(EXIT_FAILURE);
    }
  }
}

void calc_aves_jks(double** trans_rans,double* aves_trans_rans,double** jks_trans_rans,double*** db_jks_trans_rans,int n,int nrv)
{
  for (int j=0; j<nrv; j++) {
    aves_trans_rans[j]=0.0;
    for (int i=0; i<n; i++)
      aves_trans_rans[j]+=trans_rans[i][j];
    aves_trans_rans[j]/=double(n);
    for (int i=0; i<n; i++) {
      jks_trans_rans[i][j]=(aves_trans_rans[j]*double(n)-trans_rans[i][j])/double(n-1);
      for (int k=0; k<n; k++)
	db_jks_trans_rans[i][k][j]=(aves_trans_rans[j]*double(n)-trans_rans[i][j]-trans_rans[k][j])/double(n-2);
    }
  }
}

void calc_aves_errs(double* aves_trans_rans,double** jks_trans_rans,double*** db_jks_trans_rans,double* aves_single,double* errs_single,double* err_errs_single,int n,int nrv,int nq)
{
  //Calculate averages of quantities calculated from random
  //variables.
  calc_quantities(aves_trans_rans,aves_single,nrv,nq);
  
  //Calculate jackknife values and double jackknife values
  //of quantities calculated from random variables.
  double** jackknifes=new double* [n];
  double*** db_jackknifes=new double** [n];
  for (int i=0; i<n; i++) {
    jackknifes[i]=new double [nq];
    calc_quantities(jks_trans_rans[i],jackknifes[i],nrv,nq);
    db_jackknifes[i]=new double* [n];
    for (int j=0; j<n; j++) {
      db_jackknifes[i][j]=new double [nq];
      calc_quantities(db_jks_trans_rans[i][j],db_jackknifes[i][j],nrv,nq);
    }
  }

  //Use jackknife values to calculate errors in quantities,
  //and errors on the errors.
  double* err_jk=new double [n];
  for (int j=0; j<nq; j++) {
    errs_single[j]=0.0;
    for (int i=0; i<n; i++) {
      double tmp=jackknifes[i][j]-aves_single[j];
      errs_single[j]+=tmp*tmp;
      err_jk[i]=0.0;
      for (int k=0; k<n; k++)
	if (k!=i) {
	  tmp=db_jackknifes[i][k][j]-jackknifes[i][j];
	  err_jk[i]+=tmp*tmp;
	}
      err_jk[i]*=double(n-1)/double(n-2);
      err_jk[i]=sqrt(err_jk[i]);
    }
    errs_single[j]*=double(n)/double(n-1);
    errs_single[j]=sqrt(errs_single[j]);
    err_errs_single[j]=0.0;
    for (int i=0; i<n; i++) {
      double tmp=err_jk[i]-errs_single[j];
      err_errs_single[j]+=tmp*tmp;
    }
    err_errs_single[j]*=double(n)/double(n-1);
    err_errs_single[j]=sqrt(err_errs_single[j]);
  }

  //Delete memory that was allocated
  for (int i=0; i<n; i++) {
    delete [] jackknifes[i];
    for (int j=0; j<n; j++)
      delete [] db_jackknifes[i][j];
    delete [] db_jackknifes[i];
  }
  delete [] jackknifes;
  delete [] db_jackknifes;
  delete [] err_jk;
}

void calc_quantities(double* ran_vars,double* quantities,int nrv,int nq)
{
  //Modify this function to calculate the quantities you want
  //as a function of the random variables.
  
  //Equal to the uniform random variable
  quantities[0]=ran_vars[0];
  //Equal to the Gaussian random variable
  quantities[1]=ran_vars[1];
  //Equal to the uniform random variable squared
  quantities[2]=ran_vars[2];
  //Functions of a single random variable
  quantities[3]=exp(ran_vars[0]);
  quantities[4]=exp(ran_vars[1]);
  quantities[5]=exp(ran_vars[2]);
  quantities[6]=ran_vars[0]*ran_vars[0];
  quantities[7]=ran_vars[1]*ran_vars[1];
  quantities[8]=ran_vars[2]*ran_vars[2];
  quantities[9]=cos(ran_vars[1]);
  quantities[10]=cos(ran_vars[1]*ran_vars[1]);
  quantities[11]=cos(ran_vars[1]*ran_vars[1]*ran_vars[1]);
  //Combining two random variables derived from the same
  //uniform random variable
  int ntmp=11;
  for (int j=0; j<=2; j++)
    for (int k=0; k<=2; k++)
      if (j!=k) {
	ntmp++;
	quantities[ntmp]=ran_vars[j]*ran_vars[k];
      }
  for (int j=0; j<=2; j++)
    for (int k=0; k<=2; k++)
      if (j!=k) {
	ntmp++;
	quantities[ntmp]=exp(ran_vars[j]*ran_vars[j]*ran_vars[k]);
      }
  //Multiplication of uncorrelated variables
  for (int j=0; j<=2; j++)
    for (int k=3; k<=5; k++) {
      ntmp++;
      quantities[ntmp]=ran_vars[j]*ran_vars[k];
    }
  //More complicated functions of two uncorrelated variables
  for (int j=0; j<=2; j++)
    for (int k=3; k<=5; k++) {
      ntmp++;
      quantities[ntmp]=cos(ran_vars[j]*ran_vars[j]*ran_vars[k]);
    }
  for (int j=0; j<=2; j++)
    for (int k=3; k<=5; k++) {
      ntmp++;
      quantities[ntmp]=ran_vars[j]/(ran_vars[k]+10.0);
    }
  //Combining two correlated variables (correlated via linear
  //combination)
  for (int displ=0; displ<4; displ++) {
    int var_num=5;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	var_num++;
	ntmp++;
	quantities[ntmp]=ran_vars[j]*ran_vars[var_num+9*displ];
      }
  }
  for (int displ=0; displ<4; displ++) {
    int var_num=5;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	var_num++;
	ntmp++;
	quantities[ntmp]=cos(ran_vars[j]*exp(ran_vars[var_num+9*displ]));
      }
  }
  //Combining two correlated variables (correlated via taking
  //functions)
  for (int displ=0; displ<2; displ++) {
    int var_num=41;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	var_num++;
	ntmp++;
	quantities[ntmp]=ran_vars[j]*ran_vars[var_num+9*displ];
      }
  }
  int check_num;
  for (int displ=0; displ<2; displ++) {
    int var_num=41;
    for (int j=0; j<=2; j++)
      for (int k=3; k<=5; k++) {
	var_num++;
	ntmp++;
	quantities[ntmp]=cos(ran_vars[j]*exp(ran_vars[var_num+9*displ]));
      }
    check_num=var_num+9*displ;
  }
  //Check that I'm not accessing ran_vars that I haven't defined.
  if (check_num!=nrv-1) {
    cout << "Error: nrv=" << nrv << " but last ran_var accessed was " << check_num << "\n";
    exit(EXIT_FAILURE);
  }
  //Check that I allocated the correct number of quantities.
  ntmp++;
  if (ntmp!=nq) {
    cout << "Error: nq=" << nq << " but tried to define " << ntmp << " quantities.\n";
    exit(EXIT_FAILURE);
  }
}

void output_results(double** aves,double** errs,double** err_errs,int N,int nq)
{
  //Show info about averages
  double ave_ave[nq];
  double sigma_ave[nq];
  cout << "Averages of Quantities\n\n";
  for (int j=0; j<nq; j++) {
    double ave=0.0, sigma=0.0;
    for (int i=0; i<N; i++)
      ave+=aves[i][j];
    ave/=double(N);
    ave_ave[j]=ave;
    for(int i=0; i<N; i++) {
      double tmp=aves[i][j]-ave;
      sigma+=tmp*tmp;
    }
    sigma/=double(N-1);
    sigma=sqrt(sigma);
    sigma_ave[j]=sigma;
    cout << "Quantity " << j << ":  ave=" << ave << "  sigma=" << sigma << "\n";
  }
  cout << "\n";

  //Show info about jackknife errors
  double ave_err[nq];
  double sigma_err[nq];
  cout << "Jackknife Errors of Quantities\n\n";
  for (int j=0; j<nq; j++) {
    double ave=0.0, sigma=0.0;
    for (int i=0; i<N; i++)
      ave+=errs[i][j];
    ave/=double(N);
    ave_err[j]=ave;
    for(int i=0; i<N; i++) {
      double tmp=errs[i][j]-ave;
      sigma+=tmp*tmp;
    }
    sigma/=double(N-1);
    sigma=sqrt(sigma);
    sigma_err[j]=sigma;
    cout << "Quantity " << j << ":  ave=" << ave << "  sigma=" << sigma << "\n";
  }
  cout << "\n";

  //Show info about errors on the errors
  double ave_err_err[nq];
  double sigma_err_err[nq];
  cout << "Errors on Errors of Quantities\n\n";
  for (int j=0; j<nq; j++) {
    double ave=0.0, sigma=0.0;
    for (int i=0; i<N; i++)
      ave+=err_errs[i][j];
    ave/=double(N);
    ave_err_err[j]=ave;
    for(int i=0; i<N; i++) {
      double tmp=err_errs[i][j]-ave;
      sigma+=tmp*tmp;
    }
    sigma/=double(N-1);
    sigma=sqrt(sigma);
    sigma_err_err[j]=sigma;
    cout << "Quantity " << j << ":  ave=" << ave << "  sigma=" << sigma << "\n";
  }
  cout << "\n";

  //Compare sigma_ave and ave_err.
  cout << "Compare jackknife error estimate and actual sigma of mean\n\n";
  for (int j=0; j<nq; j++)
    cout << "Jackknife error/Actual error " << j << ":  val=" << ave_err[j]/sigma_ave[j] << "  err=" << sigma_err[j]/sigma_ave[j] << "\n";
  cout << "\n";

  //Compare sigma_err and ave_err_err.
  cout << "Compare jackknife estimate of error on the error and\n";
  cout << "sigma of the jackknife error\n\n";
  for (int j=0; j<nq; j++)
    cout << "Jackknife error on error/Actual error on jackknife error " << j << ":  val=" << ave_err_err[j]/sigma_err[j] << "  err=" << sigma_err_err[j]/sigma_err[j] << "\n";
  
  //Output histograms of averages and jackknife errors.
  for (int j=0; j<nq; j++) {
    stringstream ss;
    ss.str("");
    ss << "histograms/q" << j << "_ave.dat";
    double tmp[N];
    for (int i=0; i<N; i++)
      tmp[i]=aves[i][j];
    histogram(tmp,N,ss.str());
    ss.str("");
    ss << "histograms/q" << j << "_err.dat";
    for (int i=0; i<N; i++)
      tmp[i]=errs[i][j];
    histogram(tmp,N,ss.str());
  }
}

void histogram(double list[],int N,string filename)
{
  double ave=0.0;
  double sigma=0.0;
  for (int i=0; i<N; i++)
    ave+=list[i];
  ave/=double(N);
  for (int i=0; i<N; i++) {
    double tmp=list[i]-ave;
    sigma+=tmp*tmp;
  }
  sigma/=double(N-1);
  sigma=sqrt(sigma);
  
  double max=list[0], min=list[0];
  for (int i=1; i<N; i++) {
    if (list[i]>max)
      max=list[i];
    if (list[i]<min)
      min=list[i];
  }
  
  double min_sig=(min-ave)/sigma, max_sig=(max-ave)/sigma;
  double min_bin=round_nearest(min_sig,0.25);
  double max_bin=round_nearest(max_sig,0.25);
  int num_bins=round_nearest(4.0*(max_bin-min_bin)+1.0,1.0);
  double hist[num_bins][2];
  for (int i=0; i<num_bins; i++) {
    hist[i][0]=(min_bin+0.25*double(i))*sigma+ave;
    hist[i][1]=0.0;
  }
  for (int i=0; i<N; i++) {
    double x_sig=(list[i]-ave)/sigma;
    double x_bin=round_nearest(x_sig,0.25);
    int bin_num=round_nearest(4.0*(x_bin-min_bin),1.0);
    hist[bin_num][1]++;
  }
  
  //Output to file
  ofstream fout(filename.c_str());
  for (int i=0; i<num_bins; i++)
    fout << hist[i][0] << " " << hist[i][1] << "\n";
  fout.close();
}

double round_nearest(double val,double interval)
{
  return interval*floor(val/interval+0.5);
}

double inv_erf(double x)
{
  int max_itr=1000;
  double y=0.0;
  double epsilon=1E-8;
  double func=erf(y);
  double dfunc=2.0*exp(-y*y)/sqrt(Pi);
  double res=fabs(func-x);
  for (int i=0; i<max_itr; i++) {
    if (res<=epsilon) break;
    y=(x-func)/dfunc+y;
    func=erf(y);
    dfunc=2.0*exp(-y*y)/sqrt(Pi);
    res=fabs(func-x);
  }
  return y;
}
