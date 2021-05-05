
///////////////DECLARATION SEQUENCE///////////////////////

#include <iostream>
#include <string.h>
#include <fstream>
#include <math.h>
#include <time.h>
#include <cstdlib>
#include <ctime>
#include <conio.h>
#include <stdlib.h>    
#include <stdio.h>

///////////////RAN2PARAMETERS///////////////////////

#define IM1   2147483563
#define IM2   2147483399
#define AM    (1.0/IM1)
#define IMM1  (IM1-1)
#define IA1   40014
#define IA2   40692
#define IQ1   53668
#define IQ2   52774
#define IR1   12211
#define IR2   3791
#define NTAB  32
#define NDIV  (1+IMM1/NTAB)
#define EPS   1.2e-7
#define RNMX  (1.0 - EPS)
#define SWAP(a,b) tempr=(a);(a)=(b);(b)=tempr
#define Pi 3.14159

///////////////PARAMETERS///////////////////////

#define T 200000 //Time of simulation in msec
#define N 100//number of neurons
#define dt 0.1 //time step for Euler Maruyama integration
#define w_o 0.21//synaptic connectivity weight (criticality) 
#define beta 25 //activation function non-linear gain
#define rho 0.15//connection probability
#define Omega 10.0 //network spatial scale in mm
#define h 0.10//activation function threshold (inflexion point)
#define c_min 0.1 //minimal conduction velocity in m/s
#define c_max 100//maximal conduction velocity in m/s
#define Q 300//counter

//////////////////////FUNCTION DEFINITIONS///////////////////////
double f(double u);//neuron response function
double f_prime(double u);//derivative of f
float ran2(long *idum);//random number generator - initial cond.

//////////////////////OBJECT DEFINITIONS///////////////////////
double x[N];//neuron position along x axis
double y[N];//neuron position along y axis
double z[N];//neuron position along z axis

double c[N][N];//conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal tract lengths matrix
int tau[N][N];//delays matrix
double Velocity[Q];//velocity matrix

double xi[T];//noise
double u[T];//dynamics

double Mean[Q];//mean value of the linearized mean activity
double Variance[Q];//variance of the linearized mean activity
double TD_FRANK_Variance[Q];//variance from the TD Frank paper: 
//"Kramers-Moyal expansion for stochastic differential equations 
//with single and multiple delays: Applications to financial physics and neurophysics"
//Physics Letters A 360, 552-562 (2007).
using namespace std;


int main()
{
	
				for (int q=0;q<Q;q++)
				{
	
					Velocity[q] = 0.01+10/((double)Q)*q;
//////////////////////DEFINING CONNECTIVITY, DELAYS AND CONDUCTION VELOCITIES///////////////////////

					for (int i=0;i<N;i++)
					{
						long seed1 = (long)21+i*187198+56*i+12*1+1+q*485;
						long seed2 = (long)21+i*56+56*i+11*1+13783*q;
						long seed3 = (long)21+i*789+56*i+8*1+1535*q+q;
		                          	  					
						x[i] = ran2(&seed1)*Omega;
						y[i] = ran2(&seed2)*Omega;
						z[i] = ran2(&seed3)*Omega;
					}
					
					for(int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							w[i][j] = 0;
							l[i][j] = fabs(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])));
							c[i][j] = Velocity[q];
							tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
							
						}
						
					}
				
//////////////////////MAXIMAL DELAY///////////////////////

					int tau_max =-20000;
		
					for (int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							if (tau[i][j]>tau_max)
							{
								tau_max = tau[i][j];
							}
						}
					}
					
					
//////////////////////NOISE///////////////////////
					
					for (int t=0;t<T;t++)
					{	
						long noisyseed1 = (long)21*t+878+243+767*t+565656*q+656;
            			long noisyseed2 = (long)69*t+11+t+66+t*t+34+656*q+q+87; 
	             
						xi[t] = sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
						
					}
					
					
//////////////////////HISTORY///////////////////////			

					for (int t=0;t<tau_max;t++)
					{			
						long noisyseed1 = (long)21*t+57584+243+767*q+8767*q;
            			long noisyseed2 = (long)69*t+11+t+565+t*t+34+656+7567*q+q; 
	             
						u[t] = sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
					}
					
					
//////////////////////EULER MARUYAMA INTEGRATION SCHEME///////////////////////
					double fixedPoint=0.00334;
					double noiseIntensity = w_o*w_o*f(fixedPoint)/(rho*(double) N);	
 		
					for (int t=tau_max;t<T;t++)
 					{
 						double sum=0;
			
 						for(int i=0;i<N;i++)
 						{
 							for(int j=0;j<N;j++)
 							{
 								sum = sum+1/((double) N*N)*u[t-tau[i][j]];
							}
						}
 			
 						u[t+1] = u[t]+dt*(-u[t]+w_o*f_prime(fixedPoint)*sum)+sqrt(noiseIntensity*dt)*xi[t];
					}
		

//////////////////////MEAN VALUE AND VARIANCE OF THE MEAN ACTIVITY///////////////////////
		
   				    double Gamma=0;
			
        			for (int i=0;i<N;i++)
        			{
        				for (int j=0;j<N;j++)
        				{
        					Gamma = Gamma+1/((double)N*N)*exp(-tau[i][j]*dt);
						}
					}
		
        			TD_FRANK_Variance[q] = noiseIntensity/(2*(1-w_o*f_prime(fixedPoint)*Gamma));     
        
    
    				double mean=0;
    				double var=0;
    	
    				for (int t=T/2;t<T;t++)
    				{
    					mean = mean+1/((double)T/2)*u[t];
					}
		
					for (int t=T/2;t<T;t++)
					{
						var = var+1/((double)T/2)*(u[t]-mean)*(u[t]-mean);	
					}
		
					Mean[q] = mean;
					Variance[q] = var;
		
					cout<<Velocity[q]<<" "<<Mean[q]<<"	"<<Variance[q]<<" "<<TD_FRANK_Variance[q]<<endl;
		
				}// q loop
		
					
			
					
/////////////////////OUTPUT ROUTINE/////////////////////////////					
					
		ofstream outfile;
      	
    	outfile.open("ADM  -  Mean and Variance.txt", ios::out);
    	for(int q=0;q<Q;q++)
    	{
	   		outfile<<Velocity[q]<<"	"<<Mean[q]<<"	"<<Variance[q]<<endl;
    	}  
    	outfile.close(); 
      	
    
    	outfile.open("ADM  -  Variance and TD-Frank Variance.txt", ios::out);
    	for(int q=0;q<Q;q++)
    	{
	   		outfile<<Velocity[q]<<"	"<<Variance[q]<<" "<<TD_FRANK_Variance[q]<<endl;
    	}  
    	outfile.close(); 
					
		cout<<"Simulations complete..."<<endl;
      
return 0;    
}


///////////////FUNCTION DEFINITIONS///////////////////////

float ran2(long *idum)
{
  int j;
  long k;
  static long idum2 = 123456789;
  static long iy = 0;
  static long iv[NTAB];
  float temp;

  if (*idum <= 0) {                             /* initialize */
    if (-(*idum) < 1)                           /* prevent idum == 0 */
      *idum = 1;
    else
      *idum = -(*idum);                         /* make idum positive */
    idum2 = (*idum);
    for (j = NTAB + 7; j >= 0; j--) {           /* load the shuffle table */
      k = (*idum) / IQ1;
      *idum = IA1 * (*idum - k*IQ1) - k*IR1;
      if (*idum < 0)
        *idum += IM1;
      if (j < NTAB)
        iv[j] = *idum;
    }
    iy = iv[0];
  }

  k = (*idum) / IQ1;
  *idum = IA1 * (*idum - k*IQ1) - k*IR1;
  if (*idum < 0)
    *idum += IM1;
  k = idum2/IQ2;
  idum2 = IA2 * (idum2 - k*IQ2) - k*IR2;
  if (idum2 < 0)
    idum2 += IM2;
  j = iy / NDIV;
  iy = iv[j] - idum2;
  iv[j] = *idum;
  if (iy < 1)
    iy += IMM1;
  if ((temp = AM * iy) > RNMX)
    return RNMX;                                /* avoid endpoint */
  else
    return temp;
}



double f(double u)
{
	double output;
	output = 1/(1+exp(-beta*(u-h)));
	return output;
	
}

double f_prime(double u)
{
	double output;
	output = 1/(1+exp(-beta*(u-h)))*1/(1+exp(-beta*(u-h)))*exp(-beta*(u-h))*beta;
	return output;
	
}






