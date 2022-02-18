///////////////DECLARATION SEQUENCE///////////////

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

///////////////RAN2PARAMETERS///////////////

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

///////////////PARAMETERS///////////////

#define T 1024*8//Time of simulation in msec
#define T_window 1000//size of time window (ms)
#define S 21//stimulated neurons
#define Q 60//counter
#define N 100//number of neurons
#define dt 0.1//time step for integration
#define Dt 0.001//Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec
#define w_o 0.21//synaptic connectivity weight (criticality)
#define I_o 0//input amplitude (a.u.)
#define temporal_window 5//Temporal window for spike coherence measures
#define beta 25//activation function non-linear gain
#define rho 0.15//connection probability
#define Omega 100.0//network spatial scale in mm
#define h 0.10//activation function threshold (inflexion point) 
#define D 0.00000//0.015// ADD NOISE TO SUPPRESS SEIZURE Noise variance
#define range_ee N
#define c_min 0.1//minimal conduction velocity in m/s
#define c_max 100//maximal conduction velocity in m/s

#define alpha_retraction 0.0001//rate constant for ADM retraction
#define alpha_formation 0.001//rate constant for ADM formation

#define t_stim_on 0//stimulation onset time
#define t_stim_off T//stimulation offset time

#define T_plastic_on T//ADM onset time (default to initial simulation starting time 0)
#define T_plastic_off T//ADM offset time (default to full simulation time T)

#define t_delay_statistics_initial 1000//time for baseline delay distribution calculations
#define t_delay_statistics_final T-100//time for stabilized delay distribution calculations
#define delay_graining 50//bins for delay distribution calculations
#define speed_graining 50//bins for velocity distribution calculations
#define max_delay 200//maximum for delay distribution calculations
#define max_speed 30//maximum for velocity distribution calculations

//stimulated nodes
#define node_1 0
#define node_2 30

//recording nodes
#define node_3 70
#define node_4 100

///////////////FUNCTION DEFINITIONS///////////////

double f(double u, double threshold);//neuron response function 
double f_prime(double u, double threshold);//neuron response function derivative
float ran2(long *idum);//random number generator - initial cond.
void shuffle_indices(int Array[], int size);//shuffle an array with int  entries
void shuffle(double Array[], int size);//shuffle an array with double  entries
void sort(double Array[], int size);//sort array 
void four1(double data[], unsigned long nn, int isign);// Fourier transform
double coherence(double X1[], double X2[], int window);// coherence

///////////////OBJECT DEFINITIONS///////////////

double u_e[N][T];//Dynamics of E units
double mean_u_e[T];//mean network activity E units
double rate_stim[T];//stimulus rate
double rate_rec[T];
double Input[N][T];//input matrix
double c[N][N];//conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal tract lengths matrix
int tau[N][N];//delays matrix
double delay[N][N];//delays 
int tau_o[N][N];//baseline delay matrix
double X_e[N][T];//spiking activity
double r_e[N][T];//rates
double Freq[T/2];//Array of output frequencies for PSD display
double PSD_network[T];//power spectral density network
double PSD_signal[T];//power spectral density signal
double MI_single_trial[Q];//mutual information single trial
double MI[S][S];//mutual information
double Var_MI[S][S];//mutual information variance
double FREQUENCY[S];//frequency vector
double AMPLITUDE[S];//amplitude vector
double MEAN_MI_ACROSS_FREQUENCIES[S];// mean mutual information across frequences
double MEAN_MI_ACROSS_FREQUENCIES_var[S];// variance of mean mutual information across frequences
double CV[S];//conduction velocities

double Delay_distribution_final[delay_graining];//array for delay distribution - stabilized
double Delay_distribution_initial[delay_graining];//array for delay distribution - baseline
double Speed_distribution_final[speed_graining];//array for velocity distribution - stabilized
double Speed_distribution_initial[speed_graining];//array for velocity distribution - baseline

using namespace std;


int main()
{
	
srand (time(NULL));


		cout<<"Importing data..."<<endl;
		ifstream sfile("LENGTHS100.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			sfile>>l[i][j];
         		
			}
                         
         }  
        
        
        ifstream ufile("CVS100.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			ufile>>c[i][j];
         		
			}
                         
      }  
          			
		double sum=0;
		for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         			if(i==j)
							{
								w[i][j]=0;
							}
							else
							{
								if (fabs(i-j)<range_ee)
								{
									long seed=(long)21+i*j*187198+56*j+12*1+1;
		                          	  
									if(ran2(&seed)<rho)
									{
											w[i][j]=w_o/rho;
											
									}
								}
								else
								{
									w[i][j]=0;
								}
								
							}
						
         		
			}
                         
         }  
									
								
				
			
		for(int s1=0;s1<S;s1++)
		{
			CV[s1] =c_min+100*s1/((double) S-1);
			
			MEAN_MI_ACROSS_FREQUENCIES[s1]=0;
			MEAN_MI_ACROSS_FREQUENCIES_var[s1]=0;
				
			for(int s2=0;s2<S;s2++)
			{
				AMPLITUDE[s2] = 0.05;
			
				FREQUENCY[s2]=1+100*s2/((double)S);
		
				cout<<s1<<"	"<<s2<<endl;
				
				for (int q=0;q<Q;q++)
				{
					MI_single_trial[q]=0;
				}
				for(int q=0;q<Q;q++)
				{
					
					///////////////reset connectivity///////////////

						for (int i=0;i<N;i++)
						{
							for (int j=0;j<N;j++)
							{
								w[i][j]=0;
					    		
							}
						}
						
						double sum=0;
						for(int i=0;i<N;i++)
				         {
				         	for (int j=0;j<N;j++)
				         	{
				         			if(i==j)
											{
												w[i][j]=0;
											}
											else
											{
												if (fabs(i-j)<range_ee)
												{
													long seed=(long)21+i*j*187198+56*j+s2+s1*q+12*1+18*q;
						                          	  
													if(ran2(&seed)<rho)
													{
															w[i][j]=w_o/rho;
															
													}
												}
												else
												{
													w[i][j]=0;
												}
												
											}
												
											c[i][j] = CV[s1];
											delay[i][j] = (double) (l[i][j]/(c[i][j]));
						         			tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
								 			tau_o[i][j]=tau[i][j];
											
											
				         		
							}
				                         
				         }  
						//reset init cond and set input
			        	for (int t=0;t<T;t++)
						{
							
							double sum=0;
							for(int n=0;n<N;n++)
			        		{
			                            long d=rand()%105;      
			                            long noisyseed1=(long)21*t+n*s2*187198+34*s1+q+56*d+12*1+1;
			                            long noisyseed2=(long)69*t*n+11*s2+t+45*d+2*q+10*s1; 
			                          
			            				 
			            				u_e[n][t]=sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));
			            				X_e[n][t]=0;
			            				
			            				if (t>t_stim_on&&t<t_stim_off)
										{
												if (n>node_1 && n<node_2)
												{
													Input[n][t] = AMPLITUDE[s2]*sin(2*Pi*FREQUENCY[s2]/100*t*dt); 
												}
												else
												{
													Input[n][t]=0;
												}
									//I=I_o;
										}
							}
							mean_u_e[t]=0;
							rate_stim[t]=0;
							rate_rec[t]=0;
							PSD_signal[t]=0;
							PSD_network[t]=0;
						
			           }
			           
///////////////EULER SCHEME///////////////				

					for (int t=0;t<T;t++)
					{
						
							
						double sum_ee; double sum_theo=0;double Gamma=0;
						for (int i=0;i<N;i++)
						{
						
							sum_ee=0;
						
							for (int j=0;j<N;j++)
							{								
									sum_ee=sum_ee+1/((double)N)*w[i][j]*X_e[j][t-tau[i][j]];
									
	 								
							}
							
						
							u_e[i][t+1]=u_e[i][t]+dt*(-1*u_e[i][t]+sum_ee+Input[i][t]);
					
										
										   // Poisson neurons
										 
										  long seed_q =rand()%101;  
						   				  long seed_e = (long) 12+i*s1+t+15+q*t*seed_q+56*s2+i+99*t+2*s1;              
			                              double p_e =ran2(&seed_e);
			                              double p_fire_e = ( 1-exp(-f(u_e[i][t+1],h)*dt)); 
			                              
			                              
			                              if (p_e<p_fire_e)//spike occurs
			                             
			                               {                   
			                                               X_e[i][t+1] =1/dt; 
			                                               
			        										
			                               }
			                               else//no spike
			                               {
			                                               X_e[i][t+1]=0; 
			                                               			                               
										   }
										   
										   
			                              
						}
								
				
						for (int i=0;i<N;i++)
						{
									mean_u_e[t+1] = mean_u_e[t+1]+1/((double)N)*u_e[i][t+1];
				
						}
										
						
					}// t loop
					

///////////////Firing Rate///////////////
								int twindow=20;
								int range_stim= (int) node_2-node_1;
								int range_rec = (int) node_4-node_3 ;
								for (int t=twindow;t<T;t++)
								{
									double rate_stim_counter=0;
									double rate_rec_counter=0;
								
									double counts_stim=0;
									double counts_rec=0;
								
							
									for (int time=0;time<twindow;time++)
									{
										for (int index=node_1;index<node_2;index++)
										{
												if(X_e[index][t-(int)twindow/2+time]>0.1)
												{
													rate_stim_counter=rate_stim_counter+1;
													counts_stim = counts_stim+1;
												}
												else{}
										}
										for (int index=node_3;index<node_4;index++)
										{
												if(X_e[index][t-(int)twindow/2+time]>0.1)
												{
													rate_rec_counter=rate_rec_counter+1;
													counts_rec = counts_rec+1;
												}
												else{}
										}
									}
									
										rate_stim[t]=rate_stim_counter/((double) twindow*Dt)*1/((double)range_stim);
										rate_rec[t]=rate_rec_counter/((double) twindow*Dt)*1/((double)range_rec);
										
								}
									
									
			///////////////MUTUAL INFORMATION///////////////
						for (int k=0;k<(int)T;k++)
						{
												
												PSD_network[k] = 	rate_rec[k];
												PSD_signal[k]=		Input[node_1+1][k];													 		 
						}
									
						unsigned long nn=T/2;
						four1(PSD_network-1, nn,1);
																 
						for (int k=0;k<T/2;k++)
						{
											PSD_network[k] = 1/((double)T*T)*(fabs(PSD_network[k])*fabs(PSD_network[k])+fabs(PSD_network[(int)T-k])*fabs(PSD_network[(int)T-k]));
						}
						
						four1(PSD_signal-1, nn,1);
																 
						for (int k=0;k<T/2;k++)
						{
											PSD_signal[k] = 1/((double)T*T)*(fabs(PSD_signal[k])*fabs(PSD_signal[k])+fabs(PSD_signal[(int)T-k])*fabs(PSD_signal[(int)T-k]));
						}
						//normalize power spectra
						double norm1=0;double norm2=0;
						for (int k=0;k<T/2;k++)
						{
											norm1 = norm1 + PSD_network[k];
											norm2 = norm2 + PSD_signal[k];
			
						}
						for (int k=0;k<T/2;k++)
						{
										PSD_network[k] = PSD_network[k] /norm1;
										PSD_signal[k]=PSD_signal[k]/norm2;
			
						}
						double convo=0;
						double autoconv1 =0;
						double autoconv2 =0;
						for (int k=5;k<T/2;k++)
						{
											convo=convo+PSD_network[k]*PSD_signal[k];
											autoconv1 = autoconv1 + PSD_network[k]*PSD_network[k];
											autoconv2 = autoconv2 + PSD_signal[k]*PSD_signal[k];
										
						}
						
						
						MI_single_trial[q]=0.5*log2(1+fabs(convo/sqrt(autoconv1*autoconv2))/(1-fabs(convo/sqrt(autoconv1*autoconv2))));
						
						}//q trial
						
						MI[s1][s2]=0;
						for (int q=0;q<Q;q++)
						{
							MI[s1][s2] = MI[s1][s2] +1/((double) Q)*MI_single_trial[q];
						}
						Var_MI[s1][s2]=0;
						for (int q=0;q<Q;q++)
						{
							Var_MI[s1][s2] = Var_MI[s1][s2] +1/((double) Q)*(MI[s1][s2]-MI_single_trial[q])*(MI[s1][s2]-MI_single_trial[q]);
						}
						
								MEAN_MI_ACROSS_FREQUENCIES[s1]=MEAN_MI_ACROSS_FREQUENCIES[s1]+1/((double)S)*MI[s1][s2];	
								MEAN_MI_ACROSS_FREQUENCIES_var[s1] = 	MEAN_MI_ACROSS_FREQUENCIES_var[s1] +1/((double)S)*Var_MI[s1][s2];
					}//s2 loop
					
					
}//s1 loop
							
		for(int k=0;k<(int) T/2;k++)
					
					{
								Freq[k] = k*1/((double) 2* (T)*Dt);
					}
		
///////////////OUTPUT///////////////	
	
			
        ofstream outfile;
        
        
     	   outfile.open("ADM -  PSD.txt", ios::out);
         for(int k=5;k<(int) 50*(2*T)*Dt;k++)
         {
         	
         		  			outfile<<Freq[k]<<"	"<<PSD_network[k]<<"	"<<PSD_signal[k]<<endl;
    
         }  
        outfile.close(); 
     	 
		//
        
          outfile.open("ADM  -  membrane potential dynamics.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<u_e[0][t]<<"	"<<u_e[2][t]<<"	"<<u_e[3][t]<<"	"<<u_e[N-1][t]<<endl;
         	
         }  
        outfile.close(); 
        
		//       
        
          outfile.open("ADM  -  Mean activity and Firing rate.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	
         		  			outfile<<t<<"	"<<rate_stim[t]<<"	"<<rate_rec[t]<<"	"<<Input[node_1+1][t]<<endl;
         	
         }  
        outfile.close(); 
        
        //
    
        outfile.open("ADM  -  Connectivity.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
  		//
        
          outfile.open("ADM  -  Mean mutual Information vs CV.txt", ios::out);
         for(int s1=0;s1<S;s1++)
         {
         			
         		  			outfile<<CV[s1]<<"	"<<MEAN_MI_ACROSS_FREQUENCIES[s1]<<"	"<<MEAN_MI_ACROSS_FREQUENCIES_var[s1]<<endl;
			                         
         }  
         outfile.close(); 
         
         //
        
          outfile.open("ADM  -  Mutual Information.txt", ios::out);
         for(int s1=0;s1<S;s1++)
         {
         	for(int s2=0;s2<S;s2++)
         	{
			 
         		  			outfile<<CV[s1]<<"	"<<FREQUENCY[s2]<<"	"<<MI[s1][s2]<<"	"<<Var_MI[s1][s2]<<endl;
         	}
			                         
         }  
        outfile.close(); 
        
        //
        
         outfile.open("ADM  -  E Spikes.txt", ios::out);
         for(int t=0;t<T;t++)
         {
         	for (int n=0;n<N;n++)
         	{
         		if(X_e[n][t]>0.1)
         		{
         			
         		  			outfile<<t<<"	"<<n<<endl;
         		}
			}
                         
         }  
        outfile.close(); 
       
        
		cout<<"Simulations complete..."<<endl;
      
return 0;    
}

///////////////FUNCTION DEFINITIONS///////////////

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



double f(double u,double threshold)
{
	double output;
	output = 1/(1+exp(-beta*(u-threshold)));
	return output;
	
}

double f_prime(double u,double threshold)
{
	double output;
	output = 1/(1+exp(-beta*(u-threshold)))*1/(1+exp(-beta*(u-threshold)))*exp(-beta*(u-threshold))*beta;
	return output;
	
}





void sort(double Array[],int size)
{
	
	for (int q=0;q<size;q++)
	{
	
		for (int i=0;i<size;i++)
		{
				double temp;
				if (Array[i]>Array[i+1])
				{
				temp = Array[i+1];
				Array[i+1]= Array[i];
				Array[i] =temp;
				
				}
		
		}
	}
	
	
}



void shuffle_indices(int Array[], int size)
{
   int temporary;
   int randomNum;
   int last;
 
   for (last = size; last > 1; last=last-1)
   {
      randomNum = (int) floor(rand() % last);
      temporary = Array[randomNum];
      Array[randomNum] = Array[last - 1];
      Array[last - 1] = temporary;
   }
}

void shuffle(double Array[], int size)
{
   double temporary;
   int randomNum;
   int last;
 
   for (last = size; last > 1; last=last-1)
   {
      randomNum = rand() % last;
      temporary = Array[randomNum];
      Array[randomNum] = Array[last - 1];
      Array[last - 1] = temporary;
   }
}


/******************************************************************************/
void four1(double data[], unsigned long nn, int isign)
/*******************************************************************************
Replaces data[1..2*nn] by its discrete Fourier transform, if isign is input as
1; or replaces data[1..2*nn] by nn times its inverse discrete Fourier transform,
if isign is input as -1.  data is a complex array of length nn or, equivalently,
a real array of length 2*nn.  nn MUST be an integer power of 2 (this is not
checked for!).
*******************************************************************************/
{
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	double tempr,tempi;

	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2) { /* This is the bit-reversal section of the routine. */
		if (j > i) {
			SWAP(data[j],data[i]); /* Exchange the two complex numbers. */
			SWAP(data[j+1],data[i+1]);
		}
		m=nn;
		while (m >= 2 && j > m) {
			j -= m;
			m >>= 1;
		}
		j += m;
	}

	mmax=2;
	while (n > mmax) { /* Outer loop executed log2 nn times. */
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax); /* Initialize the trigonometric recurrence. */
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) { /* Here are the two nested inner loops. */
			for (i=m;i<=n;i+=istep) {
				j=i+mmax; /* This is the Danielson-Lanczos formula. */
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr; /* Trigonometric recurrence. */
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}
   
double coherence(double X1[], double X2[], int window)
{
		int Bins = T/((int) window);
		double Y1[Bins];
		double Y2[Bins];
		//create binned spike train representations
		for (int k=0;k<Bins;k++)
		{
			Y1[k]=0;
			Y2[k]=0;
			for (int s=0;s<window;s++)
			{
				if (X1[k*window+s]>0)
				{
					Y1[k]=1;
				}
				if (X2[k*window+s]>0)
				{
					Y2[k]=1;
				}
				
			}
					
		}
		//compute correlation
		double sumY1Y2=0;
		double sumY1=0;
		double sumY2=0;
		for (int k=0;k<Bins;k++)
		{
			sumY1Y2 = sumY1Y2+Y1[k]*Y2[k];
			sumY1 = sumY1+Y1[k];
			sumY2 = sumY2+Y2[k];
		}
		
		double output;
		if (sumY1>0&&sumY2>0)
		{
				output = sumY1Y2/sqrt(sumY1*sumY2);
		}
		else{ output =0;}		 
		return output;
		
		

}




