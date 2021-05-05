
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

#define T 1000*1000//Time of simulation in msec
#define T_start 2000 //start of simulation (ms)
#define T_window 1000//size of time window (ms
#define Delta_T_rate 50 //time window for rate calculations (ms)
#define Sample 10//number of sampled non-zero connections
#define N 100//number of neurons
#define dt 0.1 //time step for Euler Maruyama integration
#define Dt 0.001// Value in sec of 1 dt e.g. Dt=0.001 means that 1 dt equals 1 msec
#define w_o 0.21//synaptic connectivity weight (criticality) 
#define temporal_window 10 //Temporal window for spike coherence measures
#define beta 25 //activation function non-linear gain
#define rho 0.15//connection probability
#define epsilon 0.3//myelination rate
#define Omega 10.0 //network spatial scale in mm
#define h 0.10//activation function threshold (inflexion point) 
#define I_o 0.1 //input amplitude (a.u.)
#define c_min 0.1 //minimal conduction velocity in m/s
#define c_max 100//maximal conduction velocity in m/s

#define alpha_retraction 0.0000001 //rate constant for ADM retraction
#define alpha_formation 0.0000001//rate constant for ADM formation 

#define t_stim_on T/4 //stimulation onset time
#define t_stim_off 2*T/4//stimulation offset time

#define T_plastic_on 0 //ADM onset time (default to initial simulation starting time 0)
#define T_plastic_off T //ADM offset time (default to full simulation time T)

#define t_delay_statistics_initial 20 // time for baseline delay distribution calculations
#define t_delay_statistics_final T-T_window-100 // time for stabilized delay distribution calculations
#define delay_graining 50// bins for  delay distribution calculations
#define speed_graining 50// bins for  velocity distribution calculations
#define max_delay 500// maximum for delay distribution calculations
#define max_speed 40// maximum for velocity distribution calculations

//////////////////////FUNCTION DEFINITIONS///////////////////////
double f(double u, double threshold);//neuron response function 
float ran2(long *idum);//random number generator - initial cond.

//////////////////////OBJECT DEFINITIONS///////////////////////
double u_e[N][T_window];//Dynamics (au)
double mean_u_e[T];//mean network activity
double mean_rate[T];//mean firing rate 
double Sigma[T];// variance in time
double x[N];//neuron position along x axis
double y[N];//neuron position along y axis
double z[N];//neuron position along z axis

double c[N][N];//conduction velocity matrix
double w[N][N];//synaptic weight matrix
double l[N][N];//axonal tract lengths matrix
int tau[N][N];//delays matrix
int tau_o[N][N];//baseline delay matrix
double X_e[N][T_window]; //spiking activity
double r_e[N][T_window];//rates
double CV_in_time[Sample][T];//conduction velocity in time for sampled non-zero connections
double Delay_in_time[Sample][T];//conduction delay in time for sampled non-zero connections

double  Delay_distribution_final[delay_graining];//array for delay distribution - stabilized
double  Delay_distribution_initial[delay_graining];//array for delay distribution - baseline
double Speed_distribution_final[speed_graining];//array for velocity distribution - stabilized
double Speed_distribution_initial[speed_graining];//array for velocity distribution - baseline
using namespace std;


int main()
{
	
srand (time(NULL));


//////////////////////DEFINING CONNECTIVITY, BASELINE DELAYS AND CONDUCTION VELOCITIES///////////////////////

					for (int i=0;i<N;i++)
					{
						long seed1=(long)21+i*187198+56*i+12*1+1;
						long seed2=(long)21+i*56+56*i+11*1+1;
						long seed3=(long)21+i*789+56*i+8*1+1;
		                          	  					
						x[i] = ran2(&seed1)*Omega;
						y[i] = ran2(&seed2)*Omega;
						z[i] = ran2(&seed3)*Omega;

					}
					double mean_w=0;
					for(int i=0;i<N;i++)
					{
						for (int j=0;j<N;j++)
						{
							l[i][j]= fabs(sqrt((x[i]-x[j])*(x[i]-x[j])+(y[i]-y[j])*(y[i]-y[j])+(z[i]-z[j])*(z[i]-z[j])));
							c[i][j]=c_min;//c_o;
							if(i==j)
							{
								w[i][j]=0;
							}
							else
							{
								
									long seed=(long)21+i*j*187198+56*j+12*1+1;
		                          	  
									if(ran2(&seed)<rho)
									{
											w[i][j]=w_o/rho;
											
									}
								
									else
									{
										w[i][j]=0;
									}
								
							}
							mean_w=mean_w+1/((double)N*N)*w[i][j];
		 					tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
		 					tau_o[i][j]=tau[i][j];
		 					
						}
						
					}
				
				
//////////////////////SELECTING SAMPLE TRAJECTORIES FOR PLOTTING PURPOSES///////////////////////				
					int sampled=0;
					int sampled_index_i[Sample];
					int sampled_index_j[Sample];
				
					for (int i=0;i<N;i++)
					{
							for (int j=0;j<N;j++)
							{
											
											if(w[i][j]>0&&sampled<Sample)
											{
												sampled_index_i[sampled] =i;
												sampled_index_j[sampled] =j;
												sampled++;
												
											}

						}				
					}
				
				
				
//////////////////////INITIAL CONDITIONS///////////////////////
					
				
		        	for (int s=0;s<T_window-1;s++)
					{
					
						for(int n=0;n<N;n++)
		        		{
		                            long d=rand()%105;      
		                            long noisyseed1=(long)21*s+n*n*187198+56*d+12*1+1;
		                            long noisyseed2=(long)69*s*n+11*n+s+45*d+2*1+1; 
		            				
         							u_e[n][s]=sqrt(-2*log(ran2(&noisyseed1)))*cos(2*3.14159*ran2(&noisyseed2));	
         							r_e[n][s] = 0;
						}
						for(int i=0;i<Sample;i++)
						{
							CV_in_time[i][s]=c[sampled_index_i[i]][sampled_index_j[i]];
							Delay_in_time[i][s]=tau[sampled_index_i[i]][sampled_index_j[i]];
						}
						
						
					
		           }
		               
					
					

					for (int s=0;s<T-T_window;s++)
					{
						cout<<s<<endl;
						
//////////////////////SLIDING TIME WINDOW///////////////////////
						for (int t=0;t<T_window-1;t++)
						{
							for (int i=0;i<N;i++)
							{
									u_e[i][t] = u_e[i][t+1]	;
									X_e[i][t] = X_e[i][t+1];
									r_e[i][t] = r_e[i][t+1];
								
								
							}	
						}
						
						
	//////////////////////EULER MARUYAMA INTEGRATION SCHEME///////////////////////
						double I=0;
						if (s>t_stim_on&&s<t_stim_off)
						{
							I=I_o;
						}
					
							
						double sum_ee; double sum_theo;
						for (int i=0;i<N;i++)
						{
						
							sum_ee=0;
							sum_theo=0;
							for (int j=0;j<N;j++)
							{								
									sum_ee=sum_ee+1/((double)N)*w[i][j]*X_e[j][T_window-1-tau[i][j]];
									sum_theo  =sum_theo + 1/((double)N)*w[i][j]*f(u_e[j][T_window-1-tau[i][j]],h);
	 								
							}

			
							u_e[i][T_window-1]=u_e[i][T_window-1]+dt*(-1*u_e[i][T_window-1]+sum_ee+I);
					
										
										   // Poisson neurons
										 
										  long seed_q =rand()%101;  
						   				  long seed_e = (long) 12+i*i+s+15+s*seed_q+56*i+i+99*s+2;              
			                              double p_e =ran2(&seed_e);
			                              double p_fire_e = ( 1-exp(-f(u_e[i][T_window-1],h)*dt)); 
			                              if (p_e<p_fire_e)//spike occurs
			                               {                   
			                                               X_e[i][T_window-1] =1/dt; 
			                                               
			        										
			                               }
			                               else//no spike
			                               {
			                                               X_e[i][T_window-1]=0; 
			                                               			                               
										   }
										   double sum=0;
										   for (int t_index=0;t_index<Delta_T_rate;t_index++)
										   {
										   			sum=sum+1/((double)Delta_T_rate*Dt)*X_e[i][T_window-1-t_index]*dt;
										   			
										   }
										   r_e[i][T_window-1]=sum;

			                              
						}

//////////////////////ADM RULE///////////////////////
						
						if(s>T_plastic_on&&s<T_plastic_off)
						{
								
						for (int i=0;i<N;i++)
						{
							for (int j=0;j<N;j++)
							{
								if (i==j)
								{
									c[i][j]=c_min;//
								}
								else
								{
									
									c[i][j]=c[i][j]+dt*((c_min-c[i][j])*alpha_retraction+epsilon*X_e[j][T_window-1]*(l[i][j]/(c[i][j]))*alpha_formation);
									
								}
								 if(c[i][j]<0)
								 {
								 	c[i][j]=c_min;
								 }
								 else
								 {
								 	
								 }
								tau[i][j] = (int) floor(l[i][j]/(c[i][j]));
								
							}
						}
						}
						
						for(int i=0;i<Sample;i++)
						{
							CV_in_time[i][T_window+s]=c[sampled_index_i[i]][sampled_index_j[i]];
							Delay_in_time[i][T_window+s]=tau[sampled_index_i[i]][sampled_index_j[i]];
						}
						

						for (int i=0;i<N;i++)
						{
									mean_u_e[T_window+ s] = mean_u_e[T_window+ s]+1/((double)N)*u_e[i][T_window-1];
									mean_rate[T_window+ s] = mean_rate[T_window+ s]+1/((double)N)*r_e[i][T_window-1];
						}
						
//////////////////////BASELINE DELAY STATISTICS///////////////////////					
					
						
						if (s==t_delay_statistics_initial)
						{
						 
													for (int ss=0;ss<delay_graining;ss++)
													{
																			Delay_distribution_initial[ss]=0;
													}
																					
													for (int tt=0;tt<delay_graining;tt++)
													{
																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																			double dttau = (max_delay)/((double) delay_graining);
																							
																			for (int i=0;i<N;i++)
																			{
																				for (int j=0;j<N;j++)
																				{
																					
																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																									double dttau = (max_delay)/((double) delay_graining);
																									double delay  = (double) tau[i][j];			
																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
																									{
																												Delay_distribution_initial[tt] = Delay_distribution_initial[tt]+w[i][j];
													
																									}	
																				}
																			}
													}
									
														
										for (int ss=0;ss<speed_graining;ss++)
										{
													Speed_distribution_initial[ss]=0;
										}
										for (int tt=0;tt<speed_graining;tt++)
										{
																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																double dcc = (max_speed)/((double) speed_graining);
																
																for (int i=0;i<N;i++)
																{
													
																		for (int j=0;j<N;j++)
																		{
																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																			double dcc = (max_speed)/((double) speed_graining);
																			double speed = (double) c[i][j];			
																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
																			{
																						Speed_distribution_initial[tt] = Speed_distribution_initial[tt]+w[i][j];
																						
																						
																			}	
																		}
															
																}
																
										}
										
						}//initial distribution loop
						
//////////////////////STABILIZED DELAY STATISTICS///////////////////////
								
						if (s==t_delay_statistics_final)
						{
								
													for (int ss=0;ss<delay_graining;ss++)
													{
																			Delay_distribution_final[ss]=0;
													}
																					
													for (int tt=0;tt<delay_graining;tt++)
													{
																			double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																			double dttau = (max_delay)/((double) delay_graining);
																							
																			for (int i=0;i<N;i++)
																			{
																				for (int j=0;j<N;j++)
																				{
																					
																									double ttau = 0 + (max_delay)/((double) delay_graining)*tt;
																									double dttau = (max_delay)/((double) delay_graining);
																									double delay  = (double) tau[i][j];			
																									if (((double)delay>ttau-dttau*0.5)&&((double)delay<ttau+dttau*0.5))
																									{
																												Delay_distribution_final[tt] = Delay_distribution_final[tt]+w[i][j];
													
																									}	
																				}
																			}
													}
									
														
										for (int ss=0;ss<speed_graining;ss++)
										{
													Speed_distribution_final[ss]=0;
										}
										for (int tt=0;tt<speed_graining;tt++)
										{
																double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																double dcc = (max_speed)/((double) speed_graining);
																
																for (int i=0;i<N;i++)
																{
													
																		for (int j=0;j<N;j++)
																		{
																			double tc = 0 + (max_speed)/((double) speed_graining)*tt;
																			double dcc = (max_speed)/((double) speed_graining);
																			double speed = (double) c[i][j];			
																			if (((double)speed>tc-dcc*0.5)&&((double)speed<tc+dcc*0.5))
																			{
																						Speed_distribution_final[tt] = Speed_distribution_final[tt]+w[i][j];
																						
																						
																			}	
																		}
															
																}
																
										}
										
								}// final distribution loop
						
						
												
						
					}// s loop




			
/////////////////////OUTPUT ROUTINE/////////////////////////////			
			
			
        ofstream outfile;
        

        outfile.open("ADM - Sampled Conduction velocities in time.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<CV_in_time[0][t]<<"	"<<CV_in_time[1][t]<<"	"<<CV_in_time[2][t]<<"	"<<CV_in_time[3][t]<<"	"<<CV_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
         outfile.open("ADM - Sampled conduction delays in time.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<Delay_in_time[0][t]<<"	"<<Delay_in_time[1][t]<<"	"<<Delay_in_time[2][t]<<"	"<<Delay_in_time[3][t]<<"	"<<Delay_in_time[4][t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        
          outfile.open("ADM -   Mean activity and Firing rate time series.txt", ios::out);
         for(int t=0;t<T-T_window;t=t+200)
         {
         	
         		  			outfile<<t<<"	"<<mean_u_e[t]<<"	"<<mean_rate[t]<<"	"<<Sigma[t]<<endl;
         	
                         
         }  
        outfile.close(); 
        
        
        
      	outfile.open("ADM  - Baseline and stabilized delay distribution.txt", ios::out);
		
		
		
				for (int k=0;k<delay_graining;k++)
				{
				 
				 	
		         					outfile<<(max_delay)/((double) delay_graining)*k<<"	"<<Delay_distribution_initial[k]<<"	"<<Delay_distribution_final[k]<<endl;
		         		 
		     	}
		     	
   
        outfile.close(); 
        
        	outfile.open("ADM  -  Baseline and stabilized Conduction Velocity distribution.txt", ios::out);
		
		
		
				for (int k=0;k<speed_graining;k++)
				{
				 
				 	
		         					outfile<<(max_speed)/((double) speed_graining)*k<<"	"<<Speed_distribution_initial[k]<<"	"<<Speed_distribution_final[k]<<endl;
		         		 
		     	}
		     	
   
        outfile.close(); 
        
       
        

        outfile.open("ADM  -  Connectivity Matrix.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<i<<"	"<<j<<"	"<<c[i][j]<<"	"<<tau[i][j]<<"	"<<w[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
      outfile.open("ADM - axonal lengths.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<l[i][j]<<endl;
         		
			}
                         
         }  
        outfile.close(); 
        
         outfile.open("ADM - conduction velocities.txt", ios::out);
         for(int i=0;i<N;i++)
         {
         	for (int j=0;j<N;j++)
         	{
         	
         			
         		  			outfile<<c[i][j]<<endl;
         		
			}
                         
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



double f(double u,double threshold)
{
	double output;
	output = 1/(1+exp(-beta*(u-threshold)));
	return output;
	
}









