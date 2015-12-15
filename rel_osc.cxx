#include <iostream>
#include <cmath>

using namespace std;

void calculation(const double y1, const double y2, double f[2]);


int main() {
      const double p_0 = 0.1;
      
      double y[2];
      double yt;
      
      const double dt = 0.1;
      
      double k1[2];
      double k2[2];
      double k3[2];
      double k4[2];
      double t=0;
      
      for(double p = p_0; p <= 5.0; p += 0.1) {
	  
	  y[0] = p; 	// p
	  y[1] = 0.0; 	// p_dot
	 
	  int i = 0;
	  t=0;
	  while(true) {
	      yt = y[1];
	      calculation(y[0], y[1], k1);
	      calculation(y[0] + 0.5*dt*k1[0], y[1] + 0.5*dt*k1[1], k2);
	      calculation(y[0] + 0.5*dt*k2[0], y[1] + 0.5*dt*k2[1], k3);
	      calculation(y[0] + dt*k3[0], y[1] + dt*k3[1], k4);
	      
	      t+=dt;
	      
	      for(int j = 0; j < 2; j++) {
		  y[j] += dt/6.0 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]);
	      }
	     
	      i++;
	      
	      
	      if (yt > 0 && y[1] < 0) { 
		  //cout << i*dt << "\t" << p << endl;
		  break; 	
	      }
	      
	  }
	  
	  double thetaLinks = 0;
	  double thetaRechts = 1;
	  double theta = 0.5;
	  
	  double b1, b2, b3, b4;
	  double error = 1;
	  
	  while(abs(error) > 1e-5) {
	      b1 = theta - 1.5*theta*theta + 2*theta*theta*theta/3.0;
	      b2 = theta*theta - 2*theta*theta*theta/3.0;
	      b3 = b2;
	      b4 = - 0.5*theta*theta + 2*theta*theta*theta/3.0;
	      
	      error = yt + dt*(b1*k1[1] + b2*k2[1] + b3*k3[1] + b4*k4[1]);
	      
	      if(error > 0) 
		thetaLinks = theta;
	      else
		thetaRechts = theta;
	      
	      theta = 0.5*(thetaLinks + thetaRechts);
	  }
	  
	  cout << p << "\t" << (i-1)*dt + theta*dt << endl;
      }
      
      return 0;
}
      
      
void calculation(const double y1, const double y2, double f[2]) {
      f[0] = y2;
      f[1] = - y1/sqrt(1+y1*y1);
}