#include "ConstantHeaderDeterministic.h"
#include "RungeKuttaFunctions.h"
#include <math.h> 
#include <stdlib.h> 
#include <time.h>
#include <fstream>
#include <ctype.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <random>
using namespace std;

int main() {

  for (unsigned int iA = 0; iA < 50; iA++){
    cout << iA << endl;
    Vol = 2.0*pow(10,-14);	NA = 6.02*pow(10,23);	Cal = NA*Vol*pow(10,-6);
    Ntot = 1*Cal; KN = 0.029*Cal;	kt = 1.03/Cal; 	a = 1.05/Cal; KI = 0.035*Cal;	IKKtot = 2*Cal;

    double V1,V2,V3,V4,V5;

    for (unsigned int itest = 40; itest < 80; itest++){
      double ETA,Amp1,Amp2,omega1,omega2,pha1;
      ETA = 0.1;
      Amp1 = 0.005*iA; Amp2 = 0.1; 
      omega1 = 1.0/(40+itest); omega2 = 1.0/50; 
      pha1 = 0;

      ostringstream newname;
      newname << "Datafiles/DoublePendulum_A"; newname << iA; newname << "_"; newname << itest; newname << ".txt";
      std::ofstream newfile (newname.str().c_str());

      N = 1500;    IRna = 4979.9; I = 35442;    IKKa = 2284.81;  IKKi = 19072.8; RT = 16.189;
      int cmax = 50; 
      int c = 0; int cinit = 0;  double injtime = 0; click = 0;     double I0 = 0; double N0 = 0.0;
      double Poin1 = RT;
      double Nin = N; double IRnain = IRna;     double Iin = I; double IKKain = IKKa;     double IKKiin = IKKi; double RTin = RT; 
      double Del0, Del1, Del2, Del3, Del4, Del5,mui,sigi;
      Del0 = 0;      
      double rf,wf; double con = 1.0;
      rf = 0.1; wf = 0;
      ////////////////////////////////  RunScript //////////////////////////////
      double Ti [cmax];      
      while (c < cmax) {
	RT = RT + dt;          it++;

	TNF = 0.5 + Amp1*sin(2*3.141592*omega1*RT + pha1);
	ETA = 0.5 + Amp2*sin(2*3.141592*omega2*RT);
	
	//////////////////////////////////// ACTUAL UPDATES ////////////////////////////////////////////////////////////////
	ka1 = f1(N,I)*dt;	  ja1 = f2(N,IRna,ETA)*dt;	  la1 = f3(N,IRna,I,IKKa)*dt;	                                    
	ma1 = f4(IKKa,IKKi,TNF)*dt;	  na1 = f5(IKKa, IKKi,TNF)*dt;                                                      //
	
	ka2 = f1(N + 0.5*ka1, I + 0.5*la1)*dt;	  ja2 = f2(N + 0.5*ka1, IRna + 0.5*ja1,ETA)*dt;	                    // 
	la2 = f3(N + 0.5*ka1, IRna + 0.5*ja1, I + 0.5*la1, IKKa + 0.5*ma1)*dt;
	ma2 = f4(IKKa +0.5*ma1, IKKi + 0.5*na1, TNF)*dt;	  na2 = f5(IKKa +0.5*ma1, IKKi + 0.5*na1, TNF)*dt;          //
	
	ka3 = f1(N + 0.5*ka2, I + 0.5*la2)*dt;	  ja3 = f2(N + 0.5*ka2, IRna + 0.5*ja2,ETA)*dt;                         //
	la3 = f3(N + 0.5*ka2, IRna + 0.5*ja2, I + 0.5*la2, IKKa + 0.5*ma2)*dt;	  
	ma3 = f4(IKKa +0.5*ma2, IKKi + 0.5*na2, TNF)*dt;	  na3 = f5(IKKa +0.5*ma2, IKKi + 0.5*na2, TNF)*dt;          //
	
	ka4 = f1(N + ka3, I + la3)*dt;	  ja4 = f2(N + ka3, IRna + ja3,ETA)*dt;	                                    //
	la4 = f3(N + ka3, IRna + ja3, I + la3, IKKa + ma3)*dt;
	ma4 = f4(IKKa + ma3, IKKi + na3, TNF)*dt;	  na4 = f5(IKKa + ma3, IKKi + na3, TNF)*dt;                         //
      
	N = N + 1./6*(ka1 + 2.*ka2 + 2.*ka3 + ka4);	  IRna = IRna + 1./6*(ja1 + 2.*ja2 + 2.*ja3 + ja4);	            //
	I = I + 1./6*(la1 + 2.*la2 + 2.*la3 + la4);
	IKKa = IKKa + 1./6*(ma1 + 2.*ma2 + 2.*ma3 + ma4);	  IKKi = IKKi + 1./6*(na1 + 2.*na2 + 2.*na3 + na4);         // 
      
	if (N < 0){	    N = 0;}	  if (IRna < 0){	    IRna = 0;}	  if (I < 0){	    I = 0;}	            //
	if (IKKa < 0){	    IKKa = 0;}	  if (IKKi < 0){	    IKKi = 0;}
	///////////////////////////////////    END    //////////////////////////////////////////////////////////////////////
	
	if (RT > 0.1*click){
	  click++;
	  newfile << N << "\t" << TNF << "\t" << RT << "\n";

	}	

	if (N > 1500.0 && N0 <= 1500.0 && IRna < 8000.0 && (RT - Poin1) > 10.0){
	  Del0 = RT - Poin1;	  Poin1 = RT;	    c++;  Ti[c] = Del0;   	  
	  // cout << Del0 << endl;	  
	}
	N0 = N;
	I0 = I;
      }

    }    
  }
}
