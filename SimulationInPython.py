from __future__ import division
import numpy as np
import pandas as pd
from tqdm import tqdm
from numba import njit, prange
import matplotlib.pyplot as plt
import time

@njit
def Single_run_numba(N):
    SIRfile = []
    Vol = 2.0*10**(-14);
    NA = 6.02*10**(23);
    Cal = NA*Vol*10**(-6);

    ########### Contants for N %%%%%%%%%%%%%%%                                                                     
    kNin = 5.4;
    Ntot = 1*Cal;
    KN = 0.029*Cal;
    kIin = 0.018;
    ########### Contants for IRNA %%%%%%%%%%%%%%%                                                                
    kt = 1.03/Cal;
    ym = 0.017;
    ########### Contants for I %%%%%%%%%%%%%%%                                                                      
    ktl = 0.24;
    a = 1.05/Cal;
    KI = 0.035*Cal;
    ########### Contants for IKKa %%%%%%%%%%%%%%                                                                  
    ka = 0.24;
    IKKtot = 2*Cal;
    ki = 0.18;
    ########## Contants for IKKi %%%%%%%%%%%%%%%                                                                   
    kp = 0.036;
    kA20 = 0.0018;
    A20 = 0.0026;
   
    RT = 0; click = 0; ts = 0.1; Tmax = 10000; dt = 0.001
    TNF = 0.5
    Amp = 0.3
    x = np.ones(5)
    dx = np.ones(5)
    xtmp = np.ones(4)
    RKV = np.ones((5,4))
    S = np.ones(4); S[0] = 0; S[1] = 0.5; S[2] = 0.5; 
    while(RT < Tmax):
        RT += dt
        TNF = 0.5 + Amp*np.sin(2*3.141592/40.0*RT)

        for i in range(4):
            xtmp = x + S[i]*dx
            dx[0] = dt*(kNin*(Ntot-xtmp[0])*KI/(KI+xtmp[2]) - kIin*xtmp[2]*xtmp[0]/(KN+xtmp[0]));
            dx[1] = dt*(kt*xtmp[0]**2 - ym*xtmp[1]);
            dx[2] = dt*(ktl*xtmp[1] - a*xtmp[3]*(Ntot-xtmp[0])*x[2]/(KI+xtmp[2]));
            dx[3] = dt*(ka*TNF*(IKKtot-xtmp[4]-xtmp[3]) - ki*xtmp[3]);
            dx[4] = dt*(ki*xtmp[3] - kp*xtmp[4]*kA20/(kA20+A20*TNF));

            RKV[:,i] = dx

        x += 1.0/6.0*(RKV[:,0] + 2*RKV[:,1] + 2*RKV[:,2] + RKV[:,3])
        
        if (RT > ts*click):                        
            SirTmp = np.zeros(7)
            SirTmp[0] = RT;
            SirTmp[1] = x[0];
            SirTmp[2] = x[1];
            SirTmp[3] = x[2];
            SirTmp[4] = x[3];
            SirTmp[5] = x[4];
            SirTmp[6] = TNF;
            SIRfile.append(SirTmp)
            click += 1

    return SIRfile



@njit(parallel = True)
def multiple_loops(N_loops,N):
    SIRfiles = []
    for i in prange(N_loops):
        SIRfile = Single_run_numba(N)
        SIRfiles.append(SIRfile)
    return SIRfiles

N_loops = 1




for test2 in range (1):
    Pers = []
    for test in range(1):
        start = time.time()
        N = 40;
        SIRfiles = multiple_loops(N_loops,N)
        
        end = time.time()
        print(f"\nElapsed (with compilation) = {end - start:.2f}",test)
        
        for i0 in range(len(SIRfiles)):
            XN = SIRfiles[i0]
            tt = np.zeros((len(XN),7))
            for j2 in range(7):
                for j in range(len(XN)):
                    tt[j][j2] = XN[j][j2]
                if (j2 == 1):
                    plt.plot(tt[:,0],tt[:,j2])
            plt.show()            
#            np.savetxt("Data/Figure2Av3_ph%s_%s.txt"%(test,i0),TT)
