clear all; close all; clc
T0 = 500; A0 = 0.5;
for test = 1:2

    RT = 0; Tmax = 10000; ts = 1; dt = 0.01; Saf = zeros(round(Tmax/ts),3); click = 0;

    gam1 = 1.0; del1 = .1; Gm1 = 0.1; Dl1 = .1;
    gam2 = 1.0; del2 = .1; Gm2 = 0.2; Dl2 = .1;
    gam3 = 1.0; del3 = .1; Gm3 = 0.3; Dl3 = .1;
    
    cpp = 1; h = 8;

    psi = 0.01; eps = 0.0001;
    m1 = 20; P1 = 20; P0 = 1.2*10^4; 
    m2 = 20; P2 = 20; 
    m3 = 20; P3 = 20; 
    CP = psi*P1^h/(P0^h + P1^h)/eps;
    if (test == 1)
        T = T0; 
        A = A0;
    elseif (test == 2)
        T = 2*T0; 
        A = A0;
    elseif (test == 3)
        T = T0; A = 2*A0;
    end
    c15 = 0;
    while RT < Tmax
        RT = RT + dt;
        P = 1 + A*cos(2*pi/T*RT);

        dm1 = gam1*P - del1*m1 + 0.1*randn;
        dP1 = Gm1*m1/(P1+P2+P3) - Dl1*P1*P3 + randn;
        dm2 = gam2*P - del2*m2 + 0.1*randn;
        dP2 = Gm2*m2/(P1+P2+P3) - Dl2*P2*P1 + randn;
        dm3 = gam3*P  - del3*m3 + 0.1*randn;
        dP3 = Gm3*m3/(P1+P2+P3) - Dl3*P3*P2 + randn;
        
        m1 = m1 + dt*dm1;
        m2 = m2 + dt*dm2;
        m3 = m3 + dt*dm3;
        P1 = abs(P1 + dt*dP1);
        P2 = abs(P2 + dt*dP2);
        P3 = abs(P3 + dt*dP3);

        
        if (RT > click*ts)
            click = click+1;
            Saf(click,1) = RT;
            Saf(click,2) = P1;
            Saf(click,5) = P2;
           Saf(click,3) = CP;
           Saf(click,4) = P;

        end
    end
    c15
    t = Saf(1:click,1);
    p1 = Saf(1:click,2);
    p2 = Saf(1:click,5);
    p = Saf(round(click/2):click,4);
    cp = Saf(end/2:end,3);
    
   plot(t,p1); hold on
   %plot(t,p2)

    %plot(CC(1:end-T)/mean(CC),CC(T+1:end)/mean(CC)); hold on

end

