clear all; close all; clc
N = 20;
Vals = zeros(N,1);
Ts = zeros(N,1);

Tmax = 5000;
LM = 100;
CT = 100;
for test = 1:N
    test
    if test == N
        T = 1000000;
    else
        T = 10*test/1.0;
    end
    RT = 0;


    lamI = 1;
    lam = lamI*3.0;
    Amp = lamI;
    c1 = 0;
    Cpos1 = [];
    c2 = 0;
    Cpos2 = [];
    c3 = 0;
    Cpos3 = [];
    c4 = 0;
    Cpos4 = [];
    nc1 = 0;
    nc2 = 0;
    nc3 = 0;
    nc4 = 0;
    nf1 = 0;
    nf2 = 0;
    nf3 = 0;
    nf4 = 0;

    cn = 0;
    cmu1 = 0;
    cmu2 = 0;
    cmu3 = 0;
    cmu4 = 0;
    while RT < Tmax
        cn = cn+1;
        c = c1+c2+c3+c4;
        %%%%%%%%%%%%%%% For G1 %%%%%%%%%%%%%%%%
        Ri1 = (lamI + Amp*sin(2*pi/T*RT))*(CT-c)/(1.0*CT);        
        if (length(Cpos1)>0)
            if Cpos1(1) == 1;
                Ri1 = 0;
            end
        end
        if (c1>0)
            if (Cpos1(c1) < LM-1)
                Re1 = lam;
            else
                Re1 = lam/1.0;
            end
            for i = 1:c1-1
                if Cpos1(i+1)-Cpos1(i) > 1
                    Re1 = Re1+lam;
                end
            end
        else
            Re1 = 0;
        end

        %%%%%%%%%%%%%%% For G2 %%%%%%%%%%%%%%%%
        Ri2 = (lamI + Amp*sin(2*pi/T*RT))*(CT-c)/(1.0*CT);        
        if (length(Cpos2)>0)
            if Cpos2(1) == 1;
                Ri2 = 0;
            end
        end
        if (c2>0)
            if (Cpos2(c2) < LM-1)
                Re2 = lam;
            else
                Re2 = lam/1.0;
            end
            for i = 1:c2-1
                if Cpos2(i+1)-Cpos2(i) > 1
                    Re2 = Re2+lam;
                end
            end
        else
            Re2 = 0;
        end

                %%%%%%%%%%%%%%% For G3 %%%%%%%%%%%%%%%%
        Ri3 = (lamI + Amp*sin(2*pi/T*RT))*(CT-c)/(1.0*CT);        
        if (length(Cpos3)>0)
            if Cpos3(1) == 1;
                Ri3 = 0;
            end
        end
        if (c3>0)
            if (Cpos3(c3) < LM-1)
                Re3 = lam;
            else
                Re3 = lam/1.0;
            end
            for i = 1:c3-1
                if Cpos3(i+1)-Cpos3(i) > 1
                    Re3 = Re3+lam;
                end
            end
        else
            Re3 = 0;
        end

          %%%%%%%%%%%%%%% For G3 %%%%%%%%%%%%%%%%
        Ri4 = (lamI + Amp*sin(2*pi/T*RT))*(CT-c)/(1.0*CT);        
        if (length(Cpos4)>0)
            if Cpos4(1) == 1;
                Ri4 = 0;
            end
        end
        if (c4>0)
            if (Cpos4(c4) < LM-1)
                Re4 = lam;
            else
                Re4 = lam/1.0;
            end
            for i = 1:c4-1
                if Cpos4(i+1)-Cpos4(i) > 1
                    Re4 = Re4+lam;
                end
            end
        else
            Re4 = 0;
        end



        Tot = Ri1+Ri2+Ri3+Re1+Re2+Re3+Ri4+Re4;
        RT = RT - log(rand)/Tot;
        A = rand;
        if (A < Ri1/Tot)
            c1 = c1+1;
            nc1 = nc1+1;
            Cpos1 = [1; Cpos1];
        elseif (A < (Ri1+Re1)/Tot)
            Rtot = Ri1/Tot;
            ce = 0;
            while Rtot < A
                ce = ce+1;
                if ce < c1
                    if (Cpos1(ce+1)-Cpos1(ce) > 1)
                        Rtot = Rtot + lam/Tot;
                    end
                else
                    Rtot = A;
                end
            end
            Cpos1(ce) = Cpos1(ce)+1;
            if (Cpos1(ce) == LM)
                Cpos1 = Cpos1(1:end-1);
                c1 = c1-1;
                nf1 = nf1+1;
            end
        elseif (A < (Ri1+Re1+Ri2)/Tot)
            c2 = c2+1;
            nc2 = nc2+1;
            Cpos2 = [1; Cpos2];
        elseif A < (Ri1+Re1+Ri2+Re2)/Tot
            Rtot = (Ri1+Re1+Ri2)/Tot;
            ce = 0;
            while Rtot < A
                ce = ce+1;
                if ce < c2
                    if (Cpos2(ce+1)-Cpos2(ce) > 1)
                        Rtot = Rtot + lam/Tot;
                    end
                else
                    Rtot = A;
                end
            end
            Cpos2(ce) = Cpos2(ce)+1;
            if (Cpos2(ce) == LM)
                Cpos2 = Cpos2(1:end-1);
                c2 = c2-1;
                nf2 = nf2+1;
            end
        elseif (A < (Ri1+Re1+Ri2+Re2+Ri3)/Tot)
            c3 = c3+1;
            nc3 = nc3+1;
            Cpos3 = [1; Cpos3];
        elseif A < (Ri1+Re1+Ri2+Re2+Ri3+Re3)/Tot
            Rtot = (Ri1+Re1+Ri2+Re2+Ri3)/Tot;
            ce = 0;
            while Rtot < A
                ce = ce+1;
                if ce < c3
                    if (Cpos3(ce+1)-Cpos3(ce) > 1)
                        Rtot = Rtot + lam/Tot;
                    end
                else
                    Rtot = A;
                end
            end
            Cpos3(ce) = Cpos3(ce)+1;
            if (Cpos3(ce) == LM)
                Cpos3 = Cpos3(1:end-1);
                c3 = c3-1;
                nf3 = nf3+1;
            end

        elseif (A < (Ri1+Re1+Ri2+Re2+Ri3+Re3+Ri4)/Tot)
            c4 = c4+1;
            nc4 = nc4+1;
            Cpos4 = [1; Cpos4];
        elseif A < (Ri1+Re1+Ri2+Re2+Ri3+Re3+Ri4+Re4)/Tot
            Rtot = (Ri1+Re1+Ri2+Re2+Ri3+Re3+Ri4)/Tot;
            ce = 0;
            while Rtot < A
                ce = ce+1;
                if ce < c4
                    if (Cpos4(ce+1)-Cpos4(ce) > 1)
                        Rtot = Rtot + lam/Tot;
                    end
                else
                    Rtot = A;
                end
            end
            Cpos4(ce) = Cpos4(ce)+1;
            if (Cpos4(ce) == LM)
                Cpos4 = Cpos4(1:end-1);
                c4 = c4-1;
                nf4 = nf4+1;
            end

        end

        cmu1 = (cmu1*(cn-1)+c1)/cn;
        cmu2 = (cmu2*(cn-1)+c2)/cn;
        cmu3 = (cmu3*(cn-1)+c3)/cn;
        cmu4 = (cmu4*(cn-1)+c4)/cn;
    end

    Vals(test) = (nf1+nf2+nf3+nf4);
    Ts(test) = T;
end

plot(Ts(1:N-1),Vals(1:N-1)/Vals(N),'LineWidth',3); goodplot


