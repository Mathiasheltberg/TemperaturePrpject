clear all; close all; clc

x = 1; y = 1; z = 1;

c1 = 2; c2 = 1; c3 = 2;

dt = 0.01; RT = 0; Tmax = 100; click = 0; ts = 0.1; SavF = zeros(round(Tmax/ts),3);
lam = 0;
while RT < Tmax
    RT = RT + dt;
    
    dx = (-c1*z-lam*x);
    dy = (c2*x-lam*y);
    dz = (c3*y-lam*z);
    x = x + dt*dx;
    y = y + dt*dy;
    z = z + dt*dz;

    if (RT > click*ts)
        click = click+1;
        SavF(click,1) = RT;
        SavF(click,2) = x;
        SavF(click,3) = z;
    end
end
plot(SavF(:,1),SavF(:,2)./exp(0.8*SavF(:,1))); goodplot

figure; 

x = linspace(0.1,3,100000)'; k = 0.01; fx = x./(1+k*x);fpx = 1./((1+k*x).^2); hold on; plot(x,fx); plot(x,fpx)
figure

b = 2.0;
b2 = 2*b;
yx = 1./(b*fx);
yy = x.^2;
yx2 = 1./(b2*fx);

idx = find(yy(1:end-1) < yx(1:end-1) & yy(2:end) > yx(2:end));
plot(x,yx); hold on;
plot(x(idx),yy(idx),'o')
plot(x,yy); goodplot



idx2 = find(yy(1:end-1) < yx2(1:end-1) & yy(2:end) > yx2(2:end));
plot(x,yx2); hold on;
plot(x(idx2),yy(idx2),'*')
plot(x(idx2),yx(idx2),'*')
DelX = (x(idx2)-x(idx));

(yy(idx2)-yy(idx))/(b2-b)

plot(x,yy); goodplot

yyp = 2*yy(idx); 
yxp = -1/(b*yy(idx)^2);
dbyx = -1./(b^2*x(idx)./(1+k*x(idx)));

Dels = [dbyx/(yyp-yxp) (x(idx2)-x(idx))/(b2-b) DelX*2*x(idx) (yy(idx2)-yy(idx)) dbyx]
plot([x(idx) x(idx2)],[yx(idx) yx(idx)+DelX*2*x(idx)],'-*')

plot(x,yx-yy,'--k')
plot(x,zeros(length(x),1),'-k')


figure
x = linspace(0.01,2,100000)';
zz = x./(.01+1*x);
yx = 1./(10*fx);
y0 = x;
yy = x.^2;
plot(x,yx-yy,'--r');hold on;
plot(x,yx-zz,'--b')
plot(x,yx-y0,'--c')
plot(x,yx,'g'); 
plot(x,yy,'r'); goodplot
plot(x,y0,'c'); goodplot
plot(x,zz,'b'); goodplot
plot(x,zeros(length(x),1),'-k')


