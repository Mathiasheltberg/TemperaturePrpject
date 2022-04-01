clear all; close all; clc

x = linspace(0,1,1000)';
M = zeros(1000,1000);
for test = 1:1000;
    a = abs(randn);
    b = 10000*abs(rand);
    l = b*x./(a+x);
    M(test,:) = l;
    plot(l); hold on;
end

plot(mean(M),'k','LineWidth',3);
goodplot