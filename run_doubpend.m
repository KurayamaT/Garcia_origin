% All code by Mariano Garcia (msg5@cornell.edu), 1998
% Some documentation of sorts and accompanying figures can
% be found in chaps 2 and 3 of my thesis, which can be downloaded from
% http://tam.cornell.edu/students/garcia/msghomepage.html
% or
% http://www.ccmr.cornell.edu/~ruinalab/
% also read the README.txt file.
% This simulation was used to produce results for the paper
% "The Simplest Walking Model: Stability, Complexity, And Scaling"
% ASME Journal of Biomechanical Engineering Vol 120, No. 2 pp. 218-288, 1998.

% Note that theta2 is the exterior angle between the
% legs, not the same as phi in the ASME paper above.
% theta2 = pi - phi

% see the accompanying figure pointfoot.cartoon2.eps


global graphics
graphics=0;
format long
clf reset
clear

% leg parameters
l = 1;   % leg length
m1 = 1;  % hip mass
m2= 0;   % foot mass
beta=m2/m1;

gam = 0.009;  % slope, in rad
g = 1;        % gravity

%initial conditions gam=0.009 long step to 1e-10
% loc  = -0.19120868949940 +/- 0.55633363687106i 
theta1 = 0.20031090049483;
theta2= pi-2*theta1;
theta1dot = -0.19983247291764;
theta2dot=-theta1dot*(1-cos(2*theta1));

%initial conditions gam=0.009 short step to 1e-10
% loc = 4.000, 0.459
%theta1 = 0.19393736960624;
%theta2= pi-2*theta1;
%theta1dot = -0.20386692731962;
%theta2dot=-theta1dot*(1-cos(2*theta1));

y0 = [theta1 theta2 theta1dot theta2dot]

tol = 1e-10;
tfinal =5;
%tic
[t,y]=int_doubpend(0,tfinal,tol,y0,beta,g,gam,l);
%toc


% plot output angles
plot(t,-y(:,1),t,pi-y(:,1)-y(:,2));
%plot(t,y(:,3),t,y(:,4));
%plot(t,y(:,1),t,y(:,2));
%plot(t,y(:,2)-pi+2.*y(:,1));
%plot(t,cos(y(:,1))-cos(pi-y(:,1)-y(:,2)));



