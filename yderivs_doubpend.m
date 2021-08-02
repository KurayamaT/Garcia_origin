function ydot = yderivs_doubpend(y,betavec);

% 要するにRunge-Kutta-Fehlberg Methodを使って、導関数を導くためのプログラム

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

% This file contains the derivatives of the simplest walker.


beta = betavec(1); % Runge-Kutta-Fehlberg Method の係数行列が入っています。
gam = betavec(2);
g = betavec(3);
l = betavec(4);

% prepare for define Garcia's equations. 
cO1 = cos(y(1));
sO1 = sin(y(1));

cO2 = cos(y(2));
sO2 = sin(y(2));

sO1g = sin(y(1)-gam);
sO12g=sin(y(1)+y(2)-gam);

% define Garcia's equation(upper Matrix p3).
% 左辺第1項
M(1,1)=1+2*beta*(1+cO2);
M(1,2)=beta*(1+cO2);
M(2,1)=1+cO2;
M(2,2)=1;

% 左辺第2項
V(1)=-beta*sO2*y(4)*(2*y(3)+y(4));
V(2)=sO2*y(3)^2;

% 左辺第3項
G(1)=-g/l*(beta*sO12g+sO1g*(1+beta));
G(2)=-g/l*sO12g;

%-V'-G'

Oddot=M\(-V'-G');


ydot(1)=y(3);
ydot(2)=y(4);
ydot(3)=Oddot(1);
ydot(4)=Oddot(2);


ydot = ydot';


