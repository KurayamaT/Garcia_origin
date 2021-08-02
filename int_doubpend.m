function [tout,yout]=int_doubpend(t0,tfinal,tol,y0,beta,g,gam,l);
% kludge to run mex file- assumes g=l=1, stores [beta gamma] in BETA
betavec=[beta,gam,g,l];
% now there are two different betas floating around

% Does forward integration for a "simplest" walker with
% 2 DOF. 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%       C.B. Moler, 3-25-87, 10-5-91, 6-3-93.
%       Copyright (c) 1984-93 by The MathWorks, Inc.
%       Modified for walking simulations
%       by Mariano Garcia (msg5@cornell.edu), 1996
%%%%%%%%%%%%%%%%%%%%%%%%%% Program begins  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialization
% The Fehlberg coefficients:
% https://ja.wikipedia.org/wiki/%E3%83%AB%E3%83%B3%E3%82%B2%EF%BC%9D%E3%82%AF%E3%83%83%E3%82%BF%E6%B3%95%E3%81%AE%E3%83%AA%E3%82%B9%E3%83%88
% http://godfoot.world.coocan.jp/Runge-Kutta-Fehlberg.htm
% https://maths.cnam.fr/IMG/pdf/RungeKuttaFehlbergProof.pdf
% http://ruina.tam.cornell.edu/research/topics/locomotion_and_robotics/
beta  = [ [    1      0      0     0      0    0]/4
          [    3      9      0     0      0    0]/32
          [ 1932  -7200   7296     0      0    0]/2197
          [ 8341 -32832  29440  -845      0    0]/4104
          [-6080  41040 -28352  9295  -5643    0]/20520 ]';
gamma = [ [902880  0  3953664  3855735  -1371249  277020]/7618050
          [ -2090  0    22528    21970    -15048  -27360]/752400 ]';
pow = 1/5; %累乗(power)
f = zeros(length(y0),6);

if nargin < 5, tol = 0.001; end % nargin = Number of function input arguments
t = t0;

hmax = (tfinal - t)/16;
h = hmax/8;
y = y0(:);
chunk = 128;
tout = zeros(chunk,1);
yout = zeros(chunk,length(y));
k=1;
tout(k) = t;
yout(k,:) = y.';
tau = tol * max(norm(y, 'inf'), 1);

clf reset;

axis([-2 2 -2 2]);

% Get coordinates for graphics
% cfx,cfy are coordinates of center of stance foot in a world fixed frame
% hip, knee, and end are coordinates of hip, swing knee, and end foot center

hipx=-l*sin(y(1));
hipy=l*cos(y(1));
endx=hipx-l*sin(y(1)+y(2));
endy=hipy-l*cos(pi-y(1)-y(2));

%plot graphics initially
walker_data=[0,0;hipx,hipy;endx,endy];
walker_bits=line('xdata',walker_data(:,1),'ydata',walker_data(:,2),...
'erasemode','xor');
drawnow;

%%%%%%%%%%%%%%%%%%%%%  Main Loop Begins  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
while (t < tfinal)

while  ((y(2)+h*y(4) <= pi-2*(y(1)+h*y(3))))...
 & t < tfinal & abs(y(1)) < pi/2 
         % Compute the slopes
	 f(:,1) = yderivs_doubpend(y,betavec);
   	    for j = 1:5
      	      f(:,j+1) = yderivs_doubpend(y+h*f*beta(:,j),betavec);
   	    end


        % Estimate the error and the acceptable error
   delta = norm(h*f*gamma(:,2),'inf');
   tau = tol*max(norm(y,'inf'),1.0);


        % Update the solution only if the error is acceptable
        ts = t;
        ys = y;
        if delta <= tau
      	   t = t + h;
	y = y + h*f*gamma(:,1);
       	   k = k+1;
              if k > length(tout)
 	          tout = [tout; zeros(chunk,1)];
		  yout = [yout; zeros(chunk,length(y))];
              end
	   tout(k) = t;
	yout(k,:) = y.';
	   %Update Graphics



	hipx=-l*sin(y(1));
	hipy=l*cos(y(1));

	endx=hipx-l*sin(y(1)+y(2));
	endy=hipy-l*cos(pi-y(1)-y(2));

	walker_data=[0,0;hipx,hipy;endx,endy];
        set(walker_bits,'xdata',walker_data(:,1),'ydata',walker_data(:,2));

	drawnow;

	end;

        % Update the step size

      if delta ~= 0.0
      h = min(hmax, 0.9*h*(tau/delta)^pow);
     end
end;

if abs(y(1)) > pi/2 
t=tfinal
end;

while  ((y(2)+h*y(4) >= pi-2*(y(1)+h*y(3))) ) ...
 & t < tfinal & abs(y(1)) < pi/2
         % Compute the slopes
         f(:,1) = yderivs_doubpend(y,betavec);
            for j = 1:5
              f(:,j+1) = yderivs_doubpend(y+h*f*beta(:,j),betavec);
            end

        % Estimate the error and the acceptable error
   delta = norm(h*f*gamma(:,2),'inf');
   tau = tol*max(norm(y,'inf'),1.0);


        % Update the solution only if the error is acceptable
        ts = t;
        ys = y;
        if delta <= tau
           t = t + h;
	y = y + h*f*gamma(:,1);
           k = k+1;
              if k > length(tout)
                  tout = [tout; zeros(chunk,1)];
                  yout = [yout; zeros(chunk,length(y))];
              end
           tout(k) = t;
        yout(k,:) = y.';
           %Update Graphics



        hipx=-l*sin(y(1));
        hipy=l*cos(y(1));

        endx=hipx-l*sin(y(1)+y(2));
        endy=hipy-l*cos(pi-y(1)-y(2));

        walker_data=[0,0;hipx,hipy;endx,endy];
        set(walker_bits,'xdata',walker_data(:,1),'ydata',walker_data(:,2));
        drawnow;

        end;

        % Update the step size

      if delta ~= 0.0
      h = min(hmax, 0.9*h*(tau/delta)^pow);
     end
end;

if abs(y(1)) > pi/2
t=tfinal
end;

    horig=h;
%%%%%%%%%%%%%% End THree Link Loop %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% now we are less than h seconds away from strike
% need to zero in on collision hypersurface
	q=0;

if t<tfinal
	while abs(pi-2*y(1)-y(2)) > tol/100
	  h=(pi-2*y(1)-y(2))/(y(4)+2*y(3));
            f(:,1) = yderivs_doubpend(y,betavec);
            for j = 1:5
              f(:,j+1) = yderivs_doubpend(y+h*f*beta(:,j),betavec);
            end

%don't worry about error computation in here, since h small enough
% assume solution is within tolerable error
	 
	  t = t + h;
	y = y + h*f*gamma(:,1);
          k = k+1;

		if k > length(tout)
		  tout = [tout; zeros(chunk,1)];
		  yout = [yout; zeros(chunk,length(y))];
		end
	  tout(k) = t;
	  yout(k,:) = y.';
	  q=q+1;

end;
end;

if t<tfinal
%Update Graphics

        hipx=-l*sin(y(1));
	hipy=l*cos(y(1));
	endx=hipx-l*sin(y(1)+y(2));
	endy=hipy-l*cos(pi-y(1)-y(2));
	walker_data=[0,0;hipx,hipy;endx,endy];
        set(walker_bits,'xdata',walker_data(:,1),'ydata',walker_data(:,2));

        drawnow;
  
h=horig;

% now we hit strike. AMB about strike point gives


y(3)=(y(3)*cos(2*y(1)))/(1+betavec(1)*(sin(2*y(1)))^2);
y(4)=-y(3)*(1+cos(y(2)));

y(1)=-y(1);
y(2)=pi-2*y(1);

y
t

k=k+1;
tout(k) = t;
yout(k,:) = y.';

end;
end;

tout = tout(1:k);
yout = yout(1:k,:);

