function raytracing(sstep,depmat,cmat,thetas,zsources,f0,numsteps,debug)

%%%%%%Variables%%%%%%%
% sstep- step size in s
% depmat- depth vector for the water column, [rval1 depval1;rval2 depval2...]
% cmat- c matrix
% thetas- launch angle vector for beams in radians, between -pi/2 and pi/2 [theta1;theta2;theta3;...]
% zsources- source depths
% f0 - source frequency
% numsteps - number of steps to calculate
% debug flag- if debug==true, then process will plot at each bottom/top
% interaction and pause so you can look at it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%define A0
A0=1;

%initialize figures (raytrace figure, intensity figure)
raytraceFig=figure;
hold on
plot(depmat(:,1),depmat(:,2))
set(gca,'YDir','Reverse')
intensityFig=figure;
hold on
%soundspeedFig=figure;
%hold on
cc=hsv(length(thetas)*length(zsources));
for l=1:length(zsources)
z0=zsources(l);

for k=1:length(thetas)

[c0,cz0]=getCVal(0,z0,cmat);
% For each theta value, calculate ray
theta0=thetas(k);

%initialize raytrace vectors
sray=[0];%s-values for ray
xiray=[cos(theta0)/c0];% xi values for ray
zetaray=[sin(theta0)/c0];% zeta values for ray
thetavec=[theta0]; % theta values for ray
rray=[0];%range vector for ray
zray=[z0];%depth vector for ray
cray=[c0];%soundspeed vector for ray
tau=[0];%tau vector for ray
q=[0];
p=[1/c0];
A=[A0];
cz=[cz0];
J0=0;
i=1;
rsign=1;

% Iterate over numsteps to calculate out the ray
while i<numsteps

r=rray(i);%r from previous step
z=zray(i);%z from previous step

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INSERT CODE HERE TO GET r, z OF RAY IN NEW STEP %%%%%%

dr=sstep*abs(cos(thetavec(i)));
dz=sstep*abs(sin(thetavec(i)));
r=rray(i)+dr;
thetavec(i);
if thetavec(i) > pi
thetavec(i) = thetavec(i) - 2*pi;
end
if thetavec(i)>0 & thetavec(i)<pi
z=zray(i)+dz;
else
z=zray(i)-dz;%moving up
end

if abs(thetavec(i))>pi/2
r=rray(i)-dr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%figure out if there is bottom or surface interaction,
%calculate new values
[interact,sray,zray,rray,zetaray,xiray,cray,tau,thetavec,c0,theta0,cz]=getBottomInteraction(r,z,...
                                                                                            depmat,sray,zray,rray,zetaray,xiray,cray,theta0,cmat,tau,thetavec,c0,cz);

%If there is an interaction, add two to 'i' and get new p, q, A
if interact

%add 2 to i to account for the two points from interaction fctn
[q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz); %First addition
[q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i+1,p,q,A,c0,cz); %Second
i=i+2; %get off of the surface
if debug %Plot the interaction if debug
plot(rray,zray)
hold on
plot(depmat(1:1000,1),depmat(1:1000,2))
hold off
%figure
%plot(sray(:,1:100),A(:,1:100))
pause
end
continue
end

% add the new rray, zray values
rray(i+1)=r; % the values you just calculated
zray(i+1)=z;

% Get new soundspeed for new location
[cray(i+1),cz(i+1)]=getCVal(rray(i+1),zray(i+1),cmat);


thetavec(i+1)=thetavec(i);%Placeholder
xiray(i+1)=xiray(i);%Placeholder
zetaray(i+1)=zetaray(i);%Placeholder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PUT YOUR CODE TO GET NEW thetavec(i+1) HERE  %%%%%%%%%%

if thetavec(i)>0 & thetavec(i) <= pi
thetavec(i+1)=acos(cray(i+1)/cray(i)*cos(thetavec(i)));
elseif thetavec(i) < 0 & thetavec(i) > -pi
thetavec(i+1)=-acos(cray(i+1)/cray(i)*cos(thetavec(i)));
else
thetavec(i)=0;
if i==1
thetavec(i+1)=.0000001;
else
thetavec(i+1)=-thetavec(i-1);
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% PUT YOUR CODE TO GET zetaray(i+1) HERE %%%%%%%%%%%%%%
zetaray(i+1)=zetaray(i)-sstep*cz(i+1)/cray(i+1)^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


sray(i+1)=sray(i)+sstep;% Not a placeholder
dtau=sstep/cray(i); % Not a placeholder
tau(i+1)=tau(i)+dtau; % Not a placeholder

%next, calculate the amplitudes,p and q:
[q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz);

i=i+1;% Iterate i
if rem(i/100,1)==0
disp('on step number ')
i
end
end

figure(raytraceFig)
plot(rray,zray,'color',cc(length(thetas)*(l-1)+k,:))

%plot(depmat(:,1),depmat(:,2))

figure(intensityFig)
plot(sray,20*log10(abs(A)),'color',cc(length(thetas)*(l-1)+k,:))

%figure(soundspeedFig)
%plot(cmat(:,2),cmat(:,1))
end
end
%figure
%plot(rray,20*log10(abs(A)))
%plot(sray,10*log10(A))

function [q,p,A]=getAmplitudes(cray,zray,rray,sstep,xiray,theta0,A0,i,p,q,A,c0,cz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% getAmplitudes function: get q, p, A values for ray input %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize:
q(i+1)=0;
p(i+1)=0;
A(i+1)=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% INSERT YOUR CODE HERE TO GET p, q, A values   %%%%%%%%%%%%%%%%%

%calculate the q and p values from equations pg 205 in book

dq=sstep*cray(i+1)*p(i);

q(i+1)=dq+q(i);

cz1=cz(i+1);
cz0=cz(i);
if zray(i+1)-zray(i) == 0
czz=0;
elseif i>1
czz=(cz(i+1)-cz(i))/(zray(i+1)-zray(i));%czz=(cz(i+1)-cz(i))/(zray(i+1)-zray(i));
%czz=1500*.00737*exp(2)*exp(-zray(i)/650)/(650^2);
%czz=(cray(i+1)-2*cray(i)+cray(i-1))/(zray(i+1)-zray(i))^2

else
czz=0;
end

%the soundspeed is only a f(z), so cnn is f(czz) only
cnn=(cray(i+1)*xiray(i+1))^2*czz;
%get new p value
dp=-sstep*cnn*q(i+1)/cray(i+1)^2;
p(i+1)=p(i)+dp;
J=rray(i+1)*q(i+1);
%calculate the amplitude-
A(i+1)=1/(4*pi)*A0*sqrt(abs(cos(theta0)*cray(i+1)/(J*c0)));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [interact,sray,zray,rray,zetaray,xiray,cray,tau,thetavec,c0,theta0,cz]=getBottomInteraction(r,z,depmat,sray,zray,rray,zetaray,xiray,cray,theta0,cmat,tau,thetavec,c0,cz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% getBottomInteraction function: return interaction %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m=length(sray);%current index

%find the depth and bottom angle for your range:
[D,bot_ang]=getDepth(r,depmat);

%initialize vars
interact=false;%no interaction initialized
if (z > D)

%there is a bottom interaction because you want to go lower than
%the bottom
disp('bottom interaction!')
interact=true;
%calculate the pathlength needed to reach the bottom exactly
dsnew=(D-zray(m))/(sin(thetavec(m)));
dtau=dsnew/cray(m);
tau(m+1)=tau(m)+dtau;

% Placeholders
rray(m+1)=rray(m);%placeholder
zray(m+1)=zray(m);%placeholder
xiray(m+1)=xiray(m);%placeholder
thetavec(m+1)=thetavec(m);%placeholder
sray(m+1)=sray(m);%placeholder
zetaray(m+1)=zetaray(m);%placeholder
cray(m+1)=cray(m);%placeholder
cz(m+1)=cray(m);%placeholder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+1  %%%%%%%%%%%%%%%%
%%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%


%get the time for the ray to go that distance

%calculate the range and depth for the bottom interaction
if abs(thetavec(m)) > pi/2
rray(m+1)=rray(m)-abs(dsnew*cray(m)/cray(m-1)*cos(thetavec(m-1)));
else
rray(m+1)=rray(m)+abs(dsnew*cray(m)/cray(m-1)*cos(thetavec(m-1)));
end
zray(m+1)=D;


%get values for the point on the bottom (m+1)
[cray(m+1),cz(m+1)]=getCVal(rray(m+1),zray(m+1),cmat);%get c at the bottom
%calculate xi, thetavec(m+1) at new point

%calculate thetavec(m+1)
if thetavec(m)-bot_ang < pi/2
thetavec(m+1)= -thetavec(m)+2*bot_ang;
rray(m+2)=rray(m+1)+dsnew*abs(cos(thetavec(m+1)));
else
disp('greater than 90')
thetavec(m+1) = 2*pi - (thetavec(m)-2*bot_ang);
thetavec(m+1)*180/pi

xiray(m+1)=-xiray(m);
rray(m+2)=rray(m+1)-abs(dsnew*cos(thetavec(m+1)));
end



%calculate zeta, s at new point
zetaray(m+1) = sin(thetavec(m+1))/cray(m+1);
sray(m+1)=sray(m)+dsnew;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Calculate for next step: going up one step.

% reset c0 and theta0- as though rays are now starting with new angles
% necessary for when there is a sloped bottom
c0=cray(m+1);
theta0=thetavec(m+1);
dtau=dsnew/cray(m);
tau(m+2)=tau(m+1)+dtau;

% Placeholder for m+2:
rray(m+2)=rray(m);%placeholder
zray(m+2)=zray(m);%placeholder
xiray(m+2)=xiray(m);%placeholder
thetavec(m+2)=thetavec(m);%placeholder
sray(m+2)=sray(m);%placeholder
zetaray(m+2)=zetaray(m);%placeholder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+2  %%%%%%%%%%%%%%%%
%%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%



%calculate r, z, c, s, tau, zeta, xi, theta at the following pointthetavec(i+1)=

zray(m+2)=zray(m+1)+dsnew*sin(thetavec(m+1));
%get c value
[cray(m+2),cz(m+2)]=getCVal(rray(m+2),zray(m+2),cmat);
%calculate m+2 values
sray(m+2)=sray(m+1)+dsnew;
dtau=dsnew/cray(m);
tau(m+2)=tau(m+1)+dtau;

if abs(thetavec(m+1)) > pi/2
thetavec(m+2)=2*pi-acos(cray(m+2)/c0*cos(theta0));
else
thetavec(m+2)=-acos(cray(m+2)/c0*cos(theta0));
end
xiray(m+2)=xiray(m+1);
zetaray(m+2)=(zetaray(m+1)-dsnew*cz(m+1)/cray(m+1)^2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



elseif (z < 0)
%there is a top interaction
disp('top interaction!')
interact=true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% INSERT YOUR CODE TO GET values at m+1  %%%%%%%%%%%%%%%%
%%%% of rray, zray, zetaray, cray, cz, xiray, sray, thetavec %%%%%%%%%
dsnew=abs(zray(m)/(sin(thetavec(m))));
%calc rray, zray
rray(m+1)=rray(m)+dsnew*cray(m)/c0*cos(theta0);

zray(m+1)=0;
%get c value and first derivative
[cray(m+1),cz(m+1)]=getCVal(rray(m+1),zray(m+1),cmat);
%get sray
sray(m+1)=sray(m)+dsnew;
sray(m+2)=sray(m)+2*dsnew;
%get theta
if abs(thetavec(m)) < pi
thetavec(m+1)=-thetavec(m);
else
thetavec(m+1) = 2*pi-thetavec(m);
end
%reset theta0, c0
c0=cray(m+1);

theta0=thetavec(m+1);
%calc zeta, xi
zetaray(m+1)=-zetaray(m);
xiray(m+1)=xiray(m);
%calc theta, rray, zray
theta_inc=acos(cray(m)/c0*cos(theta0));
rray(m+2)=rray(m+1)+dsnew*cray(m+1)/c0*cos(theta0);

zray(m+2)=zray(m+1)+dsnew*cray(m)/c0*sin(theta0);
%get c
[cray(m+2),cz(m+2)]=getCVal(rray(m+2),zray(m+2),cmat);
%get zeta, xi
zetaray(m+2)=zetaray(m+1)-dsnew*cz(m+2)/cray(m+2)^2;
xiray(m+2)=xiray(m+1);
%get theta
thetavec(m+2)=thetavec(m+1);
dtau=dsnew/cray(m);
%get taus
tau(m+1)=tau(m)+dtau;
tau(m+2)=tau(m+1)+dtau;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function [D,ang]=getDepth(r,depmat)


D=-100;ang=0;
% r- range value in m
% depmat- matrix of depth v. r

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TO GET DEPTH D and bottom angle ang at cur r, depmat  %%%%%%%%%%

%return a depth value for a given range value (r)
D=-100;

for l=1:size(depmat,1)-1
if depmat(l,1)==r
D=depmat(l,2);
ang=atan((depmat(l+1,2)-depmat(l,2))/(depmat(l+1,1)-depmat(l,1)));

elseif (depmat(l,1)<r & depmat(l+1,1)>r)
D=(depmat(l,2)+depmat(l+1,2))/2;

ang=atan((depmat(l+1,2)-depmat(l,2))/(depmat(l+1,1)-depmat(l,1)));

end

end

if D==-100
if r>=max(depmat(:,1))
D=depmat(end,2);
deplen=size(depmat,1);
ang=0;
%         if lastr < max(depmat(:,1))
%             disp('you are beyond the available range- assuming flat from last defined...')
%
%         end
else

ang=0;
D=depmat(1,2);
%         if lastr > 0
%             disp('Ray going backwards past 0! Assuming first depth correct.')
%
%         end

%error('desired r exceeds matrix range')
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [c,cz]=getCVal(r,z,cmat)
%return a c value in m/s for a given r,z value:
%search grid
c=-100; % soundspeed
cz=-100; % slope

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% INSERT YOUR CODE HERE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% TO GET soundspeed and soundspeed gradient for r,z  %%%%%%%%%%%%

for l=1:size(cmat,1)
if cmat(l,1)==z
c=cmat(l,2);
if l>1
cz=(cmat(l,2)-cmat(l-1,2))/(cmat(l,1)-cmat(l-1,1));
else
cz=(cmat(l,2)-cmat(l+1,2))/(cmat(l,1)-cmat(l+1,1));
end
elseif (l>1 & cmat(l,1)>z & cmat(l-1,1)<z)
%linearly interpolate between the values to get the answer
z;
cmat(l-1,1);
distm1=(z-cmat(l-1,1))/(cmat(l,1)-cmat(l-1,1));
distm2=(cmat(l,1)-z)/(cmat(l,1)-cmat(l-1,1));
c=(cmat(l-1,2)*distm1+cmat(l,2)*distm2 );
cz = (cmat(l,2)-cmat(l-1,2))/(cmat(l,1)-cmat(l-1,1));
elseif (z > max(cmat(:,1)))
c=cmat(end,2);
cz = 0;
%disp('using highest depth value for c... you exceeded dimensions of matrix')
end
end

if (z < min(cmat(:,1)))
%if between the last noted values and the surface, use linear
%extrapolation
cz=(cmat(1,2)-cmat(2,2))/(cmat(1,1)-cmat(2,1));
c=cmat(1,2)+cz*(cmat(1,1)-cmat(2,1));

end


if c==-100
%no c value was found
%return error: requested depth exceeds input range

error('requested depth exceeds input range...no c value can be returned')
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


