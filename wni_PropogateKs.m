%calculate w, P matrices
function[Pvec]= wni_PropogateKs(cvec,kzms,zints,omega,rhow,rhob,cb,sourceind,kr,krmax)
%number of layers
warning('off','all');
numlays =  length(zints); %number of interfaces
Vr = zeros(2,numlays); %homogeneous initialization
Vs=zeros(2,numlays);
Vr(:,numlays)=[1;-rhob*omega^2/(i*kzms(end))];

Sw=4*pi/(rhow*omega^2);
Vs(:,sourceind) = [Sw/(2*pi);0];%source initialization


RN=1;
RS=1;
%propogate up from the bottom
for m=1:numlays-1


% kmz1 = kzms(numlays-m);
kmz1 = kzms(numlays-m+1);

z0=zints(numlays-m);
z1=zints(numlays-m+1);

if m==1
rho1 = rhob;
else
rho1 = rhow;
end
%This set has odd (Though not crazy) poles and better looking TL
Cm0 = [-i*kmz1*exp(-i*kmz1*z0),i*kmz1*exp(i*kmz1*z0);-rho1*omega^2*exp(-i*kmz1*z0),-rho1*omega^2*exp(i*kmz1*z0)];
Cm1 = [-i*kmz1*exp(-i*kmz1*z1),i*kmz1*exp(i*kmz1*z1);-rho1*omega^2*exp(-i*kmz1*z1),-rho1*omega^2*exp(i*kmz1*z1)];

P = Cm0*inv(Cm1);
RN=RN*P;
Vr(:,numlays-m)=RN*Vr(:,numlays);

%add in source
if (numlays-m+1)<=sourceind

RS=RS*P;
Vs(:,numlays-m) = RS*Vs(:,sourceind);

end

end
wn= -(Vs(2,1)/Vr(2,1));
Vr_norm = wn(end).*Vr;
Vr=Vr_norm;
for m=numlays-sourceind+1:numlays-1

%you are above the source: calculate up to the surface for each
Vs=zeros(2,numlays);
Vs(:,numlays-m+1) = [Sw/(2*pi);0];
RS=1;

for n=m:numlays-1
kmz1 = kzms(numlays-n+1);
%kmz1 = kzms(numlays-n+1);

z0=zints(numlays-n);
z1=zints(numlays-n+1);
rho1=rhow;
Cm0 = [-i*kmz1*exp(-i*kmz1*z0),i*kmz1*exp(i*kmz1*z0);-rho1*omega^2*exp(-i*kmz1*z0),-rho1*omega^2*exp(i*kmz1*z0)];
Cm1 = [-i*kmz1*exp(-i*kmz1*z1),i*kmz1*exp(i*kmz1*z1);-rho1*omega^2*exp(-i*kmz1*z1),-rho1*omega^2*exp(i*kmz1*z1)];

P = Cm0*inv(Cm1);
RS=RS*P;
%add in source

Vs(:,numlays-n) = RS*Vs(:,numlays-m+1);


end
wn= -(Vs(2,1)/Vr(2,1));

Vr_norm(:,numlays-m) = Vr(:,sourceind)*wn;


end


Pvec=Vr_norm(2,:)/50;
