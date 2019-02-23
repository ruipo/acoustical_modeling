function [pr,r]=WNI(f,cmatw,numlay_water,sourcelay,rmin,rmax,rhob,cb,samp_mult,max_coeff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Run wavenumber integration with the following inputs %%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f- frequency
% cmatw- water column soundspeed [z0,c0;z1,c1...]
% numlay_water- number of layers to break the water part into
% sourcelay- layer containing the source
% rmin- minimum r
% rmax- maximum r
% rhob- bottom density
% cb - bottom soundspeed
% samp_mult - multiplier on baseline dkr to use in WNI (larger -> larger
                                                        % step size)
% max_coeff - maximum multiplier of kr to use in WNI
%
% Assumes that top is vacuum and bottom is a fluid halfspace with cb, rhob.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%% Initialization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rhow=1000; %Set density of water
omega=2*pi*f; %frequency

[cvec,zints] = getCLayers(cmatw,numlay_water); % Get soundspeed in each layer
cvec(end+1) = cb; % Set the last soundspeed to the bottom soundspeed
klays = 2*pi*f./cvec; % Get k values in each layer
disp(zints)
disp(size(zints))
disp(size(cvec))

krmax=omega/min(cmatw(:,2)); % Maximum value of kr in water column
delta=3/(rmax*log10(exp(1))); % Set delta value

dkr=pi/(2*(rmax-rmin)); % Set dkr baseline

krs=0:dkr*samp_mult:krmax*max_coeff; % Set range of krs for integration

krs=krs*(1-i*delta); % Want to integrate across complex wavenumber with branch cut below real axis.

wvec=zeros(length(zints),length(krs)); % Initialize matrix of ws
Pmat=zeros(length(krs),length(zints)); % Initialize matrix of pressures

%row corresponds to z depth, column corresponds to kr
sourceind=sourcelay; % The source index is the source layer

kmat=[]; % Initialize kmat
row=1;% row=1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             INSERT YOUR CODE HERE                %%%%%%%%%%%%%%
%%%%%%%%% Loop across kr values, populate P values in Pmat %%%%%%%%%%%%%%
%%%%%%%%%     Feel free to add functions to do this.       %%%%%%%%%%%%%%

for ind=1:length(krs)
kzms = sqrt(klays.^2 - krs(ind).^2);
kr=krs(ind);
Pvec = PropogateKs(cvec,kzms,zints,omega,rhow,rhob,cb,sourceind,krs(ind),krmax);
if isnan(Pvec(1))
disp('got isnan')
else
Pmat(row,:)=Pvec;
row=row+1;
end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude v. horizontal wavenumber at 46m
% NOTE: IF YOU CHANGE NUMBER OF LAYERS SO IT IS NOT 100, CHANGE INDEX

plot_idx=46; % Index for plotting

figure(1)
hold on
plot(krs,abs(Pmat(:,plot_idx)))

% Plot kernel
figure(2)
hold on
h=pcolor(abs(krs),zints',20*log10(abs(Pmat')))
set(h,'LineStyle','None')

% run ffp to do actual integration
[pr,r]=ffp(Pmat',krs,zints);
           
           % Calculate transmission loss
           TL = 20*log10(abs(pr));
           
           % Plot transmission loss at 46m
           figure(3)
           hold on
           plot(r,20*log10(abs(pr(plot_idx,:))))
           axis([0,5000,-60,-10])
           
           % Plot transmission loss color plot
           figure(4)
           hold on
           surface('XData',r,'YData',zints,'ZData',TL,'CDATA',TL,'LineStyle','None')
           title('Plot of Transmission Loss v range and depth')
           xlabel('Range (m)')
           ylabel('Depth (m)')
           
           
           function [c]=getCVal(z,cmat)
           %return a c value in m/s for a given r,z value:
           %search grid
           diffvec = abs(cmat(:,1)-z);
           [zm,ind]=min(diffvec);
           c=cmat(ind,2);
           
           function [cvec,zints]=getCLayers(cmatw,numlay_water);
           maxd=max(cmatw(:,1))
           %divide up the cprofile up into layers
           zints =0:maxd/(numlay_water-1):maxd;
           cvec(1) = 340; % Air
           for m=2:numlay_water
           z_cur = (zints(m) + zints(m-1))/2;
           cvec(m) = getCVal(z_cur,cmatw);
           end
