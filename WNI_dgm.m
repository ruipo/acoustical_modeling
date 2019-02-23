function [pr,r,Pmat]=WNI_dgm(f,cmatw,numlay_water,sourcelay,rmin,rmax,rhob,cb,samp_mult,max_coeff)
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

[cvec,z_bounds,z_cents] = getCLayers(cmatw,numlay_water); % Get soundspeed in each layer
cvec(end+1) = cb; % Set the last soundspeed to the bottom soundspeed
klays = 2*pi*f./cvec; % Get k values in each layer

krmax=omega/min(cmatw(:,2)); % Maximum value of kr in water column
delta=3/(rmax*log10(exp(1))); % Set delta value

dkr=pi/(2*(rmax-rmin)); % Set dkr baseline

krs=0:dkr*samp_mult:krmax*max_coeff; % Set range of krs for integration

krs=krs*(1-1i*delta); % Want to integrate across complex wavenumber with branch cut below real axis. 

Pmat=zeros(length(krs),length(z_cents)); % Initialize matrix of pressures

%row corresponds to z depth, column corresponds to kr
sourceind=sourcelay; % The source index is the source layer


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%             INSERT YOUR CODE HERE                %%%%%%%%%%%%%%
%%%%%%%%%    Across kr values, populate P values in Pmat   %%%%%%%%%%%%%%
%%%%%%%%%     Feel free to add functions to do this.       %%%%%%%%%%%%%%
for ndx_kr = 1:length(krs)
    % initialize DGM
    dgm = sparse(zeros(2*numlay_water+1));
    % calculate kz at each layer for the current kr value
    kz = zeros(numlay_water,1);
    for ndx_kz = 1:numlay_water
        kz(ndx_kz) = get_kz(krs(ndx_kr),klays(ndx_kz));
    end
  
    % Step through layer by layer to populate the DGM from the local
    % coefficients matrices

    % Assume a Pressure Release boundary for the surface: Boundary condition: simga_z = 0;
    p_fl1 = fluid_layer(omega,kz(1),rhow);
    dgm(1,1:2) = p_fl1(2,1:2);
    %dgm(1,2) = dgm(1,2)*exp(1i*kz(1)*z_bounds(2));
    for pointer_ndx = 2:numlay_water
        dgm((2*pointer_ndx-2):(2*pointer_ndx-1),(2*pointer_ndx-3):...
            (2*pointer_ndx-2)) = -fluid_layer(omega,kz(pointer_ndx-1),rhow)...
            *ee(kz(pointer_ndx-1),z_bounds(pointer_ndx));
        dgm((2*pointer_ndx-2):(2*pointer_ndx-1),(2*pointer_ndx-1:2*pointer_ndx))...
            = fluid_layer(omega,kz(pointer_ndx),rhow)...
            *ee(kz(pointer_ndx),z_bounds(pointer_ndx));
    end

    % Assume a fluid halfspace for the bottom condition
    dgm((end-1):end,(end-2):(end-1)) = -fluid_layer(omega,kz(numlay_water),rhow)*...
        ee(kz(numlay_water),z_bounds(end));
    kzb = get_kz(krs(ndx_kr),2*pi*f/cb);
    p_bhs = fluid_layer(omega,kzb,rhob)*ee(kzb,z_bounds(end));

    dgm((end-1):end,end) = p_bhs(1:2,1);    


    Sw = -4*pi/(50*rhow*omega^2);
    rhs = sparse(zeros(length(dgm),1));
    rhs_stamp = [Sw*exp(1i*kz(sourceind)*(z_cents(sourceind)-z_bounds(sourceind-1)))/(4*pi);...
        Sw*exp(1i*kz(sourceind)*(z_cents(sourceind)-z_bounds(sourceind-1)))*rhow*omega^2/(1i*4*pi*kz(sourceind));...
        Sw*exp(1i*kz(sourceind)*(z_bounds(sourceind+1)-z_cents(sourceind)))/(4*pi);...
        -Sw*exp(1i*kz(sourceind)*(z_bounds(sourceind+1)-z_cents(sourceind)))*rhow*omega^2/(1i*4*pi*kz(sourceind))];

    rhs(2*(sourceind-1):2*(sourceind-1)+3) = rhs_stamp;
    dgm(isnan(dgm)) = 0;
    dgm(~isfinite(dgm)) = 0;
    amp = dgm\rhs;


    for z_ndx = 1:length(z_cents)
    
        Pmat(ndx_kr,z_ndx) = rhow*omega^2*(amp(2*z_ndx-1)*exp(1i*kz(z_ndx)*z_cents(z_ndx))+...
        amp(2*z_ndx)*exp(-1i*kz(z_ndx)*z_cents(z_ndx)));
    end
    Pmat(ndx_kr,sourceind) = Pmat(ndx_kr,sourceind)+rhow*omega^2*Sw/(4*pi*1i*kz(sourceind));

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Plot magnitude v. horizontal wavenumber at 46m
% NOTE: IF YOU CHANGE NUMBER OF LAYERS SO IT IS NOT 100, CHANGE plot_idx
% 
% plot_idx=46; % Index for plotting
% 
% figure(1)
% hold on
% plot(krs,abs(Pmat(:,plot_idx)))
% xlabel('Horizontal Wavenumber (m^-1)')
% ylabel('Magnitude')
% % Plot kernel
% figure(2)
% hold on
% h=pcolor(abs(krs),flipud(z_cents),20*log10(abs(Pmat')));
% set(h,'LineStyle','None')
% colormap('jet')
% colorbar
% xlabel('Horizontal Wavenumber (m^-1)')
% ylabel('Depth (m)')

% run ffp to do actual integration
[pr,r,r_tot]=ffp(Pmat',krs,z_cents);

% Calculate transmission loss
TL = 20*log10(abs(pr));

TL_tot = zeros(size(TL));

TL_tot(:,size(TL,2)/2+1:end) = TL(:,1:size(TL,2)/2);
TL_tot(:,1:size(TL,2)/2) = fliplr(TL(:,1:size(TL,2)/2));

TL_tot = 0.5*(TL_tot + TL);

% Plot transmission loss at 46m
% figure(3)
% hold on
% plot(r,20*log10(abs(pr(plot_idx,:))))
% axis([0,5000,-60,-10])
% xlabel('Range (m)')
% ylabel('Loss (dB)')

% Plot transmission loss color plot
figure(4)
hold on
surface('XData',r_tot,'YData',z_cents,'ZData',flipud(TL),'CDATA',flipud(TL),'LineStyle','None')
title('Plot of Transmission Loss v range and depth')
xlabel('Range (m)')
ylabel('Depth (m)')
colormap('jet')
colorbar

    function [c]=getCVal(z,cmat)
        %return a c value in m/s for a given r,z value:
        %search grid
        diffvec = abs(cmat(:,1)-z);
        [zm,ind]=min(diffvec);
        c=cmat(ind,2);
    end

    function [cvec,z_bounds,z_cent]=getCLayers(cmatw,numlay_water);
        maxd=max(cmatw(:,1));
        %divide up the cprofile up into layers
        z_bounds =0:maxd/(numlay_water):maxd;
        z_cent = zeros(numlay_water,1); 
        for m=1:numlay_water
            z_cur = (z_bounds(m) + z_bounds(m+1))/2;
            z_cent(m) = z_cur;
            cvec(m) = getCVal(z_cur,cmatw);
        end
    end
% Define a function that creates the submatrix for a homogenous fluid layer (2x2)        
    function p_fl = fluid_layer(omega,kz,rho)
        
        % w = j*k_z*A_p-j*k_z*A_m
        % sigma_zz = -rho*omega^2*A_p-rho*omega^2*A_m
        p_fl = [1i*kz, -1i*kz;-rho*omega^2,-rho*omega^2];
    end

    function kz = get_kz(kr,k)
        if real(kr) <=k
            kz = sqrt(k^2-kr^2);
        else
            kz = 1i*sqrt(kr^2-k^2);
        end
    end

    function emat = ee(kz,z)
        emat = diag([exp(1i*kz*z),exp(-1i*kz*z)]);
    end

     
end