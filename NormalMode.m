function [TL,rmat,zmat]=NormalMode(rs,zs,f,f_vel,rrec,zrec,cmatw,rhow,cb,rhob,inc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inputs
% rs     source range
% zs     source depth
% f      frequency
% f_vec  freqs to use in calculating u, v
% rrec   reciever range vector, [r1,r2...]
% zrec   reciever depth vector, [z1,z2...]
% cmatw  soundspeed profile for the water column [z1,c1;z2,c2...]
% rhow   density of the water column (constant)
% cb     soundspeed in the bottom (scalar)
% rhob   density of the bottom
% inc    increment for layer size (for calculating along)
% 
% outputs:
% TL     transmission loss gridded sized r x z
% rmat   redundant grid of range (for plotting)
% zmat   redundant grid of depth (for plotting)

% the top interface is always assumed to be air, i.e. pressure release
% this code assumes two different fluid layers, rhow and rhob, with no
% boundary interactions between them

%% Initialize parameters

% Get the max water depth
maxd=max(cmatw(:,1));

% Make depth matrix for water
z=0:inc:maxd;

% Make depth and height matrices for water
rho=z*0+rhow;

% Calculate cmatrix for layers, given z of layers
cw=z*0+1500;
for q=1:length(z)
    cw(q)=getCVal(z(q),cmatw);   
end

% Add in bottom
zbot = (z(end)+inc):inc:z(end)*3; % go to a depth 3*cur depth
rhobot=0*zbot+rhob;
cbot=0*zbot+cb;
z=[z,zbot]; nz = length(z);
h(1:length(z))=inc;
rho=[rho,rhobot];
cw=[cw,cbot];

%% --------------------- Mode Shapes, krm values -------------------------
% ---------------------- Figures A, B ------------------------------------
psi_norm=zeros(length(z),101); %Initialize psi_norm matrix 
krm = [];
%actual psi_norm will have size length(z) x length(krm), initialize like this
%for plotting template

disp('calculating mode shapes and krm values...')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% INSERT YOUR CODE HERE TO GET psi_norm, krm at freq f %%%%%%%%%%%%

% Hint: Make sure you normalize psi and sort the eigenvalues/vectors to get
% the modes in the right order. 

% set up finite difference matrix COA pg 364
%eig value method

omega = 2*pi*f;
E = 1./(h(2:end).*rho(2:end));
D = (-2+h.^2.*(omega^2./cw.^2))./(h.*rho);
A = diag(D) + diag(E,1) + diag(E,-1);
[eigV,eigD] = eig(-A);
kappa = diag(eigD);
krm = sqrt(-kappa.*rhow./inc);
psi_norm = -eigV;

% use positive modes
num_modes = min(sum(krm>0),length(krm));
krm = krm(1:num_modes);

nbytes = fprintf('krm: 0 of %d \n', length(krm));
%COA pg 367
for m = 1:length(krm)
    while nbytes > 0
        fprintf('\b');
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('krm: %d of %d \n', m, length(krm));
    
    psi_norm(:,m) = -eigV(:,m);
    
%     w = ones(length(z),1);
%     for k = 1:10 % significantly large!
%         psi_norm(:,m) = (A-kappa(m)*eye(length(z)))\w;
%         psi_norm(:,m) = psi_norm(:,m)/norm(psi_norm(:,m));
%     end
end

% Figure A
% plot krm v. mode number
figure
plot(1:length(krm),krm, 'linewidth', 2);
title('Plot of krm v. mode number')
xlabel('mode number')
ylabel('krm (1/m)')
axis tight

% Figure B
% plot the first, second, and 30th normalized modes
figure
count = 0;
for mn = [1 2 30]
    count = count + 1;
    subplot(1,3,count);
    plot_psi = psi_norm(:,mn)./max(psi_norm(:,mn)); %bounds at -1,1
    plot(plot_psi,z','linewidth',1.2)
    hold on
    plot([-2 2],[maxd maxd],'k:')
    plot([0 0],[3*maxd 0],'k--')
    hold off
    title(['Mode No. ' num2str(mn)])
    xlim([-1.1 1.1]); ylim([0 2*maxd]);
    view(0,-90)
end
subplot(1,3,1)
ylabel('Depth(m)');
   
%% ------------------- Transmission Loss -------------------------------
%  ------------------- Figures C, D ------------------------------------
disp('calculating Transmission Loss...')
TL = zeros(length(z),length(rrec));
zmat=TL;% Populate this as you populate TL
rmat=TL;% Populate this as you populate TL

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INSERT YOUR CODE HERE TO GET TRANSMISSION LOSS %%%%%%%%%%%%%%
zs_ind = getZInd(zs,z);
% assume point source, pg 340
for j = 1:length(z)
    inner = 0;
    for m = 1:num_modes
        inner = inner + psi_norm(zs_ind,m)*psi_norm(j,m)*exp(1i*krm(m)*abs(rrec-rs))./sqrt(krm(m));
    end
    TL(j,:) = 20*log10(abs(1/rho(zs_ind)*sqrt(2*pi./abs(rs-rrec)).*inner));
    zmat(j,:) = z(j)*ones(1,length(rrec));
    rmat(j,:) = rrec;
end

% Figure C
% transmission loss vs range and depth
figure
surface('XData',rmat/1000,'YData',zmat,'ZData',TL,'CDATA',TL,'LineStyle','None')
set(gca,'Ydir','reverse')
title('Transmission Loss v range and depth')
xlabel('Range (km)')
ylabel('Depth (m)')
colormap 'jet';
% Statistics for colorbar - ECB
medTL = median(TL(:));
caxis([medTL-5 medTL+25])

% Figure D
% transmission loss @ source depth
figure
plot(rrec/1000,TL(zs_ind,:),'linewidth',1.2)
title('TL vs range @ z = zs')
xlabel('Range (km)')
ylabel('Transmission Loss (dB)')

%% -------------------------- Phase and Group Velocities ----------------
%  -------------------------- Figure E ----------------------------------

disp('calculating phase and group velocities...')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INSERT YOUR CODE HERE TO GET PHASE AND GROUP VELOCITIES %%%%%%%%

% Run through process to find eigenvalues each time for different freqs
vn=zeros(3,length(f_vel)); %phase velocity
un=zeros(3,length(f_vel)); %group velocity
first3krm = vn';


nbytes = fprintf('freq: 0 of %d \n',length(f_vel));
for df = 1:length(f_vel)
    while nbytes > 0
        fprintf('\b');
        nbytes = nbytes - 1;
    end
    nbytes = fprintf('freq: %d of %d \n', df, length(f_vel));
    omega = 2*pi*f_vel(df);
    
    for kk = 1:nz
        if kk > 1
            A(kk,kk-1) = 1/(inc*rho(kk));
        end
        A(kk,kk) = (-2 + inc.^2*(omega.^2./cw(kk).^2))/(inc*rho(kk));
        if kk < nz
            A(kk,kk+1) = 1/(inc*rho(kk));
        end
    end
    kappa = eigs(A,nz);
    kappa = sort(kappa,'descend');
    krm = sqrt(kappa.*rho.'./inc);
    krm = krm(1:3);
    
    first3krm(df,:) = krm;
    vn(:,df) = omega./krm;
end

for mm = 1:3
    un(mm,:) = gradient(2*pi*f_vel, first3krm(:,mm));
end

% Figure E
% plot un and vn for the first 3 modes
figure
% clean vn, un
ndx = find(imag(vn(:))); vn(ndx) = NaN;
ndx = find(imag(un(:))); un(ndx) = NaN;
plot(f_vel,vn,'linewidth',2)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(f_vel,un,'--','linewidth',2)
hold off
xlabel('Frequency (Hz)');
ylabel('Velocity (m/s)');
title('Phase (v) and group (u) velocities vs frequency');

lgd = legend('v_1','v_2','v_3','u_1','u_2','u_3');
lgd.NumColumns = 2;

fprintf('Done! \n')
end



%% Helper Function : getZInd
function ind = getZInd(zdes,z)
diffvec=z - zdes;
[~,minind]=min(abs(diffvec));
ind = minind;
end

%% Helper Function: getCVal
function [c]=getCVal(z,cmat)
%return a c value in m/s for a given r,z value:
%search grid
diffvec = abs(cmat(:,1)-z);
[~,ind]=min(diffvec);
c=cmat(ind,2);
end