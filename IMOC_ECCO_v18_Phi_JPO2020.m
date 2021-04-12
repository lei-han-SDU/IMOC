

%% LH created on 2019-6-1 15:04:03£º 
% LH re-created on 2019-8-5 21:15:58 With capability to handle topography


reload = 0;

if reload; clear; reload=1; end

set_fig_font;

p = genpath('gcmfaces/'); addpath(p);

clear global;
% predefine grid
grid_load; 


%% Load data


global mygrid myenv;

% Load time variable
nc=netcdf.open('release3/nctiles_monthly/oceTAUX/oceTAUX.0001.nc',0);
vv = netcdf.inqVarID(nc,'tim');
time_ecco_raw=netcdf.getVar(nc,vv);
time_ecco = datetime(1992,1,0)+time_ecco_raw;

ntime = length(time_ecco);


%% Load ECCO data (most time-consuming part, caution re-updating)

if reload
    tic
    myenv.nctilesdir=fullfile('release3',filesep,'nctiles_monthly',filesep);
    listVars={'UVELMASS','VVELMASS','WVELMASS','THETA'};

    for vvv=1:length(listVars)
        vv=listVars{vvv};
        tmp1=read_nctiles([myenv.nctilesdir vv '/' vv],vv);
        %tmp1=mean(tmp1,4); %tmp1(mygrid.mskC==0)=NaN;
        eval([vv '=tmp1;']);
    end
    
    clear tmp1;

    disp('Finish loading UVEL and VVEL')
    toc
   
    %clear VVELMASS UVELMASS WVELMASS THETA

end



%%
% Extract lat for [V,W] in YZ plan
lat_v = convert2array(mygrid.YG);
lat_v_1d = lat_v(100,:);
lat_w = convert2array(mygrid.YC);
lat_w_1d = lat_w(100,:);

lon_v = convert2array(mygrid.XC);

jm = length(lat_v_1d);
km = length(mygrid.RC);
im = jm;

% Grid sizes
DYG2d = convert2array(mygrid.DYG);
DXG2d = convert2array(mygrid.DXG);
DYF2d = convert2array(mygrid.DYF);
DXF2d = convert2array(mygrid.DXF);
DRF   = mygrid.DRF;
DYF1d = DYF2d(100,:);
DYG1d = DYG2d(100,:);
% For Poisson LH 2019-6-7 23:24:28
DYC2d = convert2array(mygrid.DYC);
DYC1d = DYC2d(100,:);
DRC   = mygrid.DRC;


%% Stretched vertical coordinate 
dep_uv = - mygrid.RC;
dep_w  = - mygrid.RF;   dep_w(end) = dep_uv(end);

dep_uv_strt = sqrt(dep_uv);
dep_w_strt  = sqrt(dep_w );

% define labels in stretched coordinate
dep_strt_cell = {'5900','5000','4000','3000','2000','1000','500','100'} ;  % Note: must in descending order!
dep_strt_num  = str2num(char(dep_strt_cell)) ;
dep_strt_tick_uv = interp1(dep_uv, dep_uv_strt, dep_strt_num);
dep_strt_tick_w  = interp1(dep_uv, dep_uv_strt, dep_strt_num);


%% Revise mask at ITF inlet. LH 2019-6-6 14:33:48

basin_name = {'ind'};
% list of available basins:
% list0={	'pac','atl','ind','arct','bering',...
%         'southChina','mexico','okhotsk','hudson','med',...
%         'java','north','japan','timor','eastChina','red','gulf',...
%         'baffin','gin','barents'};
% 
% define basin mask
basin_msk = v4_basin(basin_name);

basin_msk_2d = convert2array(basin_msk);
basin_msk_2d(153:165,139:152) = 0;  % revise ITF inlet


%% East and west boundary LH 2019-6-6 14:36:25
% Manually identify the Open boundaries in the west and east 

i_west = ones(jm,1);
i_east = ones(jm,1);

i_west(195:196) = 98;
i_west(183:185) = 90;
i_west( 62:123) = 59;

i_east(176:178) = 137;  % Malacca
i_east(155)     = 145;  % Sunda
i_east(139:151) = 153;  % ITF+Lombok
% f4
i_east( 67:111) = 186;  % Must be 186
i_east(116:117) = 186;

%% Extract snapshot and Convert to double array LH 2019-8-4 15:29:15
nt = 61+6;
time_nt_str = datestr(time_ecco(nt));

vcur3d = convert2array(VVELMASS(:,:,:,nt));
ucur3d = convert2array(UVELMASS(:,:,:,nt));
wcur3d = convert2array(WVELMASS(:,:,:,nt));
theta3d = convert2array(THETA(:,:,:,nt));

%% Extract u_east and u_west  LH 2019-6-6 14:38:28

u_east = zeros(jm,km);
u_west = zeros(jm,km);

for j=1:jm
  u_west(j,:) = ucur3d(i_west(j),j,:);
  u_east(j,:) = ucur3d(i_east(j),j,:);
  if j<=125  % area f4: zonal vel = vcur
    u_east(j,:) = vcur3d(i_east(j),j,:);
  end    
  if i_west(j)==1; u_west(j,:) = 0; end
  if i_east(j)==1; u_east(j,:) = 0; end
end

% Make sure boundary value is 0 instead of nan otherwise du will be
% contaminated
u_east(isnan(u_east))=0;
u_west(isnan(u_west))=0;

du_WE = u_west - u_east;

figure;
subplot(3,1,1)
%h = pcolor(lat_w_1d,-dep_uv_strt,du');colorbar;%cmocean('curl','pivot',0); 
h = pcolor(1:360,-dep_uv_strt,u_west');colorbar;%cmocean('curl','pivot',0); 
caxis([-.2 .2]);colormap('redblue')
set(h, 'edgecolor','none'); title('u:west')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
set(gca,'color',111*[1 1 1]/255)

subplot(3,1,2)
%h = pcolor(lat_w_1d,-dep_uv_strt,du');colorbar;%cmocean('curl','pivot',0); 
h = pcolor(1:360,-dep_uv_strt,u_east');colorbar;%cmocean('curl','pivot',0); 
caxis([-.2 .2]);colormap('redblue')
set(h, 'edgecolor','none'); title('u:east')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
set(gca,'color',111*[1 1 1]/255)

subplot(3,1,3)
%h = pcolor(lat_w_1d,-dep_uv_strt,du');colorbar;%cmocean('curl','pivot',0); 
h = pcolor(1:360,-dep_uv_strt,du_WE');colorbar;%cmocean('curl','pivot',0); 
caxis([-.2 .2]);colormap('redblue')
set(h, 'edgecolor','none'); title('u:west-u:east')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
set(gca,'color',111*[1 1 1]/255)



%% Calculate the divergence vdx and wdx  LH 2019-6-6 14:43:20

S_vw = zeros(jm,km);
S_v = S_vw;
S_w = S_vw;

vdx3d_msk = zeros(im,jm,km);
wdx3d_msk = zeros(im,jm,km);
theta_msk = zeros(im,jm,km);
for i=1:im
    theta_msk(i,:,:) = theta3d(i,:,:).*basin_msk_2d(i,:);
    if i<=180  % for data in .f1 and .f2 tiles
      vdx3d_msk(i,:,:) = vcur3d(i,:,:).*DXG2d(i,:).*basin_msk_2d(i,:);
      wdx3d_msk(i,:,:) = wcur3d(i,:,:).*DXF2d(i,:).*basin_msk_2d(i,:);
    else  % for data in .f4 tile
      vdx3d_msk(i,2:end,:) = -ucur3d(i,1:end-1,:).*DYG2d(i,1:end-1).*basin_msk_2d(i,1:end-1);
      %vdx3d_msk(i,:,:) = -ucur3d(i,:,:).*DYG2d(i,:).*basin_msk_2d(i,:);
      wdx3d_msk(i,:,:) =  wcur3d(i,:,:).*DYF2d(i,:).*basin_msk_2d(i,:);
    end
end

% Make Vertical plane mask 
yz_msk_c = squeeze(nansum(abs(theta_msk),1));
yz_msk_c(yz_msk_c==0) = nan;
yz_msk_c(~isnan(yz_msk_c)) = 1;


%%

% zonally integrated V and W
Vdx = squeeze(nansum(vdx3d_msk,1));
Wdx = squeeze(nansum(wdx3d_msk,1));

% divergence of velocity field
for j=1:jm-1
    for k=1:km-1
        S_v(j,k) = (Vdx(j+1,k)-Vdx(j,k)) / DYF1d(j);
        S_w(j,k) = (Wdx(j,k)-Wdx(j,k+1)) / DRF(k);
    end
end
S_vw = S_v + S_w;


% Correlation of S_v and S_w, & rms(S_vw)
x1 = S_v  ;
x2 = S_w  ;
x3 = S_vw ;
x1(isnan(x3))=[];
x2(isnan(x3))=[];
x3(isnan(x3))=[];

cc=corrcoef(x1,x2);

disp(['RMS of Svw = ',num2str(rms(x3))])
disp(['Correlation of Sv and Sw = ',num2str(cc(1,2))])

%% Check du_WE == S_vw ? LH 2019-6-6 14:50:34


figure;
subplot(3,1,1)
h = pcolor(lat_w_1d,-dep_uv_strt,du_WE');colorbar;%cmocean('curl','pivot',0); 
%h = pcolor(1:360,-dep_uv_strt,du_WE');colorbar;%cmocean('curl','pivot',0); 
caxis([-.1 .1]);colormap('redblue')
set(h, 'edgecolor','none'); title('u(west)-u(east)')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
xlim([-70 20])
set(gca,'color',111*[1 1 1]/255)

subplot(3,1,2)
h = pcolor(lat_w_1d,-dep_uv_strt,S_vw');colorbar;%cmocean('curl','pivot',0); 
%h = pcolor(1:360,-dep_uv_strt,S_vw');colorbar;%cmocean('curl','pivot',0); 
caxis([-.1 .1])
set(h, 'edgecolor','none'); title('\partial(Vdx)/\partialy+\partial(Wdx)/\partialz')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
xlim([-70 20])
set(gca,'color',111*[1 1 1]/255)

subplot(3,1,3)
h = pcolor(lat_w_1d,-dep_uv_strt,S_vw'-du_WE');colorbar;%cmocean('curl','pivot',0); 
%h = pcolor(1:360,-dep_uv_strt,S_vw');colorbar;%cmocean('curl','pivot',0); 
caxis([-.001 .001])
set(h, 'edgecolor','none'); title('\partial(Vdx)/\partialy+\partial(Wdx)/\partialz - [u(west)-u(east)]')
set(gca,'ytick',-dep_strt_tick_w); set(gca,'yticklabel',dep_strt_cell)
%set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
xlim([-70 20])
set(gca,'color',111*[1 1 1]/255)



% 1) Total transports
du_sum = 0;
svw_sum = 0;
ue_sum = 0;
uw_sum = 0;

for j=1:270
    for k=1:km
        if ~isnan(S_vw(j,k))
        svw_sum=svw_sum+S_vw(j,k)*DYG1d(j)*DRF(k);
        end
        if ~isnan(du_WE(j,k))
        du_sum=du_sum+du_WE(j,k)*DYG1d(j)*DRF(k);
        end
        if ~isnan(u_east(j,k))
        ue_sum=ue_sum+u_east(j,k)*DYG1d(j)*DRF(k);
        end
        if ~isnan(u_west(j,k))
        uw_sum=uw_sum+u_west(j,k)*DYG1d(j)*DRF(k);
        end
    end
end

du_sum=du_sum*1e-6;
ue_sum=ue_sum*1e-6;
uw_sum=uw_sum*1e-6;
svw_sum=svw_sum*1e-6;

disp(['du Ue Uw Svw'])
[du_sum ue_sum uw_sum svw_sum]


% 2) Correlation of S_vw and du_WE
x1=S_vw(:);
x2=du_WE(:);
x3=S_vw(:)-du_WE(:);
x1(isnan(x3))=[];
x2(isnan(x3))=[];
x3(isnan(x3))=[];

cc=corrcoef(x1,x2);
disp(['========== Date: ',datestr(time_ecco(nt)),' ============='])
disp(['RMS of (Svw-du) = ',num2str(rms(x3),3),'.  Correlation of Svw and du = ',num2str(cc(1,2))])





%% Compute meridional transport across IO lat-section

vdxdz3d_msk = zeros(im,jm,km);

for k=1:km
    vdxdz3d_msk(:,:,k) = vdx3d_msk(:,:,k) * DRF(k);
end

vol_v_yz = squeeze(nansum(vdxdz3d_msk,1))*1e-6;
vol_v_y  = squeeze(nansum(vol_v_yz,2));   % Integrated meridional transport in Sv





%% Streamfunction psi_up and psi_dn  LH: 2019-6-3 21:25:15

% Calculate streamfunction as AMOC definition: 
% two ways: 1) cumsum downward from top; 2) cumsum upward from bottom
psi_MOC_yz_dn0 =  cumsum(vol_v_yz,2,'omitnan');
psi_MOC_yz_up0 = -cumsum(vol_v_yz,2,'reverse','omitnan');

% Add the starting zero for both Psi
psi_MOC_yz_dn = [zeros(jm,1), psi_MOC_yz_dn0];
psi_MOC_yz_up = [psi_MOC_yz_up0, zeros(jm,1)];

psi_MOC_yz_dn(isnan(lat_v_1d),:) = [];
psi_MOC_yz_up(isnan(lat_v_1d),:) = [];
lat_v_nan = lat_v_1d;
lat_v_nan(isnan(lat_v_1d)) = [];


% Select domain of plotting Streamfunction
j0=1;   
j1=197;  




%% SOR solving Poisson equation

% make mask -- removed by LH on 2019-8-3 09:58:21


% Define sub-domain to implement SOR (Principle: make sure at least 5 land points included)
j0_SOR  = 61;    % ...61:land  62:water 
j1_SOR  = 197;   %  196:water 197:land...
jm_SOR0  = j1_SOR-j0_SOR+1;
km_SOR0  = km;  

S_SOR0 = zeros(jm_SOR0,km_SOR0);

% Define source term
source_term = 'S_v_w';   % '\deltau'=du_WE    'S_v_w'=S_vw
%source_term = '\Deltau';   

if strcmp(source_term,'\Deltau')
    S_SOR0 = du_WE(j0_SOR:j1_SOR,:);
elseif strcmp(source_term,'S_v_w')
    S_SOR0 = S_vw(j0_SOR:j1_SOR,:);
end


% LH 2019-6-18 15:42:35: Interp conservedly 
% Step-1) Svw interpolated in Y
dy0 = 1e4;
[S0,y0,ym,yf] = interp1d_conserve(S_SOR0(:,2),DYF1d(j0_SOR:j1_SOR),dy0);
jmi = length(S0);
Si_y = zeros(jmi,km_SOR0);
for k = 1:km_SOR0
    [S0,y0] = interp1d_conserve(S_SOR0(:,k),DYF1d(j0_SOR:j1_SOR),dy0);
    Si_y(:,k) = S0;
end
% Step-2) Svw interpolated in Z
dz0 = 10;
[S1,z0,zm,zf] = interp1d_conserve(Si_y(10,:),DRF,dz0);
kmi = length(S1);
Si_yz = zeros(jmi,kmi);
for j = 1:jmi
    [S1,z0] = interp1d_conserve(Si_y(j,:),DRF,dz0);
    Si_yz(j,:) = S1;
end


% For surface OBC
Wdx_SOR_surf = Wdx(j0_SOR:j1_SOR,1);
% Wdx Interpolation LH 2019-6-19 19:22:04
Wdx_SOR_surfi = interp1d_conserve(Wdx_SOR_surf,DYF1d(j0_SOR:j1_SOR),dy0);


% Make SORi mask
yz_msk_c_SOR = yz_msk_c(j0_SOR:j1_SOR,:) ;
yz_msk_c_SORi = zeros(jmi,kmi);

for j=1:jmi-1
    j0 = sum(y0(j)>=yf);    
    kc = nansum(yz_msk_c_SOR(j0,:));
    if kc>0; yz_msk_c_SORi(j,z0<=zf(kc+1)-dz0*0.5) = 1; end
end

%yz_msk_w_SORi(yz_msk_w_SORi==0)=nan;
%figure;pcolor(ym,zm,yz_msk_w_SOR');ylim([0 5670])
%figure;pcolor(y0(1:jmi)+0.5*dy0,z0(1:end-1),yz_msk_w_SORi');shading interp;ylim([0 5670])


%% Define new grid to be the dense interpolated

km_SOR = kmi;
jm_SOR = jmi;


S_SOR = Si_yz;
Phi = zeros(jm_SOR,km_SOR+1);



% stretch grid size
asp_y = 1/dy0;
asp_z = 1/dz0;
DYC_SOR = diff(y0);
DRC_SOR = diff(z0);
DYC_SOR = DYC_SOR * asp_y;
DRC_SOR = DRC_SOR * asp_z;  DRC_SOR(end+1) = DRC_SOR(end);
dh2 = DRC_SOR(1)*DYC_SOR(1);

% SOR iteration

omega = 1.995;  % optimum choice : 1.995
%n_ite = 4000;  % 1000 not stable, but 2000 is enough


tic
Phi0 = Phi;
Phi1 = Phi;
phi_err = [];
max_err = 1e-8;
diff_rms_err = 100;
n_ite = 1;
min_ite_steps = 2000;

yz_msk_Phi = [zeros(jm_SOR,1)  yz_msk_c_SORi]; % Add one layer for surface

N_water = zeros(jm_SOR,km_SOR);
for j=2:jm_SOR-1
    for k=2:km_SOR
        N_water(j,k) = yz_msk_Phi(j-1,k) + yz_msk_Phi(j+1,k) ...
                     + yz_msk_Phi(j,k-1) + yz_msk_Phi(j,k+1) ;
    end
end

phi_err = [];
while diff_rms_err > max_err || n_ite < min_ite_steps

    for k = 2:km_SOR
        for j = 2:jm_SOR-1
            
            if N_water(j,k)==0 || yz_msk_Phi(j,k)==0; continue; end
            
            tmp0 = Phi0(j-1,k)*yz_msk_Phi(j-1,k) + Phi0(j+1,k)*yz_msk_Phi(j+1,k) ...
                 + Phi0(j,k-1)*yz_msk_Phi(j,k-1) + Phi0(j,k+1)*yz_msk_Phi(j,k+1) ;
            if 1 && k==2
                tmp0 = Phi0(j-1,k)*yz_msk_Phi(j-1,k) + Phi0(j+1,k)*yz_msk_Phi(j+1,k) ...
                     + Phi0(j,k+1)*yz_msk_Phi(j,k+1) + Wdx_SOR_surfi(j)*DRC_SOR(1)*asp_z ;
            end
            
            tmp1 = tmp0 - dh2*S_SOR(j,k-1);

            Phi0(j,k) = (1-omega)*Phi0(j,k) + omega*tmp1/N_water(j,k);
        end
        
    end
    
    phi_err(n_ite) = rms(Phi1(:)-Phi0(:));
    if n_ite>1 
        diff_rms_err = abs(phi_err(n_ite)-phi_err(n_ite-1));
    end
    
    Phi1 = Phi0;

    n_ite=n_ite+1;
    
    %[n,phi_err(n)]
    
end

phi_err_all(nt) = phi_err(n_ite-1);
ite_num_all(nt) = n_ite-1;

Phi = Phi0;

disp(['omega=',num2str(omega),' Max_err=',num2str(max_err),' Min_steps=',num2str(min_ite_steps),'  Iteration steps=',num2str(n_ite)])

toc
        

%%

% Derive the Potential [v,w]
Vdx_pot = zeros(jm_SOR,km_SOR);
Wdx_pot = zeros(jm_SOR,km_SOR);

% mask for [V,W] in interpolated coordinate
yz_msk_Phi_v = zeros(jm_SOR,km_SOR);
yz_msk_Phi_w = zeros(jm_SOR,km_SOR);
for j=2:jm_SOR
    for k=1:km_SOR
        yz_msk_Phi_v(j,k) = yz_msk_Phi(j,k+1) * yz_msk_Phi(j-1,k+1);
        yz_msk_Phi_w(j,k) = yz_msk_Phi(j,k+1) * yz_msk_Phi(j,k);
    end
end
yz_msk_Phi_w(:,1)=yz_msk_Phi_w(:,2);


for k = 1:km_SOR
    for j = 2:jm_SOR
        Vdx_pot(j,k) =  (Phi(j,k+1)-Phi(j-1,k+1)) / DYC_SOR(j) / asp_y ...
            * yz_msk_Phi(j,k+1) * yz_msk_Phi(j-1,k+1);
    end
end
for k = 2:km_SOR
    for j = 1:jm_SOR
        Wdx_pot(j,k-1) = -(Phi(j,k)-Phi(j,k-1)) / DRC_SOR(k) / asp_z  ...
            * yz_msk_Phi_w(j,k-1);
    end
end
Wdx_pot(:,1) = yz_msk_Phi_w(:,1).*Wdx_SOR_surfi'; % Set surface w to be w-total as BC requires: LH 2019-8-23 14:18:09
Wdx_pot(:,km_SOR+1)=0;  % bottom floor

%
% Interpolate Vdx_tot LH 2019-6-18 19:22
% vertical interp requires conservation, so need to apply interp_conserve.func
Vi = interp1d_conserve(Vdx(100,:),DRF,dz0);
kmi_v = length(Vi);
Vi_z = zeros(jm_SOR0+1,kmi_v);
Wi_z = zeros(jm_SOR0+1,kmi_v);

Wdx_p1 = Wdx; Wdx_p1(:,km+1)=0;
DRC_10m = DRC; DRC_10m(1)=10;
for j = j0_SOR:j1_SOR+1
    Vi = interp1d_conserve(Vdx(j,:),DRF,dz0);
    Vi_z(j-j0_SOR+1,:) = Vi;
    Wi = interp1d_conserve(Wdx_p1(j,:),[DRC_10m; 2*(kmi*dz0-sum(DRC))],dz0); % to extend DRC to deeper than z0
    %Wi = interp1d_conserve(Wdx(j,:),DRF,dz0); % to extend DRC to deeper than z0
    Wi_z(j-j0_SOR+1,:) = Wi(1:km_SOR);
end

Vi = interp1(yf,Vi_z(:,5),y0); 
jmi_v = length(Vi);
Vdx_SOR = zeros(jmi_v,kmi_v);
Wdx_SOR = zeros(jmi_v,kmi_v);
for k = 1:kmi_v
    Vdx_SOR(:,k) = interp1(yf,Vi_z(:,k),y0); 
    Wdx_SOR(:,k) = interp1(yf,Wi_z(:,k),y0); 
end
Vdx_SOR(end,:) = [];
Wdx_SOR(end,:) = [];
Wdx_SOR(:,km_SOR+1)=0;  % bottom floor

% Both interp1_conserve
Vi = interp1d_conserve(Vi_z(1:end-1,5),DYF1d(j0_SOR:j1_SOR),dy0); 
Wi = interp1d_conserve(Wi_z(1:end-1,5),DYF1d(j0_SOR:j1_SOR),dy0); 
jmi_v = length(Vi);
Vdx_SOR3 = zeros(jmi_v,kmi_v);
Wdx_SOR3 = zeros(jmi_v,kmi_v);
for k = 1:kmi_v
    temp_v3 = interp1d_conserve(Vi_z(1:end,k),DYC1d(j0_SOR:j1_SOR+1),dy0); 
    Vdx_SOR3(:,k) = temp_v3(1:jmi); 
    Wdx_SOR3(:,k) = interp1d_conserve(Wi_z(1:end-1,k),DYF1d(j0_SOR:j1_SOR),dy0); 
end
%Vdx_SOR3(end,:) = [];
%Wdx_SOR3(end,:) = [];
Wdx_SOR3(:,km_SOR+1)=0;  % bottom floor


% 2D interpolation by interp2 directly
[zvi, yvi] = meshgrid(z0+5,y0);
Vdx_interp = Vdx(j0_SOR:j1_SOR+1,:);
Vdx_interp = [Vdx_interp , Vdx_interp(:,km)];
Vdx_SOR2 = interp2([reshape(zm,1,km) 6134.5],yf',Vdx_interp,zvi,yvi);
Vdx_SOR2(:,end)=[];
Vdx_SOR2(end,:)=[];

[zwi, ywi] = meshgrid(z0,y0);
Wdx_interp = Wdx(j0_SOR:j1_SOR+1,:);
Wdx_interp = [Wdx_interp , zeros(j1_SOR+2-j0_SOR,1)];
Wdx_SOR2 = interp2(zf,yf',Wdx_interp,zwi,ywi);
Wdx_SOR2(:,end)=[];
Wdx_SOR2(end,:)=[];
Wdx_SOR2(:,km_SOR+1)=0;  % bottom floor


% Calculate the residual -- which is streamfunction
Vdx_str = Vdx_SOR2 - Vdx_pot;
Wdx_str = Wdx_SOR2 - Wdx_pot;


Vdxdz_tot = zeros(jm,km);
for k=1:km
    Vdxdz_tot(:,k) = Vdx(:,k) * DRF(k) * 1e-6;
end
vol_v_tot  = squeeze(nansum(Vdxdz_tot(j0_SOR:j1_SOR,:),2));


Vdxdz_tot1 = zeros(jmi,kmi);
for k=1:kmi
    Vdxdz_tot1(:,k) = Vdx_SOR2(:,k) * dz0 * 1e-6;
end
vol_v_tot1  = squeeze(nansum(Vdxdz_tot1,2));


% calculate volume transport and stream-Psi

Vdxdz_pot = zeros(jm_SOR,km_SOR);

for k=1:km_SOR
    Vdxdz_pot(:,k) = Vdx_pot(:,k) * dz0 * 1e-6;
end
vol_y_pot  = squeeze(nansum(Vdxdz_pot,2));

%figure;plot(lat_v_nan,vol_y_pot(1:270))

Vdxdz_str = zeros(jm_SOR,km_SOR);
for k=1:km_SOR
    Vdxdz_str(:,k) = Vdx_str(:,k) * dz0 * 1e-6;
end
vol_v_str  = squeeze(nansum(Vdxdz_str,2));


% y0 corresponded Latitude
R_earth = 6.3781e6;   % in meters
pi_180 = pi/180;
dlat_y0 = dy0/R_earth/pi_180;
lat_SORi_v = lat_v_1d(j0_SOR) + dlat_y0*[0:jm_SOR-1];
lat_SORi_w = lat_v_1d(j0_SOR) + dlat_y0*[0:jm_SOR-1] + dlat_y0*0.5;

figure;
asp_char = [num2str(asp_y*10^-floor(log10(asp_y))),'x10^-^',num2str(-floor(log10(asp_y)))];
%plot(yf(1:jm_SOR0),vol_v_tot, y0(1:end-1),vol_y_pot,  'linewidth',3) ;
plot(lat_SORi_v,vol_v_tot1, lat_SORi_v,vol_y_pot,'-', lat_SORi_v,vol_v_str,  'linewidth',3) ;
title(['Meridional transport [rms(residual)=',num2str(rms(vol_v_str),3),'Sv]',newline,...
       '(Poisson parameters: \omega=',num2str(omega),', N-ite=',num2str(n_ite),...
       ', asp-y=',asp_char,', Source term:',source_term,')',newline,...
       ' [ ECCO(v4r3) Date:',datestr(time_ecco(nt)),' ]   dy0=',num2str(dy0*1e-3,3),'km',newline,...
       'max-error=10^-^',num2str(abs(log10(max_err))),' min-iteration-steps=',num2str(min_ite_steps)])
hold on;
plot(lat_v_1d(j0_SOR:j1_SOR),vol_v_tot,'k-')

legend('Total','Potential component','Residual','Original total',...
               'location','southeast','AutoUpdate','off');
ylabel('Sv'); xlabel('Latitude')
plot([lat_SORi_v(1) lat_SORi_v(end)],[0 0],'k--');
ylim([-140 40]); %xlim([-70 30])
        
%close

%close

%% Calculate streamfunction as AMOC definition: 
% two ways: 1) cumsum downward from top; 2) cumsum upward from bottom
psi_MOC_yz_dn0 =  cumsum(Vdxdz_str,2,'omitnan');
psi_MOC_yz_up0 = -cumsum(Vdxdz_str,2,'reverse','omitnan');

% Add the starting zero for both Psi
psi_MOC_yz_dn = [zeros(jm_SOR,1), psi_MOC_yz_dn0];
psi_MOC_yz_up = [psi_MOC_yz_up0, zeros(jm_SOR,1)];

%psi_MOC_yz_dn(isnan(lat_v_1d),:) = [];
%psi_MOC_yz_up(isnan(lat_v_1d),:) = [];
lat_v_nan = lat_v_1d;
lat_v_nan(isnan(lat_v_1d)) = [];

% Make mask for interpolated grids:  LH 2019-6-18 23:04:28
yz_msk_v_SORi_psi = [yz_msk_Phi_v zeros(jmi,1)];
yz_msk_v_SORi_psi_nan = yz_msk_v_SORi_psi;
yz_msk_v_SORi_psi_nan(yz_msk_v_SORi_psi_nan==0)=nan;


% Stretch surface layer
z0_strt = z0.^(0.5); 

% define labels in stretched coordinate
z0_strt_cell = {'6000','5000','4000','3000','2000','1000','500','100'} ;  % Note: must in descending order!
z0_strt_num  = str2num(char(z0_strt_cell)) ;
z0_strt_tick = interp1(z0, z0_strt, z0_strt_num);

figure;
    tmp_dep = -z0_strt;
    subplot(2,1,1)
      xx=psi_MOC_yz_dn.*yz_msk_v_SORi_psi_nan;
      contourf(lat_SORi_v,tmp_dep,xx',[-30:2:25]); 
      %contourf(xx',[-21:2:13]); 
      colorbar; 
      %contourf(j0_SOR:j1_SOR,tmp_dep,psi_MOC_yz_dn',16); colorbar; 
      set(gca,'ytick',-z0_strt_tick); set(gca,'yticklabel',z0_strt_cell)
      %set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
      ylabel('Depth/m')
      caxis([-30 15])
      cmocean('balance','pivot',0); 
      title(['\psi integrated top down',newline,' [ ECCO(v4r3) Date:',datestr(time_ecco(nt)),' ]'])
      hold on;
      %plot([0 0],[0 tmp_dep(end)],'w--', 'linewidth', 3 );
      
    subplot(2,1,2)
      tmp_dep = -z0_strt;
      xx=psi_MOC_yz_up.*yz_msk_v_SORi_psi_nan;
      contourf(lat_SORi_v,tmp_dep,xx',[-30:2:25]); 
      %contourf(y0(1:jmi),z0,xx',16); colorbar; 
      %contourf(xx',16); 
      colorbar; 
      set(gca,'ytick',-z0_strt_tick); set(gca,'yticklabel',z0_strt_cell)
      ylabel('Depth/m');  xlabel('Latitude')
      %set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
      caxis([-30 15])
      cmocean('balance','pivot',0); 
      title(['\psi integrated bottom up ',newline,...
          '(Poisson parameters: \omega=',num2str(omega),', N-ite=',num2str(n_ite),...
       ', asp-y=',asp_char,', Source term:',source_term,')'])
      %xlim([-35 25])
      hold on;
      %plot([0 0],[0 tmp_dep(end)],'w--', 'linewidth', 3 );
      
%close
      


%% Calculate the divergence and curl of the three velocities

%% 2019-8-22 17:31:42
% Convert the original Vdx_pot from the finer original grid to the coarser ECCO grid

Vdx_tot_i = Vdx_SOR3;
Wdx_tot_i = Wdx_SOR3;

% Linear interp in y
Vpot_ecco_y = zeros(jm_SOR0,kmi);
Vtot_ecco_y = zeros(jm_SOR0,kmi);
Wpot_ecco_y = zeros(jm_SOR0,kmi+1);
Wtot_ecco_y = zeros(jm_SOR0,kmi+1);
tic
for k=1:kmi
   Vpot_ecco_y(:,k) = interp1(y0(1:jmi),Vdx_pot(:,k),yf(1:jm_SOR0));
   Vtot_ecco_y(:,k) = interp1(y0(1:jmi),Vdx_tot_i(:,k),yf(1:jm_SOR0));
   Wpot_ecco_y(:,k) = interp1d_coarser_uneven(Wdx_pot(:,k),dy0,DYF1d(j0_SOR:j1_SOR));  % DYF(j0_SOR:j1_SOR)=diff(yf)
   Wtot_ecco_y(:,k) = interp1d_coarser_uneven(Wdx_tot_i(:,k),dy0,DYF1d(j0_SOR:j1_SOR));  % DYF(j0_SOR:j1_SOR)=diff(yf)
end
toc

Vpot_ecco_yz = zeros(jm_SOR0,km_SOR0);
Vtot_ecco_yz = zeros(jm_SOR0,km_SOR0);
Wpot_ecco_yz = zeros(jm_SOR0,km_SOR0+1);
Wtot_ecco_yz = zeros(jm_SOR0,km_SOR0+1);
tic
for j=1:jm_SOR0
   Vpot_ecco_yz(j,:) = interp1d_coarser_uneven(Vpot_ecco_y(j,:),dz0,DRF);
   Vtot_ecco_yz(j,:) = interp1d_coarser_uneven(Vtot_ecco_y(j,:),dz0,DRF);
   Wpot_ecco_yz(j,:) = interp1(z0,Wpot_ecco_y(j,:),zf);
   Wtot_ecco_yz(j,:) = interp1(z0,Wtot_ecco_y(j,:),zf);
end
toc

Vtot_SOR = Vdx(j0_SOR:j1_SOR,:);
Wtot_SOR = Wdx(j0_SOR:j1_SOR,:); Wtot_SOR(:,km+1)=0;

% volume transport along latitudes
Vtot_SOR_vol = nansum(Vtot_SOR.*DRF'*1e-6,2);
Vtot_ecco_yz_vol = nansum(Vtot_ecco_yz.*DRF'*1e-6,2);
Vpot_SOR_vol = nansum(Vdx_pot.*dz0*1e-6,2);
Vpot_ecco_yz_vol = nansum(Vpot_ecco_yz.*DRF'*1e-6,2);


% calculate Vstr_ecco
Vstr_ecco = Vtot_SOR - Vpot_ecco_yz;
Wstr_ecco = Wtot_SOR - Wpot_ecco_yz;


%% Calculte the divergence field

S_vw_pot = zeros(jm_SOR0,km_SOR0);
S_v_pot = S_vw_pot;
S_w_pot = S_vw_pot;
S_vw_str = zeros(jm_SOR0,km_SOR0);
S_v_str = S_vw_str;
S_w_str = S_vw_str;

% divergence of velocity field
for j=1:jm_SOR0-1
    for k=1:km_SOR0
        S_v_pot(j,k) = (Vpot_ecco_yz(j+1,k)-Vpot_ecco_yz(j,k)) / DYF1d(j+j0_SOR-1);
        S_w_pot(j,k) = (Wpot_ecco_yz(j,k)-Wpot_ecco_yz(j,k+1)) / DRF(k);
        S_v_str(j,k) = (Vstr_ecco(j+1,k)-Vstr_ecco(j,k)) / DYF1d(j+j0_SOR-1);
        S_w_str(j,k) = (Wstr_ecco(j,k)-Wstr_ecco(j,k+1)) / DRF(k);
    end
end
S_vw_pot = S_v_pot + S_w_pot;
S_vw_str = S_v_str + S_w_str;
S_vw_tot = S_vw(j0_SOR:j1_SOR,:);




%% Calculte the CURL field

curl_vw_tot = zeros(jm_SOR0,km_SOR0);
curl_v_tot = curl_vw_tot;
curl_w_tot = curl_vw_tot;
curl_vw_pot = zeros(jm_SOR0,km_SOR0);
curl_v_pot = curl_vw_pot;
curl_w_pot = curl_vw_pot;
curl_vw_str = zeros(jm_SOR0,km_SOR0);
curl_v_str = curl_vw_str;
curl_w_str = curl_vw_str;

% curl of velocity field
for j=2:jm_SOR0
    for k=2:km_SOR0
        curl_v_tot(j,k) = (Vtot_SOR(j,k)-Vtot_SOR(j,k-1)) / DRC(k);
        curl_w_tot(j,k) = (Wtot_SOR(j,k)-Wtot_SOR(j-1,k)) / DYC1d(j+j0_SOR-1);

        curl_v_pot(j,k) = (Vpot_ecco_yz(j,k)-Vpot_ecco_yz(j,k-1)) / DRC(k);
        curl_w_pot(j,k) = (Wpot_ecco_yz(j,k)-Wpot_ecco_yz(j-1,k)) / DYC1d(j+j0_SOR-1);

        curl_v_str(j,k) = (Vstr_ecco(j,k)-Vstr_ecco(j,k-1)) / DRC(k);
        curl_w_str(j,k) = (Wstr_ecco(j,k)-Wstr_ecco(j-1,k)) / DYC1d(j+j0_SOR-1);
    end
end
curl_vw_tot = curl_v_tot + curl_w_tot;
curl_vw_pot = curl_v_pot + curl_w_pot;
curl_vw_str = curl_v_str + curl_w_str;




%% Figure 5 paper Fig.5

% Stretch surface layer
zm_strt = zm.^(0.5); 
% define labels in stretched coordinate
zm_strt_cell = {'5000','4000','3000','2000','1000','500','100','0'} ;  % Note: must in descending order!
zm_strt_cell_num = {'5000','4000','3000','2000','1000','500','100','10'} ;  % Note: must in descending order!
zm_strt_num  = str2num(char(zm_strt_cell_num)) ;
zm_strt_tick = interp1(zm, zm_strt, zm_strt_num);

figure;
subplot(321)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,S_vw_tot');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,S_vw_tot');  
colorbar; 
caxis([-0.1 0.1])
%cmocean('balance','pivot',0);
colormap('redblue')
title('Div original')
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
set(gca,'layer','top')
ylabel('Depth/m')
grid on;

subplot(323)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,S_vw_pot');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,S_vw_pot');  
colorbar; 
caxis([-0.1 0.1])
%cmocean('balance','pivot',0);
colormap('redblue')
title('Div(\phi)')
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
set(gca,'layer','top')
ylabel('Depth/m')
grid on;

subplot(325)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,S_vw_str');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,S_vw_str');  
colorbar; 
caxis([-0.1 0.1])
%cmocean('balance','pivot',0);
colormap('redblue')
title('Div(\psi)')
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
set(gca,'layer','top')
ylabel('Depth/m')
grid on;


% Stretch surface layer
zf_strt = zf(1:km).^(0.5); 
% define labels in stretched coordinate
zf_strt_cell = {'500','400','300','200','100','0'} ;  % Note: must in descending order!
zf_strt_cell_num = {'500','400','300','200','100','10'} ;  % Note: must in descending order!
zf_strt_num  = str2num(char(zf_strt_cell_num)) ;
zf_strt_tick = interp1(zf(1:km), zf_strt, zf_strt_num);



subplot(322)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,curl_vw_tot');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,curl_vw_tot');  
colorbar; 
caxis([-1 1]*1e4)
%cmocean('balance','pivot',0);
colormap('redblue')
title('Curl original')
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
set(gca,'layer','top')
ylim([-sqrt(550) -sqrt(-1+str2num(zf_strt_cell_num{end}))])
ylabel('Depth/m')
grid on;

subplot(324)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,curl_vw_pot');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,curl_vw_pot');  
colorbar; 
caxis([-1 1]*1e4)
%cmocean('balance','pivot',0);
colormap('redblue')
title('Curl(\phi)')
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
set(gca,'layer','top')
ylim([-sqrt(550) -sqrt(-1+str2num(zf_strt_cell_num{end}))])
ylabel('Depth/m')
grid on;


subplot(326)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,curl_vw_str');  set(h,'edgecolor','none'); 
%h=contourf(lat_v_1d(j0_SOR:j1_SOR),zm,curl_vw_str');  
colorbar; 
caxis([-1 1]*1e4)
%cmocean('balance','pivot',0);
colormap('redblue')
title('Curl(\psi)')
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
set(gca,'layer','top')
%close 
ylim([-sqrt(550) -sqrt(-1+str2num(zf_strt_cell_num{end}))])
ylabel('Depth/m')
grid on;



%%


%% Calculate streamfunction with ECCO grid Vstr: 


Vdxdz_str = zeros(jm_SOR0,km_SOR0);
for k=1:km_SOR0
    Vdxdz_str(:,k) = Vstr_ecco(:,k) * DRF(k) * 1e-6;
end

% two ways: 1) cumsum downward from top; 2) cumsum upward from bottom
psi_MOC_yz_dn0 =  cumsum(Vdxdz_str,2,'omitnan');
psi_MOC_yz_up0 = -cumsum(Vdxdz_str,2,'reverse','omitnan');

% Add the starting zero for both Psi
psi_MOC_yz_dn = [zeros(jm_SOR0,1), psi_MOC_yz_dn0];
psi_MOC_yz_up = [psi_MOC_yz_up0, zeros(jm_SOR0,1)];

%psi_MOC_yz_dn(isnan(lat_v_1d),:) = [];
%psi_MOC_yz_up(isnan(lat_v_1d),:) = [];
lat_v_nan = lat_v_1d;
lat_v_nan(isnan(lat_v_1d)) = [];

% Make mask for interpolated grids:  LH 2019-6-18 23:04:28
yz_msk_v_SORi_psi = [yz_msk_Phi_v zeros(jmi,1)];
yz_msk_v_SORi_psi_nan = yz_msk_v_SORi_psi;
yz_msk_v_SORi_psi_nan(yz_msk_v_SORi_psi_nan==0)=nan;


% Stretch surface layer
zf_strt = zf.^(0.5); 

% define labels in stretched coordinate
zf_strt_cell = {'6000','5000','4000','3000','2000','1000','500','100'} ;  % Note: must in descending order!
zf_strt_num  = str2num(char(zf_strt_cell)) ;
zf_strt_tick = interp1(zf, zf_strt, zf_strt_num);

figure;
    tmp_dep = -zf_strt;
    subplot(2,1,1)
      xx=psi_MOC_yz_dn;%.*yz_msk_c_SOR;
      contourf(lat_v_1d(j0_SOR:j1_SOR),tmp_dep,xx',[-30:2:25]); 
      %contourf(xx',[-21:2:13]); 
      colorbar; 
      %contourf(j0_SOR:j1_SOR,tmp_dep,psi_MOC_yz_dn',16); colorbar; 
      set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
      %set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
      ylabel('Depth/m')
      caxis([-30 15])
      cmocean('balance','pivot',0); 
      title(['\psi integrated top down',newline,' [ ECCO(v4r3) Date:',datestr(time_ecco(nt)),' ]'])
      hold on;
      %plot([0 0],[0 tmp_dep(end)],'w--', 'linewidth', 3 );
      
    subplot(2,1,2)
      tmp_dep = -zf_strt;
      xx=psi_MOC_yz_up;%.*yz_msk_c_SOR;
      contourf(lat_v_1d(j0_SOR:j1_SOR),tmp_dep,xx',[-30:2:25]); 
      %contourf(y0(1:jmi),zf,xx',16); colorbar; 
      %contourf(xx',16); 
      colorbar; 
      set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
      ylabel('Depth/m');  xlabel('Latitude')
      %set(gca,'xtick',lat_label_num); set(gca,'xticklabel',lat_label_cell)
      caxis([-30 15])
      cmocean('balance','pivot',0); 
      title(['\psi integrated bottom up ',newline,...
          '(Poisson parameters: \omega=',num2str(omega),', N-ite=',num2str(n_ite),...
       ', asp-y=',asp_char,', Source term:',source_term,')'])
      %xlim([-35 25])
      hold on;
      %plot([0 0],[0 tmp_dep(end)],'w--', 'linewidth', 3 );
      
%close
      




%% V W components 3x2 plots

% Stretch surface layer
zm = -mygrid.RC;
zm_strt = zm.^(0.5); 
% define labels in stretched coordinate
zm_strt_cell = {'5900','5000','4000','3000','2000','1000','500','100'} ;  % Note: must in descending order!
zm_strt_num  = str2num(char(zm_strt_cell)) ;
zm_strt_tick = interp1(zm, zm_strt, zm_strt_num);

figure;
subplot(321)
%h=pcolor(lat_v_1d(j0_SOR:j1_SOR),zm,Vstr_ecco_yz');set(h,'edgecolor','none'); colorbar;
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,Vtot_SOR');set(h,'edgecolor','none'); colorbar;
caxis([-5 5]*1e4);   
grid;
colormap('redblue')
%ylim([0 1000]); 
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
set(gca,'layer','top');
title('V-tot')

subplot(323)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,Vpot_ecco_yz');set(h,'edgecolor','none'); colorbar;
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
grid ; set(gca,'layer','top')
caxis([-5 5]*1e4);   
colormap('redblue')
title('V-pot ')

subplot(325)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zm_strt,Vstr_ecco');set(h,'edgecolor','none'); colorbar;
set(gca,'ytick',-zm_strt_tick); set(gca,'yticklabel',zm_strt_cell)
grid ; set(gca,'layer','top')
caxis([-5 5]*1e4);   
colormap('redblue')
title('V-str ECCO')




subplot(322)
%h=pcolor(lat_v_1d(j0_SOR:j1_SOR),zm,Vstr_ecco_yz');set(h,'edgecolor','none'); colorbar;
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,Wtot_SOR');set(h,'edgecolor','none'); colorbar;
caxis([-5 5]*1e1);   
colormap('redblue')
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
grid ; set(gca,'layer','top');
title('W-tot')

subplot(324)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,Wpot_ecco_yz');set(h,'edgecolor','none'); colorbar;
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
grid ; set(gca,'layer','top')
caxis([-5 5]*1e1);   
colormap('redblue')
title('W-pot ')


subplot(326)
h=pcolor(lat_v_1d(j0_SOR:j1_SOR),-zf_strt,Wstr_ecco');set(h,'edgecolor','none'); colorbar;
set(gca,'ytick',-zf_strt_tick); set(gca,'yticklabel',zf_strt_cell)
grid ; set(gca,'layer','top')
caxis([-5 5]*1e1);   
colormap('redblue')
title(['W-str ECCO'])





