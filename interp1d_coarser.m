function [S1,zf] = interp1d_coarser(S0,dz0,dz1)

% must dz1>=dz0! (coarser)
% This function deals with re-interpolation1d-conserved with
% interp1d-conserved uniform grid data, to coarser grid, to be applied in
% the sigma-coordinate. LH created on 2019-8-3 22:35:26
% Input data S0 is vector under uniform grid with spacing dz0, dz1 is
% coarser grid spacing


im0 = length(S0);
len0 = dz0*im0;
z0 = 0:dz0:len0;

z1 = 0:dz1:len0;
if z1(end)<len0; z1 = [z1 z1(end)+dz1]; end

im1 = length(z1);
S1 = zeros(im1-1,1);

for i=1:im1-1
    z1_0 = z1(i);
    z1_1 = z1(i+1);
    j0 = sum(z0<=z1_0);  % the index of i-th z1 grid in z0
    j1 = sum(z0<=z1_1);
    
    dz_j0 = z0(j0+1)-z1_0;
    dz_j1 = z1_1 - z0(j1);
    
    if j1<=im0
        sum_S0_dz1 = S0(j0)*dz_j0 + S0(j1)*dz_j1;
    else
        sum_S0_dz1 = S0(j0)*dz_j0;
    end
    if j1-j0>1
        sum_S0_dz1 = sum_S0_dz1 + sum(S0(j0+1:j1-1)*dz0);
    end
    S1(i) = sum_S0_dz1/dz1;
end

zf = z1;




