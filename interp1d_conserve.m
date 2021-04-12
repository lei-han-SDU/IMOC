

function [S0,y0,ym,y] = interp1d_conserve(S_SOR,DYC_SOR,dy0)

% DYC_SOR(i): integral domain of S_SOR(i)
% y0(i): absolute corner distance of grid S0(i) from starting point y0(1)=0;
% ym(i): absolute distance of median point of grid S_SOR(i)
% y(i): absolute corner distance of grid S_SOR(i)


    S = S_SOR;%  (2:end-1,5);
    dy = DYC_SOR; %(2:end-1);

    y = [0;reshape(cumsum(dy),length(dy),1)];
    jy = length(y);

    ym = zeros(1,jy-1);
    for j=1:jy-1
      ym(j) = y(j)+0.5*dy(j);
    end


    %dy0 = 1e2;
    y0 = 0:dy0:y(jy);
    if y0(end)<y(jy); y0 = 0:dy0:y(jy)+dy0; end  % make sure y0 covers y LH 2019-7-28
        
    jy0 = length(y0);

    y0i = zeros(jy0,1);
    for j=1:jy0
        y0i(j) = sum(y-y0(j)<=0);
    end

    i_crs = [diff(y0i);0];
    i_crs(end-1)=0;


    S0 = zeros(1,jy0-1);
    for j=1:jy0-1
    %for j=1:7
        if i_crs(j)==0
            S0(j) = S(y0i(j));
        elseif i_crs(j)==1
            yi = y0i(j);
            if yi==jy-1; continue; end
            y_edge = y(yi+1);
            S1 = S(yi)*(y_edge-y0(j));
            S2 = S(yi+1)*(y0(j+1)-y_edge);
            S0(j) = (S1+S2)/dy0;
        end
    end

    %figure;plot(ym,S,y0(1:end-1),S0)
    
    %[sum(S.*dy') sum(S.*dy')-sum(S0.*dy0)]
    %sum(S.*dy')-sum(S0.*dy0)


end




