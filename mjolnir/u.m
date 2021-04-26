function u(Mh, Rho, Pressure, lon, lat, num, tsp, ps0, nv, sigmaref)

format long
 
% Set the reference pressure.
Pref = ps0*sigmaref;
d_sig = size(sigmaref);

% Set the latitude-longitude grid.
res_deg = 0.005;
[loni,lati] = meshgrid(0:res_deg:2*pi, -pi/2:res_deg:pi/2);
d_lon = size(loni);

%%%%%%%%%%%%%%%%%%
% Zonal Momentum %
%%%%%%%%%%%%%%%%%%

% Initialize arrays
ZonalMii = zeros(nv);
ZonalMi  = zeros(num, d_sig(2));
ZonalM   = zeros(d_lon(1),d_lon(2),d_sig(2), tsp);

% Compute zonal winds
for t = 1:tsp
    for i = 1:num
        sigma = Pressure(i,:,t);
        for lev = 1:nv
            ZonalMii(lev) = (Mh(1,i,lev,t) * (-sin(lon(i))) +...
                             Mh(2,i,lev,t) * ( cos(lon(i))) +...
                             Mh(3,i,lev,t) * (0))/Rho(i,lev,t);        
        end
        % Interpolate atmospheric column to the reference pressure.
        aux = interp1(sigma, ZonalMii, Pref,'PCHIP');
        for lev = 1:d_sig(2)
            ZonalMi(i,lev) = aux(lev); 
        end
    end
    % Convert icosahedral grid into lon-lat grid
    for lev = 1:d_sig(2)
        ZonalM(:,:,lev,t) = griddata(lon,lat, ZonalMi(:,lev), loni, lati,'nearest');
    end
end
    
clear zonalMii ZonalMi;

% Initialize arrays
d_z = size(ZonalM);
ZonalMl  = zeros(d_z(1), d_sig(2), tsp);
ZonalMlt = zeros(d_z(1), d_sig(2));

% Averaging in time and longitude.
if(tsp>1) 
    for j = 1:d_z(1)
        for lev = 1:d_sig(2)
            for t = 1:d_z(4)
                ZonalMl(j,lev,t) = mean(ZonalM(j,:,lev,t));
            end
        end
    end
    clear ZonalM;
    for lev = 1:d_sig(2)
        for j = 1:d_z(1)
            ZonalMlt(j,lev)  = mean(ZonalMl(j,lev,:));
        end
    end
    clear ZonalMl;
else
    for j = 1:d_z(1)
        for lev = 1:d_sig(2)
            ZonalMlt(j,lev) = mean(ZonalM(j,:,lev));
        end
    end
    clear ZonalM;  
end

%%%%%%%%%%%%%%%%%
% Create figure %
%%%%%%%%%%%%%%%%%

% Latitude
latp = -pi/2:res_deg:pi/2;

% Contour plot.
[C,h]=contourf(latp*180/pi, Pref/100, ZonalMlt', 40,'LineColor','none');
h = colorbar;
grid on
colormap(jet);
caxis([-10 30]);
set(h,'ytick',[-10:5:30]);
set(gca,'YTick',[0:100:1000])
set(gca,'XTick',[-80:20:80]) 
set(get(gca,'XLabel'),'String','Latitude (deg)','FontSize',20);
set(gca,'YDir','reverse');
set(get(gca,'YLabel'),'String','Pressure (mba)','FontSize',20);
set(gca,'FontSize',20);
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);
ylim([1 1000]);
xlim([-87 87]);

%%%%%%%%%
% Print %
%%%%%%%%%

print('-djpeg100','-zbuffer',['figures/u_ref.jpg']);
% print('-depsc2','-zbuffer','-r150',['figures/u_ref.eps']);

end
