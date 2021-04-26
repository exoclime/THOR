function temperature(Rho, Pressure, lon, lat, num, tsp, Rd, ps0, nv, sigmaref)

format long
 
% Set the reference pressure.
Pref = ps0*sigmaref;
d_sig = size(sigmaref);

% Set the latitude-longitude grid.
res_deg = 0.005;
[loni,lati] = meshgrid(0:res_deg:2*pi, -pi/2:res_deg:pi/2);
d_lon = size(loni);

%%%%%%%%%%%%%%%
% Temperature %
%%%%%%%%%%%%%%%

% Initialize arrays
Temperatureii = zeros(nv);
Temperaturei  = zeros(num, d_sig(2));
Temperature   = zeros(d_lon(1),d_lon(2),d_sig(2), tsp);

% Compute temperatures
for t = 1:tsp
    for i = 1:num
        sigma = Pressure(i,:,t);
        for lev = 1:nv
            Temperatureii(lev) = Pressure(i,lev,t)/(Rd*Rho(i,lev,t));                
        end
        % Interpolate atmospheric column to the reference pressure.
        aux = interp1(sigma, Temperatureii, Pref,'PCHIP');
        for lev = 1:d_sig(2)
            Temperaturei(i,lev) = aux(lev); 
        end
    end
    % Convert icosahedral grid into lon-lat grid
    for lev = 1:d_sig(2)
        Temperature(:,:,lev,t) = griddata(lon, lat, Temperaturei(:,lev), loni, lati,'nearest');
    end
end
    
clear Temperatureii Temperaturei;

% Initialize arrays
d_z = size(Temperature);
Temperaturel  = zeros(d_z(1), d_sig(2), tsp);
Temperaturelt = zeros(d_z(1), d_sig(2));

% Averaging in time and longitude.
if(tsp>1) 
    for j = 1:d_z(1)
        for lev = 1:d_sig(2)
            for t = 1:d_z(4)
                Temperaturel(j,lev,t) = mean(Temperature(j,:,lev,t));
            end
        end
    end
    clear Temperature;
    for lev = 1:d_sig(2)
        for j = 1:d_z(1)
            Temperaturelt(j,lev)  = mean(Temperaturel(j,lev,:));
        end
    end
    clear Temperaturel;
else
    for j = 1:d_z(1)
        for lev = 1:d_sig(2)
            Temperaturelt(j,lev) = mean(Temperature(j,:,lev));
        end
    end
    clear Temperature;  
end

%%%%%%%%%%%%%%%%%
% Create figure %
%%%%%%%%%%%%%%%%%

% Latitude
latp = -pi/2:res_deg:pi/2;

% Contour plot.
[C,h]=contourf(latp*180/pi, Pref/100, Temperaturelt', 40,'LineColor','none');
h = colorbar;
grid on
colormap(jet);
caxis([180 310]);
set(h,'ytick',[190:10:300]);
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

print('-djpeg100','-zbuffer',['figures/temperature_ref.jpg']);
% print('-depsc2','-zbuffer','-r150',['figures/temperature_ref.eps']);

end
