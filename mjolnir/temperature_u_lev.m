function temperature_u_lev(Mh, Rho, Pressure, lon, lat, num, tsp, Rd, nv, Pref)

format long

% Set the latitude-longitude grid-
res_deg = 0.005;
[loni,lati] = meshgrid(0:res_deg:2*pi, -pi/2:res_deg:pi/2);
d_lon = size(loni);

%%%%%%%%%%%%%%%%%%%%%%%
% Winds & temperature %
%%%%%%%%%%%%%%%%%%%%%%%

% Initialize arrays
Temp = zeros(nv,1);
Pr   = zeros(nv,1);
Mx   = zeros(nv,1);
My   = zeros(nv,1);
Mz   = zeros(nv,1);
Tempi= zeros(num,1);
Rhot = zeros(num,1);
Mxf  = zeros(num,1);
Myf  = zeros(num,1);
Mzf  = zeros(num,1);
Ui   = zeros(num,1);
Vi   = zeros(num,1);
Temperatureii = zeros(d_lon(1), d_lon(2), tsp);
Uii  = zeros(d_lon(1), d_lon(2), tsp);
Vii  = zeros(d_lon(1), d_lon(2), tsp);

% Compute winds and temperatures 
for t = 1:tsp 
    for i = 1:num
        for lev = 1:nv
             Pr(lev)  = Pressure(i,lev,t);
             Temp(lev)= Pressure(i,lev,t)/(Rd*Rho(i,lev,t));
             Mx(lev) = Mh(1,i,lev,t);
             My(lev) = Mh(2,i,lev,t);
             Mz(lev) = Mh(3,i,lev,t);
        end  
        % Interpolate in pressure.
        Rhot(i) = interp1(Pr, Rho(i,:,t), Pref,'linear');
        Mxf(i)  = interp1(Pr, Mx, Pref,'linear');
        Myf(i)  = interp1(Pr, My, Pref,'linear');
        Mzf(i)  = interp1(Pr, Mz, Pref,'linear');
        Tempi(i)= interp1(Pr, Temp      , Pref,'linear');   
    end
    
    % Convert icosahedral grid into lon-lat grid  
    for i = 1:num            
            Ui(i) =(Mxf(i) * (-sin(lon(i))) +...
                    Myf(i) * ( cos(lon(i))) +...
                    Mzf(i) * (0))/Rhot(i); 
            Vi(i) =(Mxf(i) * (-sin(lat(i))*cos(lon(i))) +...
                    Myf(i) * (-sin(lat(i))*sin(lon(i))) +...
                    Mzf(i) *  cos(lat(i)))/Rhot(i);          
    end
%     Temperatureii(:,:,t) = griddata(lon, lat, Tempi, loni, lati, 'nearest');
%     Uii(:,:,t)           = griddata(lon, lat, Ui   , loni, lati, 'nearest');
%     Vii(:,:,t)           = griddata(lon, lat, Vi   , loni, lati, 'nearest');
    Temperatureii(:,:,t) = griddata(lon, lat, Tempi, loni, lati, 'cubic');
    Uii(:,:,t)           = griddata(lon, lat, Ui   , loni, lati, 'cubic');
    Vii(:,:,t)           = griddata(lon, lat, Vi   , loni, lati, 'cubic');   
    
end

% Initialize arrays
d_z = size(Temperatureii);
Temperature  = zeros(d_z(1), d_z(2));
Uiii = zeros(d_z(1), d_z(2));
Viii = zeros(d_z(1), d_z(2));

% Averaging in time and longitude.
if(tsp>1) 
    for i = 1:d_z(1)
        for j = 1:d_z(2) 
            Temperature(i,j) = mean(Temperatureii(i,j,:));
            Uiii(i,j)        = mean(Uii(i,j,:))          ;        
            Viii(i,j)        = mean(Vii(i,j,:))          ;
        end
    end
else
    for i = 1:d_z(1)
        for j = 1:d_z(2) 
            Temperature(i,j) = Temperatureii(i,j);
            Uiii(i,j)        = Uii(i,j)          ;        
            Viii(i,j)        = Vii(i,j)          ;
        end
    end
end
% Clear memory
clear Tempi Uii Vii Temperaturei Ui Vi Rhot;

% Wind arrays
lonqi = 0:res_deg:2*pi    ;
latqi = -pi/2:res_deg:pi/2;
d_z  = size(Uiii)      ;
spacing = 80;
count = 1              ;
for i = 1:spacing:d_z(1)
    for j = 1:spacing:d_z(2)
        U(count)    = Uiii(i,j);
        V(count)    = Viii(i,j);
        lonq(count) = lonqi(j)*180/pi;
        latq(count) = latqi(i)*180/pi;   
        count = count + 1;
    end
end
clear Uiii Viii;

%%%%%%%%%%%%%%%%%
% Create figure %
%%%%%%%%%%%%%%%%%

% Latitude and Longitude
latp = -pi/2:res_deg:pi/2;
lonp = 0    :res_deg:2*pi;

hold on
[C,h]=contourf(lonp*180/pi,latp*180/pi,Temperature,50,'LineColor','none');
h = colorbar;
colormap(jet);
grid on
set(gca,'XTick',[0:60:360])
set(gca,'YTick',[-80:20:80]) 
set(get(gca,'XLabel'),'String','Longitude (deg)','FontSize',20);
set(get(gca,'YLabel'),'String','Latitude (deg)','FontSize',20);
set(gca,'FontSize',20);
set(get(gca,'XLabel'),'FontSize',20);
set(get(gca,'YLabel'),'FontSize',20);

quiver(lonq, latq, U, V, 2, 'k');
xlim([0 360]);
ylim([-87 87]);
hold off;

%%%%%%%%%
% Print %
%%%%%%%%%

print('-djpeg100','-zbuffer',['figures/temperature-uv_lev.jpg']);
% print('-depsc2','-zbuffer','-r150',['figures/temperature-uv_lev.eps']);

end
