%------ code by Joao Mendonca ---------------------------------------------

function mjolnir

format long;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Options
%
% nview plot
% 1     Averaged (time and longitude) zonal winds.
%       2D Map Latitude Vs Pressure.
%
% 2     Averaged (time and longitude) temperatures.
%       2D Map Latitude Vs Pressure.
%
% 3     Averaged (time) temperature and wind field.
%       2D Map Longitude Vs Latitude.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nview    = 3  ; % type of plot
ntsi     = 10 ; % initial file id number
nts      = 10 ; % last file id number

simulation_ID = 'Earth';

%%%%%%%%%%
% Planet %
%%%%%%%%%%

fileh5 = ['../results/esp_output_' simulation_ID '.h5'];
% h5disp(fileh5)
A            = h5read(fileh5,'/A'           );
Rd           = h5read(fileh5,'/Rd'          );
Omega        = h5read(fileh5,'/Omega'       );
P_Ref        = h5read(fileh5,'/P_Ref'       );
Top_altitude = h5read(fileh5,'/Top_altitude');
Cp           = h5read(fileh5,'/Cp'          );

%%%%%%%%
% Grid %
%%%%%%%%

fileh5 = ['../results/esp_output_grid_' simulation_ID '.h5'];
% h5disp(fileh5)
Altitude  = h5read(fileh5,'/Altitude' );
Altitudeh = h5read(fileh5,'/Altitudeh');
areasT    = h5read(fileh5,'/areasT');
lonlat    = h5read(fileh5,'/lonlat');
point_num = h5read(fileh5,'/point_num');
nv        = h5read(fileh5,'/nv');
nvi= nv+1;

lon = zeros(point_num,1);
lat = zeros(point_num,1);
for i = 1:point_num
   lon(i) = lonlat((i-1)*2 + 1);
   if(lon(i) < 0)
       lon(i) = lonlat((i-1)*2 + 1) + 2*pi;
   end
   lat(i) = lonlat((i-1)*2 + 2);
end

%%%%%%%%%%%%%%%
% Diagnostics %
% %%%%%%%%%%%%%


% Initialize arrays
Rho      = zeros(point_num,nv,ntsi-nts+1);
Pressure = zeros(point_num,nv,ntsi-nts+1);
Mh       = zeros(3,point_num,nv,ntsi-nts+1);
Wh       = zeros(point_num,nvi,ntsi-nts+1);

% Read model results
for t = ntsi:nts
    ts= num2str(t);
    fileh5 = ['../results/esp_output_' simulation_ID '_' ts '.h5'];
    % h5disp(fileh5)
    Rhoi     = h5read(fileh5,'/Rho' );
    Pressurei= h5read(fileh5,'/Pressure' );
    Mhi      = h5read(fileh5,'/Mh' );
    Whi      = h5read(fileh5,'/Wh' );

    for i = 1:point_num
        for j = 1:nv
            Rho(i,j,t-ntsi+1) = Rhoi((i-1)*nv + j);
            Pressure(i,j,t-ntsi+1) = Pressurei((i-1)*nv + j);
            for k = 1:3
                Mh(k,i,j,t-ntsi+1) = Mhi((i-1)*nv*3 + (j-1)*3 + k);
            end
        end
    end

    for i = 1:point_num
        for j = 1:nv+1
            Wh(i,j,t-ntsi+1) = Whi((i-1)*nvi + j);
        end
    end
end

%%%%%%%%%
% Plots %
%%%%%%%%%

% Sigma values for the plotting.
sigmaref  = [1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3 0.2 0.1 0.05];

if(nview == 1)
    % Averaged Zonal winds (latitude vs pressure)
    u(Mh, Rho, Pressure, lon, lat, point_num, nts-ntsi+1, P_Ref, nv, sigmaref);
elseif(nview == 2)
    % Averaged temperature (latitude vs pressure)
    temperature(Rho, Pressure, lon, lat, point_num, nts-ntsi+1, Rd, P_Ref, nv, sigmaref);
elseif(nview == 3)
    % Averaged temperature and wind field (longitude vs latitude)
    % PR_LV - Pressure level (Pa)
    PR_LV = 25000.0;
    temperature_u_lev(Mh, Rho, Pressure, lon, lat, point_num, nts-ntsi+1, Rd, nv, PR_LV);
end

end
