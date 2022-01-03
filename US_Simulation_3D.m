close all, clc
clearvars;

% =========================================================================
%  MOUSE SETUP
% =========================================================================


close all, clc
clearvars;
%% Create the computational GRID
% x = Nx*dx

PML_size = 10;

Nx = 432 - 2*PML_size;            % number of grid points in the x direction
Ny = 432 - 2*PML_size;            % number of grid points in the y direction
Nz = 800 - 2*PML_size;            % number of grid points in the z direction

x = 70e-3;
y = 70e-3;
z = 200e-3;
d = 0.2e-3;
dx = d;        % grid point spacing in the x direction [m]
dy = d;        % grid point spacing in the y direction [m]
dz = d;        % grid point spacing in the z direction [m]

kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%% Define the properties of the propagation MEDIUM
% Water
c0 = 1480;
rho0 = 1000;

[medium,setup] = Acoustic_prop(c0,rho0,Nx,Ny,Nz,d);


%% Set the TIME STEP
c_ref = c0;
c_max = 2800;
t_end = 300e-6;          % [s]
CFL = 0.3;
kgrid.t_array = makeTime(kgrid, c_max, CFL, t_end);

%% Create the TRANSDUCER

% Source mask
[transducer,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\1 mm models\Trasduttore.nrrd');
% Pre - processing
%
% transducer = cat(1,transducer,zeros(364 - size(transducer,1), size(transducer,2), size(transducer,3)));
% transducer = cat(2,transducer,zeros(size(transducer,1), 364 - size(transducer,2), size(transducer,3)));
% transducer = cat(3,transducer, zeros(size(transducer,1), size(transducer,2), 1004 - size(transducer,3)));
transducer = single(transducer);
%
source.p_mask = zeros(size(transducer));
index = find(transducer);
source.p_mask(index) = 1;
% voxelPlot(source.p_mask), title('Transducer')

% Source input signal
source_freq = 500e3;   % [Hz]
source_mag = 1;      % [Pa]
source.p = toneBurst(10e6,source_freq,10); 
% source.p = source_mag * sin(2 * pi * source_freq * kgrid.t_array(1:700));
figure,
plot(kgrid.t_array, [source.p , zeros(1,length(kgrid.t_array) - length(source.p))])
title('Input pressure signal')
xlabel('Time - [s]'), ylabel('Pressure - [Normalized]')

%% Define the SENSOR mask
sensor.mask = ones(size(transducer));
% sensor.mask(round(size(transducer,1)/2),:,:) = 1;
% % voxelPlot(sensor.mask), title('Sensor Mask')
sensor.record = {'p_final','p_max','p_rms'};

%% SIMULATION
% voxelPlot(mask + source.p_mask), title('Overall Setup Structure')
% input arguments
input_args = {'PlotSim', false,'PMLInside',false,...
     'PlotPML',false, 'PMLSize',PML_size};
% run the simulation
filename = 'C:\Users\Gabri\Desktop\Mouse_Sim.h5';
sensor_data = kspaceFirstOrder3D(kgrid, medium, source, sensor, input_args{:},'SaveToDisk',filename);

sensor_data.p_max = reshape(sensor_data.p_max, [Nx, Nz]);
sensor_data.p_rms = reshape(sensor_data.p_rms, [Nx, Nz]);

%% Plot the BEAM PATTERN of the Max Pressure
x_beam = kgrid.x_vec;
y_beam = [kgrid.y_vec(1)*2:dy:kgrid.y_vec(end)*2+dy]';

figure;
ax = gca;
imagesc(y_beam* 1e3, (x_beam - min(x_beam)) * 1e3  , sensor_data.p_max);
xlabel(['x-position [mm]']);
ylabel('z-position [mm]');
title('Total Beam Pattern Using MAX Of Recorded Pressure');
colormap(jet)
c = colorbar;
ax.YDir = 'normal';
ylabel(c, 'Pressure [Normalized]');
axis image;

