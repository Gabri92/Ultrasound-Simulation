% Script per la lettura dell'output della simulazione

clc, clear all, close all

%% Take the file from PC

% % [transducer,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\0.2 mm models\Trasduttore.nrrd');
% % [mouse,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\0.2 mm models\Mouse_Segmentation.nrrd');
% % [platform,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\0.2 mm models\Piattaforma.nrrd');
% % index = find(platform);
% % platform(index) = 12;
% % clear indexc
% % index = find(transducer);
% % transducer(index) = 12;
% % clear index
% % setup = transducer + mouse + platform;
% % setup = cat(1,setup,zeros(364 - size(setup,1), size(setup,2), size(setup,3)));
% % setup = cat(2,setup,zeros(size(setup,1), 364 - size(setup,2), size(setup,3)));
% % setup = cat(3,setup, zeros(size(setup,1), size(setup,2), 1004 - size(setup,3)));
% % setup = squeeze(setup(round(size(setup,1)/2),:,:));
% % setup = single(setup)';

load("C:\Users\Gabri\Desktop\Simulazioni 3D Risultati\750 KHz\simGrid_750KHz.mat");
struct = simGrid_750KHz;
clear simGrid_750KHz

%--------------------------------------------------------------------------
piattaforma = zeros(size(struct));
index = find(struct == 8);
piattaforma(index) = 1;
index = find(struct == 100);
piattaforma(index) = 1;
index = find(struct == 101);
piattaforma(index) = 1;
index = find(struct == 102);
piattaforma(index) = 1;
piattaforma = cat(1,piattaforma,zeros(412 - size(piattaforma,1),...
    size(piattaforma,2), size(piattaforma,3)));
piattaforma = cat(2,piattaforma,zeros(size(piattaforma,1),...
    412 - size(piattaforma,2), size(piattaforma,3)));
piattaforma = cat(3,piattaforma, zeros(size(piattaforma,1),...
    size(piattaforma,2), 780 - size(piattaforma,3)));
%--------------------------------------------------------------------------
filepath = "C:\Users\Gabri\Desktop\";
% filename = "C:\Users\Gabri\Desktop\0.2mm_withAir.h5";   %% In acqua
% filename = "C:\Users\Gabri\Desktop\Simulazioni topo\Simulazioni 0.2mm\With air\Mouse_Sim_Res_wAir.h5"; %% Con aria
filename = "C:\Users\Gabri\Desktop\Mouse_Sim_Res_750KHz.h5"; %% Con aria


%% create the computational grid  - Automatizzare esportando la griglia
PML_size = 10;
Nx = 432 - 2*PML_size;            % number of grid points in the x direction
Ny = 432 - 2*PML_size;            % number of grid points in the y direction
Nz = 1125 - 2*PML_size;            % number of grid points in the z direction
d = 0.15e-3;
dx = d;        % grid point spacing in the x direction [m]
dy = d;        % grid point spacing in the y direction [m]
dz = d;        % grid point spacing in the z direction [m]
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%% Re-order the data
sensor_data.p_max = h5read(filename,"/p_max");
sensor_data.p_max = reshape(sensor_data.p_max, [Nx Ny Nz]);
sensor_data.p_max = sensor_data.p_max(Nx/2,:,:);
data = (squeeze(sensor_data.p_max))';
data = data./max(max(data));
% temp = sensor_data.p_max(1:Nz,1:Nx);
% data = temp;
% clear sensor_data temp
% sensor_data.p_max = sensor_data.p_max./max(max(sensor_data.p_max));
% setup = setup./10;
% % index = find(setup == 0);
% % setup(index) = data(index);
% % clear index
% % index = find(setup == 12);
% % data(index) = max(max(data))+0.01;
% % clear index
% temp = zeros(160,size(setup,1));
% temp = setup(600:760,:);
% clear setup
% setup = temp;
% clear temp
% setup = imrotate(setup,-90);
% imagesc(setup)
% colormap([[parula];[0 0 0]])

%% PLOT
d = 0.15;
x_beam = kgrid.x_vec;
z_beam = [kgrid.z_vec(1):dz:kgrid.z_vec(end)+dz]';
% setup_interval = [600:760]*d;

figure
% ax = axes('XLim',[-Nx/2 Nx/2]*d,'YLim',[0 Nz]*d);
im_beam = imshow(data(:,54:353),'XData',[-Nx/2+54 Nx/2-(412 - 353)]*d,'YData',[0 Nz]*d);
ax = gca;
% ax.XLim = im_beam.XData
% ax.YLim = im_beam.YData
xlabel('x-position [mm]');
ylabel('z-position [mm]');
title('Max Recorded Pressure');
colormap([[jet];[0 0 0]])
c = colorbar;
ax.YDir = 'normal';
ylabel(c, 'Pressure [Normalized]');
hold on
% im_setup = imagesc(x_beam*1e3, setup_interval*1e3, setup)
% set(im_setup,'AlphaData',0.3);
axis image;
axis on

point = [Nx/2 Ny/2 700];
normal = [0 -1 0];
[platform2D,x,y,z] = obliqueslice(piattaforma,point,normal);
figure
ax2 = axes('XLim',[0 Nx]*d,'YLim',[0 Ny]*d,'ZLim',[0 Nz]*d);
surf(x,y,z,platform2D,'EdgeColor','None','HandleVisibility','off','Parent',ax2);

[~,bw] = bwboundaries(platform2D(:,54:353));
% bw = imrotate(bw,180);

figure(1)
hold on
handler = visboundaries(ax,bw,'Color','black','LineWidth',0.1)
handler.Children(1).XData = handler.Children(1).XData*d - 30.5;
handler.Children(1).YData = handler.Children(1).YData*d;
handler.Children(2).XData = handler.Children(2).XData*d - 30.5;
handler.Children(2).YData = handler.Children(2).YData*d; 
% set(gca,'XLim',[-Nx/2 + 54 Nx/2 - (412 - 353)] * d,'YLim',[0 Nz] * d)


% %% Figura intera
% 
% h = figure
% im_setup = imagesc(x_beam*1e3, (z_beam-min(z_beam)) * 1e3, setup)
% ax = gca;
% xlabel('x-position [mm]');
% ylabel('z-position [mm]');
% title('Max Recorded Pressure');
% ax.YDir = 'normal';
% % ylabel(c, 'Pressure [Normalized]');
% hold on
% im_beam = imagesc(x_beam*1e3, (z_beam-min(z_beam)) * 1e3, data);
% set(im_beam,'AlphaData',0.7);
% set(h,'color','none','visible','off');
% colormap([[parula];[0 0 0]])
% c = colorbar;
% axis image;
% 
% %% Zoom
% 
% z_abs = (z_beam - min(z_beam));
% index = find(z_abs >= 0.12 & z_abs <= 0.16); 
% z_zoomed = z_abs(index)
% 
% % figure, 
% % im_beam = imagesc(x_beam * 1e3, z_zoomed * 1e3, data(600:760,:))
% % ax = gca;
% % xlabel('x-position [mm]');
% % ylabel('z-position [mm]');
% % title('Max Recorded Pressure');
% % colormap([[jet];[0 0 0]])
% % % c = colorbar;
% % ax.YDir = 'normal';
% % % ylabel(c, 'Pressure [Normalized]');
% % % hold on
% % % im_beam = imagesc(x_beam*1e3, z_zoomed * 1e3 , data(600:760,:));
% % % set(im_beam,'AlphaData',0.7);
% % axis image;