clc, clear all, close all

%% -- Definizione della griglia (che deve essere la stessa usata in simulazione)

PML_size = 10;

Nx = 432 - 2*PML_size;            % number of grid points in the x direction
Ny = 432 - 2*PML_size;            % number of grid points in the y direction
Nz = 800 - 2*PML_size;            % number of grid points in the z direction
d = 0.2e-3;
dx = d;
dy = d;
dz = d;
kgrid = kWaveGrid(Nx, dx, Ny, dy, Nz, dz);

%% Carico il modello, utile per plottare colon e trasduttore in 3D

load simGrid_500KHz.mat
grid = simGrid_500KHz;
clear simGrid_500KHz

%--------------------------------------------------------------------------
% Allineo le dimensioni del modello con quello della griglia(che deve essere
% una potenza di 2,3 o al più 5), escluso il PML. Le etichette della
% segmentazione dipende dall'ordine dato alle varie parti su 3D slicer
%
% Colon
colon = zeros(size(grid));
index = find(grid == 9);
colon(index) = 1;
colon = cat(1,colon,zeros(412 - size(colon,1), size(colon,2), size(colon,3)));
colon = cat(2,colon,zeros(size(colon,1), 412 - size(colon,2), size(colon,3)));
colon = cat(3,colon, zeros(size(colon,1), size(colon,2), 780 - size(colon,3)));
% Transducer 
transducer = zeros(size(grid));
index = find(grid == 101);
transducer(index) = 1;
transducer = cat(1,transducer,zeros(412 - size(transducer,1),...
    size(transducer,2), size(transducer,3)));
transducer = cat(2,transducer,zeros(size(transducer,1),...
    412 - size(transducer,2), size(transducer,3)));
transducer = cat(3,transducer, zeros(size(transducer,1),...
    size(transducer,2), 780 - size(transducer,3)));
clear index
%--------------------------------------------------------------------------

%% Carico i dati della simulazione 3D

%--------------------------------------------------------------------------
% Dopo averli caricati li riordino con 'reshape' e li normalizzo
filename = "C:\Users\Gabri\Desktop\Simulazioni 3D Risultati\500 KHz\Mouse_Sim_Res_Pos2.h5";
sensor_data.p_max = h5read(filename,"/p_max");
data = reshape(sensor_data.p_max, [Nx Ny Nz]);
data = data./max(max(max(data)));
%--------------------------------------------------------------------------

%% Figure 1
d = 0.2; 

%--------------------------------------------------------------------------
% Nella prima figura visualizzo il plot 3D. Con la funzione 'axes' imposto 
% il range da visualizzare nel plot. La funzione 'meshgrid' restituisce una
% griglia 3D in base ai vettori dati in ingressi. Questa griglia verrà
% utilizzata, una volta riportata al valore in mm (moltiplicando per 'd'),
% per plottare il colon ed il trasduttore con le due funzioni successive.
figure(1)
ax2 = axes('XLim',[-Nx/2  Nx/2-1]*d,'YLim',[-Ny/2  Ny/2-1]*d,'ZLim',[0 Nz-1]*d);
[X,Y,Z] = meshgrid(-Nx/2: Nx/2 -1, -Ny/2 : Ny/2 -1, 0 : Nz-1);
patch(isosurface(X*d, Y*d, Z*d, colon),'FaceColor', 'white', 'EdgeColor', 'black','Parent',ax2);
patch(isosurface(X*d, Y*d, Z*d, transducer),'FaceColor', 'white', 'EdgeColor', 'black','Parent',ax2);
hold on

% La funzione 'obliqueslice' prende dalla griglia 3D un piano passante per
% 'point' e con vettore uscente perpendicolare al piano definito in
% 'normal'. 
% 'axis equal' gestisce le proporzioni tra gli assi e 'view' gestisce la
% posizione della telecamera nello spazio
point = [0 Ny/2 0];  
normal = [0 1 0];
[B,x,y,z] = obliqueslice(data,point,normal);
surf((x-Nx/2)*d,(y-Ny/2)*d,(z+1)*d,B,'EdgeColor','None','HandleVisibility','off','Parent',ax2);
point = [Nx/2 Nx/2 643];  
normal = [0 0 1];
[B,x,y,z] = obliqueslice(data,point,normal);
surf((x-Nx/2)*d,(y-Ny/2)*d,z*d,B,'EdgeColor','None','HandleVisibility','off','Parent',ax2);
axis equal
view([-62.6 13.8])
plot3(point(1),point(2),point(3),'or','MarkerFaceColor','r');
plot3(normal(1),normal(2),normal(3),'ob','MarkerFaceColor','b');
xlabel('x - grid voxels')
ylabel('y - grid voxels')
zlabel('z - grid voxels')
%--------------------------------------------------------------------------

figure(2),
h = imshow(imrotate(B,180),'XData',[-Nx/2 Nx/2]*d,'YData',[-Ny/2 Ny/2]*d)
axis on

clear x y z
figure('Name','rectum 3D'),
ax3 = axes('XLim',[0 Nx],'YLim',[0 Ny],'ZLim',[0 Nz]);
colon = single(colon); % Come int16 non la prende. E' necessaria una conversione
[rect,x,y,z] = obliqueslice(colon,point,normal);
surf(x,y,z,rect,'EdgeColor','None','HandleVisibility','off','Parent',ax3);
axis equal
view([-62.6 13.8])
% max_data = max(max(max(data)));
% B(end,end) = max_data;
% ax3 = axes('XLim',[0 368],'YLim',[0 364]);
figure('Name','rectum 2D'),
imshow(rect)
colormap([[1 1 1];[0 0 0]])

% bw = imrotate(activecontour(rect, rect, 5, 'edge'),180);
[~,bw] = bwboundaries(rect);
bw = imrotate(bw,180);

%--------------------------------------------------------------------------
% Nuova colormap
color = [linspace(0,1,256)',linspace(0,1,256)',linspace(0,1,256)'];
%--------------------------------------------------------------------------

figure(2)
hold on
handler = visboundaries(bw,'Color','black','LineWidth',0.5)
handler.Children(1).XData = (handler.Children(1).XData - Nx/2)*d;
handler.Children(1).YData = (handler.Children(1).YData - Ny/2)*d;
handler.Children(2).XData = (handler.Children(2).XData - Nx/2)*d;
handler.Children(2).YData = (handler.Children(2).YData - Ny/2)*d;
colormap(jet)
xlabel('x-position [mm]');
ylabel('y-position [mm]');
c = colorbar
c.Location = 'southoutside';