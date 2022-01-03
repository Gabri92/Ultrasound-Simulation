function [medium,setup] = Acoustic_prop(c0,rho0,Nx,Ny,Nz,d)

% clc, clear all, close all
% 
% c0 = 1480;
% rho0 = 1000;
% Nx = 351;
% Ny = Nx;
% Nz = 1000;
% d = 0.2e-3;
%% Import Models
[mouse,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\1 mm models\Mouse_Segmentation.nrrd');
[platform,~] = nrrdread('C:\Users\Gabri\Desktop\Novembre 2020\Lavoro x modello su topo\Risultati 3D slicer\1 mm models\Piattaforma.nrrd');

platform = single(platform);
mouse = single(mouse);

%% Map and creation of the environment
index = find(platform);
platform(index) = 10;
clear index
%------------------------------
% Water/Air = 0
% body = 1;
% backbone = 2;
% back muscles = 3;
% bladder = 4;
% cecum = 5;
% rectum = 6;
% colon = 7;
% platform = 10;
%-----------------------------

setup = mouse + platform;
% setup = cat(1,setup,zeros(364 - size(setup,1), size(setup,2), size(setup,3)));
% setup = cat(2,setup,zeros(size(setup,1), 364 - size(setup,2), size(setup,3)));
% setup = cat(3,setup, zeros(size(setup,1), size(setup,2), 1004 - size(setup,3)));
%
% index = find(setup);
% setup(index) = 1;
% clear index
%
% index = find(mouse);
% mouse(index) = 1;

% figure, volshow(setup)

%% Assignment of acoustical properties

medium.sound_speed = single(zeros(size(setup)));
medium.density = single(zeros(size(setup)));    

height = 135e-3;
height_grid = round(height/d);
% Water
medium.sound_speed(:,:,:) = 1480;	
medium.density(:,:,:) = 1000;
% Air
medium.sound_speed(:,:,height_grid:end) = 343;	
medium.density(:,:,height_grid:end) = 1.20;


% Polycarbonate
index = find(setup == 10);
medium.sound_speed(index) = 2270;
medium.density(index) = 1200;
clear index

% Body
index = find(setup == 1);
medium.sound_speed(index) = 1480;
medium.density(index) = 1000;
clear index

% BackBone
index = find(setup == 2);
medium.sound_speed(index) = 2800;
medium.density(index) = 1900;
clear index

% Back muscles
index = find(setup == 3);
medium.sound_speed(index) = 1590;
medium.density(index) = 1065;
clear index

% Bladder
index = find(setup == 4);
medium.sound_speed(index) = 1480;
medium.density(index) = 1000;
clear index

% Cecum
index = find(setup == 5);
medium.sound_speed(index) = 1480;
medium.density(index) = 1000;
clear index

% Rectum
index = find(setup == 6);
medium.sound_speed(index) = 1480;
medium.density(index) = 1000;
clear index

% Colon
index = find(setup == 7);
medium.sound_speed(index) = 1480;
medium.density(index) = 1000;
clear index

medium.sound_speed_ref = c0;

% figure, volshow(medium.sound_speed)
% figure, volshow(medium.density)
% voxelPlot(medium.sound_speed)
% voxelPlot(medium.density)
