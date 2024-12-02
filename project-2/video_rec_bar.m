%% 2.5 - 2d Video of Scattering off Rectangular Barrier

close all;
clear; clc;
format long;

% Simulation maximum time 
tmax = 0.05;
% Discretization level
level = 8;
% Delta t by Delta x ratio
lambda = 0.05;

% idtype = 0   ->  Exact family (sine wave)
% idtype = 1   ->  Boosted Gaussian
idtype = 1;
%x0      = idpar(1);      y0 = idpar(2);    
%delta_x = idpar(3); delta_y = idpar(4); 
%p_x     = idpar(5);     p_y = idpar(6);   
idpar = [0.5, 0.3, 0.1, 0.1, 0.0, 30];

% vtype = 0   ->  No potential
% vtype = 1   ->  Rectangular barrier or well
% vtype = 2   ->  Double Slit
vtype = 1;
x_min = 0.2;   x_max = 0.8;    
y_min = 0.7;   y_max = 0.8; 
%Vc    = vpar(5); 
vpar = [x_min, x_max, y_min, y_max, 1e8];

% Compute solution 
[x y t psi psire psiim psimod v] = ...
    sch_2d_adi(tmax, level, lambda, idtype, idpar, vtype, vpar);

% Dimensions of matrix 
[nt, nx, ny] = size(psimod);

% Create a VideoWriter object
video = VideoWriter('../../output/problem2/rec_bar.avi');
video.FrameRate = 30; 
open(video);

% Create a figure and define where the rectangle indicating where
% the potential barrier is will be drawn
figure;
rect_xmin = x_min * (nx - 1) + 1;
rect_ymin = y_min * (ny - 1) + 1;
rect_xlen = (x_max - x_min) * (nx - 1);
rect_ylen = (y_max - y_min) * (ny - 1);

% Loop over time steps
for n = 1:nt
    % reshape ψ to create a 2d matrix at this timestep
    psi_n = reshape(psimod(n,:,:), nx, ny);
    
    % Create filled contour plot
    contourf(psi_n, 20, 'LineStyle', 'none');
    colormap("default");
    xlabel('x');
    ylabel('y');
    title({'2d Schrödinger Equation Simulation'
           '|ψ| Scattering off a Rectangular Barrier' 
           ['tmax = ', num2str(tmax), ', level = ', num2str(level), ...
            ', lambda = ', num2str(lambda), ', idpar = [', ...
            num2str(idpar(1)), ' ', num2str(idpar(2)), ' ', ...
            num2str(idpar(3)), ' ', num2str(idpar(4)), ' ', ...
            num2str(idpar(5)), ' ', num2str(idpar(6)), ']']
          ['Time Step n = ', num2str(n)]});
    ax = gca;
    ax.FontSize = 12;
    
    % Set axis limits for consistency
    axis([1 nx 1 ny]);
    % Set color axis limits to match data range
    clim([min(psimod(:)) max(psimod(:))]);

    % Scale the axes tick labels to range from 0 to 1
    xticks(linspace(1, nx, 11));
    yticks(linspace(1, ny, 11));
    xticklabels(linspace(0, 1, 11));
    yticklabels(linspace(0, 1, 11));

    % Draw rectangle where the potential is 
    rectangle('Position', [rect_xmin, rect_ymin, rect_xlen, rect_ylen], ...
              'EdgeColor', 'black', 'LineWidth', 1);
    
    % Write to video file
    frame = getframe(gcf);
    writeVideo(video, frame);
end

% Close the video file and figure
close(video);