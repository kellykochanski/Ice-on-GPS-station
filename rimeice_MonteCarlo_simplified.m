function[ ] = rimeice_MonteCarlo(lat, long, ht, name, plot_on, path_out, navfile)
%% Rime plume monte-carlo
% Kelly Kochanski 2015, cleaned up 2018

% Simulates a variety of ice plumes on a hemispherical GPS station
% located at a position defined by:
%  lat      : latitude  (degrees)
%  long     : longitude (degrees)
%  ht       : elevation above sea level (m)
%  name     : name of site, string (used in filenames, keep short)
% Parameters controlling input and output:
%  plot_on  : binary. if True (1), plot delay in signal due to rime plume.
%  path_out : path for saving output data, string (set to 0 if not saving)
%              matlab must have write privilges in this folder.
%  navfile  : broadcast ephemeris file with GPS satellite tracks
%           : (defaults to 'brdc0760.12n')

%% DEFAULT ARGUMENTS
if nargin < 6;      navfile  = 'brdc0760.12n';  end
if nargin < 5;      disp('defaults'); plot_on  = 1;               end
if nargin < 4;      path_out = '';              end
if path_out == 0;   saving_on = 0;              
                    else saving_on = 1;         end

%% PARAMETERS FOR MONTE-CARLO

%% THE SHAPE OF AN ICE PLUME
% Function defining the average shape, or boundaries, of an ice plume,
% where the boundary is a distance from the center of the sphere at each
%  azimuthal angle   : a
%  elevation angle   : e
% This shape will be randomly perturbed in N monte carlo trials: 
N = 100;
% by spherical harmonics with size:
noiselevel = 0.2;

r_s     = 1; % radius of station, m
r_r     = 0; % initialize float
plumeboundaries = @(a, e)...
    ( ( - pi/2 <= a).*(pi/2 >= a) ).*...
    ( min([r_r, ...
           r_s./sin(e), ...
           r_s./cos(e)./abs(sin(a))])...
      - r_s );
  
%% Get satellite tracks
% Vector of satellite azimuths and elevations as seen from this site.
azzd = GPS_tracks(0, lat, long, ht, navfile);

%% Run a whole bunch of monte-carlo simulations
tic;

% SET UP OUTPUT
column_titles = 'lat, long, ht, r_r, rimedir, phi, dx, dy, dz, dclock, datm';
if saving_on;
    % If saving, write line-by-line to a csv file
    output_filename = strcat( path_out, 'rimeice-mc-output-', name);
    file = fopen(output_filename, 'w');
    strcat('Saving to : ', output_filename)
    dlmwrite(output_filename, column_titles, ...
        'delimiter', ',');
else
    % Else print output to display
    disp(column_titles)
end
    
% For a bunch of possible sizes (r_r), directions (rimedir) and widths
% (phi) of rime plumes, do N trials for each:

for r_r     =  r_s + exp(-4:0.5:5);
for rimedir = -pi + pi/50:pi/5:pi;
for phi     = -pi/4:pi/4:pi/4;
tmp_output = zeros(11, N);
parfor i    = 1:N; % parallizable with parfor
    output_vec      = run_one_trial(r_r, r_s, rimedir, phi,...
    noiselevel, plumeboundaries, azzd, plot_on);
    tmp_output(:,i) = [lat, long, ht, output_vec];
end % done with N trials (done with parallelized section)
for i = 1:N
    line = tmp_output(:,i);
    if saving_on; write_line_of_csv(output_filename, line);
    else          write_line_output(line); end
end
end % done with plume width phi
end % done with plume bearing rimedir
end % done with plume size r_r

toc;
end % end function rimeice_MonteCarlo

function [output_vector] = run_one_trial(r_r, r_s, rimedir, phi,...
    noiselevel, plumeboundaries, azzd, plot_on)
    % Create a bunch of random coefficients for spherical harmonics -
    %  these define an irregular, random rime plume
    coeffs = r_r*2*sqrt(pi).*[1, noiselevel*rand(1,15)-0.5*noiselevel];

    % Cut out variations with wavelengths less than that of L1        
    % Note: we'll be missing some spherical harmonic terms if r_r>36cm
    if      r_r <= 0.0605;  coeffs = coeffs(1);
    elseif  r_r <= 0.121;   coeffs = coeffs(1:4);
    elseif  r_r <= 0.242;   coeffs = coeffs(1:9);
    end

    % Determine the delay of a wave passing straight through the plume
    %  as a function of azimuth/elevation angle (direction)
    delayfunction = @(a, e) ...
        plumeboundaries(a, e).*1030.*sphericalHarmonic(a,e,coeffs);

    if plot_on;  
        disp('Plotting ice-caused delay as a function of satellite angle relative to GPS station')
        plot_delay_function(delayfunction, rimedir, r_s); 
    end

    %% Find shift in apparent position of station
    [ dx, dy, dz, dclock, datm ] = ...
        rime_effect(rimedir, azzd, delayfunction );

    % How do we average the observations from 30 satellites?
    % (Naive linear average here)
    Dx = mean(dx); Dy = mean(dy); Dz = mean(dz); 
    Dclock = mean(dclock); Datm = mean(datm);

    % Put independent vars + results into an array...
    output_vector = [
        r_r, rimedir, phi, ...
        Dx, Dy, Dz, Dclock, Datm];
    
    fprintf('Finished trial: rimedir=%d r_r=%d phi=%d \n', ...
        rimedir, r_r, phi)
end

function [ dx, dy, dz, dclock, datm ] = ...
    rime_effect(rimedir, azzd, delayfunction )
%rime_effect is the output of simulateRime.m at ONE point in time
%  rimedir  : direction of rime plume, radians
%  azzd     : produced by GPS_tracks.m - a set of satellite tracks.
%       for satellite s at position j:
%                 azzd(1,j,s) = azimuth angle,  radians, -pi < a < pi
%                                               0 is south
%                                               pi/2 is west
%                 azzd(2,j,s) = zenith angle,   radians,   0 < e < pi/2
%                                               0 is the apex (zenith)
%                                               pi/2 is the horizon=
%  delayfunction : a function handle returning ice delay as a function of 
%                 azimuth and zenith angles.

% Get number of satellites and number of time steps from azzd array
azzdsize = size(azzd);
n_sat    = azzdsize(3);
nt       = azzdsize(2);

% filter out any measurements which contain NaN values, or elevations
% greater than pi/2.
%azzd(isnan(azzd)) = pi/2;

% Arrays to store apparent displacements for each satellite
dx      = zeros(1,n_sat); 
dy      = zeros(1,n_sat); 
dz      = zeros(1,n_sat);
dclock  = zeros(1,n_sat); 
datm    = zeros(1,n_sat);

for sat = 1:n_sat % For each satellite track
    
    %  Don't allow elevation angles within 10 degrees of horizon
    %  (this operation also removes any NaN values)
    
    good_idx  = find(azzd(2,:,sat) < 1.396);
    
    zeniths    = azzd(2, good_idx, sat);
    azimuths_s = azzd(1, good_idx, sat);     % bearings from south
    azimuths_n = mod(azimuths_s + pi, 2*pi); %          from north
    %disp([max(max(azimuths_s)), min(min(azimuths_s))]);
    
    % Ice delay at each satellite position
    delays  = delayfunction(azimuths_n, zeniths);
    
    % Least-squares estimate of position (dx, dy, dz), clock error (dclock)
    % and atmospheric delay parameter (datm) due to ice
    A = [-sin(zeniths); ...
        cos(zeniths).*cos(azimuths_s); ...
        cos(zeniths).*sin(azimuths_s); ...
        ones(1,length(good_idx)); ...
        1./(sin(zeniths) + 10^-9)];
    
    % Least-squares solution to A*displacements = delays
    [L,U]          = lu(A'); 
    displacements  = U \ (L\delays');
    
    % Results of least squares solution
    dz(sat)     = displacements(1); 
    dy(sat)     = displacements(2); 
    dx(sat)     = displacements(3); 
    dclock(sat) = displacements(4); 
    datm(sat)   = displacements(5);
    
end % end for
end % end function

function [] = write_line_of_csv(output_filename, output_array)
% Function for saving output as one line in a csv file.
% Defined separately so matlab can save within parfor loop.
        dlmwrite(output_filename, output_array, 'delimiter', ',', '-append');
end

function [] = write_line_output(output_array)
% Writes one line to display
        disp(output_array);
end

function[ ] = plot_delay_function(delayfunction, rimedir, r_s)
% Plot delay function in spherical coordinates (azimuth a, elevation e)

avec = -pi:pi/100:pi;
evec = 0:pi/100:pi/2;

colormap('hot')
% Get delay at all angles around the hemisphere
delayed = zeros(length(avec), length(evec));
for k = 1:length(avec); 
    for j = 1:length(evec);
        atmp = mod(avec(k) - rimedir + pi, 2*pi) - pi;
        delayed(k,j) = delayfunction(atmp, evec(j));
    end; 
end;

[A, E]      = meshgrid(avec, evec);
[X, Y]      = azimuth_zenith_2_x_y(A, E, r_s);
 contourf(X, Y, delayed', 300, 'edgecolor', 'none')
 axis('square')
 title(fprintf('delay in picoseconds'));
 c = colorbar; 
 c.Label.String = 'Delay (ns)';
% pause();
end % End plotting function

function [X, Y] = azimuth_zenith_2_x_y(A, Z, r)
% Convert arrays of azimuth and zenith coodinates to x-y coordinates
%  within a circle of radius r. azimuth=0 on y axis and progresses cw
 distance_from_center   = r.*(1 - cos(Z)); % = r*( 1 - sin(pi/2 - Z))
 X          = distance_from_center.*sin(A);
 Y          = distance_from_center.*cos(A);
end
