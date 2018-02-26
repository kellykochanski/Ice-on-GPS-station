%% Rime plume monte-carlo

% Kelly Kochanski

% Simulates a variety of ice plumes on a hemispherical GPS station
% located at a position defined by:
%  lat:  latitude  (degrees)
%  long: longitude (degrees)
%  elev: elevation above sea level (m)

function[ ] = rimeice_MonteCarlo(lat, long, ht, name)

%% PARAMETERS FOR THIS RUN
% Does this latitude/longitude/elev belong to a named location?
if nargin < 4
    navfile = 'brdc0760.12n';
end

output_filename = 'C:/Users/Kelly/Documents/Documents/GPS_project/monte-carlo-output-20170206-AB25.csv'

% Returns a vector of satellite azimuths and elevations, 
%  as seen from this site. Takes in lat and long in radians.
azzd = HW01_2014_simplified(lat./180.*pi, long./180.*pi, ht)

%% Input parameters
etaStd = 1/3; %standard deviation of distribution of feathering length eta
r_s = 0.2;

%% Parameters for this test - varied independently
gamma = 3/2;
r_r   = r_s + 0.1;
lat = 62.92931;	long = -156.02339; ht = 961; % DO NOT SET HERE.
                                            % Set in HW01_2014.m
%% The shape of this ice plume
% Defaults so it'll compile
rimedir = 0; phi = 0;
% Function
plumeboundaries = @(a, e)...
    ( ( - pi/2 <= a).*(pi/2 >= a) ).*...
    ( min([r_r, ...
           r_s./sin(e), ...
           r_s./cos(e)./abs(sin(a))])...
      - r_s );

%% Run a whole bunch of monte-carlo simulations :)
N = 50;
rimesizes = r_s + exp(-4:0.5:5);
rimedirstep = pi/50; %radians
%phisteps = -pi/4:pi/4:pi/4;

tic;
dlmwrite(output_filename, zeros(1,12), 'delimiter', ','); 
% For a bunch of plausible sizes of rime plume
for r_r = rimesizes;
    % For each possible direction of ice plume
    for rimedir = -pi + rimedirstep :rimedirstep:pi;
        % For each reasonable direction of phi
    %    for phi = phisteps;
            % Do a bunch of trials
            output = zeros(0, 12);
            parfor i = 1:N;
                noiselevel = 0.2;
                xi = r_r*2*sqrt(pi).*[1, noiselevel*rand(1,15)-0.5*noiselevel];
                % cut out variations with wavelengths less than that of L1
                if r_r <= 0.0605; xi = xi(1);
                elseif r_r <= 0.121; xi = xi(1:4);
                elseif r_r <= 0.242; xi = xi(1:9);
                end
                % we'll be missing some spherical harmonic terms if r_r >
                % 36cm
                
                delayfunction = @(a, e) ...
                    plumeboundaries(a, e).*1030.*sphericalHarmonic(a,e,xi);
                
                %% USE THIS BLOCK TO PLOT DELAYFUNCTION(a,e)
%                 avec = -pi:pi/100:pi;
%                 evec = 0:pi/100:pi/2;
%                 hold all
%                 cmap = colormap;
%                 cmap(1, :) = [1 1 1];
%                 for k = 1:length(avec); 
%                     for j = 1:length(evec);
%                         atmp = mod(avec(k) - rimedir + pi, 2*pi) - pi;
%                         delayed(k,j) = delayfunction(atmp, evec(j));
%                         %colornum = max(1, min(64, round(3200*num)));
%                         %plot(avec(k)/pi, evec(j)/pi, '.', 'markeredgecolor', cmap(colornum, :))
%                         %disp('one done')
%                     end; 
%                 end;
%                 [A, E] = meshgrid(avec, evec);
%                 X = r_s.*cos(E).*cos(A); Y = r_s.*cos(E).*sin(A);
%                 Z = sqrt(max(0, r_s.^2-X.^2 - Y.^2));
%                 colormap('hot')
%                 caxis([0 15])
%                 title('delay in picoseconds')
%                 surf(X, Y, Z, delayed', ...
%                     'edgecolor', 'none',...
%                     'FaceColor','interp')
%                 %camlight right
%                 daspect([1 1 1]); axis off;
%                 c = colorbar; c.Label.String = 'Delay (ns)';
%                 pause();
                %% End of plotting block
                %% Solve for position
                             
                [ dx, dy, dz, dclock, datm ] = ...
                    rime_effect_improved(rimedir, azzd, delayfunction );

                % How do we average the observations from 30 satellites?
                Dx = mean(dx); Dy = mean(dy); Dz = mean(dz); 
                Dclock = mean(dclock); Datm = mean(datm);

                % Put independent vars + results into an array...
                outvec = [lat, long, ht, ...
                    gamma, r_r, rimedir, phi, ...
                    Dx, Dy, Dz, Dclock, Datm];
                output = [output; outvec];
            end
            dlmwrite(output_filename, output, 'delimiter', ',', '-append');
      %  end
    end
end

toc;
