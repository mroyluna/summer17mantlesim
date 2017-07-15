clf
% Node points
x = 0:100;
y = 0:100;

% Grid of all points
[X,Y] = meshgrid(x, y);

% Coefficients for Velocity field
K = 10;
C = 10;

% Velocity vector calcs (wavy example flow)
Vx = K + C * sin(x);
Vy = (K/2) + (C/2) * cos(y);

% Grid of velocity vectors
[VX,VY] = meshgrid(Vx, Vy);

% Plots the velocity streamslice
streamslice(X, Y, VX, VY);

% Starting points for the particles
Sx = 0;
Sy = 0:20:100;

% Grid of starting points
[SX,SY] = meshgrid(Sx, Sy);

% stream2 computes 2D streamline data
vertices = stream2(X, Y, VX ,VY, SX, SY);

% Interpolates stream line velocities
iverts = interpstreamspeed(X, Y, VX, VY, vertices, 1);
sl = streamline(iverts);
set(sl, 'Marker','o')

% % Adds animated stream particles
% streamparticles(iverts,1,'Animate',100,...
%     'ParticleAlignment','on',...
%     'MarkerFaceColor','red',...
%     'Marker','o');