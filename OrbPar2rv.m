
%% FUNCTION: 'OrbPar2rv'

% INPUT: ORBITAL PARAMETERS
% 1) a = semimajor axis [km]
% 2) e = eccentricity
% 3) i = inclination [rad]
% 4) OMEGA = ascending node anomaly [rad]
% 5) omega = pericenter anomaly [rad]
% 6) theta = position in the perifecal reference frame [rad]
% 7) mu_planet = is the product between the universal gravity constant G and
%    the mass of planet m_planet  [km^3/s^2]


% OUTPUT: POSITIONS AND VELOCITY IN PLANETOCENTRIC EQUATORIAL REFERENCE FRAME
% 1) R_EQ = column vector containing the positions in an planetocentric equatorial reference frame [km]
% 2) V_EQ = column vector containing the velocities in an planetocentric equatorial reference frame [km/s] 


%%

function [R_EQ,V_EQ] = OrbPar2rv (a,e,i,OMEGA,omega,theta,mu_planet)


if nargin > 6
    mu = mu_planet;
else
    mu = 398000;         % mu_earth [km^3/s^2]
end
   

% compute r
p = a*(1-e^2);
r = p/(1+e*cos(theta));

% Position vector in the perifocal reference frame along x and y components
R_PF = [r*cos(theta);r*sin(theta);0];

% Radial and transverse velocities
v_rad = -sqrt(mu/p)*sin(theta);
v_tr = sqrt(mu/p)*(e+cos(theta));

% Velocity in the perifocal reference frame along x and y components
V_PF = [v_rad; v_tr; 0];


% Rotation matrix about Z axis [OMEGA] 
OMEGA_ROT = [cos(OMEGA) sin(OMEGA) 0; -sin(OMEGA) cos(OMEGA) 0; 0 0 1];

% Rotation matrix about X' axis [I]  
I = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)];

% Rotation matrix about Z' axis [omega] 
omega_ROT = [cos(omega) sin(omega) 0; -sin(omega) cos(omega) 0; 0 0 1];

% COMPLESSIVE ROTATION MATRIX
T = omega_ROT*I*OMEGA_ROT;


% From perifocal to planetocentric equatorial reference frame
% Position
R_EQ = T^-1*R_PF;

% Velocity
V_EQ=T^-1*V_PF;

