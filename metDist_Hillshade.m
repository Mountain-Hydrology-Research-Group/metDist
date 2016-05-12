function [HS,SHADE] = metDist_Hillshade(S_ZA,S_AZ,X,Y,Z)
% Hillshade and topographic calculation following DHSVM.
%
% SYNTAX:
%	[HS,SHADE] = MOQ_Hillshade(S_ZA,S_AZ,X,Y,Z)
%
% INPUTS:
%	S_ZA	= Solar Zenith Angle (from AVG_EL.m)
%	S_AZ	= Solar Azimuth Angle (from solargeometry function)
%	X		= DEM x coordinate (UTM?)
%	Y 		= DEM y coordinate (UTM?)
%	Z		= DEM elevation (m), note: z must by a length(X) by length(Y) matrix
%
% OUTPUTS:
%	
%

HS = NaN(size(Z));
SHADE = NaN(size(Z));
ds = mode(diff(X(1,:)));
distance = 0:ds:5000;
% Phi in radians and degrees, taking into account Matlab's "North"
dphiR = 2*pi/72;
dphiD = 360/72;
phiR = pi/2-dphiR:-dphiR:-3*pi/2;
phiD = [dphiR:dphiR:2*pi].*180./pi;
phi_sun = phiR(round(S_AZ/dphiD));

% Convert from geographic azimuth->mathematic unit->radians
S_AZRad = 360-S_AZ+90;
if S_AZRad >= 360
	S_AZRad = S_AZRad - 360;
end
S_AZRad = S_AZRad*pi/180;
S_ZARad = S_ZA*pi/180;

for r = 2:size(Z,1)-1									% Row
	if mod(r,10) == 0
		disp(['Working on row #',num2str(r)])
	end
	for c = 2:size(Z,2)-1								% Column
		%% Topographic shading
		d = ((X-X(r,c)).^2 + (Y-Y(r,c)).^2).^(1./2);% Distance to each grid
		d(r,c) = 0;									% Force distance to current grid to be zero
		slope = (Z-Z(r,c))./d;						% Calculate the slope to each grid cell
		theta = atand(slope);						% Angle for each slope		
		xlin = distance.*cos(phi_sun)+X(r,c);		% Line along direction phi
		ylin = distance.*sin(phi_sun)+Y(r,c);		% Line along direction phi
		theta_proj = interp2(X,Y,theta,xlin,ylin,'nearest',NaN);
		SHADE(r,c) = nanmax(theta_proj) < 90-S_ZA;	% 1 if exposed, 0 if shaded
			
		%% Hillshade
		% Approximate local slope and azimuth using moving 3x3 window
		slopex = ( (Z(r-1,c+1)+2*Z(r,c+1)+Z(r+1,c+1) ) - (Z(r-1,c-1)+2*Z(r,c-1)+Z(r+1,c-1)) )/(8*ds);
		slopey = ( (Z(r+1,c-1)+2*Z(r+1,c)+Z(r+1,c+1) ) - (Z(r-1,c-1)+2*Z(r-1,c)+Z(r-1,c+1)) )/(8*ds);
		slope_apprx = atan((slopex^2 + slopey^2)^(1/2));
		if slopex ~= 0
			aspect = atan2(slopey,-slopex);
			if aspect < 0
				aspect = aspect + 2 * pi;
			end
			aspect = aspect;
		elseif slopex == 0
			if slopey > 0
				aspect = pi/2;
			elseif slopey <= 0
				aspect = 2*pi-pi/2;
			end
		end
		HS(r,c) = ((cos(S_ZARad)*cos(slope_apprx)) + sin(S_ZARad)*sin(slope_apprx)*cos(S_AZRad-aspect));
	end
end
