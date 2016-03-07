function [SV] = MOQ_Skyview(X,Y,Z)
% Hillshade calculation following DHSVM (fraction of sky visible, for diffuse calculation).
%
%
%
%

SV = zeros(size(Z));
ds = mode(diff(X(1,:)));
dphi = 2*pi/36;
distance = 0:ds:5000;
phi = pi/2:-dphi:-3*pi/2+dphi;

for r = 2:size(Z,1)-1								% Row
	if mod(r,10) == 0
		disp(['Working on row #',num2str(r)])
	end
	for c = 2:size(Z,2)-1							% Column
		d = ((X-X(r,c)).^2 + (Y-Y(r,c)).^2).^(1./2);% Distance to each grid
		d(r,c) = 0;									% Force distance to current grid to be zero
		slope = (Z-Z(r,c))./d;						% Calculate the slope to each grid cell
		theta = atand(slope);						% Angle for each slope
		% For discreet bin find angle to the horizon
		for k = 1:36
			xlin = distance.*cos(phi(k))+X(r,c);	% Line along direction phi
			ylin = distance.*sin(phi(k))+Y(r,c);	% Line along direction phi
			theta_proj = interp2(X,Y,theta,xlin,ylin,'*nearest',NaN);
			SV(r,c) = SV(r,c) + dphi*sind(max(nanmax(theta_proj),0));
		end
		SV(r,c) = SV(r,c)/(2*pi);					% Normalize to hemisphere steridians
	end
end
