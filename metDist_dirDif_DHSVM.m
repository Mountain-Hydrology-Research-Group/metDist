function [SDR_dir,SDR_dif] = metDist_dirDif_DHSVM(EL,SWdwn)
% This function will create a time series of direct-diffuse partitioning using
% DHSVM's algorithm 
%
% SYNTAX:
%	[SDR_dir,SDR_dif] = DirDifPartition_DHSVM(EL,SWdwn)
%
% INPUTS:
%	EL	= Lx1 vector of solar elevation angles (angle above horizon)
%
% OUTPUTS
%	Direct	= Lx1 vector of clear-sky direct component 
%	Diffuse	= Lx1 vector of clear-sky diffuse component
%


%%%%%%%%%%%%
%% Checks %%
%%%%%%%%%%%%
if numel(EL) ~= numel(SWdwn)
	error('The number of elements in each entry must be the same length')
end

%%%%%%%%%%%%%%%%%
%% "Constants" %%
%%%%%%%%%%%%%%%%%

S0 = sind(EL).*1367;						% Solar Constant (W/m^2) 
CI = SWdwn./S0;								% "Clearness index"
SDR_dif = NaN(size(SWdwn));
SDR_dir = NaN(size(SWdwn));

Fdif = .943 + .734.*CI - 4.9.*CI.^2 + 1.796.*CI.^3 + 2.058.*CI.^4;
Fdif(CI > .8) = .13;
Fdif(Fdif > 1) = 1;
Fdif(Fdif < .13) = .13;
Fdir = 1 - Fdif;

SDR_dir = Fdir.*SWdwn;
SDR_dif = Fdif.*SWdwn;
