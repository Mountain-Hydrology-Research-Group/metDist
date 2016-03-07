%% PURPOSE
% Find shading for applicaitons in the MntObsQC project. This file is where
% much of the testing of the shading algorithm was originally performed. NOTE:
% The logical argument for shading is an AND logical and not an OR logical.
% This is wrong and will need to be fixed in any future use of this script.
clear all
close all
clc

%% Flags
% Yosemite
YOS_DEM_FLAG = 0;								% Process DEM data (Dana Meadows)
YOS_Hori_FLAG = 0;								% Find the horizon for Yosemite CDWR sites
YOS_TransInterp_Shade = 0;						% Try to interpolate to shaded points
YOS_DistributeSW_FLAG = 1;						% Distribute shortwave to the surrounding terrain at Dana Meadows
	Skyview_FLAG = 0;
	Hillshade_FLAG = 0;
	Fill_Flag = 1;
% Washington
SNQ_DEM_FLAG = 0;								% Process DEM data for SNQ
SNQ_Hori_FLAG = 0;								% Find the horizon @ SNQ

% SBSA (Senator Beck Study Area)
SBSA_DEM_FLAG = 0;								% Process DEM data for SBSA sites
SBSA_Hori_FLAG = 0;								% Find the horizon @ SBSA sites

% General Flags
PRINT_FLAG = 0;									% Print figures

%% Directories
% General directories
HOMEdir = 'C:\Users\Karl\gdrive\';
PRINTdir = 'SnowHydrology\proj\MntObsQC\Graphics\Shading';
PROJdir = 'SnowHydrology\proj\MntObsQC\Shading';
DATdir = 'SnowHydrology\proj\MntObsQC\DATA';	% Fully QC'ed data here

% Yosemite directories
YOSDATdir = 'GroundObs\Yosemite.CDWR\';
YOSDEMdir = 'C:\DATA\Spatial.Maps\DEM\Yosemite';

% SNQ directories
SNQDEMdir = 'C:\DATA\Spatial.Maps\DEM\Washington';
SNQDATdir = 'GroundObs\Snoqualmie\';

% SBSA directories
SBSADEMdir = 'C:\DATA\Spatial.Maps\DEM\SBSA';
SBSADATdir = 'GroundObs\SenatorBeck\';


%% Load DEM (Yosemite)
if YOS_DEM_FLAG
	% Read raster (converted by ArcGIS from NED ArcGIS compatible format to
	% ascii raster). Original data in back up of Mark's folder.
	cd(YOSDEMdir)
	Z = dlmread('Yosemite.DEM30m.txt');
	Z = flipud(Z);
	
	% ASCII raster header
	% ncols         6000
	% nrows         3720
	% xllcorner     -120.66666666645
	% yllcorner     37.499999999306
	% cellsize      0.00027777777779647
	% NODATA_value  -9999
	ncols = 6000;
	nrows = 3720;
	xll = -120.66666666645;
	yll = 37.499999999306;
	ds = 0.00027777777779647;
	
	% Create X and Y vectors of coordinates @ center of pixel
	X = linspace(xll+ds/2,xll+ds*(ncols)+ds/2,ncols);
	Y = linspace(yll+ds/2,yll+ds*(nrows)+ds/2,nrows);

	% Convert to UTM (span multiple zones, takes some adjustments, hopefully what I've done isnt' stupid)
	[~,y,zone] = deg2utm(Y,ones(size(Y))*X(1));
	[x,~,zone] = deg2utm(ones(size(X))*Y(1),X);
	% Stitch together zones as continuously as possible
	xjump = x(2401)-24.5608;
	xpre = x(2400);
	x(2401:end) = x(2401:end)-xjump+xpre;
	[xdan ydan] = deg2utm(37.897,-119.257);
	xdan = xdan-xjump+xpre;
	
	cd([HOMEdir,PROJdir])
	save('YOS_DEM_30m.mat','X','Y','Z','x','y')
end

%% Find horizon
if YOS_Hori_FLAG
		
	% Point observations
	cd([HOMEdir,YOSDATdir])
	s = [{'DAN'},{'TUM'},{'TES'}];							% Yosemite sites
	load DAN.QC.mat
	load TES.QC.mat
	load TUM.QC.mat
	MET.DAN = DAN;
	MET.TES = TES;
	MET.TUM = TUM;
	
	% DEM
	cd([HOMEdir,PROJdir])
	load YOS_DEM_30m.mat
	
	for n = 1:3
		% Find grid cell containing site in lat long
		[xm ym] = meshgrid(X,Y);							% Mesh grid (degrees)
		d = ( (xm-MET.(s{n}).long).^2 + (ym-MET.(s{n}).lat).^2).^(1./2);
		[r,c] = find(min(min(d)) == d);
		[xm ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)
		xsite = xm(r,c);
		ysite = ym(r,c);

		% Briefly test that I found the right point
		figure(1)
		imagesc(Z)
		hold all
		plot(c,r,'k*')
		set(gca,'YDir','normal')

		% Calculate the slope to each grid cell
		d = ( (xm-xm(r,c)).^2 + (ym-ym(r,c)).^2).^(1./2);
		d(r,c) = 0;
		slope = (Z-Z(r,c))./d;
		% Angle for each slope
		theta = atand(slope);

		% For some discreet phi bin find the angle to the horizon
		dphi = pi/36;
		ds = 25;
		distance = 0:ds:10000;
		phi = pi/2:-dphi:-3*pi/2+dphi;
		h.(s{n}) = NaN(size(phi));
		for k = 1:length(phi)
			xlin = distance.*cos(phi(k))+xsite;
			ylin = distance.*sin(phi(k))+ysite;
			theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);
			h.(s{n})(k) = nanmax(theta_proj);		
		end
		h.(s{n})(h.(s{n}) < 0) = 0;
		h_AZ = [0:dphi:2*pi-dphi].*180./pi;

		% Hinkelman plots w/ horizon
		Hink_hand = SW_Obs_QC_PlottingTool_SolGeo(MET.(s{n}),MET.(s{n}).SWFLAG);
		figure(Hink_hand(1))
		hold all
		plot(h_AZ,90-h.(s{n}))
		figure(Hink_hand(2))
		hold all
		plot(h_AZ,90-h.(s{n}))
		
		% Print both figures
		if PRINT_FLAG
			cd([HOMEdir,PRINTdir])
			figure(Hink_hand(1))
			print('-depsc',['Trans_Horiz.',s{n},'.eps'])
			figure(Hink_hand(2))
			print('-depsc',['Trans_STD_Horiz.',s{n},'.eps'])
		end
		close all
	end
	clear Z X Y x y xm ym
end

%% YOS Transmissivity interpolation
if YOS_TransInterp_Shade
	
	% Statement of problem:
	% 1) Interpolating from times that are shaded to times that are shaded
	% isn't useful. However, we are clearly missing instances of
	% shading due to vegetation.
	% 2) Interpolating from times that don't pass QC to times that do is
	% equally troublesome. However, taking this into account makes the
	% interpolation much more difficult.
	
	%% Index for shaded and unshaded periods - DEM horizon alone
%	SWgood = MET.DAN.SWdwn;
%	ind_exp = [];
%	temp_h_AZ = 0:.05:360;
%	for k = 1:length(temp_h_AZ)-1
%			temp_h_EL = interp1(h_AZ,h.DAN,temp_h_AZ);
%			temp = find(MET.DAN.AZ >= temp_h_AZ(k) ...
%				& MET.DAN.AZ <= temp_h_AZ(k+1) ...
%				& MET.DAN.EL > temp_h_EL(k)+5);
%			ind_exp = [ind_exp; temp];
%	end
%	temp_flag = ~ismember(1:length(MET.DAN.AZ),ind_exp)';
%	temp_flag = temp_flag | MET.DAN.SWFLAG ~= 0;

	%% load point observations
	cd([HOMEdir,YOSDATdir])
	s = [{'DAN'},{'TUM'},{'TES'}];							% Yosemite sites
	load DAN.QC.mat
	load TES.QC.mat
	load TUM.QC.mat
	MET.DAN = DAN;
	MET.TES = TES;
	MET.TUM = TUM;

	MET.DAN.RH = MET.DAN.RH./100;
	MET.DAN.RH(MET.DAN.RH > 1) = 1;
	
	%% Empirically defining shading
	[hand,Tr_bar,Tr_std] = SW_Obs_QC_PlottingTool_SolGeo(MET.DAN,MET.DAN.SWFLAG);
	dEL = .5;											% Elevation bin width
	ELdiscrete = [0:dEL:90];							% Discrete elevation bin edges
	ELplot = [dEL/2:dEL:90-dEL/2];						% Middle of each EL bin
	dAZ = .5;											% Azimuth bin width
	AZdiscrete = [0:dAZ:360];							% Discrete azimuth bin edges
	AZplot = [dAZ/2:dAZ:360-dAZ/2];						% Middle of each AZ bin
	[ind_EL,ind_AZ] = find(Tr_bar > .3 & Tr_std > .1);
	ind = find(Tr_bar > .3 & Tr_std > .1);

	%% Figures
	if PRINT_FLAG
		% Mean Transmittance
		figure
		scatter(AZplot(ind_AZ),90-ELplot(ind_EL),20,Tr_bar(ind),'filled')
		% Formating
		CLim = [0 1];
		set(gca,'CLim',CLim)
		COL = colorbar;
		cmap = cbrewer('seq','Purples',11);
		cmap = cmap(2:end,:);
		colormap(cmap)
		set(gca,'YDir','reverse')
		grid on
		xlabel('Azimuth')
		ylabel('SZA')
		ylabel(COL,'Transmissivity')
		axis([0 355 5 90])
		SetFigureProperties([1.5,1],gcf)

		% Standard deviation of transmittance
		figure
		scatter(AZplot(ind_AZ),90-ELplot(ind_EL),20,Tr_std(ind),'filled')
		% Formating
		CLim = [0 .3];
		set(gca,'CLim',CLim)
		COL = colorbar;
		cmap = cbrewer('seq','Purples',11);
		cmap = cmap(2:end,:);
		colormap(cmap)
		set(gca,'YDir','reverse')
		grid on
		xlabel('Azimuth')
		ylabel('SZA')
		ylabel(COL,'Transmissivity')
		axis([0 355 5 90])
		SetFigureProperties([1.5,1],gcf)
	end

	%% 2D space search
	% Find points in time series within discrete solar geometry grid.
	ind_exp = [];										% Pre-allocate

	% Limit the 2D search space to only the area
	dEL = .5;											% Elevation bin width
	dAZ = .5;											% Azimuth bin width
	ELmax = round(max(MET.DAN.EL)./dEL).*dEL;							 
	ELmin = 0;
	AZmax = round(max(MET.DAN.AZ)./dAZ).*dAZ;
	AZmin = round(min(MET.DAN.AZ)./dAZ).*dAZ;
	
	% Discretized solar geometry space
	for k = 1:length(ind_EL)
		temp = find(MET.DAN.AZ > AZplot(ind_AZ(k)) - dAZ/2 ...
			& MET.DAN.AZ <= AZplot(ind_AZ(k)) + dAZ/2 ...
			& MET.DAN.EL > ELplot(ind_EL(k)) - dEL/2 ... 
			& MET.DAN.EL <= ELplot(ind_EL(k))+ dEL/2);
		ind_exp = [ind_exp; temp];
	end
		
 	%% Interpolate transmissivity
 	% Times not classified as exposed (not a member of exposed index) 
	ind_interp2 = ~ismember(1:length(MET.DAN.AZ),ind_exp)';
	% Times that do not pass generic QC are also excluded
	ind_interp2 = ind_interp2 | (MET.DAN.SWFLAG ~= 0 & MET.DAN.SWFLAG ~= 3);
	SWgood = MET.DAN.SWdwn;
	SWgood(ind_interp2) = NaN;
	RH = MET.DAN.RH./100;
	RH(RH > 1) = 1;
	T = MET.DAN.T;
	EL = MET.DAN.EL;
	SOLDIST = MET.DAN.SOLDIST;
	[SWproc,PFLAG] = Trans_interp(SWgood,EL,SOLDIST,RH,T,MET.DAN.t);
	cd([HOMEdir,PROJdir])
	exc_ind = ind_interp2;
	save('DAN_SWproc_TransInterp.v1.mat','SWproc','PFLAG','exc_ind')
end

% Distribute shortwave to the surrounding terrain at Dana Meadows
if YOS_DistributeSW_FLAG
	%% load data
	% Point Observations
	cd([HOMEdir,DATdir])
	load('MntObs_v4.mat')
	DAN = MET.DAN;
	clear MET
	
	%% Transmissivity interpolation using new QC version of DAN
	FILL_FLAG = 1;
	if FILL_FLAG
		% Pre-allocate qc flags matrix
		DAN.SWFLAG = zeros(length(DAN.SWQCFLAG),8);
		% Column designations:
		% 1) Minimum global SWdwn (i.e., IR leaking)
		% 2) Maximum night time
		% 3) Climate dependent maximum - ignoring
		% 4) "Physical" maximum
		% 5) Unphysically low values (below Rayleigh clear-sky scattering limit)
		% 6) Snow on dome
		% 7) Shading (requires 3 years of data)
		% 8) Failed any test
		
		% Assign each QC test to the above columns
		DAN.SWFLAG(:,1:2) = DAN.SWQCFLAG(:,1:2);
		DAN.SWFLAG(:,4:7) = [DAN.SWQCFLAG(:,4:6),DAN.SWQCFLAG(:,8)];
		% General fail flag: 1 = failed, 0 = passed
		DAN.SWFLAG(:,8) = logical(sum(DAN.SWFLAG(:,1:7),2));
		
		% Interpolate to these indices (fail tests)
		ind_interp_to = DAN.SWFLAG(:,8) ~=0;
		SWgood = DAN.SWdwn;
		SWgood(ind_interp_to) = NaN;
		RH = DAN.RH./100;
		RH(RH > 1) = 1;
		T = DAN.T;
		EL = DAN.EL;
		SOLDIST = DAN.SOLDIST;
		[SWproc,PFLAG] = Trans_interp(SWgood,EL,SOLDIST,RH,T,DAN.t);
		cd([HOMEdir,PROJdir])
		save('DAN_SWproc_TransInterp.v2.mat','SWproc')		
	end
	
	%% Distributed
	EL = DAN.EL;
	AZ = DAN.AZ;
	SWdwn = DAN.SWdwn;
	% Processed SW @ DAN
	cd([HOMEdir,PROJdir])
	load DAN_SWproc_TransInterp.v2.mat
	
	% DEM
	load YOS_DEM_30m.mat
	[xm ym] = meshgrid(x,y);					% Mesh grid (UTM - meters)
	% Delineate terrain surrounding Dana
	LL = [37.667,-119.433];
	UR = [38,-119.2];
	[Xm,Ym] = meshgrid(X,Y);
	indx = (X >= LL(2) & X <= UR(2));
	indy = (Y >= LL (1) & Y <= UR(1));
% 	ind = Xm >= LL(2) & Xm <= UR(2) & Ym >= LL (1) & Ym <= UR(1);
	xm = xm(indy,indx);
	ym = ym(indy,indx);
	% Degredate to 150m
	xdeg = xm(1,1):150:xm(1,end);
	ydeg = ym(1,1):150:ym(end,1);
	[xmdeg,ymdeg] = meshgrid(xdeg,ydeg);
	Zdeg = interp2(xm,ym,Z(indy,indx),xmdeg,ymdeg);
% 	Z = Z(indy,indx);
	
	%% Hillshade
	if Skyview_FLAG
		SV = MOQ_Skyview(xmdeg,ymdeg,Zdeg);
		cd([HOMEdir,PROJdir])
		save('Skyview_150m.Dana.mat','SV')
	end
	
	%% Hillshade
	if Hillshade_FLAG
		% Hillshade on March 15th
		for k = 5:20
			d1 = datenum(2005,3,15,k-1,59,0);
			d2 = datenum(2005,3,15,k,1,0);
			dind = FindDate(DAN.t,d1,d2);
			SZAlocal = 90-EL(dind);
			AZlocal = AZ(dind);
			[HS,SHADE] = MOQ_Hillshade(SZAlocal,AZlocal,xmdeg,ymdeg,Zdeg);
			disp('')
			oname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat']
			save(oname,'HS','SHADE')
		end
	end
	
	%% Distribute
	% Non-QC'ed obs
	[Sd,Sf] = DirDifPartition_DHSVM(EL,SWdwn);
	
	% QCed and processed obs
	[Sd_proc,Sf_proc] = DirDifPartition_DHSVM(EL,SWproc);
	
	% Skyview fraction (found earlier)
	cd([HOMEdir,PROJdir])
	load Skyview_150m.Dana.mat
	
	% Cumulative difference in W/m-2
	xdim = size(SV,1);
	ydim = size(SV,2);
	SWdist_f = NaN(xdim,ydim,12);
	SWdist_d = NaN(xdim,ydim,12);
	SWdist_f_proc = NaN(xdim,ydim,12);
	SWdist_d_proc = NaN(xdim,ydim,12);
	for k = 7:18
			d1 = datenum(2005,3,15,k-1,59,0);
			d2 = datenum(2005,3,15,k,1,0);
			dind = FindDate(DAN.t,d1,d2);
			
			cd([HOMEdir,PROJdir])
			inname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat'];
			load(inname)
			
			% Raw SW
			SWdist_f(:,:,k-6) = SV.*Sf(dind);
			SWdist_d(:,:,k-6) = HS.*SHADE.*Sd(dind);
			% Proc SW
			SWdist_f_proc(:,:,k-6) = SV.*Sf_proc(dind);
			SWdist_d_proc(:,:,k-6) = HS.*SHADE.*Sd_proc(dind);
	end
	% Sum over whole day - divide by 24
	SWdist_full = sum(SWdist_f+SWdist_d,3);
	SWdist_full_proc = sum(SWdist_f_proc+SWdist_d_proc,3);
	
	%% Figure of cumulative difference (daily average - Wm-2)
	if PRINT_FLAG
		CLim_SWdwn = [0 round(max(max(SWdist_full))./24)];
		CLim_SWdiff = [-30 30];
		
		figure
		subaxis(2,2,1)

		imagesc(SWdist_full./24)
		set(gca,'YDir','normal')
		set(gca,'CLim',CLim_SWdwn)
		cmap = cbrewer('seq','YlOrRd',50);
		cmap = cmap(5:end,:);
		colormap(cmap)
		cbfreeze(cmap)
		freezeColors

		subaxis(2,2,2)
		imagesc(SWdist_full_proc./24)
		set(gca,'YDir','normal')
		set(gca,'CLim',CLim_SWdwn)
		COL = colorbar;
		cmap = cbrewer('seq','YlOrRd',50);
		cmap = cmap(5:end,:);
		colormap(cmap)
		cbfreeze(cmap)
		freezeColors

		subaxis(2,2,3)
		imagesc( (SWdist_full-SWdist_full_proc)./24 )
		set(gca,'YDir','normal')
		set(gca,'CLim',CLim_SWdiff)
		cmap = cbrewer('div','RdBu',50);
		colormap(cmap)
		cbfreeze(cmap)
		freezeColors

		subaxis(2,2,4)
		imagesc(Zdeg)
		set(gca,'YDir','normal')
		colorbar
		colormap jet;
		cbfreeze
		freezeColors

		SetFigureProperties([2.5 2],gcf)
		cd([HOMEdir,PRINTdir])
		oname = ['DistSW.DAN.Mar15_05.CumSum.eps'];
		print('-depsc',oname)
		close
	end
	
	%% Individual time steps
% 	if PRINT_FLAG
% 		for k = 7:18
% 			d1 = datenum(2005,3,15,k-1,59,0);
% 			d2 = datenum(2005,3,15,k,1,0);
% 			dind = FindDate(DAN.t,d1,d2);
% 			SWdist_f = SV.*Sf(dind);
% 			SWdist_d = HS.*SHADE.*Sd(dind);
% 			
% 			SWdist_f_proc = SV.*Sf_proc(dind);
% 			SWdist_d_proc = HS.*SHADE.*Sd_proc(dind);
% 			
% 			CLim_SWdwn = [0 round(max(max(SWdist_f_proc+SWdist_d_proc)))];
% 			CLim_SWdiff = [-80 80];
% 
% 			cd([HOMEdir,PROJdir])
% 			inname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat'];
% 			load(inname)
% 			
% 			figure
% 			subaxis(2,2,1)
% 
% 			imagesc(SWdist_f+SWdist_d)
% 			set(gca,'YDir','normal')
% 			set(gca,'CLim',CLim_SWdwn)
% 			cmap = cbrewer('seq','YlOrRd',50);
% 			cmap = cmap(5:end,:);
% 			colormap(cmap)
% 			cbfreeze(cmap)
% 			freezeColors
% 			
% 			subaxis(2,2,2)
% 			imagesc(SWdist_f_proc+SWdist_d_proc)
% 			set(gca,'YDir','normal')
% 			set(gca,'CLim',CLim_SWdwn)
% 			COL = colorbar;
% 			cmap = cbrewer('seq','YlOrRd',50);
% 			cmap = cmap(5:end,:);
% 			colormap(cmap)
% 			cbfreeze(cmap)
% 			freezeColors
% 			
% 			subaxis(2,2,3)
% 			imagesc( (SWdist_f+SWdist_d) - (SWdist_f_proc+SWdist_d_proc))
% 			set(gca,'YDir','normal')
% 			set(gca,'CLim',CLim_SWdiff)
% 			cmap = cbrewer('div','RdBu',50);
% 			colormap(cmap)
% 			cbfreeze(cmap)
% 			freezeColors
% 			
% 			subaxis(2,2,4)
% 			imagesc(Zdeg)
% 			set(gca,'YDir','normal')
% 			colorbar
% 			colormap jet;
% 			cbfreeze
% 			freezeColors
% 			
% 			SetFigureProperties([2.5 2],gcf)
% 			cd([HOMEdir,PRINTdir])
% 			oname = ['DistSW.DAN.Mar15_05.',num2str(k),'.epsc'];
% 			print('-depsc',oname)
% 			close
% 		end		
% 	end
	
end
	
%% SNQ DEM
if SNQ_Hori_FLAG
	cd(SNQDEMdir)
	Z = dlmread('snq_dem_30m.txt');
	Z = flipud(Z);
	
% 	% ASCII raster header
% 	ncols         7212
% 	nrows         3612
% 	xllcorner     -123.0016666667
% 	yllcorner     46.99833333333
% 	cellsize      0.000277777777778
% 	NODATA_value  -9999
	ncols = 7212;
	nrows = 3612;
	xll = -123.0016666667;
	yll = 46.99833333333;
	ds = 0.000277777777778;

	% Create X and Y vectors of coordinates @ center of pixel
	X = linspace(xll,xll+ds*(ncols-1),ncols);
	Y = linspace(yll,yll+ds*(nrows-1),nrows);

	% Convert to UTM (span multiple zones, takes some adjustments, hopefully what I've done isnt' stupid)
	[~,y,~] = deg2utm(Y,ones(size(Y))*X(1));
	[x,~,~] = deg2utm(ones(size(X))*Y(1),X);
	
	cd([HOMEdir,PROJdir])
	save('SNQ_DEM_30m.mat','X','Y','Z','x','y')
end

%% SNQ Horizon
if SNQ_Hori_FLAG
	cd([HOMEdir,SNQDATdir])
	load SNQ14.CLEAN.mat
	SNQ.SWFLAG = SW_Obs_QC(SNQ,0);						% Basic QC
	
	% Observed horizon (from Nic's notes)
	obs_phi = [90 105 120 135 150 165 180 195 210 225 240 255 270] + 16.5;
	obs_theta_trees = [34 33 30 42 35 27 10 0 0 0 0 0 0];
	obs_theta_trees = atand(tand(obs_theta_trees) - 7.5/20);
	obs_theta_horiz = [0 0 0 0 0 0 0 12 10 13 14 10 16];
	
	% Find grid cell containing site in lat/long
	[xm ym] = meshgrid(X,Y);							% Mesh grid (degrees)
	d = ( (xm-SNQ.long).^2 + (ym-SNQ.lat).^2).^(1./2);
	[r,c] = find(min(min(d)) == d);
	[xm ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)
	xsite = xm(r,c);
	ysite = ym(r,c);
	
	% Briefly test that I found the right point
% 	imagesc(Z)
% 	hold all
% 	plot(c,r,'k*')
% 	set(gca,'YDir','normal')

	% Calculate the slope to each grid cell
	d = ( (xm-xm(r,c)).^2 + (ym-ym(r,c)).^2).^(1./2);
	d(r,c) = 0;
	slope = (Z-Z(r,c))./d;
	
	% Covnert slope to angle above the horizon
	theta = atand(slope);

	% For some discreet phi bin find the angle to the horizon
	dphi = pi/36;
	ds = 25;
	distance = 0:ds:10000;
	phi = pi/2:-dphi:-3*pi/2+dphi;
	h = NaN(size(phi));
	for k = 1:length(phi)
		xlin = distance.*cos(phi(k))+xsite;
		ylin = distance.*sin(phi(k))+ysite;
		theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);
		h(k) = nanmax(theta_proj);		
	end
	h(h < 0) = 0;

	% Hinkelman plot w/ Horizon
	Hink_hand = SW_Obs_QC_PlottingTool_SolGeo(SNQ,SNQ.SWFLAG);
	figure(Hink_hand(1))
	hold all
	% Horizon (DEM)
	plot([0:5:355],90-h)
	% Horizon (Nic's obs)
	plot(obs_phi,90-obs_theta_horiz,'k*-')
	plot(obs_phi,90-obs_theta_trees,'g*-')
	figure(Hink_hand(2))
	hold all
	% Horizon (DEM)
	plot([0:5:355],90-h)
	% Horizon (Nic's obs)
	plot(obs_phi,90-obs_theta_horiz,'k*-')
	plot(obs_phi,90-obs_theta_trees,'g*-')
	
	% Print both figures
	if PRINT_FLAG
		cd([HOMEdir,PRINTdir])
		figure(Hink_hand(1))
		print('-depsc',['Trans_Horiz.',s{n},'.eps'])
		figure(Hink_hand(2))
		print('-depsc',['Trans_STD_Horiz.',s{n},'.eps'])
	end
	close all
	clear Z X Y x y xm ym SNQ
end

%% Load DEM (Senator Beck Study Sites)
if SBSA_DEM_FLAG
	% Read raster (converted by ArcGIS to ascii raster). Original data in unzipped GIS files from SBSA.
	cd(SBSADEMdir)
	Z = dlmread('sbsa_dem_30m.txt');
	Z = flipud(Z);
	Z(Z == -9999) = NaN;
	
	% ASCII raster header
% 	ncols         1139
% 	nrows         1419
% 	xllcorner     258066.68592128
% 	yllcorner     4195191.3427875
% 	cellsize      10
% 	NODATA_value  -9999
	ncols = 1139;
	nrows = 1419;
	xll = 258066.68592128;
	yll = 4195191;
	ds = 10;
	
	% Create X and Y vectors of coordinates @ center of pixel
	x = linspace(xll+ds/2,xll+ds*(ncols)+ds/2,ncols);
	y = linspace(yll+ds/2,yll+ds*(nrows)+ds/2,nrows);

	% Convert to UTM (span multiple zones, takes some adjustments, hopefully what I've done isnt' stupid)
	zone = 13;
	[Y,~] = utm2ll(ones(size(y))*x(1),y,zone);
	[~,X] = utm2ll(x,ones(size(x))*y(1),zone);
	
	cd([HOMEdir,PROJdir])
	save('SBSA_DEM_30m.mat','X','Y','Z','x','y')
end


%% Find horizon
if SBSA_Hori_FLAG
		
	% DEM
	cd([HOMEdir,PROJdir])
	load SBSA_DEM_30m.mat
	
	% Point observations
	cd([HOMEdir,SBSADATdir])
	s = [{'SWA'},{'SEN'}];									% Yosemite sites
	load SEN.CLEAN.mat
	SEN.SWFLAG = SW_Obs_QC_v2(SEN,0);						% Basic QC
	SEN.Hinst = 10;											% Instruments 10m above ground
	load SWA.CLEAN.mat
	SWA.SWFLAG = SW_Obs_QC_v2(SWA,0);						% Basic QC
	SWA.Hinst = 6;											% Instrument 6m above ground
	MET.SEN = SEN;
	MET.SWA = SWA;
	
	for n = 1:2
		% Find grid cell containing site in lat long
		[xm ym] = meshgrid(X,Y);							% Mesh grid (degrees)
		d = ( (xm-MET.(s{n}).long).^2 + (ym-MET.(s{n}).lat).^2).^(1./2);
		[r,c] = find(min(min(d)) == d);
		[xm ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)
		xsite = xm(r,c);
		ysite = ym(r,c);

		% Briefly test that I found the right point
		figure(1)
		imagesc(Z)
		hold all
		plot(c,r,'k*')
		set(gca,'YDir','normal')

		% Calculate the slope to each grid cell
		d = ( (xm-xm(r,c)).^2 + (ym-ym(r,c)).^2).^(1./2);
		d(r,c) = 0;
		Z(r,c) = Z(r,c) + MET.(s{n}).Hinst;
		slope = (Z-Z(r,c))./d;
		% Angle for each slope
		theta = atand(slope);

		% For some discreet phi bin find the angle to the horizon
		dphi = pi/36;
		ds = 25;
		distance = 0:ds:10000;
		phi = pi/2:-dphi:-3*pi/2+dphi;
		h.(s{n}) = NaN(size(phi));
		for k = 1:length(phi)
			xlin = distance.*cos(phi(k))+xsite;
			ylin = distance.*sin(phi(k))+ysite;
			theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);
			h.(s{n})(k) = nanmax(theta_proj);		
		end
		h.(s{n})(h.(s{n}) < 0) = 0;
		h_AZ = [0:dphi:2*pi-dphi].*180./pi;

		% Hinkelman plots w/ horizon
		Hink_hand = SW_Obs_QC_PlottingTool_SolGeo(MET.(s{n}),~logical(MET.(s{n}).SWFLAG(:,7)));
		figure(Hink_hand(1))
		hold all
		plot(h_AZ,90-h.(s{n}))
		figure(Hink_hand(2))
		hold all
		plot(h_AZ,90-h.(s{n}))
		
		% Print both figures
		if PRINT_FLAG
			cd([HOMEdir,PRINTdir])
			figure(Hink_hand(1))
			print('-depsc',['Trans_Horiz.',s{n},'.eps'])
			figure(Hink_hand(2))
			print('-depsc',['Trans_STD_Horiz.',s{n},'.eps'])
		end
% 		close all
	end	
	clear Z X Y x y xm ym SEN SWA
end
