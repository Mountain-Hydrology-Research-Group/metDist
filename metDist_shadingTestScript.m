%% PURPOSE
% Let's make a matlab script that we can test metDist functions with. Using
% Tuolumne as our test area.
clear all
close all
clc

%% Flags
% Yosemite
flagDEM = 1;						% Process DEM data (Dana Meadows)
flagDana = 0;                       % Plot horizon on Hinkelman plot of Dana SW
flagDist = 1;
flagPrint = 1;
% flagSkyview = 0;
% flagHillshade = 0;

%% Directories
% General directories
dirHOME = '/home/lapok/gdrive/';
dirPrint = [dirHOME,'SnowHydrology/proj/metDist/Graphics'];
dirProj = [dirHOME,'SnowHydrology/proj/metDist'];

% Yosemite directories
dirGroundObs = [dirHOME,'SnowHydrology/proj/yosemiteData'];
dirDEM = '/home/lapok/Spatial.Maps/DEM/Yosemite';

%% Load DEM (Yosemite)
if flagDEM
	% Read raster (converted by ArcGIS from NED ArcGIS compatible format to
	% ascii raster). Original data in back up of Mark's folder.
	cd(dirDEM)
	Z = dlmread('Yosemite.DEM30m.txt');
	Z = flipud(Z);
	
	% ASCII raster headerpwd
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
	
	cd(dirProj)
	save('YOS_DEM_30m.mat','X','Y','Z','x','y')
end

%% Skyview @ Dana
if flagDana
	cd(dirDEM)
	load YOS_DEM_30m.mat
	
    lat = 37.8970;
    lon = -119.2570;

    % Find grid cell containing site in lat long
    [xm, ym] = meshgrid(X,Y);							% Mesh grid (degrees)
    d = ( (xm-lon).^2 + (ym-lat).^2).^(1./2);
    [r,c] = find(min(min(d)) == d);
    [xm, ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)
    xsite = xm(r,c);
    ysite = ym(r,c);

    % Briefly test that I found the right point
    figure(1)
    imagesc(Z)
    hold all
    plot(c,r,'k*')
    set(gca,'YDir','normal')

    %% Slope
    
    % Calculate the slope to each grid cell
    d = ( (xm-xm(r,c)).^2 + (ym-ym(r,c)).^2).^(1./2);
    d(r,c) = 0;
    slope = (Z-Z(r,c))./d;
    % Angle for each slope
    theta = atan(slope);                        % Angle for each slope

    %% Horizon @ Dana
    ds = 25;
    distance = 0:ds:10000;
    
    % Discrete azimuth angles to work over 
    nphi = 36; % Number of azimuth slices to work over 
    dphi = 2*pi/nphi; % Width of azimuth bins
    phi = 0+dphi/2:dphi:2*pi-dphi/2; % Phi bin angles
    h = NaN(size(phi));
    for k = 1:length(phi)
        xlin = distance.*cos(phi(k))+xsite;
        ylin = distance.*sin(phi(k))+ysite;
        theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);
        h(k) = max(nanmax(theta_proj),0);		
    end
    h(h < 0) = 0;
    
    %% Skyview
    SV = 0;
    for k = 1:nphi
       SV = SV +dphi*cos(h(k));
    end
    SV = SV/(2*pi);
    
end

% Distribute shortwave to the surrounding terrain at Dana Meadows
if flagDist
	%% load ground observation data
	cd(dirGroundObs)
	load('DanaMeadows.met_data_DANA_for_modeling.UTC_8.mat')
	
    %% Mesh grid of utm x and y
    cd(dirDEM)
	load YOS_DEM_30m.mat
    [xm, ym] = meshgrid(x,y);							% Mesh grid (UTM - meters)

    
	%% Hillshade
	cd(dirProj)
    if Skyview_FLAG
		SV = MOQ_Skyview(xm,ym,Z);
		cd(dirProj)
		save('Skyview_yosemite30m.mat','SV')
        clear SV
	end
	
	%% Hillshade\
    cd(dirProj)
	if Hillshade_FLAG
		% Hillshade on March 15th
		for k = 5:20
			d1 = datenum(2005,3,15,k-1,59,0);
			d2 = datenum(2005,3,15,k,1,0);
			dind = FindDate(DAN.t,d1,d2);
			SZAlocal = 90-DAN.EL(dind);
			AZlocal = DAN.AZ(dind);
			[HS,SHADE] = metDist_Hillshade(SZAlocal,AZlocal,xm,ym,Z);
			disp('')
			oname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat'];
			save(oname,'HS','SHADE')
            clear HS SHADE
		end
	end
	
	%% Distribute
	% Non-QC'ed obs
	[Sd,Sf] = DirDifPartition_DHSVM(DAN.EL,DAN.SWdwnGapFill);
	
	% QCed and processed obs
% 	[Sd_proc,Sf_proc] = DirDifPartition_DHSVM(DAN.EL,DAN.SWdwn);
	
	% Skyview fraction (found earlier)
	cd(dirProj)
	load Skyview_yosemite30m.mat
	
	% Cumulative difference in W/m-2
	xdim = size(SV,1);
	ydim = size(SV,2);
	SWdist_f = NaN(xdim,ydim,12);
	SWdist_d = NaN(xdim,ydim,12);
% 	SWdist_f_proc = NaN(xdim,ydim,12);
% 	SWdist_d_proc = NaN(xdim,ydim,12);
	for k = 7:18
			d1 = datenum(2005,3,15,k-1,59,0);
			d2 = datenum(2005,3,15,k,1,0);
			dind = FindDate(DAN.t,d1,d2);
			
			cd(dirProj)
			inname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat'];
			load(inname)
			
			% Raw SW
			SWdist_f(:,:,k-6) = SV.*Sf(dind);
			SWdist_d(:,:,k-6) = HS.*SHADE.*Sd(dind);
			% Proc SW
% 			SWdist_f_proc(:,:,k-6) = SV.*Sf_proc(dind);
% 			SWdist_d_proc(:,:,k-6) = HS.*SHADE.*Sd_proc(dind);
	end
	% Sum over whole day - divide by 24
	SWdist_full = sum(SWdist_f+SWdist_d,3);
	SWdist_full_proc = sum(SWdist_f_proc+SWdist_d_proc,3);
	
	%% Figure of cumulative difference (daily average - Wm-2)
	if flagPrint
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
end
