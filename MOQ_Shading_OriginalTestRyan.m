%% PURPOSE
% Find shading for applicaitons in the MntObsQC project. This file is where
% much of the testing of the shading algorithm was originally performed. NOTE:
% The logical argument for shading is an AND logical and not an OR logical.
% This is wrong and will need to be fixed in any future use of this script.
clear
close all
clc

%% Solar Radiation Distribution
% Sky view is invariant w/time, calculate once
    % This is used with the diffuse solar radiation and is fed into the
    % hillshade algorithm
% Hillshade varies through year, need for each hour of day
    % We need to check expression against Frew and Dozier 1990
    % Calculate once per month or every two weeks (following DHSVM)
    % this is multiplied by the direct beam solar radiation.
% Local slope effect: invariant in time but needs to be calculated for each
    % timestep because it depends on the weather. If its cloudy we're going to
    % have more difuse radiation and the effect will be smaller because
    % this is going to make a difference for the direct beam solar
    % radiation
% Direct/Diffuse: each time step depends on the cloudiness.


addpath(genpath('/Users/wcurrier/Desktop/Matlab_Functions/'))
%% Load DEM (Dungeness)
	% Read in geotiff raster (converted by ArcGIS from NED ArcGIS compatible format to
	% ascii raster). Original data in back up of Mark's folder.
    
	Z = dlmread('olympdem30mdd2.txt');
	Z = flipud(Z);
    
    xll = -123.11545470935;
    yll = 47.836458757832;
    ds  = 0.00009387402;
    ncols = 660;
    nrows = 682;
    
	Z(Z == -9999) = NaN;

%     xll = 491376.402;
%     yll = 5298140.697;
%     ds  = 8.5482343;
%     ncols = 541;
%     nrwos = 831;
	
    
    % Create X and Y vectors of coordinates @ center of pixel
	X = linspace(xll+ds/2,xll+ds*(ncols)+ds/2,ncols);
	Y = linspace(yll+ds/2,yll+ds*(nrows)+ds/2,nrows);
    
    [~,y,~] = deg2utm(Y,ones(size(Y))*X(1));
	[x,~,zone] = deg2utm(ones(size(X))*Y(1),X);
	
	save('OlympDung_DEM_30m.mat','X','Y','Z','x','y')

%% Find horizon
	% Point observations
    MET.s.long=-123.078798452019; % Decimal Degrees location of Dungeness SNOTEL
    MET.s.lat=47.8722382562592;   % Decimal Degrees location of Dungeness SNOTEL
    
%     MET.s.long=5301854.463;       % UTM Zone 10 location of Dungeness SNOTEL
%     MET.s.lat=494017.372;         % UTM Zone 10 location of Dungeness SNOTEL

	% DEM
	load 'OlympDung_DEM_30m.mat'

		% Find grid cell containing site in lat long
		[xm ym] = meshgrid(X,Y);                               % Mesh grid (degrees)
		d = ( (xm-MET.s.long).^2 + (ym-MET.s.lat).^2).^(1./2); % find the distance between all points and the site you're looking at
		[r,c] = find(min(min(d)) == d);                        % find the minimum distance or center of closest grid cell
		[xm ym] = meshgrid(x,y);                               % Mesh grid (UTM - meters)
		xsite = xm(r,c);                                       % get coordinates of center
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
		h.s = NaN(size(phi));
		for k = 1:length(phi)
			xlin = distance.*cos(phi(k))+xsite;
			ylin = distance.*sin(phi(k))+ysite;
			theta_proj = interp2(xm,ym,theta,xlin,ylin,'*nearest',NaN);
			h.n(k) = nanmax(theta_proj);		
		end
% 		h.(s{n})(h.(s{n}) < 0) = 0;
		h_AZ = [0:dphi:2*pi-dphi].*180./pi;
        
        
        METtmp=MET;
        	%% Distributed
            load DungMet.mat
            EL = MET.El;
            AZ = MET.AZ;
            SWdwn = MET.SWdwn;
            METtmp.t=MET.t;
            
            MET=METtmp;
	% DEM
	load OlympDung_DEM_30m.mat
	[xm ym] = meshgrid(x,y);					% Mesh grid (UTM - meters)
	% Delineate terrain surrounding Dana
	[Xm,Ym] = meshgrid(X,Y);
	% Degredate to 150m
	xdeg = xm(1,1):150:xm(1,end);
	ydeg = ym(1,1):150:ym(end,1);
	[xmdeg,ymdeg] = meshgrid(xdeg,ydeg);
	Zdeg = interp2(xm,ym,Z,xmdeg,ymdeg);
    
    % Within degraded coordinates find closest grid cell to xsite
    [~, indexX] = min(abs(xdeg-xsite));[~, indexY] = min(abs(ydeg-ysite));
    xsitedeg = xdeg(indexX);ysitedeg=ydeg(indexY);
    
% 	Z = Z(indy,indx);
        
     %% Hillshade
		SV = metDist_Skyview(xmdeg,ymdeg,Zdeg);
		save('Skyview_150m.Dung.mat','SV')
% load Skyview_150m.Dung.mat
	[EL, AZ, SOLDIST, HA] = SolarGeometry_v2(MET.t,MET.s.lat,MET.s.long,8);

		% Hillshade on March 15th
		for k = 5:20
			d2 = datenum(2014,3,15,k,0,0);
			dind = find(MET.t(:,7)==d2);
			SZAlocal = 90-EL(dind);
			AZlocal = AZ(dind);
			[HS,SHADE] = metDist_Hillshade(SZAlocal,AZlocal,xmdeg,ymdeg,Zdeg);
			disp('')
			oname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat']
			save(oname,'HS','SHADE')
        end
        
        [Sd,Sf] = DirDifPartition_DHSVM(EL,SWdwn);

        
        	% Cumulative difference in W/m-2
	xdim = size(SV,1);
	ydim = size(SV,2);
	SWdist_f = NaN(xdim,ydim,12);
	SWdist_d = NaN(xdim,ydim,12);
	for k = 7:18
			d2 = datenum(2014,3,15,k,0,0);
			dind = find(MET.t(:,7)==d2);
			
			inname = ['Hillshade.DAN.Mar15_05.',num2str(k),'.mat'];
			load(inname)
			
			% Raw SW
			SWdist_f(:,:,k-6) = SV.*Sf(dind);
			SWdist_d(:,:,k-6) = HS.*SHADE.*Sd(dind);

	end
	% Sum over whole day - divide by 24
	SWdist_full = sum(SWdist_f+SWdist_d,3);
    %%
        
    figure
    
        subplot(1,3,1),imagesc(Z),colorbar,hold on
        title('DEM','Fontsize',20)
        plot(c,r,'k*')
		set(gca,'YDir','normal')

        subplot(1,3,2),imagesc(SV),colorbar,hold on
        title('Sky View','Fontsize',20)
        plot(indexX,indexY,'k*')
        set(gca,'YDir','normal')

        subplot(1,3,3),imagesc(SWdist_full),colorbar
        title('Distributed SW Full, March 15','Fontsize',20),hold on
        plot(indexX,indexY,'k*')
        set(gca,'YDir','normal')

    figure
    for ii=1:12
        subplot(1,3,1),imagesc(Z),colorbar,hold on
        title('DEM','Fontsize',20)
        plot(c,r,'k*')
		set(gca,'YDir','normal')
        
        subplot(1,3,2),imagesc(SWdist_d(:,:,ii)),colorbar
        title('Direct SW, March 15','Fontsize',20),hold on
        plot(indexX,indexY,'k*')
        set(gca,'YDir','normal')
        
        subplot(1,3,3),imagesc(SWdist_f(:,:,ii)),colorbar,hold on
        title('Difuse SW','Fontsize',20)
        cmap=flipud(cbrewer('div','RdBu',50));
        colormap(cmap)
        plot(indexX,indexY,'k*')
        set(gca,'YDir','normal')
        pause()
    end
