function out = ittikochmod(imgi,imgcol,cw,iw,ow)
% Computes a saliency map for an image using a modified version of the
% standard Itti, Koch & Niebur algorithm (Itti, L., Koch, C. & Niebur, E. 
% (1998) A model of saliency-based visual attention for rapid scene 
% analysis. IEEE Transactions on Pattern Analysis and Machine Intelligence, 
% 20, 1254-1259). This function is based heavily on code available in the 
% implementation available from J. Harel, A Saliency Implementation in 
% MATLAB: http://www.klab.caltech.edu/~harel/share/gbvs.php.
%
% This function takes as inputs an NxM luminance (or intensity) image IMGI, 
% a NxMx4 'colour' image IMGCOL, and weighting factors for the colour,
% luminance and orietntation feature maps (CW, IW and OW, respectively),
% and outputs a struct OUT with, among others, the following fields:
% master_map: The overall saliency map
% top_level_feat_maps{1}: The colour feature map
% top_level_feat_maps{2}: The luminance feature map
% top_level_feat_maps{3}: The orientation feature map


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Define parameter values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
param.colorWeight = cw;
param.intensityWeight = iw;
param.orientationWeight = ow;
param.ittiCenterLevels = [2 3 4];
param.ittiDeltaLevels = [3 4];
param.maxcomputelevel = 8;
param.gaborangles = [0 45 90 135];
for i = 1:length(param.gaborangles)
	theta = param.gaborangles(i);
    param.gaborFilters{i}.g0 = makegaborfilters(theta,0,7/3,1,-1,pi,0);
    param.gaborFilters{i}.g90 = makegaborfilters(theta,90,7/3,1,-1,pi,0);
end
param.salmapsize = [75 100];


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1 : Compute raw feature maps from image
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
imgLum = {};
imgLum{1} = subsample(imgi);
imgCol = {};
for ci = 1:size(imgcol,3)
	imgCol{ci}{1} = subsample(imgcol(:,:,ci));
end
for i = 2:param.maxcomputelevel
	imgLum{i} = subsample(imgLum{i-1});
    for ci = 1:size(imgcol,3)
        imgCol{ci}{i} = subsample(imgCol{ci}{i-1});
    end
end

% Luminance (intensity)
obj = {};
obj.info.weight = param.intensityWeight;
obj.info.numtypes = 1;
obj.info.descriptions{1} = 'Intensity';
obj.description = 'intensity';
obj.maps = {};
obj.maps.val = {};
for ti = 1:obj.info.numtypes
    obj.maps.val{ti} = {};
    for lev = 2:param.maxcomputelevel
        mymap = imgLum{lev};
        obj.maps.origval{ti}{lev} = mymap;
        map = imresize(mymap,param.salmapsize,'bicubic');
        obj.maps.val{ti}{lev} = map;
    end
end
rawfeatmaps.intensity = obj;

% Colour
obj = {};
obj.info.weight = param.colorWeight;
obj.info.numtypes = 6;
obj.info.descriptions{1} = 'L-M';
obj.info.descriptions{2} = 'L-S';
obj.info.descriptions{3} = 'L-U';
obj.info.descriptions{4} = 'M-S';
obj.info.descriptions{5} = 'M-U';
obj.info.descriptions{6} = 'S-U';
obj.description = 'color';
obj.maps = {};
obj.maps.val = {};
for ti = 1:obj.info.numtypes
    obj.maps.val{ti} = {};
    for lev = 2:param.maxcomputelevel
        imgU = imgCol{4}{lev}; %U image
        imgS = imgCol{3}{lev}; %S image
        imgM = imgCol{2}{lev}; %M image
        imgL = imgCol{1}{lev}; %L image
        switch ti
            case 1
                mymap = safedivide(abs(imgL-imgM),imgLum{lev});
                mymap(mymap<0) = 0;
            case 2
                mymap = safedivide(abs(imgL-imgS),imgLum{lev});
                mymap(mymap<0) = 0;
            case 3
                mymap = safedivide(abs(imgL-imgU),imgLum{lev});
                mymap(mymap<0) = 0;
            case 4
                mymap = safedivide(abs(imgM-imgS),imgLum{lev});
                mymap(mymap<0) = 0;
            case 5
                mymap = safedivide(abs(imgM-imgU),imgLum{lev});
                mymap(mymap<0) = 0;
            case 6
                mymap = safedivide(abs(imgS-imgU),imgLum{lev});
                mymap(mymap<0) = 0;
        end
        obj.maps.origval{ti}{lev} = mymap;
        map = imresize(mymap,param.salmapsize,'bicubic');
        obj.maps.val{ti}{lev} = map;
    end
end
rawfeatmaps.color = obj;

% Orientation
obj = {};
obj.info.weight = param.orientationWeight;
obj.info.numtypes = length(param.gaborFilters);
for i = 1:length(param.gaborFilters)
	obj.info.descriptions{i} = sprintf('Gabor Orientation %g',param.gaborangles(i));
end
obj.description = 'orientation';
obj.maps = {};
obj.maps.val = {};
for ti = 1:obj.info.numtypes
    obj.maps.val{ti} = {};
    for lev = 2:param.maxcomputelevel
        gaborFilters = param.gaborFilters;
        j = ti;
        f0 = myconv2(imgLum{lev},gaborFilters{j}.g0);
        f90 = myconv2(imgLum{lev},gaborFilters{j}.g90);
        mymap = abs(f0) + abs(f90);
        mymap = attenuateborders(mymap,13);
        obj.maps.origval{ti}{lev} = mymap;
        map = imresize(mymap,param.salmapsize,'bicubic');
        obj.maps.val{ti}{lev} = map;
    end
end
rawfeatmaps.orientation = obj;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2 : Compute activation maps from feature maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
mapnames = fieldnames(rawfeatmaps);
mapweights = zeros(1,length(mapnames));
map_types = {};
allmaps = {};
i = 0;
for fmapi = 1:length(mapnames)
    mapsobj = eval(['rawfeatmaps.' mapnames{fmapi} ';']);
    numtypes = mapsobj.info.numtypes;
    mapweights(fmapi) = mapsobj.info.weight;
    map_types{fmapi} = mapsobj.description;
    for typei = 1:numtypes
        for centerLevel = param.ittiCenterLevels
            for deltaLevel = param.ittiDeltaLevels
                i = i + 1;
                center_ = mapsobj.maps.origval{typei}{centerLevel};
                sz_ = size(center_);
                surround_ = imresize( mapsobj.maps.origval{typei}{centerLevel+deltaLevel}, sz_ , 'bicubic' );
                allmaps{i}.map = (center_ - surround_).^2;
                allmaps{i}.maptype = [fmapi centerLevel deltaLevel];
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3 : Normalise activation maps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
norm_maps = {};
for i = 1:length(allmaps)
	norm_maps{i}.map = maxnormalise(mat2gray(imresize(allmaps{i}.map,param.salmapsize,'bicubic')));
	norm_maps{i}.maptype = allmaps{i}.maptype;
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
% STEP 4 : Average across maps within each feature channel
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%for j = 1:length(norm_maps)
    %norm_maps{j}.map = norm_maps{j}.map*param.map_weight(j);
%end
comb_norm_maps = {};
cmaps = {};
for i = 1:length(mapnames)
    cmaps{i}=0; 
end
Nfmap = cmaps;
for j = 1:length(norm_maps)
    map = norm_maps{j}.map;
    fmapi = norm_maps{j}.maptype(1);
    Nfmap{fmapi} = Nfmap{fmapi} + 1;
    cmaps{fmapi} = cmaps{fmapi} + map;
end
for fmapi = 1:length(mapnames)
	cmaps{fmapi} = maxnormalise(cmaps{fmapi});
	comb_norm_maps{fmapi} = cmaps{fmapi};
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5 : Sum across feature channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
master_idx = length(mapnames) + 1;
comb_norm_maps{master_idx} = 0;
for fmapi = 1:length(mapnames)
	comb_norm_maps{master_idx} = comb_norm_maps{master_idx} + cmaps{fmapi} * mapweights(fmapi);
end
master_map = comb_norm_maps{master_idx};
master_map = attenuateborders(master_map,4);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
master_map = mat2gray(imresize(master_map,size(imgi)));
feat_maps = {};
for i = 1 : length(mapnames)
  feat_maps{i} = mat2gray(imresize(comb_norm_maps{i},size(imgi)));
end
intermed_maps = {};
for i = 1 : length(allmaps)
 allmaps{i}.map = mat2gray(allmaps{i}.map);
 norm_maps{i}.map = mat2gray(norm_maps{i}.map);
end
intermed_maps.featureActivationMaps = allmaps;
intermed_maps.normalizedActivationMaps = norm_maps;

out = {};
out.master_map = master_map;
out.top_level_feat_maps = feat_maps;
out.map_types = map_types;
out.intermed_maps = intermed_maps;
out.rawfeatmaps = rawfeatmaps;
