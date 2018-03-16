function example()

% A step-by-step example of how to use these functions to calculate the
% relative salience of a target. This example assumes that the Image 
% Processing Toolbox is installed.

% Compile the mex files localmaxima.cc and subsample.cc, if necessary (you 
% should only need to do this once)
%mex localmaxima.cc 
%mex subsample.cc

% Load example input images, including a 600x800 luminance image and a 
% 600x800x4 LMSU 'colour' image, and the vertices of a polygon delimiting
% the target
exampleData = load('example.mat');
luminanceImage = exampleData.luminanceImage;
colourImage = exampleData.colourImage;
mothVertices = exampleData.mothVertices;

% Define the weights of the three feature channels
colourWeight = 1.0;
luminanceWeight = 1.0;
orientationWeight = 1.0;

% Run the modified Itti-Koch model
out = ittikochmod(luminanceImage,colourImage,colourWeight, ...
                  luminanceWeight,orientationWeight);

% Extract the overall saliency map and the luminance feature map
saliencyMap = out.master_map;
luminanceMap = out.top_level_feat_maps{2};

% Calculate the relative salience of the target from the overall saliency
% map
nBins = 100; % Number of histogram bins to use
targetSalience = targetsalience(saliencyMap,mothVertices,nBins);

% Plot the luminance input image
figure;
subplot(2,2,1), imshow(luminanceImage); 
title('Luminance Image');

% Plot the L, M and S channels of the colour input image
subplot(2,2,2), imshow(colourImage(:,:,1:3)); 
title('LMS Colour Image');

% Plot the overall saliency map, and highlight the target region
subplot(2,2,3), imshow(gray2ind(saliencyMap,256),hot(256));
title('Saliency Map');
hold on
line([mothVertices(:,2)' mothVertices(1,2)], ...
     [mothVertices(:,1)' mothVertices(1,1)], ...
     'LineStyle','-','LineWidth',1,'Color','w');
     % Add text giving the relative salience if the target
     text(10,25,sprintf('Target Salience = %1.3f',targetSalience),'Color','w');
hold off

% Plot the luminance feature map, which in this case contributes heavily to
% overall salience
subplot(2,2,4), imshow(gray2ind(luminanceMap,256),hot(256)); 
title('Luminance Feature Map');
