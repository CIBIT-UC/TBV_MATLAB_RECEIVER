%% Method 3 - All Raw Voxels
function [ROImeansM3,timeM3] = method3(n_rois,tbvNetInt,coordsOfVoxelsOfROI,xDim,yDim,timePoint)
tic

voxelsValsM3 = tbvNetInt.tGetRawValueOfAllVoxelsAtTime( timePoint );

ROImeansM3 = zeros(1,n_rois);

for i=1:n_rois
    
    ROImeansM3(i) = mean( voxelsValsM3(coordsOfVoxelsOfROI{i}(:,3)*xDim*yDim ...
        + coordsOfVoxelsOfROI{i}(:,2)*xDim + coordsOfVoxelsOfROI{i}(:,1) + 1) );
    
end

timeM3 = toc;
end