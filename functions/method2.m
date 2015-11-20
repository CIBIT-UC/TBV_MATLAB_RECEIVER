%% Method 2 - All Voxels
function [ROImeansM2,timeM2] = method2(n_rois,tbvNetInt,coordsOfVoxelsOfROI,xDim,yDim,timePoint)
tic

voxelsValsM2 = tbvNetInt.tGetValueOfAllVoxelsAtTime( timePoint );
% voxelsMatrix = reshape(voxelsVals, [zDim yDim xDim]);

ROImeansM2 = zeros(1,n_rois);

for i=1:n_rois
    
    ROImeansM2(i) = mean( voxelsValsM2(coordsOfVoxelsOfROI{i}(:,3)*xDim*yDim ...
        + coordsOfVoxelsOfROI{i}(:,2)*xDim + coordsOfVoxelsOfROI{i}(:,1) + 1) );
    
end

timeM2 = toc;

end