%% Method 4 - Voxels 1by1
function [ROImeansM4,timeM4] = method4(n_rois,tbvNetInt,coordsOfVoxelsOfROI,timePoint)
tic

ROImeansM4 = zeros(1,n_rois);

for i=1:n_rois
    
    temp = zeros(1,size(coordsOfVoxelsOfROI{i},1));
    for j=1:size(coordsOfVoxelsOfROI{i},1)
        
        temp(j) = tbvNetInt.tGetValueOfVoxelAtTime(coordsOfVoxelsOfROI{i}(j,1) ...
            ,coordsOfVoxelsOfROI{i}(j,2),coordsOfVoxelsOfROI{i}(j,3),timePoint);
    end
    
    ROImeansM4(i) = mean( temp );
    
end

timeM4 = toc;

end