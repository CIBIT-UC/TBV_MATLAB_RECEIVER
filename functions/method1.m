%% Method 1 - TBV
function [ROImeansM1,timeM1] = method1(n_rois,tbvNetInt,timePoint)
tic

ROImeansM1 = zeros(1,n_rois);

for i=1:n_rois
    ROImeansM1(i) = tbvNetInt.tGetMeanOfROIAtTimePoint( i, timePoint );
end

timeM1 = toc;

end