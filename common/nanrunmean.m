function data_run=nanrunmean(data,n_points)
%
% function data_run = nanrunmean(data,n_points)
%
% function that, given a vertical profile of a variable, smooths it with a 
% running mean (neglecting nans)
%
% input parameters:
%   data: the two or one column profile with values of:
%       2 columns-> depths OtherVariable
%       1 column -> OtherVariable 
%   n_points: the number of points used in the running mean (OPTIONAL:default=1)
%  
% output parameter:
%   data_run= the new profile smoothed with running mean.

if nargin<2
    n_points=1;
end
points_ext=(n_points-1)/2;
siz = size(data);
nd=siz(1);
if siz(2) == 2
    data_sort=sortrows(data,1);
    data_run=zeros(size(data_sort));
    data_run(:,1)=data_sort(:,1);
else
    data_sort = [data data];
    data_run=zeros(size(data_sort));
end

for i = 1:points_ext                                               %SMOOTHING of datas with running mean
    step=i-1;
    data_run(i,2)=nanmean(data_sort(i-step:i+step,2));
end
for i = nd-points_ext+1:nd
    step=nd-i;
    data_run(i,2)=nanmean(data_sort(i-step:i+step,2));
end
for i = points_ext+1:nd-points_ext
    step=points_ext;
    data_run(i,2)=nanmean(data_sort(i-step:i+step,2));
end

if siz(2) ==1
    data_run(:,1) = [];
end
