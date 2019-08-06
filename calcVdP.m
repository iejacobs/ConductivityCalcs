function [ condata ] = calcVdP( folder, thickness, IRevTolerance, filelist )
%CALCVDP calculates VdP conductivity assuming isotropic square sample
%   Outputs data as a struct containing raw data and calculated
%   conductivity, along with whether sample passed consistency checks and
%   some other information
%
%   Thickness is given in nm, absolute thickness error can be specified in
%   the thickness argument in the form [thickness, thicknessErr]
%
%   Conductivity is given in S/cm
%
%   Errors are given as relative errors unless otherwise stated
%
%   IRevTolerance specifies the acceptable error upon reversal of the
%   current direction (default is 0.01 = 1%)
%
%   filelist is a 4 element string array giving the filenames for the four
%   measurement data files. If no extension is given '.txt' is tried
%   automatically. Default values are ["25" "26" "35" "36"]
%
%   v1.1 Ian Jacobs, Nov 2018


%% Data import

%Set filelist to whatever the four files you save are named, or rename
%your files to fit the default format. The first and fourth files are the
%horizontal resistance measurements, while the second and third files are
%the vertical resistance measurements.

if nargin < 4
    filelist = ["25" "26" "35" "36"];
end
if nargin < 3
    IRevTolerance = .01;
end

%import data
for i=1:1:4
    try
        filepath = strcat(folder,filelist(i));
        [condata.Vdrain(:,i),condata.Idrain(:,i),condata.Vp1(:,i),condata.Vp2(:,i)] = ...
        importCond(filepath);
        condata.filepath = filepath;
    catch
        try
            filepath = strcat(folder,'/',filelist(i));
            [condata.Vdrain(:,i),condata.Idrain(:,i),condata.Vp1(:,i),condata.Vp2(:,i)] = ...
            importCond(filepath);
            condata.filepath = filepath;
        catch
            try
                filepath = strcat(folder,'/',filelist(i),'.txt');
                [condata.Vdrain(:,i),condata.Idrain(:,i),condata.Vp1(:,i),condata.Vp2(:,i)] = ...
                importCond(filepath);
                condata.filepath = filepath;
            catch
                filepath = strcat(folder,filelist(i),'.txt');
                [condata.Vdrain(:,i),condata.Idrain(:,i),condata.Vp1(:,i),condata.Vp2(:,i)] = ...
                importCond(filepath);
                condata.filepath = filepath;
            end
        end
    end

    condata.Rfull(:,i) = abs(condata.Vp1(:,i)-condata.Vp2(:,i))./abs(condata.Idrain(:,i));
end

condata.Rsample = condata.Rfull;

%% Get sample thickness and thicknessErr

%first component of thickness must always be thickness in NM
condata.thickness = thickness(1);

%use thicknessErr if given in function call
if length(thickness) > 1
    %if thickness error is specified in 2nd component of thickness
    %argument, set thicknessErr to this value
    condata.thicknessAbsErr = thickness(2);
elseif length(thickness) == 1
    %if no thickness error is given, assume the error is 10%
    condata.thicknessAbsErr = condata.thickness.*0.1;
else
    ME = MException('CalcVDP:invalidThickness', ...
        'Thickness argument must either be of the form [thickness] or [thickness, thicknessErr]',str);
    throw ME
end

%Define relative thickness error
condata.thicknessErr = condata.thicknessAbsErr./condata.thickness;


%% Current reversal consistency check and removal of bad datapoints
% Shows percent error in calculated resistance upon switching current
% direction
condata.checkIrev = abs(condata.Rsample - flipud(condata.Rsample))./condata.Rsample;

%identifies data points which have error greater than 1% (typically values
%near Vsd = 0) and removes them from the dataset
failind = condata.checkIrev > IRevTolerance;
[failrows, ~] = find(failind);
failrows = unique(failrows);
condata.Rsample(failrows,:) = [];

%if more than 10 rows fail at >1%, the dataset doesn't pass
if length(failrows) > 10
    condata.passedIrevCheck = false;
else
    condata.passedIrevCheck = true;
end

%% Reciprocity check
condata.Redgemean = mean(condata.Rsample);

%check that vertical and horizontal edges pass reciprocity checks
condata.checkRecip(1) = abs(1 - condata.Redgemean(1)/condata.Redgemean(4));
condata.checkRecip(2) = abs(1 - condata.Redgemean(2)/condata.Redgemean(3));

if condata.checkRecip > .05
    condata.passedRecipCheck = false;
else
    condata.passedRecipCheck = true;
end

%% Define VdP equation and guess its solution

%define horizontal and vertical resistivities
condata.RhorizMean = mean(condata.Redgemean([1 4]));
condata.RvertMean = mean(condata.Redgemean([2 3]));

%also define the full horizontal and vertical resistance data (for error
%calcs later)
condata.Rhoriz = condata.Rsample(:,[1 4]);
condata.Rvert = condata.Rsample(:,[2 3]);

%define van der pauw equation
condata.f = @(Rs) exp((-pi()*condata.RhorizMean)./Rs) + exp((-pi()*condata.RvertMean)./Rs) - 1;

%for a square geometry the solution to the VdP equation is Rsheet =
%Rvert*pi/log(2) ~= Rvert*4.53
condata.Rguess = mean([condata.RhorizMean condata.RvertMean]).*pi/log(2);


%% Solve VdP equation by Newton's Method (very slow but does work)
%NOT UP TO DATE WITH CURRENT VERSION OF CODE
% tol = 1E-3;
%
% condata.dfdRs = @(Rs) (pi./(Rs.^2)) .* (condata.RhorizMean .* exp((-pi().*condata.RhorizMean)./Rs)) + ...
%    (condata.RvertMean .* exp((-pi().*condata.RvertMean)./Rs));
% 
% condata.Rsheet = condata.Rguess;
% 
% while abs(condata.f(condata.Rsheet)) >= tol
%     condata.Rsheet = condata.Rsheet - condata.f(condata.Rsheet)./condata.dfdRs(condata.Rsheet);
% end

%% Solve VdP equation by interpolation (fast)

%Sets how accurately the data is interpolated. This is probably way more
%accurate that necessary but still runs fast
npoints = 1E7;

%define vector of reasonable resistance values
r = linspace((condata.Rguess./10),(condata.Rguess.*10),npoints);

%find value of r that sets VdP equation to 0 by interpolation
condata.Rsheet = interp1(condata.f(r),r,0);

%find closest value in r vector
[~, r0] = min(abs(r - condata.Rsheet));

%estimate of the error introduced by interpolation. This is probably a significant
%overestimation, since as long as the rate of change in f is slow with R
%then interpolation will still be accurate even if the r vector isn't
%granular enough, but this is only used to check that the error is small
%relative to the stdev in the data.
condata.interpolationError = abs(condata.Rsheet - r(r0+1))./condata.Rsheet;

%% Define error bars

%reshape data into column vectors
Rhoriz = reshape(condata.Rhoriz,[size(condata.Rhoriz,2).*size(condata.Rhoriz,1),1]);
Rvert = reshape(condata.Rvert,[size(condata.Rvert,2).*size(condata.Rvert,1),1]);

%take standard deviation
condata.RhorizStd = std(Rhoriz);
condata.RvertStd = std(Rvert);

condata.RhorizErr = condata.RhorizStd/condata.RhorizMean;
condata.RvertErr = condata.RvertStd/condata.RvertMean;

condata.RsheetErr = sqrt(sumsqr([condata.RhorizErr condata.RvertErr]));

%define a new vdP function with the edge means + stdevs
% condata.f_err = @(Rs) exp((-pi()*(condata.RhorizMean+condata.RhorizStd))./Rs) + ...
%     exp((-pi()*(condata.RvertMean+condata.RvertStd))./Rs) - 1;

%calculate the relative error (normalized to Rsheet)
% condata.RsheetErr = (interp1(condata.f_err(r),r,0) - condata.Rsheet)./condata.Rsheet;

%add in error due to finite size of VdP contacts (= contact diameter / distance between contacts)
condata.RsheetErr = sqrt(sumsqr(condata.RsheetErr + 0.1));

if condata.RsheetErr > 10*condata.interpolationError
    condata.passedInterpCheck = true;
else
    condata.passedInterpCheck = false;
end

%% Determine if data passed all consistency checks

condata.passedAllChecks = condata.passedInterpCheck ...
    & condata.passedRecipCheck & condata.passedIrevCheck;

%% Normalize to thickness and define conductivity

%convert thickness from nm to cm
condata.thicknessCM = thickness(1).*(10^-7);

condata.conductivity = 1./condata.Rsheet./condata.thicknessCM;

condata.conductivityRelErr = sqrt(sumsqr([condata.RsheetErr, condata.thicknessErr]));
condata.conductivityAbsErr = condata.conductivity.*condata.conductivityRelErr;


end

