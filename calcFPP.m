function [ condata ] = calcFPP( filepath, thickness, leng, width, varargin )
%calcFPP calculates four point probe conductivity for hall-bar devices
%   Output data is a struct containing raw data and calculated 
%   conductivity, along with whether sample passed consistency checks and 
%   some other information (raw data, etc.)
%
%   Thickness is given in nm, absolute thickness error can be specified in
%   the thickness argument in the form [thickness, thicknessErr]
%
%   Length (leng) and width are given in microns. Absolute dimensional
%   errors can be specified in the same was as above, i.e [leng, lengErr]
%
%   Conductivity is given in S/cm
%
%   
%
%   v1.0 Ian Jacobs, April 2019

%% Variable arguments

IRevTolerance = 0.01;

if ~isempty(varargin)
    IRevTolerance = varargin{1};
end
    


%% Data import

FE = MException("CalcFPP:invalidFilename", "Invalid file path. Try again");

try
    [condata.Vdrain,condata.Idrain,condata.Vp1,condata.Vp2] = ...
    importCond(filepath);
    condata.filename = filepath;
catch
    try
        filepath = strcat(filepath,'.txt');
        [condata.Vdrain,condata.Idrain,condata.Vp1,condata.Vp2] = ...
        importCond(filepath);
        condata.filename = filepath;
    catch
        throw(FE)
    end
    throw(FE)
end


%% Get dimension and dimension errors

% Thickness
%first component of thickness must always be thickness in NM
condata.thickness = thickness(1);
%convert thickness from nm to cm
condata.thicknessCM = condata.thickness.*(10^-7);

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
        'Thickness argument must either be of the form [thickness] or [thickness, thicknessErr]');
    throw ME
end

%Define relative thickness error
condata.thicknessErr = condata.thicknessAbsErr./condata.thickness;



% Length
%first component of length must always be length in MICRONS
condata.length = leng(1);
%convert thickness from um to cm
condata.lengthCM = condata.length.*(10^-4);

%use thicknessErr if given in function call
if length(leng) > 1
    %if length error is specified in 2nd component of thickness
    %argument, set thicknessErr to this value
    condata.lengthAbsErr = leng(2);
elseif length(leng) == 1
    %if no length error is given, assume the error is 1%
    condata.lengthAbsErr = condata.length.*0.01;
else
    ME = MException('CalcFPP:invalidLength', ...
        'Length argument must either be of the form [length] or [length, lengthErr]');
    throw ME
end

%Define relative length error
condata.lengthErr = condata.lengthAbsErr./condata.length;



%first component of width must always be width in MICRONS
condata.width = width(1);
%convert thickness from um to cm
condata.widthCM = condata.width.*(10^-4);

%use thicknessErr if given in function call
if length(width) > 1
    %if thickness error is specified in 2nd component of thickness
    %argument, set thicknessErr to this value
    condata.widthAbsErr = width(2);
elseif length(width) == 1
    %if no thickness error is given, assume the error is 10%
    condata.widthAbsErr = condata.width.*0.01;
else
    ME = MException('CalcVDP:invalidWidth', ...
        'Width argument must either be of the form [width] or [width, widthErr]');
    throw ME
end

%Define relative thickness error
condata.widthErr = condata.widthAbsErr./condata.width;


%Total dimensional error is RMS error of dimensional errors
condata.dimErr = sqrt(sum([condata.thicknessErr, condata.lengthErr, condata.widthErr].^2));


%% Calculate sample resistance

%Calculate resistance of each IV datapoint
condata.Rfull = abs(condata.Vp1-condata.Vp2)./abs(condata.Idrain);

%Current reversal consistency check and removal of bad datapoints
%Shows percent error in calculated resistance upon switching current direction
condata.checkIrev = abs(condata.Rfull - flipud(condata.Rfull))./condata.Rfull;

%identifies data points which have error greater than IRevTolarance
%(default 1%) and removes them from the dataset. These are typically values
%at/near Vsd = 0, which are dominated by noise and would skew the average
%resistance value.
failind = condata.checkIrev > IRevTolerance;
[failrows, ~] = find(failind);
failrows = unique(failrows);
condata.Rfull(failrows,:) = [];

%if more than 10 rows fail the check, something is wrong
if length(failrows) > 10
    condata.passedIrevCheck = false;
else
    condata.passedIrevCheck = true;
end


%Calculate average resistance of sample
condata.R = mean(condata.Rfull);

%Standard deviation of IV data
condata.RStd = std(condata.Rfull);

%Relative error of sample resistance
condata.RErr = condata.RStd./condata.R;


%% Determine if data passed all consistency checks

condata.passedAllChecks = condata.passedIrevCheck;

%% Normalize to thickness and define conductivity

condata.conductivity = 1./condata.R .* (condata.lengthCM./(condata.thicknessCM.*condata.widthCM));

condata.conductivityRelErr = sqrt(sum([condata.RErr, condata.dimErr].^2));
condata.conductivityAbsErr = condata.conductivity.*condata.conductivityRelErr;

%conductivity = condata.conductivity;
%errorAbs = condata.conductivityAbsErr;
end

