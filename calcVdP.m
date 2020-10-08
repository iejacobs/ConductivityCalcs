function [ condata ] = calcVdP( folder, thickness, tolerance, filelist, findCorrectOrder, silentMode)
%CALCVDP calculates VdP conductivity assuming isotropic square sample
%   Outputs data as a struct containing raw data and calculated
%   conductivity, along with results of internal consistency checks and I-V
%   data. To obtain conductivty results from the condata struct, use dot
%   notation: i.e. 
%       condata.conductivity
%   A full list of fields in condata are given below. 
%   
%   ARGUMENTS
%
%   FOLDER:     folder containing the raw data files.
%
%   THICKNESS:  film thickness in nm. Absolute thickness error can be
%   specified in the thickness argument in the form 
%   [thickness, thicknessErr]
%
%   TOLERANCE:  the tolerance for internal consistency checks (hysteresis,
%   current reversal, and reciprocity). Default value is 0.03 = 3%. If
%   internal consistency checks give errors which exceed this tolerance, a
%   warning will be shown and the corresponding check flag, e.g.
%   condata.passedHysteresisCheck, will be set to 0.

%   FILELIST:  the filenames for the four I-V sweeps. these should be
%   symmetric hysteresis sweeps (sweeping from -V to +V, then +V to -V)
%   for Vds applied along each edge of the sample. The files should be
%   ordered as: [horiz. vert. horz. vert.]
%   If no extension is given '.txt' is tried automatically. 
%   If unspecified, calcVdP will get the list of files in the folder.
%   
%   FINDCORRECTORDER:   calcVdP will try each possible ordering of I-V
%   sweeps to find which one minimizes the reciprocity error--that is, the
%   ordering in which I-V sweeps of corresponding edges appear to match. 
%   If a file list is specified in the FILELIST argument, but does give 
%   the lowest reciprocity error, a warning is displayed.
%
%   SILENTMODE: suppresses warning messages
%
%   %OUTPUT
%   
%   CONDATA: a struct containing the results.
%   List of fields:
%
%     dataChecks: struct containing I-V data used in consistency checks
%     dataExclusion: struct containing information on excluded data points
%     dataRaw: raw data as imported
%     dataResistanceCalcs: resistance data used in calculations
%
%     errorBadRowsRatio: ratio of excluded data rows to total data rows
%     errorHysteresis: inconsistency in data between hysteresis sweeps
%     errorInterpolation: error introduced by interpolating solution to VdP equation
%     errorIrev: inconsistency in data upon current reversal
%     errorRecip: average inconsistency in resistance of symmetric edges
%     errorRecipHorz: inconsistency in resistance of horizontal edges
%     errorRecipVert: inconsistency in resistance of vertical edges
%
%     passedAllChecks: flag for overall pass/fail of all consistency checks
%
%     passedDataIntegrityCheck: pass/fail of data integrity check. failure
%     means more than 20% of datapoints were excluded from calculation.
%
%     passedHysteresisCheck: pass/fail of hysteresis check. failure
%     indicates errorHysteresis > TOLERANCE
%
%     passedInterpCheck: pass/fail of interpolation check. failure
%     indicates the interpolation used in solving the vdp equation is
%     introducing a measurable error, i.e. errorInterpolation > 1% of the
%     sheet resistance error
%
%     passedIrevCheck: pass/fail of current reversal check. failure
%     indicates errorIrev > TOLERANCE
%
%     passedRecipCheck: pass/fail of reciprocity check. failure
%     indicates errorRecip > TOLERANCE
%
%     thickness: film thickness from THICKNESS argument, in nm
%     thicknessAbsErr: absolute film thickness error, either as specified
%     in THICKNESS argument or assumed as 10%
%     thicknessErr: relative error in film thickness
%
%     conductivity: calculated VdP conductivity, in S/cm
%     conductivityRelErr: conductivity relative error
%     conductivityAbsErr: conductivity absolute error, in S/cm
%
%
%
%   v2.0 Ian Jacobs, Oct 2019

%% Input parsing
%defaults to non-silent mode, which displays error messages
if nargin < 6
    silentMode = false;
end

if nargin < 5
    findCorrectOrder = true;
end

%default tolerance is 3%
if nargin < 3
    tolerance = .03;
end
condata.dataChecks.tolerance = tolerance;

%if filenames are not given, we obtain them from the folder
if nargin < 4
    filelist = getFileList(folder);
end  

%if there are not exactly 4 files (or the number of files specified is not
%4) then throw an exception

if length(filelist) ~= 4
    fileNumException = MException('CalcVDP:wrongNumFiles','Incorrect number of files given');
    throw(fileNumException)
end

%% Find which files correspond to matching edges

if findCorrectOrder
    %if we don't know which files correspond to which edges, we need to
    %check each combination to see which gives the smallest error
    %(which is likely the correct combo). We need files in the order 
    %[1 2 1 2] so we need to perform 3 calcs, assuming an initial
    %ordering:
    fileorder = [1 2 3 4; 2 1 3 4; 3 2 1 4];
  
    %loop through each combination and run calcVdP, find reciprocity error
    %in each case
    for i=1:3
        tempdata = calcVdP(folder,thickness,tolerance,filelist(fileorder(i,:)),false,true);
        tempErr(i) = tempdata.errorRecip;
    end
    
    %combination with smallest reciprocity error is correct
    [~,index] = min(tempErr);
    
    %rerun with correct file order in non-silent mode (so error messages
    %aren't suppressed) and then return
    condata = calcVdP(folder,thickness,tolerance,filelist(fileorder(index,:)),false,false);
    
    %if the filelist was specified in incorrect order, display warning
    %indicating that order was corrected
    if nargin > 3 && index ~= 1 && ~silentMode
        warning(strcat("Data file order appears to be incorrect. Specified order gives reciprocity error of ", ...
            num2str(tempErr(1)), " while file order ", num2str(fileorder(index)), ...
            " gives reciprocity error ", num2str(tempErr(index)), ...
            ". If you are sure your file order is correct, re-run with findCorrectOrder set to false: calcVdP(folder,thickness,tolerance,filelist,false)"))
    end
    return
end


%% Data import
%Set filelist to whatever the four files you save are named, or rename
%your files to fit the default format. The first and fourth files are the
%horizontal resistance measurements, while the second and third files are
%the vertical resistance measurements.

for i=1:1:4
    try
        filepath = strcat(folder,filelist(i));
        [condata.dataRaw.Vdrain(:,i),condata.dataRaw.Idrain(:,i),condata.dataRaw.Vp1(:,i),condata.dataRaw.Vp2(:,i)] = ...
        importCond(filepath);
        condata.dataRaw.filepath = filepath;
    catch
        try
            filepath = strcat(folder,'/',filelist(i));
            [condata.dataRaw.Vdrain(:,i),condata.dataRaw.Idrain(:,i),condata.dataRaw.Vp1(:,i),condata.dataRaw.Vp2(:,i)] = ...
            importCond(filepath);
            condata.dataRaw.filepath = filepath;
        catch
            try
                filepath = strcat(folder,'/',filelist(i),'.txt');
                [condata.dataRaw.Vdrain(:,i),condata.dataRaw.Idrain(:,i),condata.dataRaw.Vp1(:,i),condata.dataRaw.Vp2(:,i)] = ...
                importCond(filepath);
                condata.dataRaw.filepath = filepath;
            catch
                filepath = strcat(folder,filelist(i),'.txt');
                [condata.dataRaw.Vdrain(:,i),condata.dataRaw.Idrain(:,i),condata.dataRaw.Vp1(:,i),condata.dataRaw.Vp2(:,i)] = ...
                importCond(filepath);
                condata.dataRaw.filepath = filepath;
            end
        end
    end

    condata.dataRaw.Rfull(:,i) = abs(condata.dataRaw.Vp1(:,i)-condata.dataRaw.Vp2(:,i))./abs(condata.dataRaw.Idrain(:,i));
end


%% Check that the number of rows is even, which must be true if the data is a hysteresis scan

numRows = size(condata.dataRaw.Rfull(:,1),1);

if mod(numRows,2) == 1
    oddDataException = MException('CalcVDP:oddData', ...
        'Dataset length must be even');
    throw(oddDataException)
end

condata.dataResistanceCalcs.Rsample = condata.dataRaw.Rfull;


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
    thicknessException = MException('CalcVDP:invalidThickness', ...
        'Thickness argument must either be of the form [thickness] or [thickness, thicknessErr]');
    throw(thicknessException)
end

%Define relative thickness error
condata.thicknessErr = condata.thicknessAbsErr./condata.thickness;

%% Hysteresis check
%Shows percent error in calculated resistance upon switching sweep direction.

%Difference in resistance at each Ids when switching sweep direction
condata.dataChecks.checkHysteresis = abs(condata.dataResistanceCalcs.Rsample - flipud(condata.dataResistanceCalcs.Rsample))./condata.dataResistanceCalcs.Rsample;

%Identifies data points which have error greater than the tolerance
% (typically values near Vsd = 0) which will be removed from the dataset
failindHyst = condata.dataChecks.checkHysteresis > tolerance;
[condata.dataExclusion.failrowsHyst, ~] = find(failindHyst);

%% Current reversal consistency check
%Shows percent error in calculated resistance upon switching current
% direction. The data is two I-V sweeps in opposite directions, so we first
% must split the data in half, then flip the order of each of these halves 
% and merge them.

%index of last point in first sweep. Will always be integer since we check
%above to ensure #points is even.
mid = numRows./2; 
%Resistance values after flipping each half of the data, giving the current
%reversed resistances.
RsampleIrev = vertcat(flipud(condata.dataResistanceCalcs.Rsample(1:mid,:)),flipud(condata.dataResistanceCalcs.Rsample((mid+1):end,:)));

%compare current-reversed values and calculate
condata.dataChecks.checkIrev = abs(condata.dataResistanceCalcs.Rsample - RsampleIrev)./condata.dataResistanceCalcs.Rsample;

%identifies data points which have error greater than the tolerance
% (typically values near Vsd = 0) which will be removed from the dataset
failindIrev = condata.dataChecks.checkIrev > tolerance;
[condata.dataExclusion.failrowsIrev, ~] = find(failindIrev);

%we will remove these bad rows below...

%% Reciprocity check

%check that vertical and horizontal edges pass reciprocity checks, that is,
%that the resistances for the two horizontal and vertical edges match
condata.dataChecks.checkRecipHorz = abs(1 - (condata.dataResistanceCalcs.Rsample(:,1) + RsampleIrev(:,1)) ...
    ./ (condata.dataResistanceCalcs.Rsample(:,3) + RsampleIrev(:,3)) );
condata.dataChecks.checkRecipVert = abs(1 - (condata.dataResistanceCalcs.Rsample(:,2) + RsampleIrev(:,2)) ...
    ./ (condata.dataResistanceCalcs.Rsample(:,4) + RsampleIrev(:,4)) );

%identifies data points which have error greater than the tolerance
% which will be removed from the dataset
failindRecipHorz = condata.dataChecks.checkRecipHorz > tolerance;
failindRecipVert = condata.dataChecks.checkRecipVert > tolerance;
[dataExclusion.failrowsRecipHorz, ~] = find(failindRecipHorz);
[dataExclusion.failrowsRecipVert, ~] = find(failindRecipVert);
condata.dataExclusion.failrowsRecip = unique(vertcat(dataExclusion.failrowsRecipHorz,dataExclusion.failrowsRecipVert));

%% Removal of bad datapoints

%we don't want to use data at Vds = 0 since the noise here will be very
%large (since current is very small), so we will find these and add them to
%the rows to be removed
[~,zrow1] = min(abs(condata.dataRaw.Idrain(1:mid,:)));
[~,zrow2] = min(abs(condata.dataRaw.Idrain((mid+1):end,:)));
condata.dataExclusion.zerorows = unique(horzcat(zrow1,zrow2+mid))';

%find unique rows
condata.dataExclusion.failrows = unique(vertcat(condata.dataExclusion.failrowsHyst, ...
    condata.dataExclusion.failrowsIrev,condata.dataExclusion.failrowsRecip,condata.dataExclusion.zerorows));
condata.errorBadRowsRatio = length(condata.dataExclusion.failrows)./numRows;
condata.dataExclusion.goodrows = 1:numRows;
condata.dataExclusion.goodrows(condata.dataExclusion.failrows) = [];

%if more than 2 or 20% of rows (whichever is greater) fail, the dataset doesn't pass
if length(condata.dataExclusion.failrows) > max([2, 0.2*numRows])
    condata.passedDataIntegrityCheck = false;
    if ~silentMode
        warning(strcat("Sample ",folder," failed data integrity check. ",...
            num2str(length(condata.dataExclusion.failrows))," data points failed to pass hysteresis, current reversal, or reciprocity checks with the current tolerance setting of ", ...
            num2str(tolerance)))
    end
else
    condata.passedDataIntegrityCheck = true;
    
    %only remove failed rows if passed the check, otherwise we may remove
    %all the rows and then get errors later
    condata.dataResistanceCalcs.Rsample(condata.dataExclusion.failrows,:) = [];
end




%% Re-check tolerances with bad datapoints removed

%calculate the error for each check using trimmed data, so they don't skew 
%our results, unless we didn't pass the data integrity check in which case
%we use the full data
if condata.passedDataIntegrityCheck
    r = condata.dataExclusion.goodrows;
else
    r = 1:numRows;
end
    
condata.errorHysteresis = mean(condata.dataChecks.checkHysteresis(r));
condata.errorIrev = mean(condata.dataChecks.checkIrev(r));
condata.errorRecipHorz = mean(condata.dataChecks.checkRecipHorz(r));
condata.errorRecipVert = mean(condata.dataChecks.checkRecipVert(r));
condata.errorRecip = mean([condata.errorRecipHorz condata.errorRecipVert]);

%if hysteresis error is greater than tolerance, check fails
if condata.errorHysteresis > tolerance
    condata.passedHysteresisCheck = false;
    if ~silentMode
        warning(strcat("Sample ",folder," failed hysteresis check. Hysteresis error is ",...
            num2str(condata.errorHysteresis),"; tolerance setting is currently ", ...
            num2str(tolerance)))
    end
else
    condata.passedHysteresisCheck = true;
end


%if Irev error is greater than tolerance, check fails
if condata.errorIrev > tolerance
    condata.passedIrevCheck = false;
    if ~silentMode
        warning(strcat("Sample ",folder," failed current reversal check. Current reversal error is ",...
            num2str(condata.errorIrev),"; tolerance setting is currently ", ...
            num2str(tolerance)))
    end
else
    condata.passedIrevCheck = true;
end

%if reciproicity error is greater than tolerance, check fails
if condata.errorRecip > tolerance
    condata.passedRecipCheck = false;
    if ~silentMode
        warning(strcat("Sample ",folder," failed reciprocity check. Reciprocity error is ",...
            num2str(condata.errorRecip),"; tolerance setting is currently ", ...
            num2str(tolerance)))
    end
else
    condata.passedRecipCheck = true;
end

%% Define VdP equation and guess its solution

%calculate average resistance on each edge
condata.dataResistanceCalcs.Redgemean = mean(condata.dataResistanceCalcs.Rsample);

%define horizontal and vertical resistivities
condata.dataResistanceCalcs.RhorizMean = mean(condata.dataResistanceCalcs.Redgemean([1 3]));
condata.dataResistanceCalcs.RvertMean = mean(condata.dataResistanceCalcs.Redgemean([2 4]));

%also define the full horizontal and vertical resistance data (for error
%calcs later)
condata.dataResistanceCalcs.Rhoriz = condata.dataResistanceCalcs.Rsample(:,[1 3]);
condata.dataResistanceCalcs.Rvert = condata.dataResistanceCalcs.Rsample(:,[2 4]);

%define van der pauw equation
f = @(Rs) exp((-pi()*condata.dataResistanceCalcs.RhorizMean)./Rs) + exp((-pi()*condata.dataResistanceCalcs.RvertMean)./Rs) - 1;

%for a square geometry the solution to the VdP equation is Rsheet =
%Rvert*pi/log(2) ~= Rvert*4.53
Rguess = mean([condata.dataResistanceCalcs.RhorizMean condata.dataResistanceCalcs.RvertMean]).*pi/log(2);

%% Solve VdP equation by interpolation (fast)

%Sets how accurately the data is interpolated. This is probably way more
%accurate that necessary but still runs fast
npoints = 1E7;
scanrange = 10;

%define vector of reasonable resistance values
r = linspace((Rguess./scanrange),(Rguess.*scanrange),npoints);

%find value of r that sets VdP equation to 0 by interpolation
condata.dataResistanceCalcs.Rsheet = interp1(f(r),r,0);

%find closest value in r vector
[~, r0] = min(abs(r - condata.dataResistanceCalcs.Rsheet));

%estimate of the error introduced by interpolation. This is a significant
%overestimation, since as long as the rate of change in f is slow with R
%then interpolation will still be accurate even if the r vector isn't
%granular enough, but this is only used to check that the error is small
%relative to the stdev in the data.
condata.errorInterpolation = abs(condata.dataResistanceCalcs.Rsheet - r(r0+1))./condata.dataResistanceCalcs.Rsheet;

%% Define error bars

%reshape data into column vectors
Rhoriz = reshape(condata.dataResistanceCalcs.Rhoriz,[size(condata.dataResistanceCalcs.Rhoriz,2).*size(condata.dataResistanceCalcs.Rhoriz,1),1]);
Rvert = reshape(condata.dataResistanceCalcs.Rvert,[size(condata.dataResistanceCalcs.Rvert,2).*size(condata.dataResistanceCalcs.Rvert,1),1]);

%take standard deviation
condata.dataResistanceCalcs.RhorizStd = std(Rhoriz);
condata.dataResistanceCalcs.RvertStd = std(Rvert);

condata.dataResistanceCalcs.RhorizErr = condata.dataResistanceCalcs.RhorizStd/condata.dataResistanceCalcs.RhorizMean;
condata.dataResistanceCalcs.RvertErr = condata.dataResistanceCalcs.RvertStd/condata.dataResistanceCalcs.RvertMean;

condata.dataResistanceCalcs.RsheetErr = sqrt( sum( [condata.dataResistanceCalcs.RhorizErr condata.dataResistanceCalcs.RvertErr].^2 ) );

%define a new vdP function with the edge means + stdevs
% condata.f_err = @(Rs) exp((-pi()*(condata.dataResistanceCalcs.RhorizMean+condata.dataResistanceCalcs.RhorizStd))./Rs) + ...
%     exp((-pi()*(condata.dataResistanceCalcs.RvertMean+condata.dataResistanceCalcs.RvertStd))./Rs) - 1;

%calculate the relative error (normalized to Rsheet)
% condata.dataResistanceCalcs.RsheetErr = (interp1(condata.f_err(r),r,0) - condata.dataResistanceCalcs.Rsheet)./condata.dataResistanceCalcs.Rsheet;

%add in error due to finite size of VdP contacts (= contact diameter / distance between contacts)
%here we assume a square geometry with contact size (d) = 1mm and edge size
%(D) = 10mm. This gives an error on the order 0.1% = 0.001;
%See Rev. Sci. lnstrum. 1989, 60, 271; Rev. Sci. lnstrum. 1989, 60, 275.
condata.dataResistanceCalcs.RsheetErr = sqrt( sum( [condata.dataResistanceCalcs.RsheetErr 0.001].^2 ) );

%if more than 1% of the total error is due to interpolation error, check
%fails.
if 0.01.*condata.dataResistanceCalcs.RsheetErr > condata.errorInterpolation
    condata.passedInterpCheck = true;
else
    condata.passedInterpCheck = false;
    if ~silentMode
    warning(strcat("Sample ",folder," failed interpolation check. Interpolation error is ",...
        num2str(condata.errorInterpolation),"; which must be less than 10% of the sheet resistance relative error, ", ...
        num2str(condata.dataResistanceCalcs.RsheetErr)))
    end
end

%% Determine if data passed all consistency checks

condata.passedAllChecks = condata.passedDataIntegrityCheck && ...
    condata.passedHysteresisCheck && condata.passedInterpCheck ...
    && condata.passedRecipCheck && condata.passedIrevCheck;

%give warning if sample fails checks
if ~condata.passedAllChecks && ~silentMode
    warning(strcat("Sample ",folder," failed at least one internal consistency check."))
end

%% Normalize to thickness and define conductivity

%convert thickness from nm to cm
thicknessCM = thickness(1).*(10^-7);

condata = orderfields(condata);

condata.conductivity = 1./condata.dataResistanceCalcs.Rsheet./thicknessCM;

condata.conductivityRelErr = sqrt(sum([condata.dataResistanceCalcs.RsheetErr, condata.thicknessErr].^2));
condata.conductivityAbsErr = condata.conductivity.*condata.conductivityRelErr;

end

