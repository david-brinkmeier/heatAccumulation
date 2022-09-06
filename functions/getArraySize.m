function [sz_recommended,size_in_Gb] = getArraySize(sz,bins,precision)
% This function calculates an allowable maximum array size for 2D arrays
% of size(2^a,2^b) given the size of a first dimension (sz) and the 
% available system RAM (on Windows systems). 
% For Unix/OSX a fixed size limit is used.
% 64 bit float, 8 bit = 1 byte, Gigabyte = 1e-9 byte
%
% sz: edge length of 2D array (e.g. timesteps)
% bins: required number of bins (e.g. orthogonal direction of array)
% due to memory constraints it may not be possible to process an array of
% size (sz,bins); this function determines a maximum recommended size for
% the second array dimension
%
% the intended computation then needs to occur in chunks of size
% (sz,sz_recommended)

% effective size is quadrupled when convolving the two arrays + OVERHEAD
sz_effective = sz*4; % 2 arrays zpad to double dimension times 2 for imaginary part

switch precision
    case 'fp64'
        scale = 2;
    case 'fp32'
        scale = 1;
    otherwise
        warning('Function %s.m cannot handle provided precision type. Defaulting to fp64.',mfilename)
        scale = 2;
end

if ispc % windows host  
    userview = memory; % get system RAM info
%   in case of fast convolution array sizes are *2*2 (2 x fp64 = fp128, 2x signal and kernel, 2x zpadding)
    sizelimit = 0.7*userview.MaxPossibleArrayBytes; % limit array to 70% of free RAM (accounting for potential FFT overhead)
else % OSX or Unix host
    sizelimit = 0.5e9; % should be OK for 16gb systems
end

% maximum array length in 2nd dimension to fit within memory constraint
% err on side of caution
sz_recommended = 2^floor(log2(sizelimit*(8/64)/(sz_effective*scale)));

if sz_recommended >= bins
    sz_recommended = bins;
elseif sz_recommended < 1
    sz_recommended = 1; % this can only happen for >> 2^30 pulses...
end

size_in_Gb = sz_recommended*sz_effective*(64/8)*1e-9;

end

