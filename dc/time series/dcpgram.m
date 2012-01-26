% Plots labelled periodogram
%
% function [] = dcpgram(data,dt,numpeaks, factor)
%       data - time series
%       dt - sampling delta_time (default = 1)
%       numpeaks - top 'x' peaks to be marked
%       factor - value to scale labelled peak TIME PERIOD by
% Returns freq,pow,cn (from dcfft) and indexes of peaks in original data set
% Calls beautify.m at the end

function [freq,pow,cn,ind] = dcpgram(data,dt,mark_peaks,factor)

    if ~exist('dt','var'), dt = 1; end;

    data = squeeze(data);
    data = data-mean(data);
    
    [freq,pow,cn] = dcfft(data,dt); 

    pgram = abs(cn).^2;

    figure
    loglog(freq,pgram,'LineWidth',1.5)
    title(['Periodogram']);
    grid on
    
    if exist('mark_peaks','var') & mark_peaks ~= 0
        [pks,ind] = findpeaks(pgram,'SORTSTR','descend');
        loglinex(freq(ind(1:mark_peaks)),factor);
        ind = ind(1:mark_peaks);
    else
        ind = NaN;
    end
    
    beautify;