% Preparation
clc; clear; close;
% ======================================================================= %
% The first step focuses on identifying clusters of seismicity
% using kernel density estimation (KDE). This technique ensures that
% the study domain reflects natural groupings of earthquakes rather than
% arbitrary boundaries. For further analysis, this step is pivotal for
% identifying a set of events that are likely to be associated with the
% same seismogenesis. Studying the Gutenberg-Richter distribution within
% coherent seismicity regions, rather than in aggregate, is crucial to
% uncover localized geophysical variations and distinct earthquake types
% like dragon-kings, which broad-scale analysis obscures.
% ======================================================================= %

% Path and functions
PATH = '.\Seis_DK\';
addpath(fullfile(PATH,'Functions/'));

% Study region border
lim = [116,124,38,42];

% Read data
CAT = readmatrix([PATH,'Data\','catalog.xlsx']);

% Defining the study region by using kernel density estimation (KDE)
[REGION] = SpatialClusteringDetection( CAT, lim, PATH );

rmpath(fullfile(PATH,'Functions/'));