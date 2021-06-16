%% PREPROCESSING CALLS
set(0,'defaulttextinterpreter','latex')                                   ;%
fslab  = 14; fsleg  = 14; fstit  = 14; fstick = 14;
ep         = true                                                         ;% elasto-plasticity is true
% video settings
video      = 1                                                            ;% off = 0, on = 1
vidname    = 'press_displ'                                                ;% video name
format     = 'avi'                                                        ;% video format
nout       = 20                                                           ;% plot every 20 it
cout       = 0                                                            ;% frame numbering
ps         = 3                                                            ;% marker size for scatterplot
% data type
typeD      = 'double'                                                     ;% arithmetic precision, single or double