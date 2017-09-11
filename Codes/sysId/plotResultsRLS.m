clear;
clc;
close all;


addpath(['.' filesep 'results']);



load RLS.mat;

plot(10*log10(e3));
hold on;


load SRLS.mat;

plot(10*log10(e3));

load DSRLS.mat;

plot(10*log10(e3));

load SOBE.mat;

plot(10*log10(e3));


figure
load RLS.mat;

plot(10*log10(misalignment));
hold on;


load SRLS.mat;

plot(10*log10(misalignment));

load DSRLS.mat;

plot(10*log10(misalignment));

load SOBE.mat;

plot(10*log10(misalignment));


