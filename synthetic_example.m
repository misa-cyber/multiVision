%% Example
% generate simulates sphere image
clear all; close all; clc;
sig=5;
addpath("synthetic","graphTools","graphTools\bezier");
[I,center,rad,ang, Xe]=generateSphImage(sig);
imshow(I);