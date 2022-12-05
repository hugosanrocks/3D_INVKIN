clear all
close all
clc

nodes = load('../dat/nodesc.in');
nodes = nodes';

file = fopen('nodes.bin','w');
fwrite(file,nodes,'single');


