clear all
close all
clc

option = input('1D (1) or 2D (2) inversion results: ')

if ( option == 1 )
  sliprate1d
elseif ( option == 2 )
  sliprate2d
end

