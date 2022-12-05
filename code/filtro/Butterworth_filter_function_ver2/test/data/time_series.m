% Octave script to display the time series
close all;
clear all;
data = load('test_butterworth.txt');
plot(data(:,1),data(:,2));
hold on;
plot(data(:,1),data(:,3),'r');
legend('original','filtered');
