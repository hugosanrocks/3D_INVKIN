clear all
close all
clc


nsta=40;
n_sta_sub=15;
ncomp=3;

timesyn=load('out/syn_S001_C1.ascii');
timeobs=load('out/obs_S001_C1.a');
tsyn(:,1) = timesyn(:,1);
tsyn(:,2) = tsyn(:,1);
tsyn(:,3) = tsyn(:,1);

arrayused = [1 2 3 4 5 8 10 15 16 17 21 22 23 24 25 26 27];

[row,col] = size(timesyn);

tobs(:,1) = timeobs(:,1);
tobs(:,2) = tobs(:,1);
tobs(:,3) = tobs(:,1);

sort=[7.092     8
17.22   10
21.76   17
27.76   3
32.19   5
32.28   4
39.81   1
49.47   2
51.45   15
58.91   16
59.49   22
59.71   26
72.97   21
75.56   23
91.59   24
93.86   27
100.1   25
13.34   7
28.29   9
31.23   11
40.19   13
49.29   12
53.22   14
54.17   6
60.74   19
60.93   18
74.8    20
84.61   28
108.5   29
134.1   30];

sort2=[7.092	8
13.34	7
17.22	10
21.76	17
27.76	3
28.29	9
31.23	11
32.19	5
32.28	4
39.81	1
40.19	13
49.29	12
49.47	2
51.45	15
53.22	14
54.17	6
58.91	16
59.49	22
59.71	26
60.74	19
60.93	18
72.97	21
74.8	20
75.56	23
84.61	28
91.59	24
93.86	27
100.1	25
108.5	29
134.1	30];


%  message = sprintf('# PLOT_MAX = ',)
%  message = sprintf('# T_START = ,  0.0')
%  message = sprintf('# T_END = 65.0')
 

kk=1;
for k=1:nsta
  for i=1:ncomp
    %read synthetics
    file=sprintf('out/syn_S%03d_C%01d.ascii',k,i);
    syn2=load(file);
    syn(:,i)=syn2(:,2);
    %read observations
    fileobs=sprintf('out/obs_S%03d_C%01d.a',k,i);
    obs2=load(fileobs);
    obs(:,i)=obs2(:,2);
  end
  %normalization factor per station
  maxval=max(max(abs([syn(:,:), obs(:,:)])));
  syn = syn./maxval;
  obs = obs./maxval;
if (k < 21)
  syn = syn + k*(2);
  obs = obs + k*(2);
else
  syn = syn + kk*(2);
  obs = obs + kk*(2);
  kk = kk + 1;
end
  %write headers
  message1 = sprintf(' # PLOT_MAX = %8.2f \n',maxval);
  message2 = sprintf(' # T_START  = 0.0 \n');
  message3 = sprintf(' # T_END    = 35.0 \n');

  %fileoutsyn=fopen('seismo_to_plot.syn','a');
  for i=1:ncomp
    file=sprintf('seismo_to_plote_%02d.syn',k);
    fileoutsyne=fopen(file,'a');
    file=sprintf('seismo_to_plote_%02d.obs',k);
    fileoutobse=fopen(file,'a');
    file=sprintf('seismo_to_plotn_%02d.syn',k);
    fileoutsynn=fopen(file,'a');
    file=sprintf('seismo_to_plotn_%02d.obs',k);
    fileoutobsn=fopen(file,'a');
    file=sprintf('seismo_to_plotv_%02d.syn',k);
    fileoutsynv=fopen(file,'a');
    file=sprintf('seismo_to_plotv_%02d.obs',k);
    fileoutobsv=fopen(file,'a');
    datasyn = [tsyn(1:end,i) syn(1:end,i)];
    dataobs = [tobs(1:end,i) obs(1:end,i)];
    if (i == 1)
    fprintf(fileoutsyne,message1)
    fprintf(fileoutsyne,message2)
    fprintf(fileoutsyne,message3)
    fprintf(fileoutobse,message1)
    fprintf(fileoutobse,message2)
    fprintf(fileoutobse,message3)
    for j=1:row
     fprintf(fileoutsyne,'%14.6e %14.6e\n',[datasyn(j,1) datasyn(j,2)].');
     fprintf(fileoutobse,'%14.6e %14.6e\n',[dataobs(j,1) dataobs(j,2)].');
    end
    elseif (i == 2)
    fprintf(fileoutsynn,message1)
    fprintf(fileoutsynn,message2)
    fprintf(fileoutsynn,message3)
    fprintf(fileoutobsn,message1)
    fprintf(fileoutobsn,message2)
    fprintf(fileoutobsn,message3)
    for j=1:row
     fprintf(fileoutsynn,'%14.6e %14.6e\n',[datasyn(j,1) datasyn(j,2)].');
     fprintf(fileoutobsn,'%14.6e %14.6e\n',[dataobs(j,1) dataobs(j,2)].');
    end
    elseif (i == 3)
    fprintf(fileoutsynv,message1)
    fprintf(fileoutsynv,message2)
    fprintf(fileoutsynv,message3)
    fprintf(fileoutobsv,message1)
    fprintf(fileoutobsv,message2)
    fprintf(fileoutobsv,message3)
    for j=1:row
     fprintf(fileoutsynv,'%14.6e %14.6e\n',[datasyn(j,1) datasyn(j,2)].');
     fprintf(fileoutobsv,'%14.6e %14.6e\n',[dataobs(j,1) dataobs(j,2)].');
    end
    end
    fclose(fileoutsyne);
    fclose(fileoutobse);
    fclose(fileoutsynn);
    fclose(fileoutobsn);
    fclose(fileoutsynv);
    fclose(fileoutobsv);
    %plot(tsyn(:,i),syn(:,i),'.r'),hold on
    %plot(tobs(:,i),obs(:,i),'.b'),hold on
    %xlim([0,max(tobs(:,3))])
    %ylim([0,max(max(obs(:,:)))])
  end
end

%scale
s1 = [25 41.5; 35 41.5];
save('-ascii','graphics/scale1.txt','s1');
s2 = [25 43; 35 43];
save('-ascii','graphics/scale2.txt','s2');



