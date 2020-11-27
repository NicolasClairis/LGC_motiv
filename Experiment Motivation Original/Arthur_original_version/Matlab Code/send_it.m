clear all
close all
clc
instrreset
% I would be sending info through instance A (integers with varying sizes each time, but only one integer at a time)
   u = udp('127.0.0.1', 15010, 'LocalPort', 2012);


    fopen(u);
   
    data(1:60*3) = 0;
    data(1+60*3:60*4.5) = linspace(0,50,60*1.5);
    data(1+60*4.5:60*5) = linspace(50,100,60*0.5);
    data(1+60*5:60*5) = 100;
    data(1+60*5:60*10) = flip(linspace(50,100,60*5));
    rep_data = repmat(data,100);
    rep_data = rep_data(1,:);
    
    
    
    for i = 1 :50000
     fwrite(u,rep_data(i),'double');
     data;
     fprintf('sent it\n');
     i
     pause(1/60)
    end
     pause(5)
    
        data = 1.5;
        tic
%     for i = 1 :100
%      fwrite(u,data,'double');
%      data;
%      fprintf('sent it\n');
%      data = data + 1;
%      i
%      if i ~= 1
%     TocVal(i) = toc-TocVal(i-1)
%      else
%          TocVal(i) = toc 
%      end
%       pause(1/60)
%      
%     end
    
    
     delete(u)
     fclose('all')
