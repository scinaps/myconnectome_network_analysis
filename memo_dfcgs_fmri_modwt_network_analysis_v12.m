files=dir('\C:\Users\Ana\Documents\MATLAB\timeseries');

% data=load(files(3).name);
% file=files(3).name;
% [samples rois]=size(data);

window=15;
step=5;
samples=518;
slides=round((samples - window)/step);

coefs=10;

vector=14:104;
vector([31 33 34 52 55 90 93 ]-13)=[ ];
rois=630;


for s=1:length(vector)
    s
   % data=load(files(3).names);

 if s ~= vector(s)  
    if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'.txt') ;
    else
        s1=strcat('sub',num2str(vector(s)),'.txt') ;
    end
 end
    data=load(s1);
    [samples rois ]=size(data);
    wcoef=zeros(rois,coefs,samples);
    
    for r=1:rois
        wcoef(r,:,:)=modwt(data(:,r)','db4');
    end
    
    %%%%%%%%%%%% estimate static graph
   dfcgs=zeros(slides,rois,rois);
   
 
    for cc=1:coefs
        
     D=squareform(pdist(squeeze(wcoef(:,cc,:)),'correlation'));
     dfcgs(1,:,:)=D;
   
    t1=-step + 1;
    t2=window - step;
    
    for sl=1:slides
        [s sl]
        t1=t1+step;
        t2=t2+step;
        
        D=squareform(pdist(squeeze(wcoef(:,cc,t1:t2)),'correlation'));
        dfcgs(sl+1,:,:)=D;
    end
    
    if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    end
    save(s1 ,'dfcgs')
    end
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CORRELATE STATIC NETWORKS WITH DYNAMIC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NETWORKS
 
 corr_stat_dyn=zeros(length(vector),coefs,slides);
 
 for s=1:length(vector)
     
  for cc=1:coefs 
    if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    end
 load(s1);
   for sl=1:slides
       [s cc sl]
           cor=corrcoef(squeeze(dfcgs(1,:,:)),squeeze(dfcgs(sl+1,:,:)));
           corr_stat_dyn(s,cc,sl)=cor(1,2);
   end
     
     if vector(s) < 100       
         s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'_','corr',num2str(cc),'_dynamic.mat') ;
      else
 s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'_','corr',num2str(cc),'_dynamic.mat') ;
 end
     
 end
 end
 
 save 'corr_static_dynamic_graphs.mat' corr_stat_dyn

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% NETWORK ANALYSIS %%%%%%%%
%%%%%%%%%% STATIC NETWORK ANALYSIS

%%%%%% BASIC NETWORK METRICS
static_net=zeros(length(vector),coefs,rois,3);

%%%%%%%%%%%%%% CLUSTERING FOR STATIC NETWORKS
clustering=zeros(length(vector),coefs,rois);
quality=zeros(length(vector),coefs);  

sl=1;

for s=1:length(vector)
    
     if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(10),'_static.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(10),'_static.mat') ;
    end

load(s1);
    for cc=1:coefs
        
   
             [s cc]
             graph=squeeze(static(cc,:,:));
             %NORMALIZE THE GRAPH
             max1=max(abs(graph(:)));
             ngraph=graph./max1;
             
            %%%%%%%%%%%% GLOBAL EFFICIENCY
%             [gl_node gl distance]=global_efficiency_wu(1./ngraph);
%             static_net(s,cc,:,1)=gl_node;
            
            %%%%%%%%%%%% BETWEENNESS CENTRALITY
%             BC=betweenness_wei(ngraph);
%             BC1=BC./sum(BC);
%             static_net(s,cc,:,2)=BC1;
            
            %%%%%%%%%%%% STRENGTH
            str=sum(ngraph);
            static_net(s,cc,:,3)=str;
            
            %%%%%%%%% MODULARITY
%             [cl q]=modularity_und(ngraph);
%             clustering(s,cc,:)=cl;
%             quality(s,cc)=q;

    end
end

  %%%%%save the matrices
    save 'static_network_analysis.mat' static_net
%     save 'static_modularity.mat' clusteringd
%     

     %%%%%%%%%% ESTIMATE COEFFICIENT OF VARIATION
     cv=zeros(coefs,rois,3);
       load('static_network_analysis.mat'); 
     load('static_modularity.mat');
     load('static_modulaity_quality'); 
    
        for cc=1:coefs
            for area=1:rois
                data=squeeze(static_net(:,cc,area,1));
                mean1=mean(data);
                std1=std(data);
                cv(cc,area,1)=mean1/std1;
                
                data=squeeze(static_net(:,cc,area,2));
                mean1=mean(data);
                std1=std(data);
                cv(cc,area,2)=mean1/std1;
                
                data=squeeze(static_net(:,cc,area,3));
                mean1=mean(data);
                std1=std(data);
                cv(cc,area,3)=mean1/std1;
                
            end
        end
    
        save 'coefficient_variation_static.mat' cv
        
  
   load('coefficient_variation_static.mat')
   figure(2) ; plot(cv(2,:,3))
  
   load('static_modularity.mat'); %clustering
   load('static_modulaity_quality'); %quality
   
   load('static_network_analysis.mat');
   
   mean1=squeeze(mean(static_net,1));
   net=squeeze(mean1(1,:,1));
   figure ; plot(mean(quality)/std(quality))
   
   %%%%%%%CV metrics extraction - for each of the 3 measures, 10x coeffs%%%%%%%%%%%%%
   cv1=squeeze(cv(8,:,3))
   cv_c= transpose (cv1)
   cvtable = [cv_c];
   
  % POWER SPECTRUM OF CORRELATION TIME SERIES BETWEEN STATIC AND DYNAMIC
  % GRAPHS
  load('corr_static_dynamic_graphs.mat'); % corr_stat_dyn
  
  [scans coefs slides]=size(corr_stat_dyn);
  
  x=squeeze(corr_stat_dyn(1,1,:));
  Fs=2.2;
  N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;


  ps=zeros(scans,coefs,length(freq));
  
  for sc=1:scans
      for cc=1:coefs
           x=squeeze(corr_stat_dyn(sc,cc,:));
  Fs=2.2;
  N = length(x);
xdft = fft(x);
xdft = xdft(1:N/2+1);
psdx = (1/(Fs*N)) * abs(xdft).^2;
psdx(2:end-1) = 2*psdx(2:end-1);
freq = 0:Fs/length(x):Fs/2;
ps(sc,cc,:)=freq;     
      end
  end
  
  
 figure(2)
 
for cc=1:coefs 
subplot(2,5,cc),plot(squeeze(mean(ps(:,cc,:),1)),10*log10(psdx))

grid on
title('Periodogram Using FFT')
xlabel('Frequency (Hz)')
ylabel('Power/Frequency (dB/Hz)')
end


   
 figure(3),title('Correlation Analysis between Static and Dynamic Graphs')
 tt=window:step:518;
 tt=(tt*Fs)/60;
 
for cc=1:coefs 
subplot(2,5,cc),plot(tt,squeeze(mean(corr_stat_dyn(:,cc,:),1)))

grid on

xlabel('Time')
ylabel('Corr')
end
   
    figure(4),title('Correlation Analysis between Static and Dynamic Graphs')
 tt=window:step:518;
 tt=(tt*Fs)/60;
 
for cc=1:coefs 
subplot(2,5,cc),plot(tt,squeeze(std(corr_stat_dyn(:,cc,:),1)))

grid on

xlabel('Time')
ylabel('Corr')
end


    figure(5),title('Correlation Analysis between Static and Dynamic Graphs')
 tt=window:step:518;
 tt=(tt*Fs)/60;
 
for cc=1:coefs 
subplot(2,5,cc),plot(tt,squeeze(corr_stat_dyn(1,cc,:)))

grid on

xlabel('Time')
ylabel('Corr')
end


   %%FREQUENCIES
   
speroctave = 1;
numoctaves = 10;
a0 = 2^(1/speroctave);
Fs = (19*60)/518; % min x secs /samples
scales = a0.^(speroctave:1/speroctave:numoctaves*speroctave)/numoctaves;


Frq = scal2frq(scales,'db4',1/Fs);

Frq = Frq(:);
scales = scales(:);
T = [scales.*(1/Fs) Frq 1./Frq];
T = array2table(T,'VariableNames',{'Scale','Frequency','Period'});

   
   %%%%%%%%% PLOT STATIC NETWORKS
   graph=zeros(length(vector),coefs,rois,rois);
   
for s=1:length(vector)
    s
   % data=load(files(3).names);

    if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(10),'_static.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(10),'_static.mat') ;
    end
    
    load(s1);
    graph(s,:,:,:)=static;
    
end
   
mean_graph=squeeze(mean(graph,1));
max1=max(mean_graph(:));
mean_graph=mean_graph./max1;

%text=importdata('russome_Node_ROI.node.txt');

text=xlsread('parcellation_custom_node.xlsx');

areas=cell(1,20);

for s=1:rois
    dat=text{s};
    area=dat(6)

figure ; 

for coef=1:10
    subplot(2,5,coef) ; imagesc(squeeze(mean_graph(coef,text(:,7),text(:,7))))
end

   
% %%%%%%%%%% LEAVE IT FOR LATER - DYNAMIC NETWORK ANALYSIS
% dynet=zeros(length(vector),coefs,slides+1,rois,3);
% 
% for s=1:length(vector)
%     for cc=1:coefs
%         
%     if vector(s) < 100
%         s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'.mat') ;
%     else
%         s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'.mat') ;
%     end
% 
% load(s1);
% 
%         for sl=1:slides+1
%             [s cc sl]
%             %%%%%%%%%%%% GLOBAL EFFICIENCY
%             [gl_node gl distance]=global_efficiency_wu(1./squeeze(dfcgs(sl,:,:)));
%             dynet(s,cc,sl,:,1)=gl_node;
%             
%             %%%%%%%%%%%% BETWEENNESS CENTRALITY
%             BC=betweenness_wei(squeeze(dfcgs(sl,:,:)));
%             BC1=BC./sum(BC);
%             dynet(s,cc,sl,:,2)=BC1;
%             
%             %%%%%%%%%%%% STRENGTH
%             str=sum(sum(squeeze(dfcgs(sl,:,:))));
%             dynet(s,cc,sl,:,3)=str;
% 
%         end
%         
% end

% save 'dynami_network_analysis.mat' dynet

%%%%%%%%%%%%%% ESTIMATE THE ICC FOR EACH OF THE NETWORK METRICS AND FOR
%%%%%%%%%%%%%% EACH ROIS

icc_net=zeros(coefs,3,rois);


    for cc=1:coefs
        
    if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'_','dynet',num2str(cc),'.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'_','dynet',num2str(cc),'.mat') ;
    end
      load(s1);
    
      for area=1:rois
         for net=1:3
            out = my_ICC(3,'k',squeeze(dynet(:,cc,:,area,net)));
            icc_net(cc,net,area)=out;
         end
      end
      
    end
    
    save 'icc_dfcg_net_metrics.mat' icc_net


load('icc_dfcg_net_metrics.mat')

        
        