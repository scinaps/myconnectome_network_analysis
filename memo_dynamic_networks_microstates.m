
scans=84;
slides=101;
rois=630;
no=(rois*(rois-1))/2;
graph=zeros(scans*slides,no);

 for cc=1:coefs
     counter=0;
     for s=1:length(vector)
     if vector(s) < 100
        s1=strcat('sub0',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    else
        s1=strcat('sub',num2str(vector(s)),'_','coeff',num2str(cc),'_dynamic.mat') ;
    end

load(s1); % dfcgs 
 
for sl=1:slides
    counter=counter+1;
count=0;
 for k=1:rois
     for l=(k+1):rois
         count=count+1;
         graph(counter,count)=dfcg(sl+1,k,l); 
     end
 end
end

     end
     
        s1=strcat('dfcg_',coeff',num2str(cc),'.mat') ;
      save(s1,'graph')
    
 end
        
   