
% this script is used to identify individual subject across task conditions based on their brain activity patterns. 
% For each pair of tasks, an .nii map is produced showing for each voxel the group-mean accuracy of invidual identification. 
clear;
close all;

projpath="/Volumes/BTCruiser/GRFP";
cd (projpath); 
dest='.\results\';

addpath(genpath("/Users/bterhunecotter/MyDrive/SDSU_LLCN/NEURO/Repositories/NIFTI_toolbox")); 
mask=load_nii("PATHTOMASKHERE");
hdr=mask.hdr; % the hdr is used to convert an img matrix into a nii map.
clear mask  

% the sub-tmaps contain a set of 4D matrices - one matrix for each task (condition). 
% the first three dimensions of the matrix correspond to the t-map from the first-level analysis. 
% the fourth dimension corresponds to subject ID. 
load('.\data\sub_tmaps.mat'); % T maps seem to work better than beta maps. 

% this searchlight_list file contains the searchlights over the whole brain . 
% the spatial dimensions of the brain mask used to make the searchlights must be identical to the t-maps
load ('.\script\searchlight_list.mat');
ID_set=L.LI;
voxel=L.voxel;
sub_N=size(sign,4);
source={asl, eng};
taskname={'asl', 'eng'};

for  row=1:3
    tic
     for  cl=1:3
          if  row~=cl
              origin=source{row};
              s_name=taskname{row};
              predictor=source{cl};            
              p_name=taskname{cl};
              
%% predict subject identity across runs
              prd_id=nan(61*73*61,sub_N);
              prd_SI=nan(61*73*61,sub_N);
              for   i=1:sub_N
                    sub_origin=origin(:,:,:,i);

% calculate RS for each voxel of interest
                    for  k=1:length(ID_set)
                         ID=ID_set{k};
                         v_id=voxel(k);
                         t_sub_origin=nan(125,1);
                         t_predictor=nan(125,sub_N);
       
% get the t value within serachlight ROI
                         for   n=1:size(ID,2)
                               a=ID(:,n); 
                               t_sub_origin(n)=sub_origin(a(1),a(2),a(3));
                               t_predictor(n,:)=predictor(a(1),a(2),a(3),:); 
                         end
                         clear ID  n a
                         
 %get the predicted subject id (the subject with highest RS with the origin subjects)   
 
%                         t_sub_origin(t_sub_origin<0)=0;   % using only voxels with positive activities
%                         t_predictor(t_predictor<0)=0;
                        if ~all(t_sub_origin == 0)
                             qq=corr(t_sub_origin,t_predictor);
                             [r_value,index]=max(qq);
                             prd_id(v_id,i)=index; 
                             prd_SI(v_id,i)=r_value; % the similarity value between conditions 
                             clear r_value index qq v_id
                        end
                    end 
 
             end 
%% get accuracy of individual identification
             prd_acc=zeros(61*73*61,sub_N);
            prd_acc_SI=zeros(61*73*61,sub_N);
             for  i=1:sub_N
                  subID_prd=prd_id(:,i);       
                  ID1=find(subID_prd==i);  % if the predicted ID matches the sujbect true ID, the identification is successful        
                  
                  SI=prd_SI(:,i);
                  ID2=find(~isnan(SI));               
                  id_good=intersect(ID1,ID2); % exclue voxels with Nan value
%                   
                  id_outlier=find(SI==1);  % exclude voxles on the edges with identical values in the t-map, which would results in 1 of the SI value;
                  id_final=setdiff(id_good,id_outlier);
                  
                  prd_acc(id_final,i)=1; % set the voexl values  to 1, if voxel pattern successfully identified the given subjecti 
                 prd_acc_SI(id_final,i)=SI(id_final);% save the associated between-condition similarity in activity  pattern
                  
                  clear SI ID1 ID2 id_good id_outlier id_final 
              end
              
% save volumes
             result=mean(prd_acc,2);
             result=reshape(result,61,73,61);
             q1=make_nii(result);
             q1.hdr=hdr;
             save_nii(q1,[dest,s_name,'_',p_name,'_acc.nii']);
             
             result2=mean(prd_acc_SI,2);
             result2=reshape(result2,61,73,61);
             q2=make_nii(result2);
             q2.hdr=hdr;
             save_nii(q2,[dest,s_name,'_',p_name,'_v.nii']);  
             toc            
             save([dest,'sub_acc_',s_name,'_',p_name],'prd_acc');
             
          end
     end
end
toc