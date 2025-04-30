function L = DPA_define_searchlight_volume(ROIMask,Mask)
% step1_defineSearchlight_volume_llf('mask_spm_2level.nii','mask_spm_2level.nii')
% SPM - RSA toolbox
% Defines a volumetric searchlight using all voxels that are > 0 in both 
% ROIMask and the Mask as centres for the searchlights and as potential voxels 
% to be included in the searchlight. 
% The Mask also defines the voxel space in which the searchlight will be
% defined. 
%SPM codes will be called.
% INPUTS
%   ROIMask: a string defining the file name of ROI
%
%   Mask:      a string defining the file name of Mask
%
% OPTIONS/VARARGIN: 
%   	sphere:     definition of searchlight sphere, either one of:
%                   [R]   scalar R means fixed radius R (in mm)
%                   [R C] use maximum radius R to find approximately C voxels
%
% OUTPUT:
%       Searchlight structure with the fields: 
%          LI:      nVox x 1 cell array with linear voxel indices （矩阵坐标）
%          voxmin:  nVox x 3 Minimal voxel coordinate in x,y,z direction (
%          i.e.,the coordinate of the center
%          voxmax:  nVox x 3 maximal voxel coordinate in x,y,z
%          direction,i.e., the coordinate farest to the center.
           
%          voxel:   nVox x 3 Matrix of I,J,K voxel coordinates or 
%                   nVox x 1 vector of linear voxel indices for the centers of search-lights 
%
%
% 2/2015 - Joern Diedrichsen & Naveed Ejaz 
% modified by LLF,2018.08.12


%% 2. Getting sphere definition, and setting up reference mask

maskpath='.\';
cd(maskpath);
ROIMask='brainmask_3mm.nii'; % the mask produced by SPM during second-level activation analysis. This mask excludes voxels with NaN value in most subjects.
Mask='brainmask_3mm.nii';
sphere=[5,125];


radius      = sphere;
fixedradius = numel(sphere)==1;
if ~fixedradius
    targetvoxelcount = sphere(2);
end


%% 3. Estimate searchlight voxel indices over the provided region
% - estimate valid voxels given region and functional masks

roi=load_nii(ROIMask);
roi=roi.img;
centeridxs_roi=find(roi>0);
[x,y,z] = ind2sub(size(roi),centeridxs_roi); 
c_centers=[x';y';z'];

mask=load_nii(Mask);
mask=mask.img;
%	- calculate voxels of interest and center voxels
centeridxs_mask = find(mask>0);  
[x,y,z] = ind2sub(size(roi),centeridxs_mask); 
c_voxels = [x';y';z'];

centers     = uint32(c_centers);
voxels      = uint32(c_voxels);
ncent       = size(centers,2); 

%   - preallocate the matrices/cell arrays
li          = cell(ncent,1);        % linear indices for voxels
n           = zeros(ncent,1);       % Number of voxels 
rs          = zeros(ncent,1);       % Searchlight radius 
voxmin      = zeros(ncent,1);       % bottom left voxel
voxmax      = zeros(ncent,1);       % top right voxel


%% 4. Estimate linear indices, voxmin/voxmax for a sphere centered at each center index
for k=1:ncent
    ds = surfing_eucldist(c_centers(:,k),c_voxels);
    if fixedradius
        a       = voxels(:,ds<radius);
        rs(k,1) = radius; 
    else 
        i       = find(ds<radius(1));
        [dss,j] = sort(ds(i));
        indx    = min(targetvoxelcount,length(i)); % In case there are not enough voxels within maximal radius
        a       = voxels(:,i(j(1:indx))); 
        rs(k,1) = dss(indx); 
    end
    n(k,1)          = size(a,2);
    voxmin(k,1:3)   = min(a,[],2)';
    voxmax(k,1:3)   = max(a,[],2)';
    li{k,1}         = sub2ind(size(mask),a(1:3,:)')';
    n(k,1)          = numel(li{k});
end

% 
% for k=1:ncent;
%     n(k,1)=numel(li{k});
% end

%% 5. Setting output
   L.LI        = li;
   L.voxmin    = voxmin;
   L.voxmax    = voxmax;
   L.voxel     = centeridxs_roi;
   save('searchlight_list.mat','L')
end



function ds=surfing_eucldist(xs,ys)
% computes the euclidian distance between points
%
% DS=SURFING_EUCLDISTANCEFROM(XS,YS) 
% INPUTS:
%   XS:   3xP matrix for P coordinates
%   YS:   3xQ matrix for Q coordinates
% OUTPUT
%   DS:   PxQ matrix, where DS(i,j) is the Euclidian distance between
%         XS(:,i) and YS(:,i)
%
% This function has also been implemented in C, which runs a lot faster.
% If MEX is set up, run "mex surfing_eucldist.c" for increased speed.
%
% NNO May 2010

  [three1,n1]=size(xs);
  if three1 ~= 3, error('xs should be 3xN'); end

     [three2,n2]=size(ys);
  if three2 ~= 3, error('ys should be 3xN'); end

% assume N1 is much less than N2; if not swap the order.
  if three1>three2
    ds=surfing(eucldist(ys,xs))';
  else
  ds=zeros(n1,n2);
    for k=1:n1
        % use array indexing for faster execution
        ds(k,:)=sqrt(sum((repmat(xs(:,k),1,n2)-ys).^2))';
    end
  end
end

