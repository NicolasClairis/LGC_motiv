V=spm_vol('ROI.nii');
[I, XYZ]=spm_read_vols(V);

a=deg2rad(34.6);

ROI_temp=zeros(size(I));

%ROI_temp(XYZ(1,:)>X_vec(1)& XYZ(1,:)<X_vec(2) & XYZ(2,:)>Y_vec(1)& XYZ(2,:)<Y_vec(2) & XYZ(3,:)>Z_vec(1)& XYZ(3,:)<Z_vec(2))= 1;

Dim=size(I);

[y1, x1, z1]=COG(I);

%Center_mm=[-20.6 -26.8 3.7];

%Unit_vxl=[(max(XYZ(1,:))-min(XYZ(1,:)))./Dim(1) (max(XYZ(2,:))-min(XYZ(2,:)))./Dim(2) (max(XYZ(3,:))-min(XYZ(3,:)))./Dim(3)];


%Center_vxl=Center_mm./Unit_vxl;

%x1=Center_vxl(1);y1=Center_vxl(2);z1=Center_vxl(3);


d = [1,0,0,x1; 0,1,0,y1; 0,0,1,z1;0,0,0,1];

%Rot_mat=[cos(a) 0 sin(a) 0; 0 1 0 0; -sin(a) 0 cos(a) 0; 0 0 0 1];
 Rot_mat=[1 0 0 0; 0 cos(a) -sin(a) 0; 0 sin(a) cos(a) 0; 0 0 0 1];
 %Rot_mat=[cos(a) -sin(a) 0 0; sin(a) cos(a) 0 0; 0 0 1 0; 0 0 0 1]; For
 %right-left DL PFC

c = [1,0,0,-x1; 0,1,0,-y1; 0,0,1,-z1;0,0,0,1];

M=c'*Rot_mat*d';

%M =Rot_mat*d';


T_rot=affine3d(M);
[cb_rot] = imwarp_same(I,T_rot);


V_temp=V;
V_temp.fname='ROI_rot.nii';


V2=spm_write_vol(V_temp,cb_rot);


