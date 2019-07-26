function se = strel3D_2(rXY,rZ)
%Structuring element is determined using r which is the XY radius.
%Thus need to scale Z direction. (ATM THIS WILL JUST GIVE A 2D ARRAY BACK -
%NOT SURE IF NEED TO CHANGE THIS)
%https://uk.mathworks.com/matlabcentral/fileexchange/47937-3d-structuring-element-sphere

[y,x,z]=meshgrid(-rXY:rXY,-rXY:rXY,-rZ:rZ); 
m=sqrt(x.^2 + y.^2 + z.^2); 
b=(m <= rXY);

se=strel('arbitrary',b);
end