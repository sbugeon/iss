function se = strel3D(o,r)
%Structuring element is determined using r which is the XY radius.
%Thus need to scale Z direction. (ATM THIS WILL JUST GIVE A 2D ARRAY BACK -
%NOT SURE IF NEED TO CHANGE THIS)
%https://uk.mathworks.com/matlabcentral/fileexchange/47937-3d-structuring-element-sphere

zRadius = floor(r*o.XYpixelsize/o.Zpixelsize);
[y,x,z]=meshgrid(-r:r,-r:r,-zRadius:zRadius); 
m=sqrt(x.^2 + y.^2 + (o.Zpixelsize*z/o.XYpixelsize).^2); 
b=(m <= r);

se=strel('arbitrary',b);
end