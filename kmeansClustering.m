% cl_img_kmeans - Image classification with k-means
%
% name of image 'Output_kmeans.jpg'
%
%
clear all
%it reads the file in it transforms it from  256 x 256 x 3 to 65536 x 3
imag=double(imread('OriginalImage.jpg')); % it transforms to double to process 
[x1, y1, z1]=size(imag); %
% 3D array transformation to 2D array  x1 x y1 x z1 --> (x1 * y1) x z1
imagx=zeros((x1*y1),z1); conta=1; % it creates output table
for i=1:x1
    for j=1:y1
        for k=1:z1
            imagx(conta,k)=imag(i,j,k);
        end
        conta=conta+1;
    end
end
%it classifies
[c,u,saida,class] = Kmeans_var(imagx, 4, 500, 0.01);
c_int=round(c); % centers must be integer numbers
imagy=uint8(zeros((x1*y1),z1)); % it crates output area  unsigned int
% updating of classes centers
for i=1:(x1*y1)
    imagy(i,:)=c_int(class(i),:);
end
% it goes back 3D array
imagz=uint8(zeros(x1,y1,z1)); conta=1;
for i=1:x1
    for j=1:y1
        for k=1:z1
            imagz(i,j,k)=imagy(conta,k);
        end
        conta=conta+1;
    end
end
% it saves the image
imwrite(imagz,'Output_kmeans.jpg','jpg','Quality',100);
% it plots the points
r=imagx(:,1);
g=imagx(:,2);
b=imagx(:,3);
figure(1); plot3(r,g,b,'.'); grid on;
box on; axis tight;
title('Points Distribution','FontSize',14);
xlabel('Red','FontSize',14);
ylabel('Green','FontSize',14);
zlabel('Blue','FontSize',14);
% plota os centros
cr=c(:,1);
cg=c(:,2);
cb=c(:,3);
figure(2); hold on; grid on;
box on; axis([0 256 0 256 0 256]);
plot3(cr,cg,cb,'r*');
plot3([0; 255],[0; 255],[0; 255],'r');
title('Center Location','FontSize',14);
xlabel('Red','FontSize',14);
ylabel('Green','FontSize',14);
zlabel('Blue','FontSize',14);