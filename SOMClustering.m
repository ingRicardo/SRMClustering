% cl_img_som - Image Classification SOM
%
% name of the file 'OriginalImage.jpg'
%
% It reads the file and transforms from 256 x 256 x 3 to 65536 x 3
imag=double(imread('OriginalImage.jpg')); % it transforms it in double to be 
% processed
[x1 y1 z1]=size(imag); % 
% 3D array transformation to 2D array  x1 x y1 x z1 --> (x1 * y1) x z1
imagx=zeros((x1.*y1),z1); conta=1; % it creates output table
for i=1:x1
    for j=1:y1
        for k=1:z1
            imagx(conta,k)=imag(i,j,k);
        end
        conta=conta+1;
    end
end
% training data
treina=[];
for i=1:8:(x1*y1)
    treina=[treina; imagx(i,:)];
end
%it creates SOM network and training starts
net=newsom([0 255; 0 255; 0 255], [1 4], 'gridtop');
net.trainParam.epochs = 50;
net = train(net,treina');
% plotsom(net.iw{1,1},net.layers{1}.distances);
c_int=uint8(round(net.iw{:,:})); % get centers and rounding
c=double(c_int);
% classify
class=vec2ind(sim(net,imagx')); % get matrix with indexes
imagy=uint8(zeros((x1.*y1),z1)); % it creates output table unsigned int
%  vector of classes substitution
for i=1:(x1*y1)
    imagy(i,:)=c_int(class(i),:);
end
% It transforms to 3D array 
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
imwrite(imagz,'output_som.jpg','jpg','Quality',100);
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
% it plots the centers
cr=c(:,1);
cg=c(:,2);
cb=c(:,3);
figure(2); hold on; grid on;
box on; axis([0 256 0 256 0 256]);
plot3(cr,cg,cb,'r*');
plot3([0; 255],[0; 255],[0; 255],'r');
title('SOM Center Location','FontSize',14);
xlabel('Red','FontSize',14);
ylabel('Green','FontSize',14);
zlabel('Blue','FontSize',14);