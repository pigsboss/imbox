%% point-source test
close all;
clc;
clear
N = 128;
nlvl = 1; %white noise level
radius = 3;
alpha = 1.5;
wt = fspecial('gauss',2*radius+1,alpha);

% f0 = zeros(N);
% f0(N/4,N/4)    =200;
% f0(N/4,N/2)    =200;
% f0(N/4,N/4*3)  =200;
% f0(N/2,N/4)    =200;
% f0(N/2,N/2)    =200;
% f0(N/2,N/4*3)  =200;
% f0(N/4*3,N/4)  =200;
% f0(N/4*3,N/2)  =200;
% f0(N/4*3,N/4*3)=200;
f1 = zeros(N);
f1(N/4,N/4)    =100;
f1(N/4,N/2)    =200;
f1(N/4,N/4*3)  =300;
f1(N/2,N/4)    =400;
f1(N/2,N/2)    =500;
f1(N/2,N/4*3)  =600;
f1(N/4*3,N/4)  =700;
f1(N/4*3,N/2)  =800;
f1(N/4*3,N/4*3)=900;

% f0 = gsmooth(f0+3*nlvl,1.4);
f1 = gsmooth(f1+3*nlvl,3.2);

% g0 = normrnd(f0,nlvl);
g1 = normrnd(f1,nlvl);

% tic;x0 = inlmeans(g0,'rot',false,'reps',1024);toc
tic;[x1,~,xz] = inlmeans(g1,'rot',false,'reps',N^2);toc
% tic;y0 = inlmeans(g0,'rot',true, 'reps',1024);toc
tic;[y1,~,yz] = inlmeans(g1,'rot',true, 'reps',N^2);toc

% disp(['residual 0 w/o rot: ', num2str(std(x0(:)-f0(:)))])
% disp(['residual 0 w/  rot: ', num2str(std(y0(:)-f0(:)))])
disp(['residual 1 w/o rot: ', num2str(std(x1(:)-f1(:)))])
disp(['residual 1 w   rot: ', num2str(std(y1(:)-f1(:)))])
%%
subplot(2,4,1),imagesc(f1);title('original');axis image
subplot(2,4,5),imagesc(g1);title('noisy');axis image
subplot(2,4,2),imagesc(x1);title('denoise, w/o rotation');axis image
subplot(2,4,6),imagesc(y1);title('denoise, w/ rotation');axis image
subplot(2,4,3),imagesc(x1-g1);title('method noise, w/o rotation');axis image
subplot(2,4,7),imagesc(y1-g1);title('method noise, w/ rotation');axis image
subplot(2,4,4),imagesc(log10(xz));title('similarity, w/o rotation');axis image
subplot(2,4,8),imagesc(log10(yz));title('similarity, w/ rotation');axis image
%% spiral-galaxy image test
close all
clear
clc
N = 128;
nlvl = 5;
radius = 3;
alpha = 1.5;
f = pinwheel(N,'bw')*100;
figure('Name','spiral-galaxy image test')
subplot(2,3,1);imagesc(f);title('original image');axis image;
g = normrnd(f,nlvl);
subplot(2,3,4);imagesc(g);title('noisy image');axis image;
x = nlmeans(g,radius,alpha,nlvl);
y = inlmeans(g,radius,alpha,nlvl);
subplot(2,3,2);imagesc(x);title('nlmeans denoise');axis image;
subplot(2,3,5);imagesc(y);title('inlmeans denoise');axis image;
subplot(2,3,3);imagesc(x-g);title('nlmeans method noise');axis image;
subplot(2,3,6);imagesc(y-g);title('inlmeans method noise');axis image;
%% elliptical-galaxy image test
close all
clear
clc
N = 512;
nlvl = 5;
radius = 3;
alpha = 1.5;
f = imread('../lib/hs-2012-24-a-xlarge_web.jpg');
xc = 607; yc = 349;
rect = [xc-256,yc-256,511,511];
f = imcrop(f,rect);
f = mean(double(f),3)/255*100;
f = imresize(f,[N N]);
figure('Name','elliptical-galaxy image test')
subplot(2,3,1);imagesc(f);title('original image');axis image;
g = normrnd(f,nlvl);
% subplot(2,3,4);imagesc(g);title('noisy image');axis image;
% x = nlmeans(g,radius,alpha,nlvl);
% y = inlmeans(g,radius,alpha,nlvl);
% subplot(2,3,2);imagesc(x);title('nlmeans denoise');axis image;
% subplot(2,3,5);imagesc(y);title('inlmeans denoise');axis image;
% subplot(2,3,3);imagesc(x-g);title('nlmeans method noise');axis image;
% subplot(2,3,6);imagesc(y-g);title('inlmeans method noise');axis image;
