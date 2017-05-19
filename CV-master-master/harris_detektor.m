function [Merkmale] = harris_detektor(Image,segment_length,k,tau,min_dist,tile_size,N)
% In dieser Funktion soll der Harris-Detektor implementiert werden, der
% Merkmalspunkte aus dem Bild extrahiert
%if varargin= 'true';
if nargin==1, segment_length=6; k=0.05;tau=0.1;min_dist=10;tile_size=[200,300];N=10; end
if nargin==2, k=0.05;tau=0.1;min_dist=10;tile_size=[200,300];N=10;end
if nargin==3, tau=0.1;min_dist=10;tile_size=[200,300];N=10;end
if nargin==4, min_dist=10;tile_size=[200,300];N=10;end
if nargin==5, tile_size=[200,300];N=10;end
if nargin==6, N=10;end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%      check if it is a gray image      %%%%%%%%%%%%%%%%%%
%Imagesize=size(Image); 
%if numel(Imagesize)>2
 %   disp('Please give a Gray image!');
%else    

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%       naive harris detector           %%%%%%%%%%%%%%%%%%
sigma=1;
[a,b]=size(Image);
min_dist=min_dist^2;     
[Fx,Fy]=sobel_xy(Image);
kx=fix(a/tile_size(1,1)); 
ky=fix(b/tile_size(1,2)); 
g=fspecial('gaussian',max(1,fix(segment_length*sigma)),sigma);
Fx2=conv2(Fx.*Fx,g,'same');
Fy2=conv2(Fy.*Fy,g,'same');
Fxy=conv2(Fx.*Fy,g,'same');
H=(Fx2.*Fy2-Fxy.^2)-k*(Fx2+Fy2).^2;
mx=ordfilt2(H,segment_length^2,ones(segment_length));%%anti salt and pepper noise
[rows,cols]=find((H==mx)&(H>tau));
Merkmale=[rows,cols];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% kill those so called 'Ecke' on the boundary of the Image   %%%%%%%%
Merkmale(end,:)=[];
Bound1=find(Merkmale(:,1)==1);
Bound2=find(Merkmale(:,2)==1);
t=length(Bound1);
s=length(Bound2);
Merkmale = sortrows(Merkmale,1);
Merkmale(1:t,:)=[];
Merkmale = sortrows(Merkmale,2);
Merkmale(1:s,:)=[];
rows = Merkmale(:,1);
cols = Merkmale(:,2);

figure(1),imshow(Image),hold on,
plot(cols,rows,'ys')
title('\fontsize{20}Using naive detector');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%    set the minimum distance        %%%%%%%%%%%%%%%%%%%%%% 
b=length(Merkmale);
Merkmalenew=zeros(b,1);
for i= 1:b
    Merkmalenew(i,1)=H(Merkmale(i,1),Merkmale(i,2));
end
[Merkmalenew,index]=sort(Merkmalenew,'descend');
for i=1:b
    Merkmale(i,1)=Merkmale(index(i,1),1);
    Merkmale(i,2)=Merkmale(index(i,1),2);
end               % sort the Merkmale according to the value of H


%for i=1:b-1
%   for j=i+1:b
%     if ((Merkmale(i,1)-Merkmale(j,1)).^2+ (Merkmale(i,2)-Merkmale(j,2)).^2)<min_dist
            
%             Merkmale(j,1)=Merkmale(i,1);
%             Merkmale(j,2)=Merkmale(i,2);  
             % find if the distance is smal ler than min_dist. if so, the
             % point with smaller H change to the point with bigger H.
          
%      end
%   end
%end

for i=1:b
D = pdist(Merkmale,'euclidean');
SD = squareform(D);
SDT = tril(SD);
SDT(SDT==0)=NaN;
%if(sum(abs(SDT(:)))==0)
[j,b]=size(SDT);
if(j<1000)
    break
end   
[m,n]=find (SDT<=min_dist);
Merkmale(m(1),:)=[];
end

%Merkmale=unique(Merkmale,'rows');   % kill the repeated point.
rown=(Merkmale(:,1));
coln=(Merkmale(:,2));

figure(2),imshow(Image),hold on,
plot(coln,rown,'ys')
title('\fontsize{20}Using detector with minimum distance');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%&&&&&           block of tile_size  with limit N             %%%%%%%%%%%%%       
c=length(Merkmale);
H1=zeros(size(H));
for i=1:c
    H1(Merkmale(i,1),Merkmale(i,2))=H(Merkmale(i,1),Merkmale(i,2));   
end
test=zeros(tile_size);
for i=1:tile_size(1,1):kx*tile_size(1,1) 
   for j=1:tile_size(1,2):ky*tile_size(1,2) 
       test(1:tile_size(1,1),1:tile_size(1,2))=H1(i:(i+tile_size(1,1)-1),j:(j+tile_size(1,2)-1));
       asum=0;
       for k=1:c
          if ((i<=Merkmale(k,1))&&(Merkmale(k,1)<=(i+tile_size(1,1)))&&((j<=Merkmale(k,2))&&(Merkmale(k,2)<=(j+tile_size(1,2)))))
             asum=asum+1;
          end
       end                          % find the sum of the Merkmal in tile_size.
       if asum > N
          [aaa]=sort(test(:));
           T=aaa(end-N,1);
           test(test<T)=0;
           H1(i:(i+tile_size(1,1)-1),j:(j+tile_size(1,2)-1))=test(1:tile_size(1,1),1:tile_size(1,2));
       end
   end
end                    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  deal with the rest    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[a,b]=size(Image);
rea=a-kx*tile_size(1,1);
reb=b-ky*tile_size(1,2);
Nrea=fix(rea/tile_size(1,1)*N);
Nreb=fix(reb/tile_size(1,2)*N);
testrea=zeros(rea,tile_size(1,2));
testreb=zeros(tile_size(1,1),reb);
for j=1:tile_size(1,2):ky*tile_size(1,2)
   testrea(1:rea,1:tile_size(1,2))=H1(a-rea+1:a,j:(j+tile_size(1,2)-1));
   asum=0;
   for k=1:c
       if(((a-rea+1<=Merkmale(k,1))&&(Merkmale(k,1)<=a))&&((j<=Merkmale(k,2))&&(Merkmale(k,2)<=(j+tile_size(1,2)))))
          asum=asum+1;
       end
   end
   if asum > Nrea
     [aaa]=sort(testrea(:));
     T=aaa(end-Nrea,1);
     testrea(testrea<T)=0;
     H1(a-rea+1:a,j:(j+tile_size(1,2)-1))=testrea(1:rea,1:tile_size(1,2));
   end      
end
for i=1:tile_size(1,1):kx*tile_size(1,1)
    testreb(1:tile_size(1,1),1:reb)=H1(i:(i+tile_size(1,1)-1),b-reb+1:b);
    asum=0;
    for k=1:c
        if((i<=Merkmale(k,1))&&(Merkmale(k,1)<=(i+tile_size(1,1)))&&((b-reb+1<=Merkmale(k,2))&&(Merkmale(k,2)<=b)))
           asum=asum+1;
        end
     end
     if asum > Nreb
        [aaa]=sort(testreb(:));
         T=aaa(end-Nrea,1);
         testreb(testreb<T)=0;
         H1(i:(i+tile_size(1,1)-1),b-reb+1:b)=testreb(1:tile_size(1,1),1:reb);
     end 
end

testc=zeros(rea,reb);
Nrest=fix(Nrea*Nreb/N);
for i=a-rea+1:a
    for j=b-reb+1:b
        testc=H(a-rea+1:a,b-reb+1:b); 
        asum=0;
        for k=1:c
           if ((a-rea+1<=Merkmale(k,1))&&(Merkmale(k,1)<=a))&&((b-reb+1<=Merkmale(k,2))&&(Merkmale(k,2)<=b))
              asum=asum+1;
           end
         end
           if asum>Nrest
             [aaa]=sort(testc(:));
             T=aaa(end-Nrest,1);
             testc(testc<T)=0;
             H1(a-rea+1:a,b-reb+1:b)=testc; 
           end
     end
end
[rowm,colm]=find((H1~=0));           % output the final merkmale.
Merkmale=[rowm,colm];
                                                            
figure(3),imshow(Image),hold on,
plot(colm,rowm,'ys')
title('\fontsize{20}Using detector with blocks of given size');
end
%end