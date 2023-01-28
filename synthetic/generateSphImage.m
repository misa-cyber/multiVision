%% Sphere image simulation
%   GENERATESPHIMAGE simulates noisy image
%   Addition of additive, multiplicative noise, spota and speckles and
%   obstacles
%
%   Level of noise added increases with the value of input sig   
%
%   M. A. Isa UoN, 2021
function [I,center,rad,ang, Xe]=generateSphImage(sig)
    pad=30;
    
    aRange=[100,230]; eRange=[0,0.42];
    
    sphAmpRange=[155,198]; spotAmpRange=[235,255]; backAmpRange=[10,95]; sectAmpRange=[75,110]; nSpotRange=[30,120];%Intensity amp ranges
    sphAmp=randRange(sphAmpRange,1); backAmp=randRange(backAmpRange,1); sphAmpDev=randRange([25,55],1); sectAmp=randRange(sectAmpRange,1);
    
    numSpots=round(randRange([3,15],1)); spotSizes=randRange([15,34],numSpots); spotAmps=randRange(spotAmpRange,numSpots); spotAng=randRange([pi,2*pi],numSpots);
    noiseSpots=round(randRange(sig*[200,400],1)); noiseSizes=[4,12]; 

    

    % generate ellipse
    a=randRange(aRange,1); e=randRange(eRange,1); imSize=round(2*(a+pad))*[1,1]; 
    center=imSize/2; b=sqrt(a^2-(e*a)^2); ang=randRange([0,2*pi],1);
    [u0,v0]=plotEllipse(center,[a,b],ang,'r',round(a*7),false); Xe=[u0,v0]; rad=[a,b];

    %
    I=zeros(imSize,'uint8'); I(:,:)=backAmp+randi(15,size(I)); 
    
    for i=1:randi(4)
        sectSize=randRange([0.8*a,1.5*a],1); sectPos=[randRange([-1.1*a,1.1*a],1),randRange([-1.1*a,1.1*a],1)];
        Xb=randomShape(center+sectPos,2.5*sectSize,4*round(sectSize));
        I=insertShape(I,"FilledPolygon",Xb,"Color",[sectAmp,sectAmp,sectAmp],"Opacity",0.9); I=rgb2gray(I); 
    end
    
    I=imnoise(I,"gaussian",0,0.1*sig^2/255^2); 
    I=imnoise(I,'speckle',0.1*sig^2/255^2); 
    I=imgaussfilt(I,0.1*sig);
    
    for i=1:noiseSpots
        spotCnt=randRange([1,imSize(1)],2)'; noiseSize=randRange(noiseSizes,1); nAmp=randRange(nSpotRange,1);
        Xi=randomShape(spotCnt,noiseSize,15); 
        I=insertShape(I,"FilledPolygon",Xi,"Opacity",0.8,"Color",[nAmp,nAmp,nAmp]);
    end
    I=im2gray(I);

    I=addEllipse(I,center,[a,b],ang,sphAmp,sphAmpDev);
    
    %add spots
    %spotStruct(numSpots)=struct();
    for i=1:numSpots
        spotRad=randRange([b/2,b-spotSizes(i)],1); spotCnt=center+ spotRad*[cos(spotAng(i)), sin(spotAng(i))];
        Xi=randomShape(spotCnt,spotSizes(i),7*round(spotSizes(i))); 
        %pgon1=polyshape(Xi(:,1),Xi(:,2));
        I=insertShape(I,"FilledPolygon",Xi,"Opacity",0.8,"Color",[spotAmps(i),spotAmps(i),spotAmps(i)]);
    
        %spotStruct(i).point=Xi;
    end
    I=rgb2gray(I);

    if randi(2)==1
        fPos=center+randRange([0.85*a,1.2*a],2)'.*[cos(-spotAng(1)),sin(-spotAng(1))];
        Xf=randomShape(fPos,1.8*sectSize,2*round(sectSize));
        I=insertShape(I,"FilledPolygon",Xf,"Color",0.4*[sectAmp,sectAmp,sectAmp],"Opacity",1); I=rgb2gray(I);
    end
     
     I=imgaussfilt(I,0.3*sig);
     I=imnoise(I,"gaussian",0,0.1*sig^2/255^2); 
     I=imnoise(I,'speckle',0.1*sig^2/255^2); 
     
     
    
    %figure;
    %imshow(I); hold on; 


    
end


function I=addEllipse(I,cnt,rad,ang,amp,ampDev)
    ampDir=pi; uDir=[cos(ampDir),sin(ampDir)]; %lighting direction
    for i=1:size(I,1)
        sI=1:size(I,2); u1_=double(i)-cnt(1); v1_=double(sI)-cnt(2);
        u1=-u1_*sin(ang)*ones(size(sI))-v1_.*cos(ang);
        v1=u1_*cos(ang)*ones(size(sI))-v1_.*sin(ang);
        
        flag=(u1/rad(1)).^2+(v1/rad(2)).^2<=1;
        I_add=amp*0.75+double(I(i,flag))*0.25; %keep 30% of backgroud to increase noise
        I(i,flag)=I_add+ampDev*dot(repmat(uDir,size(sI(flag),2),1)',[u1_*ones(size(sI(flag)));v1_(flag)])/rad(1);
    end
end

function y=randRange(x,n)
%generate n random numbers between x(1) and x(2)
    y=x(1)+(x(2)-x(1))*rand(n,1);
end