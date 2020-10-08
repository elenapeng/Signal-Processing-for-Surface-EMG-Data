clc
clear
close all

file=input ('input your data file sampled @ 1000 Hz: ','s');
MVC1=input ('input your data file of MVC1 sampled @ 1000 Hz: ','s');
MVC2=input ('input your data file of MVC2 sampled @ 1000 Hz: ','s');
MVCrest=input ('input your data file resting during MVC sampled @ 1000 Hz: ','s');
onsetRMSwindow=input ('input RMSwindow required for onset time calculation: ','s');
onsetRMSsmoothingwindow=input ('input smoothing of RMS window required for onset time calculation: ','s');
threshhold=input ('input threshhold, threshold = mean of baseline+(number to be entered)*standard deviation of baseline: ','s');
onsetRMSw=str2double(onsetRMSwindow);
onsetsmoothw=str2double(onsetRMSsmoothingwindow);
threshold=str2double(threshhold);
AmpRMSw=1;

activechannels=[1 2 3 4 5];
activemuscles=[1 2 3 4];

A = importdata(file, '\t');
EMGheaders = (A.colheaders); %These are EMG Channel headers;
EMGcols = (A.data); %Corresponding Raw EMG in Channels;

A2 = importdata(MVC1, '\t');
EMGheaders2 = (A2.colheaders); %These are EMG Channel headers;
EMGcols2 = (A2.data); %Corresponding Raw EMG in Channels;

A3 = importdata(MVC2, '\t');
EMGheaders3 = (A3.colheaders); %These are EMG Channel headers;
EMGcols3 = (A3.data); %Corresponding Raw EMG in Channels;

A4 = importdata(MVCrest, '\t');
EMGheaders4 = (A4.colheaders); %These are EMG Channel headers;
EMGcols4 = (A4.data); %Corresponding Raw EMG in Channels;

datafull2=EMGcols2(2:size(EMGcols2,1),:);
datafull3=EMGcols3(2:size(EMGcols3,1),:);
datafull4=EMGcols4(2:size(EMGcols4,1),:);

MVCRF{1,1}=datafull2(:,2);
MVCRF{1,2}=datafull3(:,2);
MVCRFrest=datafull4(:,2);

MVCVMO{1,1}=datafull2(:,3);
MVCVMO{1,2}=datafull3(:,3);
MVCVMOrest=datafull4(:,3);

MVCVL{1,1}=datafull2(:,4);
MVCVL{1,2}=datafull3(:,4);
MVCVLrest=datafull4(:,4);

MVCVM{1,1}=datafull2(:,5);
MVCVM{1,2}=datafull3(:,5);
MVCVMrest=datafull4(:,5);

MVCdataF{1,1}=datafull2(:,6);
MVCdataF{1,2}=datafull3(:,6);


datafull=EMGcols(2:size(EMGcols,1),:);
Time=datafull(:,1);
RF=datafull(:,2);
VMO=datafull(:,3);
VL=datafull(:,4);
VM=datafull(:,5);
dataF=datafull(:,6);

% 1=RF
% 2=VMO
% 3=VL
% 4=VM

% 5=Force

[b,a] = butter(4, 10/(1000/2),'high');
hp(:,2) = filtfilt(b,a,VMO);
hp(:,3) = filtfilt(b,a,VL);
hp(:,4) = filtfilt(b,a,VM);
hp(:,1) = filtfilt(b,a,RF);

for i=1:size(dataF,1)
    hp(i,5)=dataF(i,1)-mean(dataF(1001:1500,1)); % Force is not filtered
end
% hp(:,5)=dataF;

for i=1:2
    hpMVC{1,i}(:,1)=filtfilt(b,a,MVCRF{1,i});
    hpMVC{1,i}(:,2)=filtfilt(b,a,MVCVMO{1,i});
    hpMVC{1,i}(:,3)=filtfilt(b,a,MVCVL{1,i});
    hpMVC{1,i}(:,4)=filtfilt(b,a,MVCVM{1,i});
    hpMVC{1,i}(:,5)=MVCdataF{1,i};
end

hpMVCrest(:,2) = filtfilt(b,a,MVCVMOrest);
hpMVCrest(:,3) = filtfilt(b,a,MVCVLrest);
hpMVCrest(:,4) = filtfilt(b,a,MVCVMrest);
hpMVCrest(:,1) = filtfilt(b,a,MVCRFrest);

pad=[0];
for i=1:(onsetRMSw+onsetsmoothw)
    pad=[pad;0];
end
pad2=[0];
for i=1:(AmpRMSw)
    pad2=[pad2;0];
end

for j=1:size(hpMVCrest,2)
    hpMVCrest2(:,j)=hpMVCrest(:,j).^2;
    restRMSmusclesMVCrest(1,j)=sqrt(sum(hpMVCrest2(751:1250,j))/500);
end


for read=1:size(activechannels,2)
    j=activechannels(1,read);
    hp3{1,j}(:,1)=[hp(:,j);pad];
    hp2(:,j)=hp3{1,j}(:,1).^2;
    for i=1:size(hp,1)-onsetRMSw-1
        onsetRMS(i,j)=sqrt(sum(hp2(i:(i+onsetRMSw-1),j))/onsetRMSw);
    end
    for i=1:size(onsetRMS,1)-onsetsmoothw-1
        onsetsRMS(i,j)=mean(onsetRMS(i:(i+onsetsmoothw-1),j));
    end
    meanoff(1,j)=mean(onsetsRMS(1001:1500,j));
    stdoff(1,j)=std(onsetsRMS(1001:1500,j));
    meanoffend(1,j)=mean(onsetsRMS((size(onsetsRMS,1)-400):(size(onsetsRMS,1)-150),j));
    stdoffend(1,j)=std(onsetsRMS((size(onsetsRMS,1)-400):(size(onsetsRMS,1)-150),j));
    onsignal{1,j}(:,1)=find(onsetsRMS(:,j)>(max(onsetsRMS(:,j))/2));
    for i=1:100000000
        if onsetsRMS(onsignal{1,j}(1,1)-i,j)<(meanoff(1,j)+(threshold*stdoff(1,j)))
            trueontime{1,j}(1,1)=onsignal{1,j}(1,1)-i+1;
            break
        end
    end
    for i=1:100000000
        if (onsignal{1,j}(size(onsignal{1,j},1),1)+i)>size(onsetsRMS,1)
            trueontime{1,j}(1,2)=size(onsetsRMS,1);
            break
        else if onsetsRMS(onsignal{1,j}(size(onsignal{1,j},1),1)+i,j)<(meanoffend(1,j)+(threshold*stdoffend(1,j)))
                trueontime{1,j}(1,2)=onsignal{1,j}(size(onsignal{1,j},1),1)+i-1;
                break
            end
        end
    end
    if (0)
    p=1;
    ontime{1,j}(1,1)=onsignal{1,j}(1,1);
    for i=1:size(onsignal{1,j},1)-1
        if onsignal{1,j}(i+1,1)-onsignal{1,j}(i,1)~=1
            ontime{1,j}(p,2)=onsignal{1,j}(i,1);
            p=p+1;
            ontime{1,j}(p,1)=onsignal{1,j}(i+1,1);
        end
    end
    ontime{1,j}(p,2)=onsignal{1,j}(size(onsignal{1,j},1),1);
    q=1;
    for i=1:size(ontime{1,j},1)
        if (ontime{1,j}(i,2)-ontime{1,j}(i,1))>1000
            trueontime{1,j}(q,:)=ontime{1,j}(i,:);
            q=q+1;
        end
    end
    end
    Amphp3{1,j}(:,1)=[hp(:,j);pad2];
    Amphp2(:,j)=Amphp3{1,j}(:,1).^2;
    for i=1:size(hp,1)-AmpRMSw-1
        AmpRMS(i,j)=sqrt(sum(Amphp2(i:(i+AmpRMSw-1),j))/AmpRMSw);
    end
    r=1;
    for i=trueontime{1,j}(1,1)+10350:(trueontime{1,j}(1,2)-999)
        stdAmpRMS{1,j}(r,1)=std(AmpRMS(i:i+999,j));
        r=r+1;
    end
    [sortedstd{1,j} isortedstd{1,j}]=sort(stdAmpRMS{1,j});
    istablestart(1,j)=isortedstd{1,j}(1,1)-1+trueontime{1,j}(1,1)+10350;
    s=1;
    for i=istablestart(1,j):(istablestart(1,j)+999)
        stabletime{1,j}(s,1)=i;
        stable{1,j}(s,1)=AmpRMS(i,j);
        s=s+1;
    end
    
end

onsetRF=trueontime{1,1}(1,1)
onsetVMO=trueontime{1,2}(1,1)
onsetVL=trueontime{1,3}(1,1)
onsetVM=trueontime{1,4}(1,1)
onsetForce=trueontime{1,5}(1,1)
% stablestartRF=istablestart(1,1)
% stablestartVMO=istablestart(1,2)
% stablestartVL=istablestart(1,3)
% stablestartVM=istablestart(1,4)
stablestartForce=istablestart(1,5)

for j=1:4
    r2=1;
    for i=501:trueontime{1,j}(1,1)-1499
        stdrestNorm{1,j}(r2,1)=std(hp2(i:i+999,j));
        r2=r2+1;
    end
    [sortedstdrestNorm{1,j} isortedstdrestNorm{1,j}]=sort(stdrestNorm{1,j});
    istablestartrestNorm(1,j)=isortedstdrestNorm{1,j}(1,1)+500;
    stablerestRMSmuscles(1,j)=sqrt(sum(hp2(istablestartrestNorm(1,j):istablestartrestNorm(1,j)+999,j))/1000);
    stableRMSmuscles(1,j)=sqrt(sum(hp2(istablestart(1,5):istablestart(1,5)+999,j))/1000);
end



% figure, plot(VMO)
% title('VMO')
% figure, plot(onsetRMS(:,2))
% title('VMO-RMS')
for i=1:5
figure, plot(onsetsRMS(:,i))
if i==1
    title('RF-smoothenedRMS')
else if i==2
        title('VMO-smoothenedRMS')
    else if i==3
            title('VL-smoothenedRMS')
        else if i==4
                title('VM-smoothenedRMS')
            else if i==5
                    title('Force-smoothenedRMS')
                end
            end
        end
    end
end
hold on
plot(trueontime{1,i}(1,1),onsetsRMS(trueontime{1,i}(1,1),i),'o')
hold on
plot(trueontime{1,i}(1,2),onsetsRMS(trueontime{1,i}(1,2),i),'o')
end
% figure, plot(AmpRMS(:,2))
% title('VMO-RMS(Amp)')
% hold on
% plot(stabletime{1,2},stable{1,2},'r')

% figure, plot(dataF)
% title('Force')
% figure, plot(onsetRMS(:,5))
% title('Force-RMS')
% figure, plot(onsetsRMS(:,5))
% title('Force-smoothedRMS')
% hold on
% plot(trueontime{1,5}(1,1),onsetsRMS(trueontime{1,5}(1,1),5),'o')
% hold on
% plot(trueontime{1,5}(1,2),onsetsRMS(trueontime{1,5}(1,2),5),'o')
% figure, plot(AmpRMS(:,5))
% title('Force-RMS(Ampcalculation)')
% hold on
% plot(stabletime{1,5},stable{1,5},'r')

for k=1:2
    for j=1:size(hpMVC{1,k},2)
    hpMVC3{1,k}{1,j}(:,1)=[hpMVC{1,k}(:,j);pad];
    hpMVC2{1,k}(:,j)=hpMVC3{1,k}{1,j}(:,1).^2;
    end
    for j=5
    for i=1:size(hpMVC{1,k},1)-onsetRMSw-1
        onsetRMSMVC{1,k}(i,j)=sqrt(sum(hpMVC2{1,k}(i:(i+onsetRMSw-1),j))/onsetRMSw);
    end
    for i=1:size(onsetRMSMVC{1,k},1)-onsetsmoothw-1
        onsetsRMSMVC{1,k}(i,j)=mean(onsetRMSMVC{1,k}(i:(i+onsetsmoothw-1),j));
    end
    meanoffMVC{1,k}(1,j)=mean(onsetsRMSMVC{1,k}(1001:1500,j));
    stdoffMVC{1,k}(1,j)=std(onsetsRMSMVC{1,k}(1001:1500,j));
    onsignalMVC{1,k}{1,j}(:,1)=find(onsetsRMSMVC{1,k}(:,j)>(max(onsetsRMSMVC{1,k}(:,j))/4));
    p=1;
    ontimeMVC{1,k}{1,j}(1,1)=onsignalMVC{1,k}{1,j}(1,1);
    for i=1:size(onsignalMVC{1,k}{1,j},1)-1
        if onsignalMVC{1,k}{1,j}(i+1,1)-onsignalMVC{1,k}{1,j}(i,1)~=1
            ontimeMVC{1,k}{1,j}(p,2)=onsignalMVC{1,k}{1,j}(i,1);
            p=p+1;
            ontimeMVC{1,k}{1,j}(p,1)=onsignalMVC{1,k}{1,j}(i+1,1);
        end
    end
    ontimeMVC{1,k}{1,j}(p,2)=onsignalMVC{1,k}{1,j}(size(onsignalMVC{1,k}{1,j},1),1);
    q1=1;
    for i=1:size(ontimeMVC{1,k}{1,j},1)
        if (ontimeMVC{1,k}{1,j}(i,2)-ontimeMVC{1,k}{1,j}(i,1))>450
            trueontimeMVC{1,k}{1,j}(q1,:)=ontimeMVC{1,k}{1,j}(i,:);
            q1=q1+1;
        end
    end
    AmphpMVC3{1,k}{1,j}(:,1)=[hpMVC{1,k}(:,j);pad2];
    AmphpMVC2{1,k}(:,j)=AmphpMVC3{1,k}{1,j}(:,1).^2;
    for i=1:size(hpMVC{1,k},1)-AmpRMSw-1
        AmpRMSMVC{1,k}(i,j)=sqrt(sum(AmphpMVC2{1,k}(i:(i+AmpRMSw-1),j))/AmpRMSw);
    end
    r=1;
    for i=trueontimeMVC{1,k}{1,j}(1,1):(trueontimeMVC{1,k}{1,j}(1,2)-499)
        stdAmpRMSMVC{1,k}{1,j}(r,1)=std(AmpRMSMVC{1,k}(i:i+499,j));
        r=r+1;
    end
    [sortedstdMVC{1,k}{1,j} isortedstdMVC{1,k}{1,j}]=sort(stdAmpRMSMVC{1,k}{1,j});
    istablestartMVC{1,k}(1,j)=isortedstdMVC{1,k}{1,j}(1,1)-1+trueontimeMVC{1,k}{1,j}(1,1);
    s=1;
    for i=istablestartMVC{1,k}(1,j):(istablestartMVC{1,k}(1,j)+499)
        stabletimeMVC{1,k}{1,j}(s,1)=i;
        stableMVC{1,k}{1,j}(s,1)=AmpRMSMVC{1,k}(i,j);
        s=s+1;
    end
    end
    for j=1:4
        stableRMSmusclesMVC{1,k}(1,j)=sqrt(sum(hpMVC2{1,k}(istablestartMVC{1,k}(1,5):istablestartMVC{1,k}(1,5)+499,j))/500);
    end
%     figure, plot(AmpRMSMVC{1,k}(:,5))
%     hold on
%     plot(stabletimeMVC{1,k}{1,5},stableMVC{1,k}{1,5},'r')
end
for k=1:2
    for j=1:4
        matstableRMSmusclesMVC2(k,j)=stableRMSmusclesMVC{1,k}(1,j);
    end
end
for j=1:4
    matstableRMSmusclesMVC(1,j)=mean(matstableRMSmusclesMVC2(:,j));
end
for read2=1:size(activemuscles,2)
    j=activemuscles(1,read2);
    Normalized(1,j)=(stableRMSmuscles(1,j)-stablerestRMSmuscles(1,j))/(matstableRMSmusclesMVC(1,j)-restRMSmusclesMVCrest(1,j));
end

NormalizedRMSstableRF=Normalized(1,1)
NormalizedRMSstableVMO=Normalized(1,2)
NormalizedRMSstableVL=Normalized(1,3)
NormalizedRMSstableVM=Normalized(1,4)