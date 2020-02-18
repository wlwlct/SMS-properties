%SPectrum added up in the code, not support records with different spectrum
%Calculate the spectrum changed with time (diff ccd).

clearvars
 solvent='F8T2N2';
 srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];


% srdir=['E:\02252019\dataset intermediates\0'];
cd (srdir)

%Threshold_box=[0,1000;1001,2000;2001,3000;3001,4000;4001,5000;5001,6000;6001,7000;7001,150000];

allnames=struct2cell(dir([ '*.mat']));
[~,len]=size(allnames);
apdintensitycombine=[];
lifetimecombine=[];
averagewavelengthcombine=[];
   E00wavelengthcombine=[];
   
   timeave=zeros(99,len);
   timeE00=zeros(len,99);
   timelifetime=zeros(99,len);
   timeE0001=zeros(len,99);
   timeintensity=zeros(len,99);
   timespectrum=zeros(100,99);
   timePosAve=zeros(len,98);
   timeNegAve=zeros(len,98);
   timePosSpectrum=zeros(100,98);
   timeNegSpectrum=zeros(100,98);
for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
try
  [~,place488]=min(abs(datasetfile.dataset.ccdt(:,1)-488));
  [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place488+1:end,3:end),[],1);
  wavelengthindex=datasetfile.dataset.ccdt(place488+maxindex,1);
  %average wavelength change with time%each second,each column
  timeave(:,len_i)=datasetfile.dataset.scatterplot.spectrum(:,1);
  %E00 change with time
  timeE00(len_i,:)=wavelengthindex;
  %lifetime change with time
     Lifindexremove=[];
   [newconti_leng,~]=size(datasetfile.dataset.newconti);
     for newconti_i=1:1:newconti_leng
       [~,co_leng]=size(datasetfile.dataset.newconti(newconti_i).co);
       for co_i=1:1:co_leng
           preparemove=datasetfile.dataset.newconti(newconti_i).co(co_i).subco(1,2:end);
           Lifindexremove=cat(2,Lifindexremove,preparemove);
       end
   end
   Lifetime=datasetfile.dataset.scatterplot.lifetime(:,2);
   Lifetime(Lifindexremove,:)=0;
  timelifetime(:,len_i)=Lifetime;
  %Intensity change with time
  timeintensity(len_i,:)=datasetfile.dataset.scatterplot.intensity(1,:);
  %spectrum change with time
  timespectrum=datasetfile.dataset.ccdt(:,3:end)+timespectrum;
  %E0001 Ratio change with time
  E00sum=sum(datasetfile.dataset.ccdt(22:30,3:end),1);
  E01sum=sum(datasetfile.dataset.ccdt(31:42,3:end),1);
  timeE0001(len_i,:)=E00sum./E01sum;
  %Disappear of Wavlength distribution
 
  trimCCD=datasetfile.dataset.ccdt(:,3:end);
  difccd=diff(trimCCD,[],2);
  posdifccd=difccd;
  posdifccd(posdifccd<0)=0;
  negdifccd=-difccd;
  negdifccd(negdifccd<0)=0;
  if len_i==1
   x=datasetfile.dataset.ccdt(:,1);
  elseif x~=datasetfile.dataset.ccdt(:,1)
      error('wavelength number are not same')
  end
  
   timePosAve(len_i,:)=sum(posdifccd.*x,1)./sum(posdifccd,1);
  timeNegAve(len_i,:)=sum(negdifccd.*x,1)./sum(negdifccd,1);
  %Disaper wavelength add up change with time
  timePosSpectrum=timePosSpectrum+posdifccd./max(posdifccd,[],1);
  timeNegSpectrum=timeNegSpectrum+negdifccd./max(negdifccd,[],1);
  
  
catch
    disp([name 'may not have length 99'])
end
end
%%
edges=500:1:700;intensityedge=min(timeintensity):(max(timeintensity)-min(timeintensity))/100:max(timeintensity);
int10000=1:100:10000;
timeavehis=zeros(99,200);timeE00his=zeros(99,200);timelifetimehis=zeros(99,243);timeintensityhis=zeros(99,100);
timeE0001his=zeros(99,20);
%histogramize
for i=1:1:99
timeavehis(i,:)=histcounts(timeave(i,:),edges);
timeE00his(i,:)=histcounts(timeE00(:,i),edges);
timelifetimehis(i,:)=histcounts(timelifetime(i,:),70:10:2500);
timeintensityhis(i,:)=histcounts(timeintensity(:,i),intensityedge);
timeint10000his(i,:)=histcounts(timeintensity(:,i),int10000);
timeE0001his(i,:)=histcounts(timeE0001(:,i),0:0.05:1);
end
for i=1:1:98
    timePosAvehis(i,:)=histcounts(timePosAve(:,i),edges);
    timeNegAvehis(i,:)=histcounts(timeNegAve(:,i),edges);
end

meanintensity=mean(timeintensity,1);
stdintensity=std(timeintensity,1);





try
    cd([srdir '/time change PN/']);
catch
    mkdir([srdir '/time change PN/']);
    cd([srdir '/time change PN/']);
end

figure
surf(edges(1,2:end),1:1:99,timeavehis,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['average wavelength change with time ' solvent])
  saveas(gcf,[solvent ' average wavelength change with time.jpg']);
  saveas(gcf,[solvent ' average wavelength change with time.fig']);
  disp('Save average wavelength change with time successfully /n');
  close all
figure
surf(edges(1,2:end),1:1:99,timeE00his,'EdgeColor','none');
colormap(jet);
view([0 0 1])
  title(['E00 change with time ' solvent])
  saveas(gcf,[solvent ' E00 change with time.jpg']);
  saveas(gcf,[solvent ' E00 change with time.fig']);
  disp('Save E00 change with time successfully /n');
  close all
figure
surf(80:10:2500,1:1:99,timelifetimehis,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['lifeitme change with time ' solvent])
  saveas(gcf,[solvent ' lifetime change with time.jpg']);
  saveas(gcf,[solvent ' lifetime change with time.fig']);
  disp('Save lifeitme change with time successfully /n');
  close all
  
figure
surf(intensityedge(1,2:end),1:1:99,timeintensityhis,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['Intensity change with time ' solvent])
  saveas(gcf,[solvent ' Intensity change with time.jpg']);
  saveas(gcf,[solvent ' Intensity change with time.fig']);
  disp('Save Intensity change with time successfully /n');
  close all  
  
figure
surf(0.05:0.05:1,1:1:99,timeE0001his,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['E0001R change with time ' solvent])
  saveas(gcf,[solvent ' E0001R change with time.jpg']);
  saveas(gcf,[solvent ' E0001R change with time.fig']);
  disp('Save E0001R change with time successfully /n');
  close all  
  
figure
surf(1:1:99,datasetfile.dataset.ccdt(:,1),timespectrum,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['Spectrum (add up) change with time ' solvent])
  saveas(gcf,[solvent ' Spectrum (add up) change with time.jpg']);
  saveas(gcf,[solvent ' Spectrum (add up) change with time.fig']);
  disp('Save Spectrum (add up) change with time successfully /n');
  close all  
 
 figure
surf(int10000(1,2:end),1:1:99,timeint10000his,'EdgeColor','none');
colormap(jet)
view([0 0 1])
  title(['Intensity 10000 change with time ' solvent])
  saveas(gcf,[solvent ' Intensity 10000 change with time.jpg']);
  saveas(gcf,[solvent ' Intensity 10000 change with time.fig']);
  disp('Save Intensity 10000 change with time successfully /n');
  close all  
  

 figure;yyaxis left; plot(meanintensity);ylabel('mean');yyaxis right;plot(stdintensity); ylabel('standard deviation')
  xlabel('time');title([solvent 'Intensity change with time'])
  saveas(gcf,[solvent ' Intensity plot change with time.jpg']);
  saveas(gcf,[solvent ' Intensity plot change with time.fig']);
  disp('Save Intensity plot change with time successfully /n');
  close all  

  figure;
  set(gcf,'position',[0,79,447,874])
  subplot(2,1,1)
  surf(edges(1,2:end),1:1:98,timePosAvehis,'EdgeColor','none');view([0 0 1]);colormap(jet)
  title('Positive Spectrum Distribution');xlabel('Wavelength');ylabel('Number')
  subplot(2,1,2)
  surf(edges(1,2:end),1:1:98,timeNegAvehis,'EdgeColor','none');view([0 0 1]);colormap(jet)
  title('negative Spectrum Distribution');xlabel('Wavelength');ylabel('Number');
  saveas(gcf,[solvent ' Positive and Negative Distribution change with time.jpg']);
  saveas(gcf,[solvent ' Positive and Negative Distribution change with time.fig']);
  disp('Save Positive and Negative distribution with time successfully /n');
  close all 
  
  figure;
  set(gcf,'position',[0,79,447,874])
  subplot(2,1,1)
  surf(1:1:98,x,timePosSpectrum,'EdgeColor','none');view([0 0 1]);colormap(jet)
  title('Add of Positive Spectrum Change with Time');xlabel('Number');ylabel('Wavelength');
  subplot(2,1,2)
  surf(1:1:98,x,timeNegSpectrum,'EdgeColor','none');view([0 0 1]);colormap(jet)
  title('Add of Negative Spectrum Change with Time');xlabel('Number');ylabel('Wavelength');
    saveas(gcf,[solvent ' Positive and Negative Add up spectrum change with time.jpg']);
  saveas(gcf,[solvent ' Positive and Negative Add up spectrum change with time.fig']);
  disp('Save Positive and Negative add upspectrum with time successfully /n');
  close all 
