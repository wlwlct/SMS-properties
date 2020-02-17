%E00 E01 ratio is related to spectrum position, 

clearvars
solvent='F8T2400nmCH apd removed without consider marker';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\02252019\dataset intermediates\0'];
cd (srdir)

%Threshold_box=[0,1000;1001,2000;2001,3000;3001,4000;4001,5000;5001,6000;6001,7000;7001,250000];

allnames=struct2cell(dir([ '*.mat']));
[~,len]=size(allnames);
apdintensitycombine=[];
lifetimecombine=[];
averagewavelengthcombine=[];
   E00wavelengthcombine=[];
   E0001ratiocombine=[];
   spectrum=[];

for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    
   apdintensitycombine=cat(2,apdintensitycombine,datasetfile.dataset.scatterplot.intensity(1,:));
   
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
   Lifetime(Lifindexremove,:)=[];
   lifetimecombine=cat(1,lifetimecombine,Lifetime);

   averagewavelengthcombine=cat(1,averagewavelengthcombine,datasetfile.dataset.scatterplot.spectrum(:,1));
  [~,place488]=min(abs(datasetfile.dataset.ccdt(:,1)-400));
    place488=1;
  [maxvalue,maxindex]=max(datasetfile.dataset.ccdt(place488:end,3:end),[],1);
  wavelengthindex=datasetfile.dataset.ccdt(maxindex,1);
  E00wavelengthcombine=cat(1,E00wavelengthcombine,wavelengthindex);
  %ratio  
  E00sum=sum(datasetfile.dataset.ccdt(27:31,3:end),1);
  E01sum=sum(datasetfile.dataset.ccdt(36:40,3:end),1);
  Ratio=E00sum./E01sum;
  E0001ratiocombine=cat(2,E0001ratiocombine,Ratio);
  %spectrum
  if len_i==1
      spectrum=sum(datasetfile.dataset.ccdt(:,3:end),2);
      x=datasetfile.dataset.ccdt(:,1);
  else
  if x==datasetfile.dataset.ccdt(:,1)
  spectrum=sum(datasetfile.dataset.ccdt(:,3:end),2)+spectrum;
  else
      disp(['Something wrong wiht the spectrum X in file ' name])
  end
  end
end
 %% 
  figure;
  histogram(apdintensitycombine,100);
  title(['APD intensity distribution in ' solvent])
  
  try
    cd([srdir '/properties/']);
catch
    mkdir([srdir '/properties/']);
    cd([srdir '/properties/']);
end
 saveas(gcf,[solvent ' APD intensity distribution.jpg']);
 saveas(gcf,[solvent ' APD intensity distribution.fig']);
disp('Save APD intensity successfully /n')

%%
figure;
  histogram(lifetimecombine,100);
  title(['Lifetime distribution in ' solvent])
  saveas(gcf,[solvent ' lifetime distribution.jpg']);
  saveas(gcf,[solvent ' lifetime distribution.fig']);
  disp('Save lifetime successfully /n')
%%
  figure;
  histogram(averagewavelengthcombine,100);
  title(['average wavelength distribution in ' solvent])
  saveas(gcf,[solvent ' average wavelength distribution.jpg']);
  saveas(gcf,[solvent ' average wavelength distribution.fig']);
  disp('Save average wavelength successfully /n');
 %% 
 figure;
   histogram(E00wavelengthcombine,100);
  title(['E00 distribution in ' solvent])
  saveas(gcf,[solvent ' E00 distribution.jpg']);
  saveas(gcf,[solvent ' E00 distribution.fig']);
  disp('Save E00 successfully /n');
  %%
  figure  
  histogram(E0001ratiocombine,100);
  title(['E00 E01 Ratio distribution in ' solvent])
  saveas(gcf,[solvent ' E0001 Ratio distribution.jpg']);
  saveas(gcf,[solvent ' E0001 Ratio distribution.fig']);
  disp('Save E0001 Ratio successfully /n');
  %%
figure;
plot(x,spectrum);
  title(['Spectrum in ' solvent])
  saveas(gcf,[solvent ' Spectrum add up.jpg']);
  saveas(gcf,[solvent ' Spectrum add up.fig']);
  disp('Save Spectrum Ratio successfully /n');
