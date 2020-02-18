clearvars
solvent='F8T2400nmCH apd removed without consider marker';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
cd (srdir)

allnames=struct2cell(dir([ '*.mat']));
[~,len]=size(allnames);
Threshold_box=[0,200;201,400;401,600;601,800;801,1000;1001,1200;1201,250000];
Threshold_leng=length(Threshold_box);

clearvars -except Threshold_box len allnames solvent srdir Threshold_leng 
for i=1:1:Threshold_leng

lifetime_combine=[];
lifetimewavelength=zeros(100,1);
Spectrum_combination=zeros(100,1);
Threshold1=Threshold_box(i,1);
Threshold2=Threshold_box(i,2);

 try
 %     openfig([srdir '/Distribution with 1000 threshold' '/' solvent ' ' num2str(Threshold) 'lifeDis']);
 %     obj=get(gca,'children');
 %     edge=get(obj,'BinEdges');
 %     N=get(obj,'BinCounts');
 %     avaliablelifetime=edge(any(N,1));
 avaliablelifetime=[1,200;201,400; 401,800;801,1800;1801,2500];
     al_leng=length(avaliablelifetime(:,1));
 lifetimewavelength = repmat(struct('sum',zeros(100,1)), al_leng, 1 );
      %lifetimewavelength=struct('sum',cell(al_leng,1));
 %     lifetimewavelength(:).sum=zeros(100,1);
 catch
 end


for len_i=1:len
    clear name
    name=char(allnames(1,len_i));
datasetfile=load([srdir '/' name]);

index=find(datasetfile.dataset.scatterplot.intensity(1,:)<Threshold2 & datasetfile.dataset.scatterplot.intensity(1,:)>Threshold1);%This is apd data
%disp('Finish load file /n')
%%
%combine all the lifetime to see distribution of lifetime
%lifetime_combine=cat(1,lifetime_combine,datasetfile.dataset.scatterplot.lifetime(index,2));
% clear datasetfile
%disp('Finish lifetime distribution once /n')
%%
%Combine all the spectrum to see the shape of the spectrum in different
% %solvents
% if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
%     error('spectrum not same')
% end

%x=datasetfile.dataset.ccdt(:,1);
%single_spectrum_sum=sum(datasetfile.dataset.ccdt(:,index+2),2);
%Spectrum_combination=Spectrum_combination+single_spectrum_sum;
%clear datasetfile
%disp('Finish add spectrum once /n')
%%
% %Combine the spectrum of same lifetime
  if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
      error('spectrum not same')
  end
  x=datasetfile.dataset.ccdt(:,1);
 
     for al_i=1:1:al_leng
 lifetime1_index=find(datasetfile.dataset.scatterplot.lifetime(index,2)>=avaliablelifetime(al_i,1));
 lifetime2_index=find(datasetfile.dataset.scatterplot.lifetime(index,2)<=avaliablelifetime(al_i,2));
 lifetime_index=intersect(lifetime1_index,lifetime2_index);
 interesting_spectrum_sum=sum(datasetfile.dataset.ccdt(:,lifetime_index+2),2);
 lifetimewavelength(al_i).sum=lifetimewavelength(al_i).sum(:,1)+interesting_spectrum_sum;
     end
 clear datasetfile
%%
%Generate after spectrum minus before spectrum, to see which part disappear
%and which apart start to show up

% if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
%     error('spectrum not same')
% end
% 
% x=datasetfile.dataset.ccdt(:,1);
% originccdt=datasetfile.dataset.ccdt(:,4:end);
% copyccdt=datasetfile.dataset.ccdt(:,3:end-1);
% difccdt=originccdt-copyccdt;
% 
% % mesh(1:1:length(difccdt(1,:)),x,difccdt);
% % colormap(jet)
% difccdt_leng=length(difccdt(1,:));
% for difccdt_i=1:1:difccdt_leng
%     figure('rend','painters','pos',[-1347 64 297 918])
%     subplot(3,1,1)
%     plot(x,originccdt(:,difccdt_i));
%     xlabel('wavelength');
%     ylabel('Intensity');
%     
%     subplot(3,1,2)
%     plot(x,copyccdt(:,difccdt_i));
%     xlabel('wavelength');
%     ylabel('Intensity');
%     
%     subplot(3,1,3)
%     plot(x,difccdt(:,difccdt_i));
%     xlabel('wavelength');
%     ylabel('Intensity');
%     title(['spectrum change from sec' num2str(difccdt_i) ' to sec ' num2str(difccdt_i+1)]);
%      saveas(gcf,[srdir '/spectrum change/' solvent ' ' num2str(difccdt_i) ' ' num2str(Threshold) 'change.jpg']);
%      saveas(gcf,[srdir '/spectrum change/' solvent ' ' num2str(difccdt_i) ' ' num2str(Threshold) 'change.fig']);
%     
% end

end


%figure
%histogram(lifetime_combine,1:1:2600);
%xlabel('Lifetime (ps)')
%ylabel('Occurance');
%title(['Lifetime distribution in ' solvent ' solution lower than threshold ' num2str(Threshold)])
%try
%    cd([srdir '/lifeDis/']);
%catch
%    mkdir([srdir '/lifeDis/']);
%    cd([srdir '/lifeDis/']);
%end
% saveas(gcf,[solvent ' below' num2str(Threshold) 'lifeDis.jpg']);
% saveas(gcf,[solvent ' below' num2str(Threshold) 'lifeDis.fig']);
%disp('Save distribution successfully /n')
% 
% 
%figure
%plot(x,Spectrum_combination)
%xlabel('wavelegnth (nm)')
%ylabel('Intensity')
%title([solvent ' spectrum combination below ' num2str(Threshold) ' threshold'])
%try
%    cd([srdir '/Spectrum/']);
%catch
%    mkdir([srdir '/Spectrum/']);
%    cd([srdir '/Spectrum/']);
%end
% saveas(gcf,[solvent ' below' num2str(Threshold) 'Spectrum.jpg']);
% saveas(gcf,[solvent ' below' num2str(Threshold) 'Spectrum.fig']);
%disp('Save spectrum successfully /n')

 for al_ii=1:1:al_leng
  figure
  plot(x,lifetimewavelength(al_ii).sum)
  xlabel('wavelegnth (nm)');
  ylabel('Intensity');
  title([solvent ' spectrum combination of lifetime ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' with ' num2str(Threshold1) ' threshold']);
  try
     cd([srdir '/100ps Spectrum corresponding to lifetime/']);
 catch
     mkdir([srdir '/100ps Spectrum corresponding to lifetime/']);
     cd([srdir '/100ps Spectrum corresponding to lifetime/']);
 end
  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' ' num2str(Threshold1) 'Spectrum.jpg']);
  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' ' num2str(Threshold1) 'Spectrum.fig']);
 end


end