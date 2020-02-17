clearvars
solvent='F8T2400nmCH apd removed without consider marker';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%  srdir=['E:\02252019\dataset intermediates\0'];
cd (srdir)

allnames=struct2cell(dir([ '*.mat']));
[~,len]=size(allnames);

clearvars -except Threshold_box len allnames solvent srdir Threshold_leng 

lifetime_combine=[];
lifetimewavelength=zeros(100,1);
Spectrum_combination=zeros(100,1);

 try
avaliablelifetime=[1,200;201,400; 401,800;801,1800;1801,2500];
     al_leng=length(avaliablelifetime(:,1));
 lifetimewavelength = repmat(struct('sum',zeros(100,1)), al_leng, 1 );
      %lifetimewavelength=struct('sum',cell(al_leng,1));
 %     lifetimewavelength(:).sum=zeros(100,1);
 catch
 end



for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
datasetfile=load([srdir '/' name]);

  if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
      error('spectrum not same')
  end
  x=datasetfile.dataset.ccdt(:,1);
 
     
     for al_i=1:1:al_leng
 lifetime1_index=find(datasetfile.dataset.scatterplot.lifetime(:,2)>=avaliablelifetime(al_i,1));
 lifetime2_index=find(datasetfile.dataset.scatterplot.lifetime(:,2)<=avaliablelifetime(al_i,2));
 lifetime_index=intersect(lifetime1_index,lifetime2_index);
 onlyccd=datasetfile.dataset.ccdt(:,3:end)./max(datasetfile.dataset.ccdt(1:end,3:end),[],1);
 interesting_spectrum_sum=sum(onlyccd(:,lifetime_index),2);
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
  plot(x,lifetimewavelength(al_ii).sum,'LineWidth',3)
  set(gca,'FontSize',16);
  set(gca,'FontWeight','Bold')
  xlabel('Wavelegnth (nm)');
  ylabel('Intensity');
  title([solvent ' spectrum combination of lifetime ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' with ' ' threshold']);
  try
     cd([srdir '/Spectrum corresponding to lifetime/']);
 catch
     mkdir([srdir '/Spectrum corresponding to lifetime/']);
     cd([srdir '/Spectrum corresponding to lifetime/']);
 end
  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' ' 'Spectrum.jpg']);
  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(al_ii,1)) 'to ' num2str(avaliablelifetime(al_ii,2)) ' ' 'Spectrum.fig']);
 end


% end
