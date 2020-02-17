clear all
srdir='C:\Users\Livi\Documents\Results\Interesting from 0930 to 1008\newset';
cd (srdir)
solvent='Chloroform'
allnames=struct2cell(dir([solvent '*.mat']));
[~,len]=size(allnames);
Threshold_box=[0,1000;1001,3000;3001,5000;5001,6000;6001,7000;7001,1500000];
Threshold_leng=length(Threshold_box);

clearvars -except Threshold_box len allnames solvent srdir Threshold_leng
for i=1:1:Threshold_leng

lifetime_combine=[];
lifetimewavelength=zeros(100,1);
Spectrum_combination=zeros(100,1);
%Threshold=Threshold_box(1,i);

% try
%     openfig([srdir '\Distribution with 1000 threshold' '\' solvent ' ' num2str(Threshold) 'lifeDis']);
%     obj=get(gca,'children');
%     edge=get(obj,'BinEdges');
%     N=get(obj,'BinCounts');
%     avaliablelifetime=edge(any(N,1));
%     al_leng=length(avaliablelifetime);
%     lifetimewavelength = repmat(struct('sum',zeros(100,1)), al_leng, 1 )
%      %lifetimewavelength=struct('sum',cell(al_leng,1));
% %     lifetimewavelength(:).sum=zeros(100,1);
% catch
% end



for len_i=1:len
    clear name
    name=char(allnames(1,len_i));
datasetfile=load([srdir '\' name]);
index_1=find(datasetfile.dataset.scatterplot.intensity(1,:)>Threshold_box(i,1));%This is apd data
index_2=find(datasetfile.dataset.scatterplot.intensity(1,:)<Threshold_box(i,2));%This is apd data
index=intersect(index_1,index_2);
%%
%combine all the lifetime to see distribution of lifetime
lifetime_combine=cat(1,lifetime_combine,datasetfile.dataset.scatterplot.lifetime(index,2));
clear datasetfile
%%
%Combine all the spectrum to see the shape of the spectrum in different
%solvents
% if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
%     error('spectrum not same')
% end
% 
% x=datasetfile.dataset.ccdt(:,1);
% single_spectrum_sum=sum(datasetfile.dataset.ccdt(:,index+2),2);
% Spectrum_combination=Spectrum_combination+single_spectrum_sum;
% clear datasetfile
%%
%Combine the spectrum of same lifetime
%  if len_i>1 && sum(x-datasetfile.dataset.ccdt(:,1))~=0
%      error('spectrum not same')
%  end
%  x=datasetfile.dataset.ccdt(:,1);
% 
%     
%     for al_i=1:1:al_leng
% lifetime_index=find(datasetfile.dataset.scatterplot.lifetime(index,2)==avaliablelifetime(1,al_i));
% interesting_spectrum_sum=sum(datasetfile.dataset.ccdt(:,lifetime_index+2),2);
% lifetimewavelength(al_i).sum=lifetimewavelength(al_i).sum(:,1)+interesting_spectrum_sum;
%     end
% clear datasetfile
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
%      saveas(gcf,[srdir '\spectrum change\' solvent ' ' num2str(difccdt_i) ' ' num2str(Threshold) 'change.jpg']);
%      saveas(gcf,[srdir '\spectrum change\' solvent ' ' num2str(difccdt_i) ' ' num2str(Threshold) 'change.fig']);
%     
% end





end


figure
histogram(lifetime_combine,1:1:2600);
xlabel('Lifetime (ps)')
ylabel('Occurance');
title(['Lifetime distribution in ' solvent ' solution wiht threshold ' num2str(Threshold)])
 saveas(gcf,[srdir '\lifeDis\' solvent ' ' num2str(Threshold) 'lifeDis.jpg']);
 saveas(gcf,[srdir '\lifeDis\' solvent ' ' num2str(Threshold) 'lifeDis.fig']);


% figure
% plot(x,Spectrum_combination)
% xlabel('wavelegnth (nm)')
% ylabel('Intensity')
% title([solvent ' spectrum combination with ' num2str(Threshold) ' threshold'])
%  saveas(gcf,[srdir '\Spectrum\' solvent ' ' num2str(Threshold) 'Spectrum.jpg']);
%  saveas(gcf,[srdir '\Spectrum\' solvent ' ' num2str(Threshold) 'Spectrum.fig']);

% for al_ii=1:1:al_leng
%  figure
%  plot(x,lifetimewavelength(al_ii).sum)
%  xlabel('wavelegnth (nm)');
%  ylabel('Intensity');
%  title([solvent ' spectrum combination of lifetime ' num2str(avaliablelifetime(1,al_ii)) ' with ' num2str(Threshold) ' threshold']);
%  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(1,al_ii)) ' ' num2str(Threshold) 'Spectrum.jpg']);
%  saveas(gcf,[solvent ' ' 'tau ' num2str(avaliablelifetime(1,al_ii)) ' ' num2str(Threshold) 'Spectrum.fig']);
% end


 end