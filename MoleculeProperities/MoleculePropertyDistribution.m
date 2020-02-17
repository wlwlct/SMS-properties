%must not include data that not with same x value.
clearvars
solvent='F8T2400nmCH apd removed without consider marker';
 srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
%srdir=['E:\02252019\dataset intermediates\0'];
cd (srdir)

%Threshold_box=[0,1000;1001,2000;2001,3000;3001,4000;4001,5000;5001,6000;6001,7000;7001,250000];

allnames=struct2cell(dir([ '*.mat']));
[~,len]=size(allnames);
Mspectrum=zeros(100,len+1);


for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    datasetfile=load([srdir '/' name]);
    disp('Finish load file /n')
    
    Mspectrum(:,1)=datasetfile.dataset.ccdt(:,1);
    Mspectrum(:,1+len_i)=sum(datasetfile.dataset.ccdt(:,3:end),2);
    clearvars datasetfile
end

Mavewav= sum(Mspectrum(1:end,2:end).*Mspectrum(1:end,1),1)./sum(Mspectrum(1:end,2:end),1);
[~,maxwavP]=max(Mspectrum(1:end,2:end),[],1);Mmaxwav=Mspectrum(maxwavP,1);


try
    cd([srdir '/Molecule Properties/']);
catch
    mkdir([srdir '/Molecule properties/']);
    cd([srdir '/Molecule properties/']);
end
edges=430:2:650;

figure;
histogram(Mavewav,edges);
set(gca,'FontSize',16,'FontWeight','Bold');
title(['Molecule Average Wavelength Distribution in ' solvent])
saveas(gcf,[solvent ' Molecule Average Wavelength Distribution.jpg']);
saveas(gcf,[solvent ' Molecule Average Wavelength Distribution.fig']);
disp('Save Molecule Average Wavelenth successfully /n')

figure;
histogram(Mmaxwav,edges);
set(gca,'FontSize',16,'FontWeight','Bold');
title(['Molecule Maximum Wavelength Distribution in ' solvent])
saveas(gcf,[solvent ' Molecule Maximum Wavelength Distribution.jpg']);
saveas(gcf,[solvent ' Molecule Maximum Wavelength Distribution.fig']);
disp('Save Molecule Maximum Wavelenth successfully /n')


