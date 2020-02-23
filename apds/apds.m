%Deal with apd files to use the shape only
clearvars
codefolder=pwd;
solvent='F8T2N2';
srdir=['/scratch/lwang74/PTU_spectrum_lifetime_bluehive/PTUdata/' solvent];
srdir=['E:\F8T2400nmCH\apd full'];
cd (srdir)

allnames=struct2cell(dir( '*.mat'));
[~,len]=size(allnames);
for len_i=1:1:len
    clear name
    name=char(allnames(1,len_i));
    apdfile=[srdir '/' name];
    disp('Finish load file /n')
    
    cd(codefolder)
    [apddata,apddataresolution]=PTUim(apdfile);
    datasource=GetDandABS(apddata,channel,'M');
    absolutetime=datasource(:,1);
    dtime=datasource(:,2);

end
