codefolder=pwd;
%For my file need to change Generalname etc.;

    Ex=actxserver('Excel.Application');
    Exwokbook=Ex.Workbooks.Add;
    Ex.Visible=1;
    Exsheets=Ex.ActiveWorkbook.Sheets;
    Firstsheet=Exsheets.get('Item',1);
    Firstsheet.Activate
    %write in folder names
    Folderrange=get(Ex.ActiveSheet,'Range','A1:A94');
    dpi = get(groot, 'ScreenPixelsPerInch');  % Get screen dpi
    picrange=get(Ex.ActiveSheet,'Range','A1:CU94');
    picrange.ColumnWidth=50;
    picrange.RowHeight=100;
    n=0;
    
    xletter=cell(99,1);
    for i=1:99
        if i<=26
            xletter{i,1}=char(64+i);
        else
            xletter{i,1}=[char(65+floor((i-27)/26)) char(65+rem(i-1,26))];
        end
    end
%%
for mol=1:94
for sec=1:99
[A.eff,A.eff_fit,~,A.numst,~]=Traceson(spectraspectrum(:,sec,mol),codefolder);
try
    figure;plot(A.eff);hold on;plot(A.eff_fit(A.numst,:))
catch
    disp('No intensity found')
end
%save the plot
n=n+1;
print(gcf, sprintf('-r%d', dpi),'-clipboard', '-dbitmap');
pause (0.5);
Ex.ActiveSheet.Range([xletter{sec,1} num2str(mol)]).PasteSpecial();
Ex.ActiveSheet.Shapes.Item(n).LockAspectRatio='msoFalse';            %Unlocking aspect ratio
Ex.ActiveSheet.Shapes.Item(n).Width=Ex.ActiveSheet.Range([xletter{sec,1} num2str(mol)]).Width;    %Adjusting width
Ex.ActiveSheet.Shapes.Item(n).Height=Ex.ActiveSheet.Range([xletter{sec,1} num2str(mol)]).Height;  %Adjusting height
Ex.ActiveSheet.Shapes.Item(n).Placement='xlMoveandSize';
close(gcf);

end
end
%%
SaveAs(Exwokbook,'All spectra.xlsx');
Close(Exwokbook)
Quit(Ex)
