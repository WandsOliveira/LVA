function [abscissa,DATA]=readuff(name)

%READUFF Converts and returns data from uff file
%[abscissa,DATA]=readuff(name)
%abscissa : x-axis
%DATA : y-axis
%name : file name WITH .uff at the end.

DATAid=fopen(name);
[~]=textscan(DATAid,'%d',1);
[~]=textscan(DATAid,'%d',1);
[~]=textscan(DATAid,'%s %*[^\n]',1);
[~]=textscan(DATAid,'%s %*[^\n]',1);
[~]=textscan(DATAid,'%s %*[^\n]',1);
[~]=textscan(DATAid,'%s',1);
[~]=textscan(DATAid,'%s',1);

functiontype=textscan(DATAid,'%f %*[^\n]',1);
abscissatype=textscan(DATAid,'%f %f %*f %f %f %*[^\n]',1);
abscissaunits=textscan(DATAid,'%f %*[^\n]',1);

[~]=textscan(DATAid,'%d %*[^\n]',1);
[~]=textscan(DATAid,'%d %*[^\n]',1);
[~]=textscan(DATAid,'%d %*[^\n]',1);

% Np=abscissatype{2};
abscissa=(abscissatype{3}:abscissatype{4}:abscissatype{4}*(abscissatype{2}-1))';

if abscissatype{1}==2||abscissatype{1}==4

    DATA=textscan(DATAid,'%f %f %f %f %f %f',abscissatype{2}/6);
    
    fclose(DATAid);
    
    DATA=cell2mat(DATA);
    nc=length(DATA(:,1));

    DATA=mat2cell(DATA,ones(nc,1),6);

    for i=1:nc

        DATA{i}=(DATA{i})';
    
    end

    DATA=cell2mat(DATA);
    DATA=DATA(1:abscissatype{2});
%     save('expdata.mat','abscissa','DATA')
    
else 
    
    DATA=textscan(DATAid,'%f %f %f %f %f %f',abscissatype{2}/3);
    
    fclose(DATAid);
    
    DATAreal=[DATA{1} DATA{3} DATA{5}];
    DATAimag=[DATA{2} DATA{4} DATA{6}];

    DATAreal=mat2cell(DATAreal,ones(length(DATA{1}),1),3);
    DATAimag=mat2cell(DATAimag,ones(length(DATA{1}),1),3);

    for i=1:length(DATA{1})

        DATAreal{i}=(DATAreal{i})';
        DATAimag{i}=(DATAimag{i})';
    
    end

    DATAreal=cell2mat(DATAreal);
    DATAreal=DATAreal(1:abscissatype{2});
    DATAimag=cell2mat(DATAimag);
    DATAimag=DATAimag(1:abscissatype{2});
    DATA=DATAreal+DATAimag*1j;

%     save('expdata.mat','abscissa','DATAcomplex')
    
end

end
