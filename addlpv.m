% Add LPV toolbox to path

format compact

cm = computer;
if cm(1) == 'M' ||  cm(1)=='G'
    addpath(pwd);
    addpath([pwd '/LPVutil']);
    addpath([pwd '/LPVutil/IQCfiles']);
    addpath([pwd '/Documentation']); 
elseif cm(1) == 'P'
    addpath(pwd);
    addpath([pwd '\LPVutil']);
    addpath([pwd '\LPVutil\IQCfiles']);
    addpath([pwd '\Documentation']);  
end
    

