function Counts = poissonThin2(Counts,~,iter)
%File: poissonThin.m
%Date: 02/4/2022
%Author: Owen Cruikshank
%Inputs:
%   -Counts, Structure
%       -Counts.o2on (range x time): integrated and background subtracted
%       o2 online combined photon counts
%       -Counts.o2off (range x time): integrated and background subtracted
%       o2 offline combined photon counts
%       -Counts.o2on_mol (range x time): integrated and background subtracted
%       o2 online molecular photon counts
%       -Counts.o2off_mol (range x time): integrated and background subtracted
%       o2 offline molecular photon counts
%       -Counts.NBins (): number of raw data bins integrated into each bin
%       -Counts.bg_o2on (1 x time): background counts for o2on
%       -Counts.bg_o2off (1 x time): background counts for o2off
%       -Counts.bg_o2on_mol (1 x time): background counts for o2on_mol
%       -Counts.bg_o2off_mol (1 x time): background counts for o2off_mol
%
%Outputs:
%   -Counts, Structure
%       -Counts.o2on (range x time): optimally smoothed, integrated, and background subtracted
%       o2 online combined photon counts
%       -Counts.o2off (range x time): optimally smoothed, integrated, and background subtracted
%       o2 offline combined photon counts
%       -Counts.o2on_mol (range x time): optimally smoothed, integrated, and background subtracted
%       o2 online molecular photon counts
%       -Counts.o2off_mol (range x time): optimally smoothed, integrated, and background subtracted
%       o2 offline molecular photon counts
    %   -Counts.Poissonthin
    %         -Counts.Poissonthin.timeWidthon (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthoff (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthon_mol (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.timeWidthoff_mol (range x 1): Gaussian filter width as number of bins in time
    %         -Counts.Poissonthin.rangeWidthon (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthoff (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthon_mol (1 x time): Gaussian filter width as number
    %         of bins in range
    %         -Counts.Poissonthin.rangeWidthoff_mol (1 x time): Gaussian filter width as number
    %         of bins in range
    %         
    %         -Counts.Poissonthin.Ezon=Ezon;
    %         -Counts.Poissonthin.Eton=Eton;
    %         -Counts.Poissonthin.Ezoff=Ezoff;
    %         -Counts.Poissonthin.Etoff=Etoff;
    %         -Counts.Poissonthin.Ezon_mol=Ezon_mol;
    %         -Counts.Poissonthin.Eton_mol=Eton_mol;
    %         -Counts.Poissonthin.Ezoff_mol=Ezoff_mol;
    %         -Counts.Poissonthin.Etoff_mol=Etoff_mol;

disp('Thining profiles')

%==== This hopefully avoids aditional errors if MATLAB runs into errors ======
terminate(pyenv);
pyenv(ExecutionMode="OutOfProcess");

%eliminate <0
tic
o2on = Counts.o2on+Counts.bg_o2on;
o2on(o2on<0)=0;
o2off = Counts.o2off+Counts.bg_o2off;
o2off(o2off<0)=0;
o2on_mol = Counts.o2on_mol+Counts.bg_o2on_mol;
o2on_mol(o2on_mol<0)=0;
o2off_mol = Counts.o2off_mol+Counts.bg_o2off_mol;
o2off_mol(o2off_mol<0)=0;
wvon = Counts.wvon+Counts.bg_wvon;
wvon(wvon<0)=0;
wvoff = Counts.wvoff+Counts.bg_wvoff;
wvoff(wvoff<0)=0;


for i = 1:iter
    Counts.fon(:,:,i) =     double(pyrunfile('PoissThin.py','F',x=int32(o2on)));
    Counts.foff(:,:,i) =    double(pyrunfile('PoissThin.py','F',x=int32(o2off   )));
    Counts.fon_mol(:,:,i) = double(pyrunfile('PoissThin.py','F',x=int32(o2on_mol)));
    Counts.foff_mol(:,:,i) = double(pyrunfile('PoissThin.py','F',x=int32(o2off_mol)));
    Counts.fwvon(:,:,i) =   double(pyrunfile('PoissThin.py','F',x=int32(wvon)));
    Counts.fwvoff(:,:,i) =  double(pyrunfile('PoissThin.py','F',x=int32(wvoff)));
end
toc

gon = Counts.o2on+Counts.bg_o2on-Counts.fon;
goff = Counts.o2off+Counts.bg_o2off-Counts.foff;
gon_mol = Counts.o2on_mol+Counts.bg_o2on_mol-Counts.fon_mol;
goff_mol = Counts.o2off_mol+Counts.bg_o2off_mol-Counts.foff_mol;
gwvon = Counts.wvon+Counts.bg_wvon-Counts.fwvon;
gwvoff = Counts.wvoff+Counts.bg_wvoff-Counts.fwvoff;

%=====Find background of thinned profiles=====
% Counts.fon_bg = (Counts.bg_o2on.*Counts.NBins)/2;% Take mean of last data points
% Counts.fon_mol_bg = (Counts.bg_o2on_mol.*Counts.NBins)/2;% Take mean of last data points
% Counts.foff_bg = (Counts.bg_o2off.*Counts.NBins)/2;% Take mean of last data points
% Counts.foff_mol_bg = (Counts.bg_o2off_mol.*Counts.NBins)/2;% Take mean of last data points

Counts.fon_bg = (Counts.bg_o2on)/2;% Take mean of last data points
Counts.fon_mol_bg = (Counts.bg_o2on_mol)/2;% Take mean of last data points
Counts.foff_bg = (Counts.bg_o2off)/2;% Take mean of last data points
Counts.foff_mol_bg = (Counts.bg_o2off_mol)/2;% Take mean of last data points
Counts.fwvon_bg = (Counts.bg_wvon)/2;
Counts.fwvoff_bg = (Counts.bg_wvoff)/2;


Counts.fon = Counts.fon-Counts.fon_bg;
Counts.foff = Counts.foff-Counts.foff_bg;
Counts.fon_mol = Counts.fon_mol-Counts.fon_mol_bg;
Counts.foff_mol = Counts.foff_mol-Counts.foff_mol_bg;
Counts.fwvon = Counts.fwvon-Counts.fwvon_bg;
Counts.fwvoff = Counts.fwvoff-Counts.fwvoff_bg;

Counts.gon = gon-Counts.fon_bg;
Counts.goff = goff-Counts.foff_bg;
Counts.gon_mol = gon_mol-Counts.fon_mol_bg;
Counts.goff_mol = goff_mol-Counts.foff_mol_bg;
Counts.gwvon = gwvon-Counts.fwvon_bg;
Counts.gwvoff = gwvoff-Counts.fwvoff_bg;

end



