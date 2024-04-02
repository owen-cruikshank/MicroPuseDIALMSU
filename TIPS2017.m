function Q= TIPS2017(molinput,isoinput,Tinput)
%File: TIPS2017.m
%Author: Owen Cruikshank
%Inputs:
%   -molinput: scalar integer HITRAN molecular identifier
%   -isoinput: scalar integer HITRAN isotope identifier
%   -Tinput:[K] (range x time) temperature input
%
%Outputs:
%   -Q: Total internal partition function
%   -Q296: Total internal partition function at HITRAN reference
%   temperature
%Description:
%Calculate total internal partition function from HITRAN 2017 files. Using
%HITRAN python function. Files are stored in QTpy folder. 
%
%
%
%%%%%%%TO AVOID PICKLING ERROR%%%%%%%%%: files must be re-saved as CR(Unix) from
%CRLF(Windows) if using a windows system.

Tinput = fillmissing(Tinput,'nearest',1);
Tinput = fillmissing(Tinput,'nearest',2);


Q = double(pyrunfile('TIPS_2017_OwenArray.py','QT',mol=int32(molinput),iso=int32(isoinput),T = py.numpy.array(Tinput)));

