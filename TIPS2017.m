function Q = TIPS2017(molinput,isoinput,Tinput)
tic
Q =zeros(size(Tinput));
% for i = 1:size(Tinput,1)
%     for j = 1:size(Tinput,2)
%         Q(i,j) = double(pyrunfile('TIPS_2017_Owen.py','QT',mol=int32(molinput),iso=int32(isoinput),T = Tinput(i,j)));
%     end
% end
toc
tic
qt =zeros(size(Q));
for i = 1:size(Tinput,1)
    qt(i,:) = double(pyrunfile('TIPS_2017_OwenArray.py','QT',mol=int32(molinput),iso=int32(isoinput),T = py.numpy.array(Tinput(i,:))));
end
toc
Q = qt;