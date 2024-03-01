function [a_0] = alpha_0(N_on,N_off,dr)

a_0 = zeros(size(N_on));
a_0(1,:,:) = log((N_on(1,:).*N_off(2,:))./(N_on(2,:).*N_off(1,:)))./2./dr;
for iii = 2:size(N_on,1)-1
    a_0(iii,:)=log((N_on(iii-1,:).*N_off(iii+1,:))./(N_on(iii+1,:).*N_off(iii-1,:)))./2./dr./2;
end
a_0(end,:,:) = log((N_on(end-1,:).*N_off(end,:))./(N_on(end,:).*N_off(end-1,:)))./2./dr;
