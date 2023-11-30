function [ res ] = TEST_SUB_mybinornd( N, p ,iter)
    %Custom fast binomial 
    [row_cnt, col_cnt] = size(N);
    res = zeros(row_cnt, col_cnt,iter);
    for ii=1:row_cnt
       for jj=1:col_cnt
           if isnan(N(ii,jj))
               res(ii,jj,:)=NaN(1,1,iter);
           else
               res(ii, jj,:) = sum(rand(1,N(ii,jj),iter)<p,2,'omitnan');
           end
       end

%        for ii=1:row_cnt
%            %if isnan(N(ii,jj))
%             %   res(ii,:)=NaN;
%            %else
%            %N(ii,isnan(N(ii,:)))=0; 
%            Nmax = max(N(ii,:),[],'omitnan');
%            if isnan(Nmax)
%                res(ii,:) = nan;
%            else
%            prop=rand(col_cnt,Nmax)<p;
% 
%            for jj = 1:col_cnt
%                if isnan(N(ii,jj))
%                    res(ii,jj) = nan;
%                else
%            res(ii, jj) = sum(prop(jj,1:N(ii,jj)),2,'omitnan');
%                end
%            end
%            end
%            %end
%        %end
%     end
    end
end