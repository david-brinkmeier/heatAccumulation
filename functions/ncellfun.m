% src: https://www.mathworks.com/matlabcentral/fileexchange/50133-execute-cellfun-on-nested-cell-array?s_tid=blogs_rc_5
%Nested cell fun 
%I often find myself writing function like cellfun(@(cell1)
%cellfun(@(cell2), fun(cell2), cell1,'un',0), outcell,'un',0). This can be
%tedious and error prone. So I wrote this generale nested cell fun which
%allows you to apply a function to a nested cell. 
%EXAMPLE
%for i=1:3
%    for j=1:2
%     mat{i}{j}=rand(10,1);
%     end
% end
%out=ncellfun(@mean, mat)
%out=ncellfun(@cell2mat,  mat,1)
%If we want to calculate the mean for each subcell, we should use to nested
%cellfun (or, of course, a double loop). With my function you can write
%out=ncellfun(@mean, mat). The nested level of cell array if determined by
%the function itself. You can set your own level by an additional parameter
%n. ncellfun(@cell2mat,  mat,1)

function out = ncellfun(fun, cell,n)
if nargin<3
    temp=cell; n=0; 
    while iscell(temp)
        n=n+1; 
        temp=temp{1};
    end    
end
str=[];
for i=1:n
    str = [str 'cellfun(@(cell' num2str(i) ')'];
end
str=[str 'fun(cell' num2str(n) ')'];

for i=n-1:-1:1
    str = [str ', cell' num2str(i) ',''un'',0)'];
end
str=[str ', cell, ''un'',0);'];
out=eval(str);
end 


