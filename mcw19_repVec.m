function mat = mcw19_frepvec(vec,dims)
%      if size(vec,1) == 1
%          vec = vec';
%      end
l = length(vec);

N = prod(dims);


for n = 1: length(dims)
    if n == 1
        s = ['[' num2str(dims(n))];
    else
        s = [s ',' num2str(dims(n))];
    end
end
s = [s ']);'];

s1 = [sprintf('mat = reshape(repmat(vec,%d,1),', N/l) s];
eval([s1])



end



