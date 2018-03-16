function tens = tensorprod(tens1,indices1,tens2,indices2)
%Takes two tensors and performs the inner product over indices in ax1 and
%indices in ax2. Negative indices are summed over

inner1 = find(indices1 < 0);
inner2 = find(indices2 < 0);

assert(numel(inner1) == numel(inner2),'Must be equal number of contractions');

outer1 = 1:numel(indices1); outer1(inner1) = [];
outer2 = 1:numel(indices2); outer2(inner2) = [];

[~,indOut] = sort([indices1(outer1),indices2(outer2)]);

[s1,ind1] = sort(indices1(inner1));
[s2,ind2] = sort(indices2(inner2));

assert(all(s1==s2),'Inner indices must match in pairs');

dim1 = size(tens1); dim2 = size(tens2);
dim1((numel(dim1)+1):numel(indices1)) = 1;
dim2((numel(dim2)+1):numel(indices2)) = 1;

assert(all(dim1(inner1(ind1)) == dim2(inner2(ind2))),'Contracted dimensions must match');

tens1 = reshape(permute(tens1,[outer1,inner1(ind1)]),prod(dim1(outer1)),prod(dim1(inner1)));
tens2 = reshape(permute(tens2,[inner2(ind2),outer2]),prod(dim2(inner2)),prod(dim2(outer2)));

p = tens1 * tens2;

tens = permute(reshape(p,[dim1(outer1),dim2(outer2)]),indOut);

end

