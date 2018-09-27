function [M] = constructMean(sampind)
    %construct mean aggregation matrix
    N = sum(sampind);
    regionnum = length(sampind);
    
    meanvector = 1./sampind;
    ivec = zeros(N,1);
    jvec = zeros(N,1);
    vvec = zeros(N,1);
    ind = 1;
    for i = 1:regionnum;
        startpoint = 1 + (i-1)*sampind(i);
        endpoint = sum(sampind(1:i));
        for j = startpoint:endpoint
            ivec(ind,1) = i;
            jvec(ind,1) = j;
            vvec(ind,1) = meanvector(i);
            ind = ind + 1;
        end
    end
    M = sparse(ivec,jvec,vvec,regionnum,N);
end