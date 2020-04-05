function [vec] = ablvec(N,n)

vec=zeros(N,1);
for nn=1:n
    x=(n-nn+1)/n;
    vec(nn)=x^2;
end
for nn=N-n+1:N
    x=(nn-(N-n))/n;
    vec(nn)=x^2;
end

vec=1-vec;