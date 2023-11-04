function Q = GS_orthonorm(A)

% Orthonormalize the columns of A to produce Q using Gram-Schmidt. 
[m,n]=size(A);
Q=zeros(m,n);
R=zeros(m,n);

for j=1:n
    v=A(:,j);
    for ii=1:j-1
        R(ii,j)=Q(:,ii)'*A(:,j);
        v = v-R(ii,j)*Q(:,ii);
    end
    R(j,j) = norm(v);
    Q(:,j) = v/R(j,j);
end

% Q is now the orthonormalized matrix

