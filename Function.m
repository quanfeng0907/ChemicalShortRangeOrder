function F = Function(P,N,c,A)
% calculation function value of nonlinear equation
%P: column vector, N: number of elements, c: elemental fraction, A: matrix
p=zeros(N);
% convert column vector to matrix
for i=1:length(P)
    temp=convert2(N,i);
    p(temp(1),temp(2))=P(i);
    p(temp(2),temp(1))=p(temp(1),temp(2));
end
for i=1:length(P)
    temp=convert2(N,i);
    F(i,1)=(c(temp(1))-sum(p(temp(1),:)))*(c(temp(2))-sum(p(temp(2),:)))-A(temp(1),temp(2))*p(temp(1),temp(2))^2;
end

end

