function f = convert2(N,N1)
% N number of elements, N1 serial number
    for i=1:N-1
        for j=i+1:N
            temp=((2*N+1)*i-i^2)/2-N+j-i;
            if temp==N1
                f(1)=i;
                f(2)=j;
                break;
            end
        end
    end  

end

