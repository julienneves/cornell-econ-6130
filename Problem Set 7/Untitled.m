reshape((eye(4)-kron(A,A'))^-1*reshape(C*C',4,1),2,2)