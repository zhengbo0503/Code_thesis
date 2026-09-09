function value = scond_columns_complex(A)
%SCOND_COLUMNS_COMPLEX - Two-norm condition number after column equilibration
%
%   Usage:
%       value = scond_columns_complex(A)
%
%   Purpose:
%       SCOND_COLUMNS_COMPLEX scales each nonzero column of A to unit
%       two-norm and computes the two-norm condition number of the scaled
%       matrix. It implements the column-scaling branch of scond.m without
%       forming the mathematically redundant left identity product D1*A.
%       Scaling factors are converted to double to retain the numerical
%       convention of scond.m.
%
%   Arguments:
%       (1) A - Real or complex floating-point or multiprecision matrix.
%
%   Outputs:
%       (1) value - Real scalar containing the scaled condition number.
%           The value is Inf if A has a zero column.
%

scaledA = A;
for j = 1:size(A,2)
    columnNorm = norm(A(:,j),2);
    if columnNorm == 0
        value = Inf;
        return
    end
    columnScale = double(columnNorm^(-1));
    scaledA(:,j) = A(:,j)*columnScale;
end
value = cond(scaledA,2);
end
