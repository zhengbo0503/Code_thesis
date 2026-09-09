function A = get_testmatrix_anymatrix( id )
if id == 1 
    A = anymatrix('regtools/blur', 10, 5, 1.4);
    A = full(A);
elseif id == 2
    A = anymatrix("gallery/kahan", 50, 1e-2);
elseif id == 3
    A = anymatrix('nessie/whiskycorr'); 
elseif id == 4
    A = anymatrix('gallery/lauchli', 500, 1e-3);
    A = A'*A;
end 
end