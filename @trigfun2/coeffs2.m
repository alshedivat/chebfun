function varargout = coeffs2( f ) 
% COEFFS2   Bivariate Fourier expansion coefficients of f. 

if ( isempty( f ) )
    varargout = { [ ] }; 
    return
end

if ( iszero( f ) ) 
    varargout = { 0 } ; 
    return
end

% CDR factorization: 
[C, D, R] = cdr( f ); 

if ( nargout <= 1 )
    % Return the matrix of coefficients
    varargout = { C.coeffs * D * R.coeffs.' }; 
elseif ( nargout <= 3 )
    varargout = { C.coeffs, D, R.coeffs };
else
    % Two output variables are not allowed.
    error('CHEBFUN:TRIGFUN2:coeffs2:outputs', ...
        'Incorrect number of outputs.');
end

end