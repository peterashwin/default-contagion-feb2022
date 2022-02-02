function x = imag_to_nan(x)
%imag_to_nan changes the values of x with imaginary parts to NANs

x(imag(x)~=0) = nan;
end

