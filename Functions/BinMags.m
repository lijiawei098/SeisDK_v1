function nmag = BinMags(mag, bin0, binsz)
    % Calculate the bin edges
    binedges = bin0 - binsz/2 : binsz : 11 + binsz/2;    
    % Calculate the bin index for each magnitude
    [~, bin] = histc(mag, binedges);    
    % Calculate the new magnitudes using bin0 and binsz
    nmag = bin0 + (bin - 1) * binsz;    
    % Round the new magnitudes to one decimal place
    nmag = round(nmag, 1);
end