;This procedure consolidates a bad pixel map after several flats and
;darks have been made, each adding to the map

pro consolidate_badpix,  filename,  threshold = threshold

if (keyword_set(threshold) eq 0) then threshold = 2
bad =  readfits(filename)
bad =  (bad + threshold - 1) < 1 > 0
writefits,  filename,  bad

end
