pro makefilenums, startnum, nfiles, hwp0FileNums, hwp225FileNums, hwp45FileNums, hwp675FileNums

allnums=indgen(nfiles)+startnum

h0i = indgen(nfiles/4)*4
hwp0FileNums=allnums[h0i]
h225i = indgen(nfiles/4)*4+1
hwp225FileNums=allnums[h225i]
h45i = indgen(nfiles/4)*4+2
hwp45FileNums=allnums[h45i]
h675i = indgen(nfiles/4)*4+3
hwp675FileNums=allnums[h675i]


end
