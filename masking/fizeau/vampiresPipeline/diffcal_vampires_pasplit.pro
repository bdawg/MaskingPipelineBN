; This is a wrapper for diffcal_vampires, which will split
; observations up by position angle. Call it in the same way as
; diffcal_vampires, with the extra PAsplit value (size of PA bins)

pro diffcal_vampires_pasplit, fileprefix, startnum, nfiles, cubeinfofile, csvoutfile, fluxvec=fluxvec,$
                      yrange=yrange, nbootstraps=nbootstraps, docals=docals, saveeps=saveeps,$
                      PAsplit=PAsplit



; First, tabulate the PAs of all data files
h0_PA = []
h225_PA = []
h45_PA = []
h675_PA = []

restore,cubeinfofile
if max(olog.tsize) GT 0 then print, 'WARNING: Calibrator star data will be IGNORED'
makefilenums,startnum,nfiles,hwp0FileNums,hwp225FileNums,hwp45FileNums,hwp675FileNums
allPAs = olog.pa[where(olog.tsize LT 0)]
;allCSz = olog.cube_sz[where(olog.tsize LT 0),2]
allFrameNums = olog.frames[where(olog.tsize LT 0)]
allUTCs = olog.UTC[where(olog.tsize LT 0)]

for i = 0,n_elements(hwp0filenums)-1 do begin
   curfilenum = hwp0filenums[i]
   curOlogInd = (where(allFrameNums eq curfilenum))[0]
   curPA = allPAs[curOlogInd]
   ;curSz = allCSz[curOlogInd]
   ;curPAvec = replicate(curPA,curSz)
   h0_PA = [h0_PA, curPA]
end

for i = 0,n_elements(hwp225filenums)-1 do begin
   curfilenum = hwp225filenums[i]
   curOlogInd = (where(allFrameNums eq curfilenum))[0]
   curPA = allPAs[curOlogInd]
   ;curSz = allCSz[curOlogInd]
   ;curPAvec = replicate(curPA,curSz)
   h225_PA = [h225_PA, curPA]
end

for i = 0,n_elements(hwp45filenums)-1 do begin
   curfilenum = hwp45filenums[i]
   curOlogInd = (where(allFrameNums eq curfilenum))[0]
   curPA = allPAs[curOlogInd]
   ;curSz = allCSz[curOlogInd]
   ;curPAvec = replicate(curPA,curSz)
   h45_PA = [h45_PA, curPA]
end

for i = 0,n_elements(hwp675filenums)-1 do begin
   curfilenum = hwp675filenums[i]
   curOlogInd = (where(allFrameNums eq curfilenum))[0]
   curPA = allPAs[curOlogInd]
   ;curSz = allCSz[curOlogInd]
   ;curPAvec = replicate(curPA,curSz)
   h675_PA = [h675_PA, curPA]
end

; Make PA positive
h0_PA[where(h0_PA LT 0)] = h0_PA[where(h0_PA LT 0)] + 360
h225_PA[where(h225_PA LT 0)] = h225_PA[where(h225_PA LT 0)] + 360
h45_PA[where(h45_PA LT 0)] = h45_PA[where(h45_PA LT 0)] + 360
h675_PA[where(h675_PA LT 0)] = h675_PA[where(h675_PA LT 0)] + 360


; Sort into PA bins
; Only use PA from hwp 0 and reuse for other HWP angles, to keep
; appropriate data together.
nMeas = n_elements(h0_PA)
paBin = lonarr(nMeas)
paBinStart = []
paBinEnd = []
newPABinIndex = [] ;Actually used for chopping up data
curBin = 0
for i = 0,nMeas-1 do begin   
   if i eq 0 then begin
      paBinStart = [paBinStart, h0_PA[i]]
      newPABinIndex = [newPABinIndex, i]
   end

   if h0_PA[i] - paBinStart[curBin] GT PASplit then begin
      paBinEnd = [paBinEnd, h0_PA[i]]
      paBinStart = [paBinStart, h0_PA[i]]
      newPABinIndex = [newPABinIndex, i]
      curBin = curBin +1
   end

   paBin[i] = curBin
end
paBinEnd = [paBinEnd, paBinStart[-1]+PASplit]
print,[strn(n_elements(paBinStart))+' PA bins']


; Make some plots
loadct,39,/sil
plot,h0_PA,ys=1
oplot,h0_PA,psym=4
for k = 0,n_elements(paBinStart)-1 do begin
   oplot,[newPABinIndex[k],newPABinIndex[k]],[-360,360],color=200,linestyle=2
end


; Now call diffcal_vampires for each set of data
for k = 0,n_elements(paBinStart)-1 do begin
   binLabel = ['PABinNum'+strn(k)+'_'+strn(paBinStart[k],len=5)+'-'+strn(paBinEnd[k],len=5)+'_']
   curStartnum = newPABinIndex[k]*4

   if k LT n_elements(paBinStart)-1 then begin  
      curNumFiles = (newPABinIndex[k+1] - newPABinIndex[k])*4
   endif else begin
      curNumFiles = nfiles - newPABinIndex[k]*4
   endelse

   diffcal_vampires, fileprefix, curStartnum, curNumFiles, cubeinfofile, csvoutfile, $
                  yrange=yrange, nbootstraps=nbootstraps, docals=docals, $
                  saveeps=saveeps,output_special=binLabel
end


end
