wdir = './results/2_GJ802/'
objName = 'source_9h_H.fits'
calName = 'calib_9h_H.fits'

obj = readfits(wdir+objName, hdrObj)
cal = readfits(wdir+calName, hdrCal)

seuilS = 1e-6
seuilC = 3e-6

;; ----------------------------------------------------------------------------
;; source
;; ----------------------------------------------------------------------------
x   = (size(obj))[1]
y   = (size(obj))[2]
nbl = (size(obj))[3]

totObj = fltarr(x,y)
han = hanning(x,y)

for i = 0, nbl - 1 do begin
  a = shift(obj[*,*,i]*han, x/2, y/2)
  b = abs(shift(fft(a, -1), x/2, y/2 ))^2
  integ = total(b)
  b = b / integ
  contr = total(b^2)
  totObj = totObj + contr * b
endfor

totObj = totObj / max(totObj)

window, /free, xsize=2*x, ysize=2*y, title='source'
tvscl, congrid(alog10(totObj > seuilS), 2*x, 2*y)

;; ----------------------------------------------------------------------------
;; le calibrateur 
;; ----------------------------------------------------------------------------
x   = (size(cal))[1]
y   = (size(cal))[2]
nbl = (size(cal))[3] 

han   = hanning(x,y)
total = fltarr(x,y)

for i = 0, nbl - 1 do begin
  a = shift(cal[*,*,i]*han, x/2, y/2)
  b = abs(shift(fft(a, -1), x/2, y/2 ))^2
  integ = total(b)
  b = b / integ
  contr = total(b^2)
  total = total + contr * b
endfor
total = total / max(total)

window, /free, xsize=2*x, ysize=2*y, title='calibrator'
tvscl, congrid(alog10(total > seuilC), 2*x, 2*y)

;; ----------------------------------------------------------------------------
;; ----------------------------------------------------------------------------
module = fltarr(x,y)

gicle = where(totObj lt seuilS)
;;gicle = where(total lt 2e-6)
totObj[gicle] = 0.0
gicle = where(total lt seuilC)
totObj[gicle] = 0.0

module = sqrt((totObj)/(total > 1e-5))
window, /free, xsize=2*x, ysize=2*y, title='visibility'
tvscl, hist_equal(congrid((module), 2*x, 2*y))

fits_write, wdir+'visibility.fits', module, hdrObj
end

