PRO plothist, arr, xhist,yhist, BIN=bin, $
	BACKGROUND=back, CHANNEL=chan, CHARSIZE=chsiz, $
	CHARTHICK=chthck, COLOR=color, 	FONT=font, LINESTYLE=linest, $
	SUBTITLE=subtit, SYMSIZ=symsiz, THICK=thick, TICKLEN=ticklen, $
	TITLE=title, XCHARSIZE=xchsiz, XMARGIN=xmargn, XMINOR=xminor, $
	XRANGE=xrange, XSTYLE=xstyle, XTICKLEN=xtickln, XTICKNAME=xticknm, $
	XTICKS=xticks, XTICKV=xtickv, XTITLE=xtitle, XTYPE=xtype, $
	YCHARSIZE=ychsiz, YMARGIN=ymargn, YMINOR=yminor, $
	YRANGE=yrange, YSTYLE=ystyle, YTICKLEN=ytickln, YTICKNAME=yticknm, $
	YTICKS=yticks, YTICKV=ytickv, YTITLE=ytitle, YTYPE=ytype, noplot=noplot
;+
; NAME:
;	PLOTHIST
; PURPOSE:
;	Plot the histogram of an array with the corresponding abcissa.
;
; CALLING SEQUENCE:
;	plothist, arr, xhist, yhist, [, BIN=bin,   ... plotting keywords]
;
; INPUTS:
;	arr - The array to plot the histogram of.   It can include negative
;		values, but non-integral values will be truncated.              
;
; OPTIONAL OUTPUTS:
;	xhist - X vector used in making the plot  
;		( = indgen( N_elements(h)) * bin + min(arr) )
;	yhist - Y vector used in making the plot  (= histogram(arr/bin))
;
; OPTIONAL INPUT KEYWORDS:
;	BIN -  The size of each bin of the histogram,  scalar (not necessarily
;		integral).  If not present (or zero), the bin size is set to 1.
;
;		Any input keyword that can be supplied to the PLOT procedure
;		can also be supplied to PLOTHIST.
; EXAMPLE:
;	Create a vector of 1000 values derived from a gaussian of mean 0,
;	and sigma of 1.    Plot the histogram of these value with a bin
;	size of 0.1
;
;	IDL> a = randomn(seed,1000)
;	IDL> plothist,a, bin = 0.1
;
; MODIFICATION HISTORY:
;	Written     W. Landsman            January, 1991
;	Modified    J. Monnier		   Aug, 1996
;		    Added NOPLOT keyword.
;-
;			Check parameters.
On_error,2

if N_params() LT 1 then begin   
	print, 'Syntax - plothist, arr, [ xhist, yhist , BIN=bin]'
	return
endif

if N_elements( arr ) LT 2 then message, $
      'ERROR - Input array must contain at least 2 elements'
arrmin = min( arr, MAX = arrmax)
if ( arrmin EQ arrmax ) then message, $
       'ERROR - Input array must contain distinct values'
;
if not keyword_set(BIN) then bin = 1. else bin = float(abs(bin))
;
;			Compute the histogram and abcissa.
;
y = fix( ( arr / bin) - (arr LT 0) ) 
yhist = histogram( y )
N_hist = N_elements( yhist )
xhist = lindgen( N_hist ) * bin + min(y*bin) + bin/2.

if not keyword_set(PSYM) then psym = 10         ;Default histogram plotting

plot_keywords, $ 

	BACKGROUND=back, CHANNEL=chan, CHARSIZE=chsiz, $
	CHARTHICK=chthck, COLOR=color, 	FONT=font, LINESTYLE=linest, $
	SUBTITLE=subtit, SYMSIZ=symsiz, THICK=thick, TICKLEN=ticklen, $
	TITLE=title, XCHARSIZE=xchsiz, XMARGIN=xmargn, XMINOR=xminor, $
	XSTYLE=xstyle, XTICKLEN=xtickln, XTICKNAME=xticknm, $
	XTICKS=xticks, XTICKV=xtickv, XTITLE=xtitle, XTYPE=xtype, $
	YCHARSIZE=ychsiz, YMARGIN=ymargn, YMINOR=yminor, $
	YRANGE=yrange, YSTYLE=ystyle, YTICKLEN=ytickln, YTICKNAME=yticknm, $
	YTICKS=yticks, YTICKV=ytickv, YTITLE=ytitle, YTYPE=ytype

if not keyword_set(XRANGE) then xrange = [ xhist(0) ,xhist(N_hist-1) ]
if (keyword_set(noplot) eq 0) then begin
plot, [xhist(0) - bin, xhist, xhist(n_hist-1)+ bin] , [0,yhist,0],  $ 
        PSYM = psym, $ 
	BACKGROUND=back, CHANNEL=chan, CHARSIZE=chsiz, $
	CHARTHICK=chthck, COLOR=color, 	FONT=font, LINESTYLE=linest, $
	SUBTITLE=subtit, SYMSIZ=symsiz, THICK=thick, TICKLEN=ticklen, $
	TITLE=title, XCHARSIZE=xchsiz, XMARGIN=xmargn, XMINOR=xminor, $
	XRANGE=xrange, XSTYLE=xstyle, XTICKLEN=xtickln, XTICKNAME=xticknm, $
	XTICKS=xticks, XTICKV=xtickv, XTITLE=xtitle, XTYPE=xtype, $
	YCHARSIZE=ychsiz, YMARGIN=ymargn, YMINOR=yminor, $
	YRANGE=yrange, YSTYLE=ystyle, YTICKLEN=ytickln, YTICKNAME=yticknm, $
	YTICKS=yticks, YTICKV=ytickv, YTITLE=ytitle, YTYPE=ytype
endif
return
end
