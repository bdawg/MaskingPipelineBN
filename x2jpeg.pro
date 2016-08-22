;Exports a window to a jpeg file

function x2jpeg,nwindow,outfile=outfile

	if (~keyword_set(outfile)) then outfile='out.jpg'
	
	prev=!d.window
	wset,nwindow
	wshow,nwindow
	
	rgb=tvrd(/true)
	write_jpeg,outfile,rgb,/true
	
	wset,prev
	
	return,1

end