;Exports a window to a gif file

function x2gif,nwindow,outfile=outfile

	if (~keyword_set(outfile)) then outfile='out.jpg'
	
	prev=!d.window
	wset,nwindow
	wshow,nwindow
	
	rgb=tvrd(/true)
	write_gif,outfile,rgb,/true
	
	wset,prev
	
	return,1

end