;This function converts a unix time to a UTC time for the current day

function unix2utc, time

jd = julday(1,1,1970,0,0,time)
caldat, jd, month,d,y,h,min,s
jd0 = julday(month,d,y,0,0,0)
utc = (jd-jd0)*60.0*60.0*24.0

return, utc

end
