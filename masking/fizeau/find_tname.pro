;######## RELEASE CLIP HERE ######## 
;This procedure finds the Keck starlist target name for a given RA and
;Dec. It is important for Palomar, and has the added benefit that the
;target names will be the same as they are in the fits header from
;Keck.
;Inputs: RA and Dec in degrees.
;Optional: A Keck starlist name. Default is root_dir + 'fizeau/nirc2/'

function find_tname,  ra, dec, starlist = starlist,  max_dist = max_dist

if not keyword_set(max_dist) then max_dist = 120.
if not keyword_set(starlist) then starlist = !ROOT_DIR + 'fizeau/nirc2/Master_Keck_Catalog'
a = READSTARLIST(starlist)
GCIRC, 1, ra/15, dec, $
                  a.raj2000/15, a.dej2000, tsep
mdist = min(tsep,  ix)
if (mdist gt max_dist) then begin
 print,  'No target within ',  strtrim(max_dist, 2),  ' arccsec!'
 stop
 return,  ''
endif

return,  strtrim(a[ix].targ, 2)

end
