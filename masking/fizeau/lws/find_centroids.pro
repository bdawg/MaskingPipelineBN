;This procedure finds the centroids for a bunch of LWS 
FUNCTION cent_func, P  
COMMON FUNC_XY, f_mask,  f_frame
RETURN, -total(shift(f_mask, p[0], p[1])*f_frame)  
END 

pro find_centroids,  cleanframes,  split_cubes,  pattern,  nod_ims = nod_ims,  nod_cubes = nod_cubes,  centroids = centroids
 COMMON FUNC_XY, f_mask,  f_frame
 pscale =  0.08
 if(pattern eq 7) then pxy=[[0,-5.],[0.,0.],[-5.,0],[-5.,-5.]] + 2.5 $
 else if (pattern eq 3) then pxy=[[-5.0, -5.0],[0.0, -5.0],[0.0, 0.0],[-5.0,0.0]]+2.5 $
 else if (pattern eq 6) then pxy=[[0.0, -5.0],[-5.0, -5.0],[-5.0,0.0],[0.0, 0.0]]+2.5 $
 else stop
 nsubarr =  (size(pxy))[2]
 pxy =  [[pxy], [0, 0]]
 pxy[0,*]=-pxy[0,*]
 pxy =  pxy/pscale + 64
 window =  exp(-(dist(128)/20.)^4)
 f_mask =  fltarr(128, 128)
 for i =  0, nsubarr-1 do $
  f_mask +=  shift(window, pxy[0, i], pxy[1, i]) 
 nframes =  (size(cleanframes))[3]
 ;Enforce nnodframes <= nframes, so that shifting works.
 if (keyword_set(nod_ims)) then nnodframes =  (size(nod_ims))[3] <  nframes
 centroids = fltarr(2, nframes)
 mn_frame = fltarr(128,128)
 print, 'Finding centroids...'
 for i =  0, nframes-1 do begin
  f_frame =  cleanframes[*, *, i]
  centroids[*, i]=  amoeba(1e-3, function_name = 'cent_func',  p0 = [0, 0],  scale = [6, 6])
  mn_frame += shift(cleanframes[*,*,i],-centroids[0,i], -centroids[1,i])
 endfor
 ;Now we have a mean frame, use this for a new optimization
 f_mask =  mn_frame
 for i =  0, nframes-1 do begin
  f_frame =  cleanframes[*, *, i]
  centroids[*, i]=  amoeba(1e-3, function_name = 'cent_func',  p0 = centroids[*, i],  scale = [3, 3])
  mn_frame += shift(cleanframes[*,*,i],-centroids[0,i], -centroids[1,i])
 endfor
 print,  'Splitting cubes...'
 split_cubes =  fltarr(64, 64, nframes,  nsubarr+1)
 if (keyword_set(nod_ims)) then nod_cubes =  fltarr(64, 64, nnodframes,  nsubarr+1)
 ;OK, first find an offset for each subarray...
 make_2d,findgen(64)-32,   findgen(64)-32,  xx,  yy
 window =   shift(exp(-(dist(64)/25.)^4), 32, 32)
 subcent =  fltarr(2, nframes)
 for i =  0, nsubarr do begin
  if (i ne nsubarr) then begin
   temp =  (shift(mn_frame,  -pxy[0, i]+32, -pxy[1, i]+32) )[0:63, 0:63]
   temp =  smooth(temp, 6)*window
   dummy =  max(temp,  m)
   m =  array_indices(temp, m)
   subcent[0, *] =  centroids[0, *] + m[0]-32
   subcent[1, *] =  centroids[1, *] + m[1]-32
  endif else subcent =  centroids
  for j =  0, nframes-1 do $
   split_cubes[*, *, j, i] =  (shift(cleanframes[*, *, j],  -subcent[0,j]-pxy[0, i]+32, -subcent[1,j]-pxy[1, i]+32) )[0:63, 0:63]
  if (keyword_set(nod_ims)) then for j =  0, nnodframes-1 do $
   nod_cubes[*, *, j, i] =  (shift(nod_ims[*, *, j],  -subcent[0,j]-pxy[0, i]+32, -subcent[1,j]-pxy[1, i]+32) )[0:63, 0:63]
 endfor
 
end
