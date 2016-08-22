;######## RELEASE CLIP HERE ######## 
;This converts the fairly exhaustive and complicated information in
;Mike's bispectrum save file into commonly used quantities
;It returns a structure:
;  {v2:fltarr(n_baselines), v2corr:fltarr(n_baselines),  $
;     v2_err:fltarr(n_baselines), cp:fltarr(n_bispect),  cp_err:fltarr(n_bispect),  cp_cov:fltarr(n_bispect, n_bispect), $
;     amp:fltarr(n_baselines), amp_err:fltarr(n_baselines),  phase:fltarr(n_baselines),  phase_err:fltarr(n_baselines), bs_amp:fltarr(n_bispect)}
;Note that phase and phase_err are just placeholders in this routine -
;subtracing cal from source closure phase increases the stability of cp2phase
;
;; EDIT: ACC added bispectral amplitude output
;
;;NB Closure phase returned is in _radians_
;
;KEYWORDS:
; i_want_a_memory_leak : This has to be set if you want to use the DWB
;  method for calibrate_v2_cp. It involves a memory leak with pointers.

function get_used_quantities,  filename,  bad_baselines,  bad_bs, root_dir, $
                               nows=nows, apply_phscorr = apply_phscorr,  $
                               correction_const = correction_const, i_want_a_memory_leak=i_want_a_memory_leak
  

  if (keyword_set(nows) eq 0) then nows =  0
  if (keyword_set(diagonal_cov) eq 0) then diagonal_cov = 0
  if (keyword_set(correction_const) eq 0) then correction_const = 0.0
  if (keyword_set(apply_phscorr) eq 0) then apply_phscorr = 0.0

  restore,  filename
  mf_filestring = root_dir +'/templates/' + mf_file
  restore, mf_filestring

  
  ;; See Re-definition in calibrate_v2 for reasons why ph_all, etc are pointers.
  quan =  {v2:fltarr(n_baselines), v2corr:fltarr(n_baselines),  $
           v2_err:fltarr(n_baselines), cp:fltarr(n_bispect),  cp_err:fltarr(n_bispect), $ 
           cp_cov:fltarr(n_bispect, n_bispect), cp_all: Ptr_New( /Allocate_Heap ), bs_all: Ptr_New( /Allocate_Heap),$
           amp:fltarr(n_baselines), amp_err:fltarr(n_baselines),  phase:fltarr(n_baselines), $
           phase_err:fltarr(n_baselines), flag:replicate(0,n_baselines), bsflag:replicate(0,n_bispect),bs_amp:fltarr(n_bispect),bs_amp_err:fltarr(n_bispect) }
                                ;First, use vis variance and/or phase slopes to find the v2 correction
  quan.v2corr = correct_vis(v2, avar, u,v, correction_const, err_avar=err_avar, err_vis=sig,  nows = nows)
  if (apply_phscorr eq 1) then quan.v2corr=quan.v2corr*phs_v2corr
                                ;Get rid of bad baselines and bs points.
  if (bad_baselines[0] ne -1) then begin
     for j = 0,n_elements(bad_baselines)-1 do begin
        v2_cov[bad_baselines[j],bad_baselines[j]] = 1.0
        v2[bad_baselines[j]] = median(v2) 
     endfor
     for j = 0,n_elements(bad_bispect)-1 do begin
        bs_var[0,bad_bispect[j]] = 1.0
        bs[bad_bispect] = median(bs)
     endfor
  endif

  ;Save V^2 and closure phase
  quan.v2 = v2
  dummy =  cov2cor(v2_cov,  sig = sig)
  quan.v2_err = sig
  ;*(quan.ph_all) = ph_all
  quan.cp = atan(bs, /phase)
  if (keyword_set(bs_all) and keyword_set(i_want_a_memory_leak)) then begin
   *(quan.cp_all) = atan( bs_all, /phase )
   *(quan.bs_all) = bs_all
  endif
  quan.cp_err = reform( atan( sqrt(bs_var[1,*])/abs(bs), /phase ) )

  ;;save bs_amp and uncertainty
  quan.bs_amp=abs(bs)
  quan.bs_amp_err=sqrt(bs_var[0,*])

 ;Find V amplitude and phase
  bs2amp, mf_filestring, v2, v2_cov, bs, bs_var, bs_cov, amp_err=amp_err, amplitude=amplitude, /quiet
  quan.amp =  amplitude
  quan.amp_err =  amp_err
  if keyword_set(cp_cov) then if (cp_cov[0] ne -1) then quan.cp_cov = cp_cov else begin
     quan.cp_cov[bscov2bs_ix[0,*],bscov2bs_ix[1,*]] = bs_cov[1,*]/abs(bs[bscov2bs_ix[0,*]]*bs[bscov2bs_ix[1,*]]) 
     quan.cp_cov =  quan.cp_cov+transpose(quan.cp_cov) 
     in =  indgen(n_bispect)
     quan.cp_cov[in, in] =  quan.cp_err^2
  endelse
  
  ;Apply correction for simplicity
  quan.v2 =  quan.v2/quan.v2corr
  quan.v2_err =  quan.v2_err/quan.v2corr

  return,  quan

end
