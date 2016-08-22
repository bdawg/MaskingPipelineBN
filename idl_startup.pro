; Setup Path Names

code_path="~/code_svn/"
susi_data_path="~/susi/data/"
susi_root_dir="~/susi/"
palomar_data_path="~/palomar/data/"
palomar_root_dir="~/palomar/"
nirc_data_path="~/nirc/data/"
nirc_root_dir="~snert/nirc/"
nirc2_data_path="~/nirc2/data/"
nirc2_root_dir="~/nirc2/"
conica_data_path="~/conica/data/"
conica_root_dir="~/conica/"
;-----------------------------
default_path =!path

!path=code_path + "masking/misc:" + !path
!path=code_path + "simulation:" + !path
!path=code_path + "astrolib/pro:" + !path
!path=!path + ":" + code_path + "susi/util"
!path=!path + ":" + code_path + "susi/analysis"
!path=!path + ":" + code_path + "masking/fizeau"
!path=!path + ":" + code_path + "masking/fizeau/nirc"
!path=!path + ":" + code_path + "masking/fizeau/nirc2"
!path=!path + ":" + code_path + "masking/fizeau/lws"
!path=!path + ":" + code_path + "masking/fizeau/pharo"
!path=!path + ":" + code_path + "masking/fizeau/conica"
!path=!path + ":" + code_path + "masking/fizeau/vampires"
!path=!path + ":" + code_path + "masking/fizeau/vampiresPipeline"
!path=!path + ":" + code_path + "masking/fizeau/trecs"
!path=!path + ":" + code_path + "masking/util"
!path=!path + ":" + code_path + "oidata"
!path=!path + ":" + code_path + "scholz_code"
!path=!path + ":" + nirc_root_dir + "Pretty/miras/analysis_scripts"
!path=!path + ":" + "/import/spiral1/snert/nirc/newmasks/LBT/special_code"
!path=!path + ":" + code_path + "masking/fizeau/lmircam"
!path=!path + ":" + code_path

!edit_input = 1500
device,retain=2
device, decomposed=0
;
; Enter Any Constants Here.
;


print, 'BNs SVN startup complete.'
