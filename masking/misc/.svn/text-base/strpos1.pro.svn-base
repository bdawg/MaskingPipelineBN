From lkramer@neosoft.com Mon Jul 17 13:38:58 PDT 1995
Article: 5024 of comp.lang.idl-pvwave
Path: agate!howland.reston.ans.net!swrinde!news.uh.edu!uuneo.neosoft.com!sam-slip-h5.NeoSoft.COM!lkramer
From: lkramer@neosoft.com (Leonard Kramer)
Newsgroups: comp.lang.idl-pvwave
Subject: Re: (NOT SO) Dumb problem with STRPOS
Date: Sun, 16 Jul 1995 16:19:41 UNDEFINED
Organization: NeoSoft Internet Services   +1 713 968 5800
Lines: 97
Distribution: world
Message-ID: <lkramer.9.0175438D@neosoft.com>
References: <3u6g67$usi@ds2.acs.ucalgary.ca>
NNTP-Posting-Host: sam-slip-h5.neosoft.com
Keywords: substring
X-Newsreader: Trumpet for Windows [Version 1.0 Rev B final beta #1]


>From: steele@aragorn.phys.ucalgary.ca (Dave Steele)
>Newsgroups: comp.lang.idl-pvwave
>Subject: Dumb problem with STRPOS
>Date: 14 Jul 1995 19:22:47 GMT
>Organization: Institute for Space Research

>Can someone explain the following behaviour to me?

>IDL> help,/st,!version                  ; Here's the situation.
>** Structure !VERSION, 3 tags, length=24:
>   ARCH            STRING    'mipseb'
>   OS              STRING    'RISC/os'
>   RELEASE         STRING    '3.1.0'
>IDL> test='c000.'                       ; define a string
>IDL> help,test                          ; check it is a string
>TEST            STRING    = 'c000.'
>IDL> print,strpos(test,'000')           ; search for a substring
>       1
>IDL> print,strpos(test,'000.')          ; and a slightly different one
>       1
>IDL> test='c0000.'                      ; add a '0'
>IDL> help,test
>TEST            STRING    = 'c0000.'
>IDL> print,strpos(test,'000')           ; search again
>       1
>IDL> print,strpos(test,'000.')          ; Why isn't this substring found???
>      -1

>I checked the manual for this version and it says nothing that would
>explain this.  BTW, it's nothing to do with octal constants.
>All suggestions gratefully received!

Dave:  

This is a definite Bug.   I am using Window 3.6.1.   Get identical response.  
You have probably already found out that:

IDL> print,test
c0000.
IDL> print,strpos(test,'0000.')	; a longer search string.
       1
IDL> print,strpos(test,'000.')	; your search string.
      -1
IDL> print,strpos(test,'00.')           ; a shorter search string.
       3

and also

IDL> print,test
c0000X
IDL> print,strpos(test,'000X')    ;  different terminal character.
      -1

This is a very scary bug.   I can't find the source for strpos.pro in the 
library so it must be part of the dynamic run executable. 
    
I re-wrote  strpos calling it strpos1.    It performs as follows:

IDL> print,test
c0000.
IDL> print,strpos1(test,'0000.')
       1
IDL> print,strpos1(test,'000.')
       2
IDL> print,strpos1(test,'00.')
       3
IDL> print,strpos1(test,'0.')
       4

and here is the text of the routine.
---------------------------cut here-------------------------
function strpos1,targ,srch,pos
   targ=string(targ)
   srch=string(srch)
   if (n_elements(pos) eq 0) then pos = 0
   t=byte(targ)
   n=n_elements(t)
   m=strlen(srch)
   p=lindgen(m)
   i=pos
   for i=pos,n do begin
       if (srch eq string(t(p+i)))then return,i
   endfor
   return ,-1
end
--------------------------end of routine ------------------
I do not have time to properly document this.  RSI really should fix 
it.

E-mail me if you have questions.  I tried to reply directly but the mail was
returned (?).  
---------------------------------------------------------------------
Leonard Kramer, Ph.D.  Physicist.
AGAR Corporation, (Process Measurement & Control)
POB 802127 Houston, TX  77280-2127
email:  lkramer@neosoft.com


