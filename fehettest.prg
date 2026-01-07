' Perform Modified-Wald test for groupwise heteroskedasticity
' @addinname Modified Wald Test
' @author Manh Hoang Ba 
' @category Panel Data
' @description testing for heteroskedasticity in fixed-effects model
' @version 1.0: to performs Modified-Wald test for large T panel
' @version 1.1: to performs Kezdi test for large N and small T panel


!debug = 1
if !debug = 0 then
	logmode +all
endif

'check workfile is panel
if @ispanel=0 then
	@uiprompt("Procedure can only be run in panel workfiles")
	stop
endif

'check that an object exists
%type = @getthistype
if %type="NONE" then
	@uiprompt("No object found, please open an Equation object")
	stop
endif
'check that _this object is an equation
if %type<>"EQUATION" then
	@uiprompt("Procedure can only be run from an Equation object")
	stop
endif

'check that equation was estimated by LS, and by list.
%method = _this.@method
if %method <> "LS" then
	@uiprompt("Equations must be estimated by Least Squares")
	stop
endif
if _this.@bylist = 0 then
	@uiprompt("Equations must be estimated by a varlist")
	stop
endif

'check equation was estimated with cross-section fixed-effects
%ops = _this.@options
if %ops <> "CX=F" or  %ops = "" then
	@uiprompt("Estimation must use cross-sectional fixed effects")
	stop
endif

'choose Modified-Wald or Kezdi
!MWald = 0
!Kezdi = 0
!result = @uidialog("caption","Test for heteroskedasticity", _
	"text", "Please select one of the two tests below:", _
	"check", !MWald, "Modified Wald test (large-T panels)", _
	"check", !Kezdi, "Kezdi test (large-N, small-T panels)")
if !result = -1 then
	stop
endif

if !MWald+ !Kezdi = 0 or  !MWald+ !Kezdi = 2 then
	@uiprompt("Please select only one test.")
	stop
endif

'=========================================
'		Large-N, small-T: Kezdi (2003) test
'=========================================
if !Kezdi = 1 then

	%smpl 	= @pagesmpl
	!n 			= _this.@regobs
	!n_g 		= _this.@ncross
	
	smpl @all
	series _id = @crossid
	!_N_m = @max(_id)
	smpl %smpl
	
	%eit 		= @getnextname("_eit")
	_this.makeresids {%eit}
	
	%eit_sq	= @getnextname("_eit_sq")
	series {%eit_sq} = {%eit}*{%eit}/(!n - !n_g)
	!s2		= @sum({%eit_sq})
	
	%Ti		= @getnextname("_Ti")
	series {%Ti} = @obsby({%eit}, @crossid)
	
	!Tmax	= @max({%Ti})
	!Tmin		= @min({%Ti})
	if !Tmax 	= !Tmin then
		!bal	= 1
	else 
		!bal	= 0
	endif
	
	'series _ti = @date - @min(@date) + 1
	%ones	= @getnextname("_ones")
	series {%ones} = 1
	
	if {%eit} = NA then
		series {%ones} = NA
	endif
	
	%ti	= @getnextname("_t")
	series {%ti} = @cumsum({%ones})
	
	
	%vars 	= @getnextname("_vars")
	%vars 	= _this.@varlist
	
	%indv 	= ""
	for !i = 2 to @wcount(%vars)
    	if @word(%vars, !i) <> "C" then
		 	%indv = %indv + " " + @word(%vars, !i)
		endif
	next
	
	!k 				= @wcount(%indv)
	%e			= @getnextname("_e")
	%Xg			= @getnextname("_Xg")
	%D0 			= @getnextname("_D0")
	%D2 			= @getnextname("_D2")
	%D3 			= @getnextname("_D3")
	
	matrix {%D0} 	= @zeros(!k, !k)
	matrix {%D2} 	= @zeros(!k, !k)
	matrix {%D3} 	= @zeros(!k, !k)
	
	if !bal = 1 then
		%D1 				= @getnextname("_D1")
		matrix {%D1} 	= @zeros(!k, !k)
	endif
	
	group {%Xg} {%indv}
	
	'estimate Omega = 1/N * sum_(i=1 to N) ei*ei'
	if !bal = 1 then
		%Omega 	= @getnextname("_Omega")
		matrix {%D1} = @zeros(!k, !k)
		matrix {%Omega} = @zeros(!Tmax, !Tmax)
	
		for !i = 1 to !_N_m		' cause crossid indexing original data
			smpl if @crossid = !i and {%eit} <> NA
			matrix _e = @convert({%eit})
			matrix {%Omega} = {%Omega} + _e*@transpose(_e)
		next
		matrix {%Omega} = {%Omega}/!N_g
	endif
	
	if !bal = 1 then
	'compute V0, V1, V2, V3
		for !i = 1 to !_N_m		' cause crossid indexing original data
			smpl if @crossid = !i and {%eit} <> NA
			!_Ti = @obs({%eit})
			matrix _M = @identity(!_Ti) - @ones(!_Ti)*@transpose(@ones(!_Ti))/!_Ti
			matrix _X = @convert({%Xg})
			matrix _e = @convert({%eit})
			matrix _te = @transpose(_e)
			matrix _MX = _M * _X
			matrix _diag_eit2 = @makediagonal(@getmaindiagonal(_e * _te))
	
			matrix {%D0} = {%D0} + @transpose(_te * _MX) * _te * _MX
			matrix {%D1} = {%D1} + @transpose(_MX) * {%Omega} * _MX
			matrix {%D2} = {%D2} + @transpose(_MX) * _diag_eit2 * _MX * !_Ti/(!_Ti - 1)
			matrix {%D3} = {%D3} + @transpose(_MX) * _MX*!s2
		next
		
		matrix {%D0} = {%D0} / !N_g
		matrix {%D1} = {%D1} / !N_g
		matrix {%D2} = {%D2} / !N_g
		matrix {%D3} = {%D3} / !N_g
	
	'compute vj
		%v0 	= @getnextname("_v0")
		%v1 	= @getnextname("_v1")
		%v2 	= @getnextname("_v2")
		%v3 	= @getnextname("_v3")
	
		matrix {%v0} = @vech({%D0})
		matrix {%v1} = @vech({%D1})
		matrix {%v2} = @vech({%D2})
		matrix {%v3} = @vech({%D3})
		
	'compute Cj
		%c1 	= @getnextname("_c1")
		%c2 	= @getnextname("_c2")
		%c3 	= @getnextname("_c3")
		
		!r_c	= !k*(!k+1)/2
		matrix {%c1} = @zeros(!r_c, !r_c)
		matrix {%c2} = @zeros(!r_c, !r_c) 
		matrix {%c3} = @zeros(!r_c, !r_c)
		
		for !i = 1 to !_N_m		' cause crossid indexing original data
			smpl if @crossid = !i and {%eit} <> NA
			!_Ti = @obs({%eit})
			matrix _M = @identity(!_Ti) - @ones(!_Ti)*@transpose(@ones(!_Ti))/!_Ti
			matrix _X = @convert({%Xg})
			matrix _e = @convert({%eit})
			matrix _te = @transpose(_e)
			matrix _MX = _M * _X
			matrix _diag_eit2 = @makediagonal(@getmaindiagonal(_e * _te))
	
			matrix __c1 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*{%Omega}*_MX)
			matrix {%c1} = {%c1} + __c1*@transpose(__c1)
	
			matrix __c3 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*_MX*!s2)
			matrix {%c3} = {%c3} + __c3*@transpose(__c3)
			
			matrix __c2 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*_diag_eit2*_MX)
			matrix {%c2} = {%c2} + __c2*@transpose(__c2)
		next
		
		matrix {%c1} = ({%c1} + @transpose({%c1})) / (2 * !N_g)
		matrix {%c2} = ({%c2} + @transpose({%c2})) / (2 * !N_g)
		matrix {%c3} = ({%c3} + @transpose({%c3})) / (2 * !N_g)
	
	'compute hj statistics
		if @issingular({%c1}) = 0 then 
			!h1_pinv = 0
			!h1	= !N_g * @transpose({%v1} - {%v0}) * @inverse({%c1}) * ({%v1} - {%v0})
		else
			!h1_pinv = 1
			!h1	= !N_g * @transpose({%v1} - {%v0}) * @pinverse({%c1}) * ({%v1} - {%v0})
		endif

		if @issingular({%c2}) = 0 then
			!h2_pinv = 0
			!h2	= !N_g * @transpose({%v2} - {%v0}) * @inverse({%c2}) * ({%v2} - {%v0})
		else
			!h2_pinv = 1
			!h2	= !N_g * @transpose({%v2} - {%v0}) * @pinverse({%c2}) * ({%v2} - {%v0})
		endif

		if @issingular({%c3}) = 0 then
			!h3_pinv = 0
			!h3	= !N_g * @transpose({%v3} - {%v0}) * @inverse({%c3}) * ({%v3} - {%v0})
		else
			!h3_pinv = 1
			!h3	= !N_g * @transpose({%v3} - {%v0}) * @pinverse({%c3}) * ({%v3} - {%v0})
		endif
		
	
	'conpute p-value
		!df 	= !r_c + 1
		!p1	= 1 - @cchisq(!h1, !df)
		!p2	= 1 - @cchisq(!h2, !df)
		!p3	= 1 - @cchisq(!h3, !df)
	
	else
	'compute D0, D1, D2, D3
		for !i = 1 to !_N_m		' cause crossid indexing original data
			smpl if @crossid = !i and {%eit} <> NA
			!_Ti = @obs({%eit})
			
			if !_Ti > 1 and !_Ti >= !Tmin then
				matrix _M = @identity(!_Ti) - @ones(!_Ti)*@transpose(@ones(!_Ti))/!_Ti
				matrix _X = @convert({%Xg})
				matrix _e = @convert({%eit})
				matrix _te = @transpose(_e)
				matrix _MX = _M * _X
				matrix _diag_eit2 = @makediagonal(@getmaindiagonal(_e * _te))
		
				matrix {%D0} = {%D0} + @transpose(_te * _MX) * _te * _MX
		'		matrix {%D1} = {%D1} + @transpose(_MX) * {%Omega} * _MX
				matrix {%D2} = {%D2} + @transpose(_MX) * _diag_eit2 * _MX * !_Ti/(!_Ti - 1)
				matrix {%D3} = {%D3} + @transpose(_MX) * _MX*!s2
			endif
		next
		
		matrix {%D0} = {%D0} / !N_g
	'	matrix {%D1} = {%D1} / !N_g
		matrix {%D2} = {%D2} / !N_g
		matrix {%D3} = {%D3} / !N_g
	
	'compute vj
		%v0 	= @getnextname("_v0")
	'	%v1 	= @getnextname("_v1")
		%v2 	= @getnextname("_v2")
		%v3 	= @getnextname("_v3")
	
		matrix {%v0} = @vech({%D0})
	'	matrix {%v1} = @vech({%D1})
		matrix {%v2} = @vech({%D2})
		matrix {%v3} = @vech({%D3})
		
	'compute Cj
	'	%c1 	= @getnextname("_c1")
		%c2 	= @getnextname("_c2")
		%c3 	= @getnextname("_c3")
		
		!r_c	= !k*(!k+1)/2
	'	matrix {%c1} = @zeros(!r_c, !r_c)
		matrix {%c2} = @zeros(!r_c, !r_c) 
		matrix {%c3} = @zeros(!r_c, !r_c)
		
		for !i = 1 to !_N_m		' cause crossid indexing original data
			smpl if @crossid = !i and {%eit} <> NA
			!_Ti = @obs({%eit})
			
			if !_Ti > 1 and !_Ti >= !Tmin then
				matrix _M = @identity(!_Ti) - @ones(!_Ti)*@transpose(@ones(!_Ti))/!_Ti
				matrix _X = @convert({%Xg})
				matrix _e = @convert({%eit})
				matrix _te = @transpose(_e)
				matrix _MX = _M * _X
				matrix _diag_eit2 = @makediagonal(@getmaindiagonal(_e * _te))
		
		'		matrix __c1 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*{%Omega}*_MX)
		'		matrix {%c1} = {%c1} + __c1*@transpose(__c1)
		
				matrix __c3 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*_MX*!s2)
				matrix {%c3} = {%c3} + __c3*@transpose(__c3)
				
				matrix __c2 = @vech(@transpose(_te*_MX)*_te*_MX - @transpose(_MX)*_diag_eit2*_MX)
				matrix {%c2} = {%c2} + __c2*@transpose(__c2)
			endif
		next
		
	'	matrix {%c1} = ({%c1} + @transpose({%c1})) / (2 * !N_g)
		matrix {%c2} = ({%c2} + @transpose({%c2})) / (2 * !N_g)
		matrix {%c3} = ({%c3} + @transpose({%c3})) / (2 * !N_g)
	
	'compute hj statistics
'		if @issingular({%c1}) = 0 then 
'			!h1	= !N_g * @transpose({%v1} - {%v0}) * @inverse({%c1}) * ({%v1} - {%v0})
'		else
'			!h1	= !N_g * @transpose({%v1} - {%v0}) * @pinverse({%c1}) * ({%v1} - {%v0})
'		endif

		if @issingular({%c2}) = 0 then
			!h2_pinv = 0
			!h2	= !N_g * @transpose({%v2} - {%v0}) * @inverse({%c2}) * ({%v2} - {%v0})
		else
			!h2_pinv = 1
			!h2	= !N_g * @transpose({%v2} - {%v0}) * @pinverse({%c2}) * ({%v2} - {%v0})
		endif

		if @issingular({%c3}) = 0 then
			!h3_pinv = 0
			!h3	= !N_g * @transpose({%v3} - {%v0}) * @inverse({%c3}) * ({%v3} - {%v0})
		else
			!h3_pinv = 1
			!h3	= !N_g * @transpose({%v3} - {%v0}) * @pinverse({%c3}) * ({%v3} - {%v0})
		endif
	
	'conpute p-value
		!df 	= !r_c + 1
	'	!p1	= 1 - @cchisq(!h1, !df)
		!p2	= 1 - @cchisq(!h2, !df)
		!p3	= 1 - @cchisq(!h3, !df)
	
	endif
	smpl %smpl
	
	' Printing result
	%tab = @getnextname("_table")
	table(11,4) {%tab}
	!rowcounter = 1
	{%tab}.title Kezdi test for heteroskedasticity
	{%tab}(!rowcounter,1) = "Kezdi test for heteroskedasticity in large-N panels"
	'{%tab}.setmerge(!rowcounter) +d
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Model: Cross-section Fixed-effects"
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Alternative hypothesis (Ha): Heteroskedasticity"
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter) +d
		
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter, 1) = "Test"
	{%tab}.setwidth(a:a) 10.0
	{%tab}.setwidth(b:b) 15.0
	{%tab}.setwidth(c:c) 8.0
	{%tab}.setwidth(d:d) 15.0
	{%tab}.setjust(!rowcounter, 1) left
	{%tab}.setjust(!rowcounter, 2) center
	{%tab}(!rowcounter, 2) = "Statistic"
	{%tab}(!rowcounter, 3) = "d.f."
	{%tab}(!rowcounter, 4) = "Probability"
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter) +d
	
	'H1: Cross-sectional homoskedasticity
	if !bal = 1 then
		!rowcounter = !rowcounter + 1
		{%tab}(!rowcounter,1) = "H1 vs. Ha"
		{%tab}.setjust(!rowcounter, 1) left
		{%tab}(!rowcounter,2) = !h1
		{%tab}.setformat(!rowcounter, 2) g.5
		{%tab}(!rowcounter,3) = !df
		{%tab}.setformat(!rowcounter, 3) f.0
		{%tab}(!rowcounter,4) = !p1
		{%tab}.setformat(!rowcounter, 4) f.6
	endif
	
	'H2: Serially uncorrelated: e_it, x_it or both
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "H2 vs. Ha"
	{%tab}.setjust(!rowcounter, 1) left
	{%tab}(!rowcounter,2) = !h2
	{%tab}.setformat(!rowcounter, 2) g.5
	{%tab}(!rowcounter,3) = !df
	{%tab}.setformat(!rowcounter, 3) f.0
	{%tab}(!rowcounter,4) = !p2
	{%tab}.setformat(!rowcounter, 4) f.6
	
	'H3: Homoskedasticity and serially uncorrelated
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "H3 vs. Ha"
	{%tab}.setjust(!rowcounter, 1) left
	{%tab}(!rowcounter,2) = !h3
	{%tab}.setformat(!rowcounter, 2) g.5
	{%tab}(!rowcounter,3) = !df
	{%tab}.setformat(!rowcounter, 3) f.0
	{%tab}(!rowcounter,4) = !p3
	{%tab}.setformat(!rowcounter, 4) f.6
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter ) +d
	
	'Hypothesis notes
	if !bal = 1 then
		!rowcounter = !rowcounter + 1
		{%tab}(!rowcounter ,1)  = "H1: Cross-sectional homoskedasticity"
	endif
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter ,1)  = "H2: Serially uncorrelated: e_it, x_it or both"
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter ,1)  = "H3: Homoskedasticity and serially uncorrelated"
	
	'Singularity notes
	if !bal = 1 then
		if !h1_pinv = 1 or !h2_pinv=1 or !h3_pinv=1 then
			!rowcounter = !rowcounter + 1
			{%tab}.setjust(!rowcounter, 1) left
			{%tab}(!rowcounter ,1)  = "Notes:"
		endif
	else
		if !h2_pinv=1 or !h3_pinv=1 then
			!rowcounter = !rowcounter + 1
			{%tab}.setjust(!rowcounter, 1) left
			{%tab}(!rowcounter ,1)  = "Notes:"
		endif
	endif 	

	if !bal = 1 then
		if !h1_pinv=1 then
			!rowcounter = !rowcounter + 1
			{%tab}(!rowcounter ,1)  = "- Covariance matrix C1 is singular. Use generalized inverse."
		endif
	endif
	if !h2_pinv=1 then
		!rowcounter = !rowcounter + 1
		{%tab}(!rowcounter ,1)  = "- Covariance matrix C2 is singular. Use generalized inverse."
	endif
	if !h3_pinv=1 then
		!rowcounter = !rowcounter + 1
		{%tab}(!rowcounter ,1)  = "- Covariance matrix C3 is singular. Use generalized inverse."
	endif
	
	'References
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter ,1)  = "References: Kezdi, G. (2003)"
			
	_this.display {%tab}
	d {%tab}
	d(noerr) _*

else

'=========================================
'		Large-T: Modified-Wald test (Greene, 2000)
'=========================================

	%smpl 	= @pagesmpl
	!n 			= _this.@regobs
	!n_g 		= _this.@ncross
	
	%eit 		= @getnextname("_eit")
	_this.makeresids {%eit}
	
	%Ti		= @getnextname("_Ti")
	series {%Ti} = @obsby({%eit}, @crossid)
	
	%eit_sq	= @getnextname("_eit_sq")
	series {%eit_sq} = {%eit}*{%eit}
	
	!s			= @mean({%eit_sq})
	
	%si		= @getnextname("_si")
	series {%si} = @meansby({%eit_sq}, @crossid)
	
	%vi		= @getnextname("_vi")
	series {%vi} = @sumsby(({%eit_sq} - {%si})^2, @crossid)/({%Ti}*({%Ti}-1))
	
	%waldi	= @getnextname("_waldi")
	series {%waldi} = @meansby(({%si}-!s)^2/{%vi}, @crossid)
	
	%ones	= @getnextname("_ones")
	series {%ones} = 1
	
	if {%eit} = NA then
		series {%ones} = NA
	endif
	
	%t	= @getnextname("_t")
	series {%t} = @cumsum({%ones})
	
	series {%waldi} = @recode({%t}>1, NA, {%waldi})
	
	!wald 	= 	@sum({%waldi})
	!prob 	= 	@chisq(!wald, !N_g)
	
	d(noerr) _*
	
	' Printing result
	%tab = @getnextname("_table")
	table(9,4) {%tab}
	!rowcounter = 1
	
	{%tab}.title Modified Wald test for groupwise heteroskedasticity
	{%tab}(!rowcounter,1) = "Modified Wald test for Groupwise Heteroskedasticity in large-T panels"
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Model: Cross-section Fixed-effects"
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Null Hypothesis: Groupwise Homoskedasticity"
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter) +d
	
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Test Statistic"
	{%tab}.setwidth(a:a) 15.0
	{%tab}.setwidth(b:b) 15.0
	{%tab}.setwidth(c:c) 8.0
	{%tab}.setwidth(d:d) 15.0
	{%tab}.setjust(!rowcounter, 1) left
	{%tab}.setjust(!rowcounter, 2) center
	{%tab}(!rowcounter,2) = "Value"
	{%tab}(!rowcounter,3) = "d.f."
	{%tab}(!rowcounter,4) = "Probability"
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter) +d
	
	!rowcounter = !rowcounter + 1
	{%tab}(!rowcounter,1) = "Chi-square"
	{%tab}.setjust(!rowcounter, 1) left
	{%tab}(!rowcounter,2) = !wald
	{%tab}.setformat(!rowcounter, 2) g.5
	{%tab}(!rowcounter,3) = !N_g
	{%tab}.setformat(!rowcounter, 3) f.0
	{%tab}(!rowcounter,4) = !prob
	{%tab}.setformat(!rowcounter, 4) f.6
	
	!rowcounter = !rowcounter + 1
	{%tab}.setlines(!rowcounter) +d
	
	!rowcounter = !rowcounter + 1	
	{%tab}(!rowcounter,1)  = "References: Greene, W. H. (2000)"
		
	_this.display {%tab}
	d {%tab}
endif


