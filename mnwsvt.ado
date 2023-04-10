/*-----------------------------*/
/*ADO file*/
/* mnwsvt  - score variance cluster test */
/* implements procedures found in this paper */
/* https://ideas.repec.org/p/qed/wpaper/1428.html */ 
/* written by Matt Webb */
/* version 1 */
/* date 03/09/23 */
/*-----------------------------*/

capture program drop mnwsvt
program define mnwsvt, rclass
  syntax varlist(fv) [if] [in], null(varname) alt(varname) [VARiables(real 1) boots(real 1) im]
  
  local numvars = wordcount("`varlist'")
  local y = word("`varlist'",1)
      
  local xvars ""

  local start = 1 + `numvars' + 1
      
  /*split indepvars into variables of interest and controls*/
  if `variables' == 1{
    local x = word("`varlist'",2)
    
    forvalues i = 3/`numvars' {
    
      local wordi = word("`varlist'",`i')
    
      local xvars = "`xvars' `wordi' "  
    }
  }
  
  if `variables' > 1 {
    local x = ""
    local numend = 1 + `variables'
    local numstart = `numend' +1
    forvalues i = 2/`numend' {
      
      local wordi = word("`varlist'",`i')
    
      local x = "`x' `wordi'"
    }
    forvalues i = `numstart'/`numvars' {
    
      local wordi = word("`varlist'",`i')
    
      local xvars = "`xvars' `wordi' "  
    } 
  }
  
  local null = "`null'"
  local alt = "`alt'"
    
  marksample temp_sample  
  markout `temp_sample' `y' `x' `xvars' `null' `alt' 
    
  /*unrestricted residuals, for scores and bootstraps*/
  tempvar temp_uhat
  qui reg `y' `x' `xvars' 
  qui predict `temp_uhat' if `temp_sample'==1, res
  
  /*values for m factors*/
    local NmK = e(df_r)
    local N = e(N)  
    mata: N = `N'
    local K = `N' - `NmK'
      
  /*fitted values for bootstrap*/
  tempvar temp_xbr
  qui predict `temp_xbr' if `temp_sample'==1, xb
  
  /*determine number of clusters under alternative*/
  tempvar temp_indexa
  qui egen `temp_indexa' = group(`alt') if `temp_sample' ==1
  qui summ `temp_indexa' 
  local G = r(max)
  mata: G = `G'
  
  /*correct for possibility of non-unique naming of fine clusters*/
  /*determine number of clusters under null*/
  tempvar temp_indexn
  qui egen `temp_indexn' = group(`null') if `temp_sample' ==1  
  tempvar temp_ints
  qui egen `temp_ints' = group(`temp_indexa' `temp_indexn')   
  qui summ `temp_ints'
  local H = r(max)
  mata: H = `H'
  
  /*put alternative cluster index into mata*/
  qui sort `temp_indexa'
  qui putmata cg = `temp_indexa' if `temp_sample'==1, omitmissing replace
  
  /*put null cluster index into mata*/
  tempvar temp_count
  qui bysort `temp_indexa' `temp_ints' : gen `temp_count' = _n
  qui putmata mg = `temp_indexa' if `temp_count'==1 & `temp_sample'==1 , omitmissing replace
  
  /*put uhat and yhat into mata*/
  qui putmata uhat = `temp_uhat' , omitmissing replace
  qui putmata xb = `temp_xbr' , omitmissing replace
    
  /*creates table of subclusters index by G*/
  mata: infoM = panelsetup(mg,1)
    
  qui putmata ch = `temp_ints' if `temp_sample'==1, omitmissing replace

  /*creates tables of obs by null cluster, or by alt cluster*/
  mata: infoH = panelsetup(ch,1)
  mata: infoG = panelsetup(cg, 1)
      
  /*calculate scores for each  XVAR*/
  local j = 0
  
  foreach xvar in `x' {
    local j = `j' + 1
    
    qui reg `xvar' `xvars'
    local varname = "temp_z`j'"
    tempvar `varname'
    qui predict ``varname'' if `temp_sample'==1, res 
    qui putmata  `varname' = ``varname''  if `temp_sample'==1, omitmissing replace
    local xnum = `j'

  } /*end of xvar*/
    
  /*m-factor under alternative*/
    local M_A = (`G'/(`G' - 1)) * ((`N'-1)/`NmK')
    mata: M_A = `M_A'

  /*m-factor under null*/
    local M_F = (`N'-1)/(`NmK')*`H'/(`H'-1) 
    mata: M_F = `M_F'
  
  /*Gk and Hk matrik*/
    mata: Gk = J(`xnum'^2,`xnum'*(`xnum'+1)/2,0)
    
  forvalues i = 1/`xnum' {  
    forvalues j = `i' /`xnum' {
      local a = (`j'-1)*`xnum' + `i'      
      local b = (`i' - 1)*`xnum' + `j'      
      local c = (`j'-1)*(`xnum'-`j'/2)+`i'    
      mata: Gk[ `a' ,`c' ] =1
      mata: Gk[ `b',`c'] =1
    } /*end of j*/
  } /*end of i */
    
  mata: Hk = invsym(Gk'*Gk)*Gk'

  /*put G into mata*/ 
    mata: maxg = `G'
    
  /*put xnum  into mata*/ 
    mata: xnum = `xnum'
      
  /*calculate tau for original sample*/
  local j = 0
  foreach xvar in `x' {
    local j = `j' + 1
    local vartwo = "temp_z`j'"
    mata: s`j' = uhat :* `vartwo'
  }
  
  /*put all scores into one matrix*/
    mata: sh = J(N,xnum,.)
    mata: st_numscalar("xnum", xnum)
    mata: st_local("xnum", strofreal(xnum))
        
    forvalues k = 1/`xnum'{
        mata: sh[.,`k'] = s`k'
    }
      
    mata{
        
      /*initialize matrices*/
         temp_sumg = J(xnum,xnum,0) 
         temp_sumh = J(xnum,xnum,0) 
         temp_num_alt = J(xnum,xnum,0)
         temp_cross = J(xnum ,xnum , 0)
         var_right = J(rows(Hk),rows(Hk),0)
         var_left = J(rows(Hk),rows(Hk),0)
         var = J(rows(Hk),rows(Hk),0)
         ALT = J(xnum ,xnum , 0)
         
         NLL = J(xnum , xnum , 0)
         theta = J(rows(Hk),1,0)
            
        /*sum scores by either alt or null clusters*/
        sh_h = panelsum(sh,infoH )
        sg = panelsum(sh,infoG)
                  
      for (g=1; g<=maxg; g++){
                    
         temp_sumh = J(xnum,xnum,0) 
         temp_var_left = J(xnum,xnum,0)
         temp_var_right = J(xnum,xnum,0) 
         
         /*alt numerator*/
          temp_sg = sg[g,.]'
          temp_num_alt = temp_sg * temp_sg'
          ALT = ALT + temp_num_alt
          
        /*which obs are in cluster g*/
        strt = infoM[g,1]
        ed = infoM[g,2]

        for (i=strt; i<=ed; i++) {
        
          /*extract relevant row, i,  from score matrix*/
          sh1 = sh_h[i,.]'
          temp_cross = sh1 * sh1'
          
          /*var left*/
            temp_sumh = temp_sumh + temp_cross
          
          /*var right*/
            temp_var_right = Hk * (temp_cross # temp_cross) * Hk'
            var_right = var_right + temp_var_right

        } /*i loop */
        
        /*var left*/
          temp_var_left = Hk * (temp_sumh # temp_sumh) * Hk'
          var_left = var_left + temp_var_left
        
        /*null numerator*/
          temp_sumg = temp_sumg + temp_sumh
          
      } /*g loop*/
      
    /*theta*/

      NLL = temp_sumg
          
      NLL = M_F:*NLL
      ALT = M_A:*ALT    

      theta = vech(ALT-NLL)
    
    /*variance*/
      var_left = 2*var_left
      var_right = 2*var_right
      var = var_left - var_right
    
    /*tau stat*/
    if (xnum == 1)  tau = theta/sqrt(var) 
      
    else if (xnum !=1) tau = theta'* invsym(var)*theta 
    
    } /*mata end*/
    
    *disp "what about now?"
  
  /*store result as local and as tauhat for bootstrap P*/
    mata: tauhat = tau
    mata: st_numscalar("tau", tau)
    mata: st_local("tau", strofreal(tau))
    
  /*calculate chi2 dof*/
    local dfchi = `variables'*(`variables'+1)/2
  
  /*calculate and display analytic P value*/
    if (`variables' == 1) {
      local MNW_P_2s = 2*min(normprob(`tau'),1-normprob(`tau'))
      local MNW_P = 1-normprob(`tau')
      
      return scalar P2s = `MNW_P_2s'

    }
    else if (`variables' >= 2) local MNW_P = chi2tail(`dfchi', `tau')     
    
    *disp "what is boots `boots'"
  /*if bootstrap p-value is requested*/
    if inrange(`boots', -0.01, 0.01) == 1 disp " " 
    else if inrange(`boots', -0.01, 0.01) != 1 {
    *if `boots'>=2 {
      
      *disp "what is boots `boots'"
      
      mata: tau_star = J(`boots',1,.)
      
      forvalues b = 1/`boots' {
          
        /*generate bootstrap samples*/
        mata: raddraw = J(N,1,.)
        
        forvalues h = 1/`H' {
          local rndnum = runiform()
          if (`rndnum'<0.5) local radweight = -1 
          if (`rndnum'>=0.5) local radweight = 1 
        
          mata: rowstart = infoH[`h',1]
          mata: rowend = infoH[`h',2]
          mata: diff = rowend - rowstart + 1
          mata: grouprad = J(diff,1,`radweight')
          mata: raddraw[rowstart::rowend,1] = grouprad
          
        }
        
        mata: ystar = xb +  uhat:*raddraw

        cap drop `ystar'
        tempvar ystar 
        
        getmata `ystar' = ystar, force
            
        qui reg `ystar' `x' `xvars' 
        
        
        tempvar temp_uhat_new
        qui predict `temp_uhat_new' if `temp_sample'==1, res
        qui putmata uhat = `temp_uhat_new' if `temp_sample'==1, omitmissing replace
        
        cap drop `temp_uhat_new'
        
        /*calculate tau for bootstrap sample*/

        local j = 0
        foreach xvar in `x' {
          local j = `j' + 1
          local vartwo = "temp_z`j'"
          mata: s`j' = uhat :* `vartwo'
        }
        
        /*put all scores into one matrix*/
        mata: sh = J(N,xnum,.)

        forvalues k = 1/`xnum'{
            mata: sh[.,`k'] = s`k'
        }
          
        mata{
            
          /*initialize matrices*/
             temp_sumg = J(xnum,xnum,0) 
             temp_num_alt = J(xnum,xnum,0)
             var_right = J(rows(Hk),rows(Hk),0)
             var_left = J(rows(Hk),rows(Hk),0)
             ALT = J(xnum ,xnum , 0)
             NLL = J(xnum , xnum , 0)
             theta = J(rows(Hk),1,0)
                
            /*sum scores by either alt or null clusters*/
            sh_h = panelsum(sh,infoH)
            sg = panelsum(sh,infoG)
                      
          for (g=1; g<=maxg; g++){
                        
             temp_sumh = J(xnum,xnum,0) 
             temp_var_left = J(xnum,xnum,0)
             temp_var_right = J(xnum,xnum,0) 
             
             /*alt numerator*/
              temp_sg = sg[g,.]'
              temp_num_alt = temp_sg * temp_sg'
              ALT = ALT + temp_num_alt
              
            /*which obs are in cluster g*/
            strt = infoM[g,1]
            ed = infoM[g,2]

            for (i=strt; i<=ed; i++) {
            
              /*extract relevant row, i,  from score matrix*/
                sh1 = sh_h[i,.]'
                temp_cross = sh1 * sh1'
              
              /*var left*/
                temp_sumh = temp_sumh + temp_cross
              
              /*var right*/
                temp_var_right = Hk * (temp_cross # temp_cross) * Hk'
                var_right = var_right + temp_var_right

            } /*i loop */
            
            /*var left*/
              temp_var_left = Hk * (temp_sumh # temp_sumh) * Hk'
              var_left = var_left + temp_var_left
            
            /*null numerator*/
              temp_sumg = temp_sumg + temp_sumh
              
          } /*g loop*/
          
        /*theta*/
          
          NLL= temp_sumg
          
          NLL = M_F:*NLL
          ALT = M_A:*ALT
          
          theta = vech(ALT-NLL)
        
        /*variance*/
          var_left = 2*var_left
          var_right = 2*var_right
          var = var_left - var_right
        
        /*tau stat*/
        if (xnum == 1)  taustar = theta/sqrt(var) 
          
        else if (xnum !=1) taustar = theta'* invsym(var)*theta 
        
        } /*end of mata*/
            
        mata: tau_star[`b',1]=taustar
        
      } /*end of bootstrap loop*/
          
    /*calculate and display bootstrap P value*/

      if (`xnum'==1) {
        mata: temp_rej_2s =  abs(tauhat[1,1]):<=abs(tau_star)
        mata: temp_rej = tauhat[1,1]:<=tau_star 
        
        /*two sided*/
        mata: temp_U = J(rows(temp_rej_2s),1,1) 
        mata: temp_sum = temp_U'*temp_rej_2s
        mata: boot_p2s = temp_sum / rows(temp_rej_2s)
        
        mata: st_numscalar("boot_p2s", boot_p2s)
        mata: st_local("boot_p2s", strofreal(boot_p2s))
        
        return scalar BP2s = `boot_p2s'
        
      }
      else if(`xnum'>1) {
        mata:  temp_rej =  tauhat[1,1]:<=tau_star 
      }
        
      mata: temp_U = J(rows(temp_rej),1,1) 
      mata: temp_sum = temp_U'*temp_rej
      mata: boot_p = temp_sum / rows(temp_rej)
      
      mata: st_numscalar("boot_p", boot_p)
      mata: st_local("boot_p", strofreal(boot_p))
      
      return scalar BP = `boot_p'
      
        
  } /*end of bootstrap "if" code*/
    
  return scalar tau = `tau'   
  return scalar xnum = `variables'
  return scalar P = `MNW_P'
    
  
  
  /*display output*/
    disp ""
    disp ""
    disp "MNWSV  - MacKinnon, Nielsen, and Webb"
    disp " "
    disp "Testing a null of `null' against an alternative of `alt'."
    disp "There are `H' fine clusters and `G' coarse clusters."
    disp "The specified coefficients are `x'."

  /*output table - depends on k>1, boots, and IM*/
  
  /*no bootstrap and two-sided test*/
  if `variables'==1 & `boots'== 1 {
    
    matrix output = J(1,3,.)
    matrix output[1,1] = `tau'
    matrix output[1,2] = `MNW_P'
    matrix output[1,3] = `MNW_P_2s'
    
    matrix colnames output = "tau" "P one-sided" " P two-sided" 
    matrix rownames output = "asymptotic" 
    matlist output, title(Score Variance Test) ///
    rowtitle(test) /// 
    cspec(& %-12s w6 | %9.6f w10 & %6.4f w12 & %6.4f w12 |) ///
    rspec(---)  
  }
  
  /*yes bootstrap and two-sided test*/
  if `variables'==1 & `boots'>= 2 {
    
    matrix output = J(2,3,.)
    matrix output[1,1] = `tau'
    matrix output[1,2] = `MNW_P'
    matrix output[1,3] = `MNW_P_2s'
    matrix output[2,1] = `tau'
    matrix output[2,2] = `boot_p'
    matrix output[2,3] = `boot_p2s'
    
    matrix colnames output = "tau" "P one-sided" " P two-sided" 
    matrix rownames output = "asymptotic" "bootstrap"
    matlist output, title(Score Variance Test) ///
    rowtitle(test) /// 
    cspec(& %-12s w6 | %9.6f w10 & %6.4f w12 & %6.4f w12 |) ///
    rspec(----) 
  }
  
  /*no bootstrap and one-sided test*/
  if `variables'>=2 & `boots'== 1 {
    
    matrix output = J(1,2,.)
    matrix output[1,1] = `tau'
    matrix output[1,2] = `MNW_P'
    
    matrix colnames output = "tau" "P one-sided" 
    matrix rownames output = "asymptotic" 
    matlist output, title(Score Variance Test) ///
    rowtitle(test) /// 
    cspec(& %-12s w6 | %9.6f w10 & %6.4f w12 |) ///
    rspec(---)  
  }
  
  /*yes bootstrap and one-sided test*/
  if `variables'>=2 & `boots'>= 2 {
    
    matrix output = J(2,2,.)
    matrix output[1,1] = `tau'
    matrix output[1,2] = `MNW_P'
    matrix output[2,1] = `tau'
    matrix output[2,2] = `boot_p'
    
    matrix colnames output = "tau" "P one-sided"  
    matrix rownames output = "asymptotic" "bootstrap"
    matlist output, title(Score Variance Test) ///
    rowtitle(test) /// 
    cspec(& %-12s w6 | %9.6f w10 & %6.4f w12|) ///
    rspec(----) 
  }
  
  
  /*----------------------------*/
  /*run IM program if requested*/
  /*----------------------------*/
  
  if "`im'" != "" {
    if `variables' >1 {
      disp as error "The IM Test only works for one variable of interest"
    }
    else if `variables'==1 {
      
      *disp "calculate IM test"
      IMTEST  `y' "`x'" "`xvars'" `null' `alt' 
      
      mata: st_matrix("imoutput", imoutput)

      matrix colnames imoutput = "variance" "P-value"
      matrix rownames imoutput = "IM test"  
      matlist imoutput, title("Ibragimov and Mueller Test") ///
      rowtitle(test) /// 
      cspec(& %-12s w6 | %10.2f w10 & %6.4f w12 |) ///
      rspec(---)  
      
    }
  }
    
end


/*IM Test program, optionally called by MNWTEST*/
capture program drop IMTEST
program define IMTEST, rclass 

  local y "`1'"
   
  local x "`2'"
   
  local CTRLVAR "`3'" 

  local null "`4'"
   
  local alt "`5'"
     
/*determine number of clusters under alternative*/ 
  marksample temp_sample  
  markout `temp_sample' `y' `x' `CTRLVAR' `alt'
  
  tempvar temp_grp
  qui egen `temp_grp' = group(`alt') if `temp_sample'==1
  qui summ `temp_grp'
  
/*store number of clusters as j in mata*/
  local j = r(max)
  mata: j = `j'
  
/*create empty matrices to store the j beta and s.e.*/
  mata: beta = J(j,1,.)
  mata: omega = J(j,1,.)
  
/*calculate the beta and s.e. per coarse cluster*/
  /*cluster s.e. under the null*/
  forvalues g=1/`j' {
    
    qui reg `y' `x' `CTRLVAR' if `temp_grp'==`g', cluster(`null')
    
    /*store the beta and s.e. estimates in the resepective matrix*/
    local beta = _b[`x']
    mata: beta[`g',1]=`beta'
    local se = _se[`x']
    mata: omega[`g',1]=`se'
  }
  
/*Calculate the IM (2012) standard error*/  
  /*it is just the variance of the j beta */
  mata: S2 = variance(beta)
  
/*matrix for all the yj estimates*/ 
mata: ybar = J(9999,1,.)

/*replication loop*/
forvalues k=1/9999{
  
  /*multiply the standard errors by standard normals*/
  mata: yj = omega:*invnormal(uniform(j,1))
  
  /*calculate the average of the yj*/
  mata avey = mean(yj[.,1])
  
  /*square the mean differences*/
  mata: sy2 = (yj:-avey):*(yj:-avey)
  
  /*sum the squares, divide by (j-1) */
  mata: ybk = (1/(j-1))* colsum(sy2[.,1])

  /*store in the matrix*/
  mata: ybar[`k',1]=ybk
}

/*calculate the p-value*/
  mata: temp_rej =  S2[1,1]:<ybar
  mata: temp_U = J(rows(temp_rej),1,1) 
  mata: temp_sum = temp_U'*temp_rej
  mata: IM_p = temp_sum / rows(temp_rej)
  
  mata: imoutput = S2, IM_p
  
end


/*------------------------*/
/*change log*/
/*------------------------*/
*1.0 - first version of mnwsvt
