mnwsvt.ado implements the score-variance test for the level of
clustering proposed in the paper

  James G. MacKinnon, Morten Ø. Nielsen, and Matthew D. Webb, "Testing
  for the appropriate level of clustering in linear regression models"

  https://ideas.repec.org/p/qed/wpaper/1428.html

The program is written by Matt Webb

mnwsvt syntax

mnwsvt yvar (xvars) (varlist) [if] [in], null(clusvar1) alt(clusvar2) 


optional arguments:

variables(#)

boots(#)

im 

yvar -- the outcome variable of interest

xvars -- the variables for which you wish to calculate the
      score-variance test. The program by default assumes that you
      are only interested in calculating the SV test for one variable.
      If you desire to calculate the SV test for more that one
      variable, use the option "variables". For instance, if you wish
      to calculate the SV test for the first two variables, then
      include "variables(2)" in the list of options

varlist -- all other control variables, factor variables are allowed

null -- the (fine) level of clustering under the null

alt -- the (coarse) level of clustering under the alternative


options:

variables(#) -- specifies how many variables after yvar you wish to
             calculate the SV test for.  By default mnwsvt assumes 1.  

boots -- specifies the number of bootstrap replications. By default, the
      program does not calculate a bootstrap P value, because doing so
      can be very expensive.

im -- calculates the Ibragimov and Mueller test statistic and P value
      using 9999 simulations.


Examples:

one variable, no bootstrap:

mnwsvt treadss1 aide_1 small_1 treadssk male nonwhite teach_nonwhite ///
totexp1 freelunch brys2 brys3 brys5 sbq2 sbq3 sbq4 hdg2 hdg3 hdg4 ///
i.schid1n, null(newid) alt(clsid)

one variable, no bootstrap, with IM test

mnwsvt treadss1 aide_1 small_1 treadssk male nonwhite teach_nonwhite ///
totexp1 freelunch brys2 brys3 brys5 sbq2 sbq3 sbq4 hdg2 hdg3 hdg4 ///
i.schid1n, null(newid) alt(clsid)

two variables, no bootstrap

mnwsvt treadss1 small_1 aide_1 treadssk male nonwhite teach_nonwhite ///
totexp1 freelunch brys2 brys3 brys5 sbq2 sbq3 sbq4 hdg2 hdg3 hdg4 ///
i.schid1n, null(newid) alt(clsid) variables(2)

two variables, bootstrap with 999 bootstrap samples

mnwsvt treadss1 small_1 aide_1 treadssk male nonwhite teach_nonwhite ///
totexp1 freelunch brys2 brys3 brys5 sbq2 sbq3 sbq4 hdg2 hdg3 hdg4 ///
i.schid1n, null(newid) alt(clsid) variables(2) boots(999)
