# Covid_AntibodyResponse_Analysis
Codes to carry out the analysis in paper "Comparing Neutralizing Antibody Activity Over Time Between Naive and Convalesced COVID-19 Vaccinated Individuals"
Comparison of Convalesced and Covid Naive NAb levels

Title: Comparing Neutralizing Antibody Activity Over Time Between Naive and Convalesced COVID-19 Vaccinated Individuals

Authors: Lisette de Pillis 1*, Rebecca Caffrey 2, Ge Chen 2, Mark
Dela 4*, Leif Eldevik 2, Hong Liu 2, Joseph P.
McConnell 2, Shahrokh Shabahang 2 and Stephen A. Varvel 3

1 Department of Mathematics, Harvey Mudd College, 301 Platt
Blvd., Claremont, 91711, CA, USA.
2 Aditxt Inc., 2569 Wyandotte Street, Suite 101, Mountain View,
94043, CA, USA.
3 Pearsanta Inc., 737 N. Fifth Street, Suite 200, Richmond, 23219,
VA, USA.
4* San Bernardino County Department of Behavioral Health, 303
E Vanderbilt Way, San Bernardino, 92415, CA, USA.
*Corresponding author(s). E-mail(s): depillis@hmc.edu;
mark.dela@dbh.sbcounty.gov;
Contributing authors: recaffrey@gmail.com; ; leldevik@live.com;
hliu@aditxt.com; joemcconnell50@gmail.com;
sshabahang@aditxt.com; svarvel@pearsanta.com;

Abstract: Longitudinal data comprised of neutralizing antibody (NAb) activity measurements from subjects who received COVID-19 vaccinations were collected between November 2020 and April 2022. To detect differences between convalesced and naive groups with respect to the evolution of NAb activity since the subject's first COVID-19 vaccine, we initially fit a linear mixed effects model to only the decay section of NAb evolution. We conclude that NAb activity, when restricted to this region, behaves differently between these two groups, with the convalesced group generally having higher neutralizing antibody levels than the naive group. We then fit a nonlinear mixed effects model over the entire NAb progression using a system of ordinary differential equations described by De Pillis et al. 2023 as our structural component to the model. This analysis not only supports the claim that over the entire progression, NAb activity behaves differently for convalesced and naive groups, but aligns with the linear analysis in confirming that NAb decay is slower in the convalesced group than the naive group. Finally, we use the estimated parameters from the nonlinear mixed effects model to predict NAb progression for each subject from their last observed measurement to 100 days past this measurement.

Structural model developed here: https://doi.org/10.1016/j.jtbi.2022.111280 
