# odontocete-lifehistory
R code and data for phylogenetic analysis of odontocete life history traits.
This repository contains R code and data for analyzing life-history traits of odontocetes using phylogenetic methods, supporting a forthcoming journal publication.

Raw data odontocete.csv

subo = suborder
f = family code
fam = family
sp = scientific species name
name = common name
group = Archetype
clus = Results from Ferguson, S.H., Higdon, J.W., Schmidt, C., Pomerleau, C. and Matthews, C.J., 2023. Investigating the relationship between body shape and life history traits in toothed whales: Can body shape predict fast-slow life histories?. Evolutionary Biology, 50(3), pp.300-317.
len = adult body length (m)
loglen = log(len)
weight = adult body weight (g)
logwt = log(wt)
wtpred = predicted weight controlled for body length
wtresid = residual weight from log(wt) ~ log(len)
gest = gestation length (mo)
loggest = log(gest)
gestpred = predicted gestation length from regression
gestresid = residual gestation from log(gest) ~ log(len)
same fro inter = interbirth interval
long = longevity (y)
newolen = neonate length (m) at birth
asm = age at sexual maturation (y)
lbathy = log(bathymetry)
Ice15sMed = Median number of months in a year with sea ice
ablat = absolute value of latitude
