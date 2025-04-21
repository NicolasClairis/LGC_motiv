# compare_correlations will compare the correlations between dmPFC-Lactate and kEp vs dmPFC-lactate and kEm and vs aIns-lactate and kEp
# In other words, always one overlapping variable (dmPFC-Lactate or kEp) and dependent groups

# install cocor package from Diedenhofen & Musch
install.packages("cocor")

 #load library
 library("cocor")

# load relevant data

r.dmpfcLac_kEp = 0.1546
r.dmpfcLac_kEm = 0.0821
r.kEp_kEm = 0.7458
nSubs = 63

# compare the correlations
cocor(~dmpfc_lac + kEp | dmpfc_lac + kEm)

# compare the correlations
cocor.dep.groups.overlap(
  r.dmpfcLac_kEp,
  r.dmpfcLac_kEm,
  r.kEp_kEm,
  nSubs,
  alternative = 'greater',
  test = "all",
  alpha = 0.05,
  conf.level = 0.95,
  null.value = 0
)