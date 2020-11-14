## deep debugging of the get_expo_risk_i function

manifest_times = manifest
report = miss_dat$report
times = miss_dat$report.times

manifested = manifest$manifested

for(i in manifested){
  cat('Person',i,': ')
  t_i = manifest_times$times[manifested == i]
  tE_i = get_expo_risk_i(i, t_i, G_all, tmax, tmin, report, times, 
                  events, recovery_times, exp_eta, details = FALSE)
  cat('Imputed exposure time:',tE_i,'\n')
}
