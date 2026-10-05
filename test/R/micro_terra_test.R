library(NicheMapR) # load the NicheMapR package
micro <- micro_terra(ystart = 2000,
                     yfinish = 2000,
                     clearsky = 0,
                     cap = 0,
                     snowmodel = 0,
                     runmoist = 0,
                     scenario = 0,
                     run.gads = 2) # the aerosol step in R; the Fortran one crashed R here
metout <- as.data.frame(micro$metout)
soil <- as.data.frame(micro$soil)
write.csv(metout, file = 'c:/git/MicroclimateMapper.jl/test/data/micro_terra/metout_monthly_terra.csv')
write.csv(soil, file = 'c:/git/MicroclimateMapper.jl/test/data/micro_terra/soil_monthly_terra.csv')
