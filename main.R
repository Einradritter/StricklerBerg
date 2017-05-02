# fluorescein data from:
# http://omlc.org/spectra/PhotochemCAD/html/010.html
# http://omlc.org/spectra/PhotochemCAD/html/037.html


# load packages -----------------------------------------------------------
library(tibble)
library(readr)
library(NISTunits)
library(RCurl)
library(bitops)
library(minpack.lm)

# functions ---------------------------------------------------------------
last <- function(x) { return( x[length(x)] ) } #get last element of vector

turn_integration_limits_into_index_interval<-function(wl_data, wl_lower_limit, wl_upper_limit) #converts wl into index vector for integration
{
  wl_min <- min(wl_data)
  wl_data_density<-length(wl_data)/(max(wl_data)-min(wl_data))
  return(seq(from = floor(wl_data_density*(wl_lower_limit-wl_min)), to = floor(wl_data_density*(wl_upper_limit-wl_min)), by = 1))
}

integrate_num<- function(x,y){sum(y[-length(y)]*abs(diff(x)))} #crude numerical data integration function

visualize_integration<-function(x, y, lower_limit, upper_limit, integral, title='') #visualize of integral with box
{
  delta<- abs(upper_limit - lower_limit)
  plot(x = x, y = y, type = 'l', main = title)
  lines(x = x, y=rep(0,length(x)), col="red") #add baseline
  lines(x = c(lower_limit,upper_limit), y=rep(integral/delta,2), col="blue") #add box
  lines(x = rep(lower_limit,2), y= c(0,integral/delta), col="blue")
  lines(x = rep(upper_limit,2), y= c(0,integral/delta), col="blue")
  lines(x = c(lower_limit,upper_limit), y=rep(0,2), col="blue")
}

convert_wavelength_to_wavenumber <- function(wavelength) {1/wavelength*10^(7)} # conversion functions
convert_wavenumber_to_wavelength <- function(wavenumber) {1/wavenumber*10^(7)}

normalize_emission <- function(emission_data) { emission_data/max(emission_data)}

# physical constants ------------------------------------------------------
speed_of_light<- kNIST2010speedOfLightInVacuum #[m/s]
Avogadro_Const<- kNIST2010AvogadroConstant #[1/mol]
# SB parameters -----------------------------------------------------------
degeneration_of_lower_state<- 1 #[-]
degeneration_of_upper_state<- 1 #[-]
refractive_index<- 1.3611 #[-]
quantum_yield<- 0.79 #[-]
results<-tibble(' ' = 'results:')
unit_factor<- 10^(-5)

# get experimental data from online source --------------------------------
file_abs<- getURL(url = 'http://omlc.org/spectra/PhotochemCAD/data/010-abs.txt')
abs <- read_delim(file = file_abs, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 22)
colnames(abs) <- c('wavelength', 'molar_extinction') #wavelength in [nm] | molar_extinction in [L/mol/cm]
abs$wavenumber<- convert_wavelength_to_wavenumber(abs$wavelength) #wavenumber in [1/cm]
abs$wavenumber_log<- log(abs$wavenumber) # logarithmic wavenumber in [1/cm]??? 

file_ems<- getURL(url = 'http://omlc.org/spectra/PhotochemCAD/data/010-ems.txt')
ems <- read_delim(file = file_ems, "\t", escape_double = FALSE, trim_ws = TRUE, skip = 15)
colnames(ems)<- c('wavelength', 'emission') #wavelength in [nm] | emission in [au]
ems$wavenumber<- convert_wavelength_to_wavenumber(ems$wavelength) # wavenumber in [1/cm]
ems$emission_normalized<- normalize_emission(ems$emission)

# visualize absorption ----------------------------------------------------
plot(x = abs$wavelength, y = abs$molar_extinction, type = 'l', main = 'absorption data')
lines(x = abs$wavelength, y=rep(0,length(abs$wavelength)), col="red") #add baseline

# set integration interval ------------------------------------------------
wl_lower_limit= 400 #[nm]
wl_upper_limit= 550 #[nm]
integration_interval<-turn_integration_limits_into_index_interval(wl_data = abs$wavelength, wl_lower_limit = wl_lower_limit, wl_upper_limit = wl_upper_limit)

# numerical integration of absorption data --------------------------------
x<- abs$wavenumber_log[integration_interval] # def x and y data
y<- abs$molar_extinction[integration_interval]
log_integral<-integrate_num(x = x, y = y) #integration
visualize_integration(x = abs$wavenumber_log, y = abs$molar_extinction, lower_limit = log(convert_wavelength_to_wavenumber(wl_lower_limit)),
                         upper_limit = log(convert_wavelength_to_wavenumber(wl_upper_limit)),integral = log_integral, title = 'logarithmic Integration of absorption data')
results$SB_abs_factor<- log_integral

# visualize emission ------------------------------------------------------
plot(ems$wavelength, ems$emission_normalized, type = 'l', main = 'emission data')
lines(x = ems$wavelength, y=rep(0,length(ems$wavelength)), col="red") #add baseline

# set integration interval ------------------------------------------------
wl_lower_limit= 450
wl_upper_limit= 745
integration_interval<-turn_integration_limits_into_index_interval(wl_data = ems$wavelength, wl_lower_limit = wl_lower_limit, wl_upper_limit = wl_upper_limit)

# numerical integration of emission data ----------------------------------
# geht besser mit der Simpson-Formel?
x<- ems$wavenumber[integration_interval]
y<- ems$emission_normalized[integration_interval]
ems_integral_nominator<- integrate_num(x = x, y = y) #nominator integral[a.u./cm]

visualize_integration(x = ems$wavenumber, y = ems$emission_normalized, lower_limit = convert_wavelength_to_wavenumber(wl_lower_limit),
                      upper_limit = convert_wavelength_to_wavenumber(wl_upper_limit),integral = ems_integral_nominator, title = 'integration of normalized emission data')

# gaussian fit ------------------------------------------------------------
#gaussian_density<- y ~ amplitude* 1/sqrt(2*pi*standard_deviation^(2))*exp(-1*(x-expected_value)^(2)/(2*standard_deviation^(2)))
#
#nlsfit<-nlsLM(formula = gaussian_density,
#              start = list(amplitude = 1, standard_deviation = 2000, expected_value = 19000)) 
#summary(nlsfit)
#plot(x = x, y = y, main = "fit", type = 'l')
#lines(x, predict(nlsfit), col= 'red')

gaussian_sum <- y ~ amplitude* 1/sqrt(2*pi*standard_deviation^(2))*exp(-1*(x-expected_value)^(2)/(2*standard_deviation^(2))) +
                    amplitude_2* 1/sqrt(2*pi*standard_deviation_2^(2))*exp(-1*(x-expected_value_2)^(2)/(2*standard_deviation_2^(2)))

nlsfit<-nlsLM(formula = gaussian_sum,
            start = list(amplitude = 1, standard_deviation = 200, expected_value = 19000, 
                         amplitude_2 = 1, standard_deviation_2 = 200, expected_value_2 = 17000),
                        control = nls.lm.control(maxiter=100))
summary(nlsfit)
plot(x = x, y = y, main = "fit", type = 'l')
lines(x, predict(nlsfit), col= 'red')
# equation (20!) ----------------------------------------------------------
y<-ems$emission_normalized[integration_interval]*x^(-3)
ems_integral_denominator<-integrate_num(x = x, y = y)
visualize_integration(x = ems$wavenumber, y = ems$emission_normalized/x^(3), lower_limit = convert_wavelength_to_wavenumber(wl_lower_limit),
                      upper_limit = convert_wavelength_to_wavenumber(wl_upper_limit),integral = ems_integral_denominator, title = 'integration of normalized emission*wavenumber-3 [equation 20]')

results$SB_em_factor<-ems_integral_nominator/ems_integral_denominator


# Sb equation (22) --------------------------------------------------------
Einstein_A<-(8*2303*pi*speed_of_light*refractive_index^(2)/Avogadro_Const*results$SB_em_factor*results$SB_abs_factor)*unit_factor
results$SB_lifetime_in_s<-1/Einstein_A #in [s]

#Versuche nun die Fluoreszenz auf 1 zu normieren.
#Emissionfactor scheint richtig zu sein