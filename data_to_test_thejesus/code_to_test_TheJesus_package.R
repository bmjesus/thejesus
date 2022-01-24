library('devtools')

devtools::document()
#check documentation
devtools::check_man()
devtools::build()
devtools::install()



1#REFLECTANCE


a<-read_ooptics("./data_to_test_thejesus/blue_tokio.txt",plot = TRUE,indices = TRUE)

a<-read_ooptics("./data_to_test_thejesus/1_1.txt",indices = TRUE)

str(a)

#ASD machine
a<-read_asd("./data_to_test_thejesus/Spectrum00018.asd.ref.txt",plot = TRUE,indices = TRUE)

a<-read_asd("./data_to_test_thejesus/Spectrum00027.asd.ref.txt",plot = TRUE,indices = TRUE)


read.table("./data_to_test_thejesus/Spectrum00027.asd.ref.txt",skip=35,sep="\t")
read.table("./data_to_test_thejesus/Spectrum00018.asd.ref.txt",skip=35,sep="\t")


#DERIVATIVES

d2(a$all_data$wl,a$all_data$reflect
   ,reflectance = TRUE,smooth_reflectance=0.04)


#INTEGRAL



#O2 files
light<-c(0,34,56,150,680)

time_intervals<-as.data.frame(rbind(c(0,4*60),c(10*60,15*60),c(20*60,25*60),c(30*60,35*60),c(40*60,45*60)))

steps<-cbind(light,time_intervals)

thejesus::light_curvesO2("./data_to_test_thejesus/O2_data.txt",steps)

light_curvesO2("./data_to_test_thejesus/O2_data.txt",steps)

thejesus::light_curvesO2("./data_to_test_thejesus/TEST3_27-02-20_Optode.txt",steps)


light_curvesO2("./data_to_test_thejesus/TEST3_27-02-20_Optode.txt",steps)



#FLUORESCENCE

d1<-read.table("./data_to_test_thejesus/etr_test.txt",header = TRUE)

fit_etr(light = d1$light,etr = d1$etr,
        model = 'JP')

fit_etr(light = d1$light,etr = d1$etr,
        model = 'P')

fit_etr(light = d1$light,etr = d1$etr,
        model = 'EP')

d2<-fit_etr(light = d1$light,etr = d1$etr,
        model = 'EP',starting_values = list(alpha = -1.4, etrmax = 4099, Eopt=2500))

fit_etr(light = d1$light,etr = d1$etr,
        model = 'JP',num_obs = 7)

fit_etr(light = d1$light,etr = d1$etr,
        model = 'JP')

fit_etr(light = d1$light,etr = d1$etr,
        model = 'JP',
        starting_values = list(alpha = 2, etrmax = 400))

d2<-read.table("npq_test.txt",header = TRUE)

fit_npq(light = d2$light,npq = d2$npq,
        starting_values = list(npqmax = 1, E50 = 40, hill=2))

fit_npq(light = d2$light,npq = d2$npq,
        starting_values = list(npqmax = 10, E50 = 90, hill=1))


d3<-read.table("eff_test.txt",header = TRUE)

fit_eff(light = d3$light,eff = d3$yield,
        model = 'EP',starting_values = list(a=1,b=1,c=3))

fit_eff(light = d3$light,eff = d3$yield,
        model = 'EP')

fit_etr(light = d3$light,etr = d3$yield*d3$light,
        model = 'EP')


fit_eff(light = (d3$light+0.0001),eff = d3$yield,
        model = 'JP',starting_values = list(alpha=0.4,Ek=200))


fit_eff(light = (d3$light+0.0001),eff = d3$yield,
        model = 'P')


0.34 * 90 * tanh(d3$light * 90^{-1}) * d3$light^{-1}


algorithm="port"


d4<-read.table("pam_test.csv",header = TRUE,sep=',')


#functions to-do

#regulated thermal energy dissipation related to NPQ
#YNPQ = F/Fm' - F/Fm

a<-ynpq(f=d4$f1,fm=d4$fm1[1],fmp=d4$fm1)
fit_npq(d4$par,a,num_obs = 8)

#"primarily constitutive losses", corresponding to the sum of non-regulated heat dissipation and fluorescence emission. Y(NO) describes the combined pathways of radiative and non-radiative deexcitation reactions,which do not lead to photochemical energy conversion and are not involving the NPQ-mechanism.
#YNO = F/Fm

yno(f=d4$f1,fm=d4$fm1[1])


#NPQ = (Fm-Fm')/Fm'

npq(fm=d4$fm1[1],fmp=d4$fm1)
npq(max(d4$fm1[1]),fmp=d4$fm1)

plot(d4$par,npq(max(d4$fm1[1]),fmp=d4$fm1))
fit_npq(d4$par,npq(max(d4$fm1[1]),fmp=d4$fm1),
        num_obs = 8)

#eff=(Fm-F)/Fm
eff(f=d4$f1,fm=d4$fm1)



#For later:
#1 - ETR calculations based in sigma values
#2 - kpi (photoinactivation)

#Read more papers...

#3 - 1st derivative, 2nd derivative
#4 - integral
#5 - light attenuation coefficients (par, color) (nls, lm)
#6 - PI models (Platt, EP, Jassby, W?) to make nice plots using user input parameters


##################################################
#MPBOM
##################################################

test_mpbom<-read.table("./data_to_test_thejesus/reflectance_test.txt",header=TRUE,sep=",")
m1<-mpbom(test_mpbom$wl,test_mpbom$reflectance,xmin=720,xmax=850)

plot(m1[[1]]$wl,m1[[1]]$abs)
m1[[2]]



##################################################
#Growth curve model
##################################################

mary<-read.table("./data_to_test_thejesus/mary",header=TRUE)

thejesus::gompertz(days=mary$days, cells = mary$cells, A = 20000000 , umax = 0.5 , l=1)

?gompertz



##################################################
#fitting NPQ models
##################################################


