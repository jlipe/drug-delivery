function [ dy ] = pk_intestine_oil_lapatinib( t,y )
%total Summary of this function goes here
%   simulation of intestinal dissolution (y(1)), absorption and elimination (y(2)) of
%   lapatinib in fed state in the absence of digestion --> output: drug conc in plasma = y(5)
%   y(3) - drug concentration in lumen, y(4) - fatty acid content in lumen
%   y(5) - drug concentration in plasma
%   y(6) - volume of lumen that drops after Tmax is reached

%constant parameters;

Dp=4.51e-10           #diffusion coefficient of solid drug in water (m2/s)
Dm=1.44e-10           #diffusion coefficient of micelles (r = 22.7 A) in water (m2/s)
V=2.5e-4                #bulk volume (m3), 250 ml;
h=2.0e-5                #static layer around particles (m), 20 um
Satw=7.0e3             #drug solubility in maleate buffer (mg/m3), 0.007 mg/ml = 0.007 mg/cm3
M0=1500                   #initial dose of drug (mg);
dens=1.318e9          #drug density (mg/m3) --> 1.318 g/cm3;
r0=1.0e-5               #initial drug particle radius (m) ---> assumed 10 um;
Chi=55.63               #solubilization of drug per mmol of surfactant (mg/mmol); (808 ug/mL - 7 ug/mL)/(12 mmol/L + 4 mmol/L - 1.6 mmoL/L)
Vmax=0.00044                 #Michaelis-Menten parameter for the absorption of fatty acids: units = mmol/cm3/s = M/s
Km=0.0022                   #Michaelis-Menten parameter for the absorption of fatty acids: units = mmol/cm3 = M

#drug uptake in oil and oil digestion;
#NOTE: some parameters need to be changed based on oil dosed and Vaq;
P=0e-8              #permeability oil/aqueous interface  (cm/s);
V_oil=54            #initial volume of oil (cm3); 2g fat for low fat breakfast and 54g fat for high fat breakfast (Koch 2009)
Satoil=0.0136            #solubility of drug in oil mg/cm3 = mg/mL;
FA0=2*V_oil*0.92/870*1000                #initial digestible fatty acids (mmol); Molarity of soybean oil = 870 g/mol
                           #NOTE: FA0 is 2*initial mmol oil dosed;
Kdig=3.6e-9             #digestion constant, in mmol/sec*cm2;
Kinh=4.3e-4             #digestion constant inhibition, 4.3e-4 1/s;
Vaq=V*10**6                  #bulk volume (cm3), 40 ml;
D0=3.86e-5                 #initial diameter of oil droplets (cm), 386 nm


Pgi=1.01e-6          #intestinal permeability (m/s), from LP in FS, pH 6.5
rgi=1.75e-2             #intestinal radius (m), 1.75 cm;

Agi=2*V/rgi            #intestinal surface area (m2), calculated as lateral surface
                               #area of a cylinder;

Vd=80e-3              #volume of distribution (m3) (don't have one for LP... they say > 2100 mL)
t_half=13.4*60*60        #drug elimination in metabolism (s), 13.4 hours for lapatinib in Humans in fasted state from Koch et al. 2009

k_travel=0.0005;         #fluids travel time to exit intestinal;
                           #Note: if V_aq is changed a lot, time to exit the
                           #small intestine might change, in this case k_travel needs
                           #to be changed accordingly

#define number of ODEs in the system;
dy=zeros(6,1)

N=V_oil/(D0^3*3.14*1/6)      #number of oil droplets;

#calculation constant parameters in dissolution equation;
par=(3*(M0^(1/3))/(h*dens*r0))

#calculation of solubility in micelles time-dependent, function of y(4);
#1e6 conversion from mg/cm3 to mg/m3;
Satm=(Chi*(14.4e-3+y(4)))*1e6

#calculation of Kmw, partition coeff micelles/buffer time dependent;
#1e6 conversion from mg/cm3 to mg/m3;
Kmw=Satm/Satw

#calculation of Koaq, partition coeff oil/aqueous;
Koaq=Satoil/((Chi*(14.4e-3+y(4)))+(Satw*10^(-6)))

#calculation of total initial oil surface area (cm2);
Aoil=N*3.14*((D0)^2)*(((FA0-(y(4)*Vaq))/FA0)^(2/3))
#unit check: N - unitless; D0 - cm; FA0 - mmol; y(4) - mmol/mL, Vaq - mL
#Aoil - cm2

kel=0.693/t_half            #elimination constant (1/s);


if t<14400:
    #t is residence time (sec)+ emptying time (sec) in fed state,
                   #3 hours+1 hour to exit the small intestine in Humans;
             #NOTE: residence time + emptying time might need to be adjusted for Rats;
    if y(1)>0:
        #y(1) is mass of drug in solid (mg);
        if (Satm-((y(3)*Kmw)/(Kmw+1)))>0:
            dy(1)= -(par*y(1)^(2/3)*(Dp*(Satw-(y(3)/(Kmw+1)))+Dm*(Satm-((y(3)*Kmw)/(Kmw+1)))))
        else:
            dy(1)=0


    #drug uptake in oil, y(2) is mass (mg)
    dy(2)=P*Aoil*y(3)*1e-6;%-(y(2)/(V_oil*Koaq)))

    if y(1)>0:
        #y(3) is drug conc. in aqueous medium in the lumen (mg/m3)
        dy(3)=(par/V*((y(1))^(2/3))*(Dp*(Satw-(y(3)/(Kmw+1)))+Dm*(Satm-((y(3)*Kmw)/(Kmw+1)))))-dy(2)/V-Agi*Pgi*(y(3)/Vd)
    else:
        dy(3)=-(Agi*Pgi*y(3)/Vd)


if FA0-(y(4)*Vaq)>0:
        #oil digestion, y(4) is concentration of FA released - concentration of FA absorbed (mmol/cm3 = M);
        dy(4)=Kdig/Vaq*Aoil-Kinh*y(4)-Vmax*y(4)/(Km+y(4))

if t<10800:   #t is residence time (sec), 3 hours in human small intestine in fed state;
             #residence time might need to be adjusted for Rats;

    #y(5) is drug conc. in plasma, mg/m3 = ug/L = ng/mL;
    dy(5)=(Agi*Pgi*y(3)/Vd)-(kel*y(5))
else:
    dy(6)=-k_travel*y(6)
    dy(5)=(y(6)*Pgi*y(3)/Vd)-(kel*y(5))



#y(3) is in mg/m3 = ng/cm3 = ng/mL
#to convert to ug/mL divide by 10^3

#To recall from the command prompt:
#options=odeset('RelTol',2.5e-13,'AbsTol',[2.5e-13 2.5e-13 2.5e-13 2.5e-13 2.5e-13 2.5e-13],'nonnegative',[1 2 3 4 5 6]);
#[T,Y]=ode45(@pk_intestine_oil_lapatinib,[1 18000],[1500 0 0 0 0 0.0228571],options);
     #where first number (250) is initial value for Dose: if Dose is changed in code,
     #first number changes accordingly;
     #last number (0.0228571) is initial value for Agi calculated in line 61 (surface area of intestine);
     #Agi depends on V=volume aqueous, if V is changed, Agi needs to be
     #changed here;
