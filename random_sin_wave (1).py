import numpy as np
import requests
import pandas as pd
import matplotlib.pyplot as plt
from scipy.signal import find_peaks
from sympy import symbols, cosh, sinh, sqrt
nb_times_extinct_A=[]
nb_times_extinct_B=[]

o=int(input("how many total reps"))
z=int(input("how many inside reps"))
for p in range(0,o):
  nb_species_A_list_time_dep=[]
  nb_species_B_list_time_dep=[]
  nb_species_A_list_time_dep_after_division=[]
  nb_species_B_list_time_dep_after_division=[]
  nb_species_A_not_rounded_list=[]
  nb_species_B_not_rounded_list=[]
  pb_quantum_tunneling_lithium_elec_1_list=[]
  kinetic_energy_electron_lithium_1=[]



  extinct_A = False
  extinct_B = False
  for d in range(1,z): #z can be represented as the total number of years (or number million years) over which we are doing this

    hbar = 1.054571817e-34
    conversion = 6.241509e18
    impedance_of_vaccum=376.73
    electron_number=1
    lithium_elec_1=0
    Vp=13
    lithium_elec_2=0
    beryllium_elec_1=0
    beryllium_elec_2=0
    total_photons=[]
    ionization_atoms=[[13.59844],[24.58738,54.41776],[5.39171, 75.64009, 122.45435], [9.32269,18.21115, 153.89620, 217.71858], [8.29803, 25.15484, 37.93064, 259.37521, 340.22580], [11.26030, 24.38332, 47.8878, 64.4939, 392.087, 489.99334], [14.53414, 29.6013, 47.44924, 77.4735, 97.8902, 552.0718, 667.046], [13.61806, 35.11730, 54.9355, 77.41353, 113.8990, 138.1197, 739.29, 871.4101], [17.42282, 34.97082, 62.7084, 87.1398, 114.2428, 157.1651, 185.186, 953.9112, 1103.1176], [21.5646, 40.96328, 63.45, 97.12, 126.21, 157.93, 207.2759, 239.0989, 1195.8286, 1362.199	], [5.13908, 47.2864, 71.6200, 98.91, 138.40, 172.18, 208.50, 264.25, 299.864, 1465.121, 1648.702], [7.64624, 15.03528, 80.1437, 109.2655, 141.27, 186.76, 225.02, 265.96, 328.06, 367.50, 1761.805, 1962.6650], [5.98577, 18.82856, 28.44765, 119.992, 153.825, 190.49, 241.76, 284.66, 330.13, 398.75, 442.00, 2085.98, 2304.1410], [8.15169, 16.34585, 33.49302, 45.14181, 166.767, 205.27, 246.5, 303.54, 351.12, 401.37, 476.36, 523.42, 2437.63, 2673.182], [10.48669, 19.7694, 30.2027, 51.4439, 65.0251, 220.421, 263.57, 309.60, 372.13, 424.4, 479.46, 560.8, 611.74, 2816.91, 3069.842]]
    m=9.1093837e-31 #mass electron
    L=1.084e-10
    intermediate=(L*sqrt(2*m*1.60218e-19))/hbar






    def generate_random_sine_wave(time):
        base_frequency = np.random.uniform(hbar, 23)
        frequency_variation = np.sin(1000*time * np.random.uniform(0.01, 10)) * 2
        amplitude = np.random.uniform(1, 4)
        phase = np.random.uniform(0, 14 * np.pi)
        # frequency fluctuations
        frequency = base_frequency + frequency_variation
        return amplitude * np.sin(2 * np.pi * 10**12 * frequency * 1000*time + phase)

    time = np.linspace(0, d, 5000*d)
    num_waves = np.random.randint(3, 100)

    combined_wave = sum(generate_random_sine_wave(1000*time) for _ in range(num_waves))
    offset = abs(min(combined_wave))
    combined_wave += offset

    peaks, _ = find_peaks(combined_wave, height=0)
    number_of_peaks = len(peaks)
    peak_amplitudes = combined_wave[peaks]
    integer_amplitudes = np.rint(peak_amplitudes)
    periods = np.diff(time[peaks])



    file_path = '/content/TSLA (1).csv'
    df = pd.read_csv(file_path)


    # save data to file
    with open('wave_data.csv', 'w') as file, open('electron_data.csv', 'w') as electron_file:
        file.write("Peak Number,Amplitude,Period (ns),Frequency (Hz),Energy (eV), Photons Flux in J/s\n")
        electron_file.write("Peak Number,Energy (eV),electron number, atom,ELECTRON\n")
        for i, (amp, period) in enumerate(zip(integer_amplitudes, periods), start=1):
            frequency = 6*10**9 / (period) # multiplied by 6 because had calculated with h not hbar
            wrapped_index = (i-1) % len(df) # make it wrap or else pb of size
            energy = 1000*conversion * hbar * frequency - float(df.iloc[wrapped_index]['Close'])/150 #famous equation, 1000 multiplication since divided the linspace by 1000 to make the computation go faster, wrapped index to go back around and getting tesla share market
            if ionization_atoms[2][0] < energy : #the first ionization energy required of lithium is smaller than the energy carried by the electron
                kinetic_energy_electron=energy-ionization_atoms[2][0] #the kinetic energy of the electron is the difference of the total energy given by the light and the energy that was required to unvound the electon

                pb_quantum_tunneling_lithium_elec_1_list.append(pb_quantum_tunneling_lithium_elec_1)
                kinetic_energy_electron_lithium_1.append(kinetic_energy_electron)
                pb_quantum_tunneling_lithithium_elec_1=1/((cosh(sqrt((Vp-kinetic_energy_electron_lithium_1[wrapped_index]))*intermediate))**2 + (1/4)*(Vp-kinetic_energy_electron_lithium_1[wrapped_index] + 1/(Vp/kinetic_energy_electron_lithium_1[wrapped_index] - 1) - 3)*(sinh(intermediate*sqrt(Vp-kinetic_energy_electron_lithium_1[wrapped_index])))**2)
                pb_quantum_tunneling_lithium_elec_1_list.append(pb_quantum_tunneling_lithium_elec_1 )
                electron_file.write(f"{i},{kinetic_energy_electron:.7f}, {ionization_atoms[2][0]},{electron_number}, lithium ionization 1, ELECTRON\n")
                electron_number=electron_number+1
                lithium_elec_1=lithium_elec_1+1
                if energy - ionization_atoms[2][0] > ionization_atoms[2][1]:
                  kinetic_energy_electron=energy-ionization_atoms[2][0]-ionization_atoms[2][1]#the total energy given minus both energy required to unbind both outer electrons and the rest is the kinetic energy of the free unbound electron
                  electron_file.write(f"{i},{kinetic_energy_electron:.7f},{ionization_atoms[2][1]},{electron_number}, lithium ionization 2, ELECTRON\n")
                  electron_number=electron_number+1
                  lithium_elec_2=lithium_elec_2+1
            if ionization_atoms[3][0]< energy:
              kinetic_energy_electron=energy-ionization_atoms[3][0]
              electron_file.write(f"{i},{kinetic_energy_electron:.7f}, {ionization_atoms[3][0]},{electron_number}, beryllium ionization 1, ELECTRON\n")
              electron_number=electron_number+1
              beryllium_elec_1=beryllium_elec_1+1
              if energy- ionization_atoms[3][0]>ionization_atoms[3][1]:
                kinetic_energy_electron=energy-ionization_atoms[3][0]-ionization_atoms[3][1]
                electron_file.write(f"{i},{kinetic_energy_electron:.7f}, {ionization_atoms[3][1]},{electron_number}, beryllium ionization 2, ELECTRON\n")
                electron_number=electron_number+1
                beryllium_elec_2=beryllium_elec_2+1
            number_of_photons=(amp**2)/(hbar*frequency*impedance_of_vaccum)
            total_photons.append(number_of_photons)
            file.write(f"{i},{int(amp)},{period},{frequency},{kinetic_energy_electron:.7f},{number_of_photons}\n")
    nb_photons=1
    for l in range(len(total_photons)):
      nb_photons=nb_photons+total_photons[l]
    print("total flux of electrons comming out is ", nb_photons, "electrons per second") #calculating the total flux of electrons that have unbinded from their atoms
    photons_per_peak=(nb_photons*10**-9*2)/(len(total_photons)-1)
    print(photons_per_peak, "is the number of photons per peak")

    for t in range(len(kinetic_energy_electron_lithium_1)): # calculating the number of electrons that quantum tunnnled
      pb_quantum_tunneling_lithium_elec_1_list=[t]

    #plt.plot(kinetic_energy_electron_lithium_1, pb_quantum_tunneling_lithium_elec_1_list)
    '''print(len(kinetic_energy_electron_lithium_1), len(pb_quantum_tunneling_lithium_elec_1_list))
    print(pb_quantum_tunneling_lithithium_elec_1_list, kinetic_energy_electron_lithium_1)'''








    #calculating the number of electrons that have each energy. We assume that every 10**15 electrons in the box(once they unbind from their atom they are all in a box), represents 1 animal of a species
    print('nb of electrons released from lithium on ionization 1 is', lithium_elec_1*photons_per_peak)
    print('nb of electrons released from lithium on ionization 2 is', lithium_elec_2*photons_per_peak)
    print("so population size of species A before quantum tunneling is ", (lithium_elec_1*photons_per_peak)/(10**(15)) + (lithium_elec_2*photons_per_peak)/(10**(15)), "which we round to", round((lithium_elec_1*photons_per_peak)/(10**(15)) + (lithium_elec_2*photons_per_peak)/(10**(15))) )
    print('nb of electrons released from beryllium on ionization 1 is', beryllium_elec_1*photons_per_peak)
    print('nb of electrons released from beryllium on ionization 2 with is', beryllium_elec_2*photons_per_peak)
    print("so population size of species B before quantum tunneling is ", (beryllium_elec_1*photons_per_peak)/(10**(15))+  (beryllium_elec_2*photons_per_peak)/(10**(15)), 'which we round to',round((beryllium_elec_1*photons_per_peak)/(10**(15))+  (beryllium_elec_2*photons_per_peak)/(10**(15))) )


    E_A_1=5.39171 #first ionization energy for lithium which is the energy assumed to ebe for all the electrons that unbinded in first ionization in lithium
    E_A_2=75.64009 #second ionization energy for lithium __
    E_B_1=9.32269 #first ionization energy for beryllium
    E_B_2=18.21115 #second ionization rate for beryllium


    #calculating the probability that the electrons with different energies that unbinded from their atoms quantum tunnel out of the box whic can be seen as a potential barrier
    pb_quantum_tunneling_lithium_elec_1= 1/((cosh(sqrt((Vp-E_A_1))*intermediate))**2 + (1/4)*(Vp/E_A_1 + 1/(Vp/E_A_1 - 1) - 3)*(sinh(intermediate*sqrt(Vp-E_A_1)))**2)
    pb_quantum_tunneling_lithium_elec_2=1 /((cosh(L* sqrt(2*m*(Vp-E_A_2)*1.60218e-19)/hbar))**2 + (1/4)*(Vp/E_A_2 + 1/(Vp/E_A_2 - 1) - 3)*(sinh(L* sqrt(2*m*(Vp-E_A_2)*1.60218e-19)/hbar))**2)
    pb_quantum_tunneling_beryllium_elec_1=1 /((cosh(intermediate*sqrt(Vp-E_B_1)))**2 + (1/4)*(Vp/E_B_1 + 1/(Vp/E_B_1 - 1) - 3)*(sinh(intermediate*sqrt(Vp-E_B_1)))**2)
    pb_quantum_tunneling_beryllium_elec_2=1 /((cosh(sqrt(2*m*(Vp-E_B_2)*1.60218e-19)/(hbar) * L))**2 + (1/4)*((Vp/E_B_2) + 1/((Vp/E_B_2) - 1) - 3)*(sinh(L* sqrt(2*m*(Vp-E_B_2)*1.60218e-19)/(hbar)))**2)

    u= np.random.uniform(-0.1, 0.1)

    #calculating the total number of  electorns, multiply by u for soime randomness/uncertainty arbitrairly up to 10%
    nb_electron_lithium_quantum_tunnel=(((pb_quantum_tunneling_lithium_elec_1*u)+pb_quantum_tunneling_lithium_elec_1)*lithium_elec_1*photons_per_peak+((pb_quantum_tunneling_lithium_elec_2*u)+pb_quantum_tunneling_lithium_elec_2)*lithium_elec_2*photons_per_peak)
    nb_electron_beryllium_quantum_tunnel=(((pb_quantum_tunneling_beryllium_elec_1*u)+pb_quantum_tunneling_beryllium_elec_1)*beryllium_elec_1*photons_per_peak+((pb_quantum_tunneling_beryllium_elec_2*u)+pb_quantum_tunneling_beryllium_elec_2)*beryllium_elec_2*photons_per_peak)
    print('pb quantum tunneling for an electron ionization 1 coming from a beryllium atom is',pb_quantum_tunneling_beryllium_elec_1, 'and probability quantum tunneling for an electron coming from a lithium ionization 1 atom is', pb_quantum_tunneling_lithium_elec_1)
    print('nb_electron lithium quantum tunnel: ',nb_electron_lithium_quantum_tunnel,'nb_electron beryllium quantum tunnel: ', nb_electron_beryllium_quantum_tunnel)

    print("pb_quantum_tunneling_beryllium_elec_1", pb_quantum_tunneling_lithium_elec_1)
    print("lithium_elec_1", lithium_elec_1)
    print("photons per peak", photons_per_peak)
    print("pb_quantum_tunneling_beryllium_elec_1", pb_quantum_tunneling_beryllium_elec_1)
    print("beryllium_elec_1", beryllium_elec_1)



    #calculating the total number of animals in eahc species by substracting nulmber that quantum tunnel
    print("so there are a total of",((lithium_elec_1*photons_per_peak +lithium_elec_2*photons_per_peak-nb_electron_lithium_quantum_tunnel)), "electrons coming from lithium atoms in the box at t=", d, "ns" )
    print("and there are a total of ", (beryllium_elec_1*photons_per_peak + beryllium_elec_2*photons_per_peak- nb_electron_beryllium_quantum_tunnel), "electrons from beryllium atoms in the box at t=",d,"ns")

    nb_species_A_not_rounded_list.append((lithium_elec_1*photons_per_peak +lithium_elec_2*photons_per_peak-nb_electron_lithium_quantum_tunnel)/(10**(15)))
    nb_species_B_not_rounded_list.append((beryllium_elec_1*photons_per_peak + beryllium_elec_2*photons_per_peak- nb_electron_beryllium_quantum_tunnel)/(10**(15)))

    nb_species_A=round((lithium_elec_1*photons_per_peak +lithium_elec_2*photons_per_peak-nb_electron_lithium_quantum_tunnel)/(10**(15)))
    nb_species_B=round((beryllium_elec_1*photons_per_peak + beryllium_elec_2*photons_per_peak- nb_electron_beryllium_quantum_tunnel)/(10**(15)))
    #since we are doing by second, every iteration say 4 ns , we divide by 4 so that average for 1 ns and add it to the previous one.
    #if there are less than 2 animals of a certain species, then that species cannot mate and thus we say that that species went extinct
    nb_species_A_list_time_dep.append(nb_species_A)
    nb_species_B_list_time_dep.append(nb_species_B)
    if d>1:
      nb_species_after_division_A=nb_species_A_not_rounded_list[d-2]/(d-1)+nb_species_after_division_A
      nb_species_after_division_B=nb_species_B_not_rounded_list[d-2]/(d-1)+nb_species_after_division_B
    else:
      nb_species_after_division_A=nb_species_A
      nb_species_after_division_B=nb_species_B

    if d>2:
      nb_species_after_division_A=nb_species_A_not_rounded_list[d-2]/(d-1)+nb_species_after_division_A-nb_species_A_not_rounded_list[d-1]/d #substract the animals who where born 2 intervals ago since they probably died.
      nb_species_after_division_B=nb_species_B_not_rounded_list[d-2]/(d-1)+nb_species_after_division_B-nb_species_B_not_rounded_list[d-1]/d


    if nb_species_after_division_A < 2:
        extinct_A = True
        #print("Species A is extinct")

    else:
        #print("species A is fine")
        extinct_B=False

    if nb_species_after_division_B < 2:
        extinct_B = True
        #print("Species B is extinct")
    else:
        #print("species B is fine")
        extinct_B=False

    if extinct_A:
          nb_species_after_division_A = -1 #once extinct goes to -1
    if extinct_B:
          nb_species_after_division_B = -1
    nb_species_A_list_time_dep_after_division.append(nb_species_after_division_A)
    nb_species_B_list_time_dep_after_division.append(nb_species_after_division_B)





    print("since population size of species A is", round(nb_species_after_division_A),"and for population size of species B", round(nb_species_after_division_B))
    print(nb_species_after_division_A, nb_species_after_division_B)



    plt.figure(figsize=(10, 6))
    plt.plot(time, combined_wave)
    plt.plot(time[peaks], combined_wave[peaks], "x")
    plt.title('Combined Random Sine Waves with Varying Frequencies')
    plt.xlabel('Time (nanoseconds)')
    plt.ylabel('Amplitude')
    plt.grid(True)
    plt.show()

  print('for species A', nb_species_A_list_time_dep_after_division)
  print('for species B', nb_species_B_list_time_dep_after_division)

  if nb_species_A_list_time_dep_after_division[-1]==-1:
    nb_times_extinct_A.append(1)
  else:
    nb_times_extinct_A.append(0)

  if nb_species_B_list_time_dep_after_division[-1]==-1:
      nb_times_extinct_B.append(1)
  else:
    nb_times_extinct_B.append(0)

print("1 is when species A got extinct", nb_times_extinct_A)
print("1 is when species B got extinct", nb_times_extinct_B)
sum_A=0
sum_B=0
for f in range(0,len(nb_times_extinct_A)):
  sum_A=sum_A+nb_times_extinct_A[f]
  sum_B=sum_B+nb_times_extinct_B[f]
print("frequency A goes extinct is", sum_A/(len(nb_times_extinct_A)), "and frequency B goes extinct is ", sum_B/len(nb_times_extinct_B))

