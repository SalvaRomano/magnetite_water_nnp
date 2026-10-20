import numpy as np
import sys
def main(filename):

   data = np.loadtxt(filename)

   data_corrected = np.zeros_like(data)
   data_corrected[:,0] = data[:,0]
   frequency = data[:,0]
   

   for col in [1,2,3,4,5,6,7]: 
       data_corrected[:,col] = data[:,col] *  mu_correction(frequency) * alpha_correction(frequency)
   #np.savetxt("corrected_"+filename,data_corrected,fmt="%.6e")
   np.savetxt("corrected_"+filename,data_corrected)
   print("")
   print("Corrected data saved to: corrected_"+filename)
   print("")

def mu_correction(omega):
   """ Compute dipole non-condon correction """
   return 1.377 + 53.03 *( (3737-omega)/6932.2 )

def alpha_correction(omega):
   """ Compute polarizability non-condon correction """
   return 1.271 + 6.287 *( (3737-omega)/6932.2 )


if __name__ == "__main__":
    filename = sys.argv[1]
    main(filename)
