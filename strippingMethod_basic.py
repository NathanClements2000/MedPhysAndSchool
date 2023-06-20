import numpy as np
import matplotlib.pyplot as plt
import xml.etree.ElementTree

'''
Python 2.7 (sorry)
Code for the stripping method
written by Chelsea Dunning (May 11, 2021)

To run this code, you will need the following:
 - Full Energy Absorption Efficiency (FEAE) curve for the detector type
    (ie. 100 120-keV x-rays incident on detector, how many of those x-rays were registered as 120 keV?)
 - Response functions of the detector to monoenergetic x-rays in increments of energy
    (ie. the detected energy spectrum in the detector from 5 keV, 6 keV, 7 keV, ..., 119 keV, 120 keV monoenergetic x-rays)

Lucky for you, I have generated a FEAE curve and the response functions of 25 mm^2 CdTe for you.
'''

def main():
    '''
    This is the main method that you will run
    :return: void
    '''
    #************************************ENTER STUFF BELOW*************************************************************

    #first, specify where the energy deposition histogram file is. This is an example, so the data here is Monte-Carlo
    #  generated (and in XML format, the helper function below puts it in a NumPy array). You will need to make a new helper
    #  function if you are reading in .mca format that's output from the real detector.
    data_directory = "/home/chelsea/topas/XFCT/monox_cdte2/stripping/broadbeam/cdte_stripping_data/"
    raw_data_file = data_directory + "example_data.xml"

    # next, specific the directory of the response functions, their file name convention, and the energy range.
    dir_response_fn = data_directory + "response_functions/"
    filename_1 = "responsefunction_" #filename in two parts because that's where energy inserted
    filename_2 = "keV5x5.npy" #if you name it better, this could just be ".npy"
    min_energy = 4.8  # keV
    max_energy = 120  # keV
    energy_bin_width = 0.2 #keV, this is what I went with
    num_primaries = 5e6  # these were the number of primaries used to generate the response functions

    #then specify the directory and file name of the FEAE curve that I've provided.
    feae_file = data_directory + "topascd2018_abseff5x5.npy"

    #*****************************THAT'S IT! below is where the magic happens******************************************

    #Let's generate a singleton of all the response functions for easy access. this will speed up the code, promise.
    MASTER_RESPONSE = generateResponseFunctionSingleton(dir_response_fn, filename_1, filename_2, min_energy, max_energy,
                                                        energy_bin_width)

    # with XFCT data, below is where you'd loop over multiple detectors and multiple position/angles of pencil beam
    #  acquisition

    #read in the example detector energy deposition histogram data
    energy, raw_spectrum = readXMLhistogram(raw_data_file, max_energy, energy_bin_width)

    #now strip the spectrum!
    stripped_spectrum = stripEnergySpectrum(raw_spectrum, MASTER_RESPONSE, min_energy, max_energy, energy_bin_width,
                                            feae_file, num_primaries)

    #let's plot raw vs. stripped spectrum
    plt.plot(energy, raw_spectrum, color="blue")
    plt.plot(energy, stripped_spectrum, color="magenta")
    #plt.bar(energy, raw_spectrum, width=energy_bin_width, color="blue")
    #plt.bar(energy, stripped_spectrum, width=energy_bin_width, color="magenta")
    plt.yscale('log')
    plt.grid(True)
    plt.ylim([1e-1, None])
    plt.xlim([0, max_energy])
    plt.xlabel("Energy [keV]")
    plt.ylabel("Counts")
    plt.title("Stripping method output")
    plt.legend(["Raw data", "stripped data"], loc=0)
    plt.show()




#*****************        HELPER FUNCTIONS      *************
def generateResponseFunctionSingleton(dir_response_fn, filename_1, filename_2, min_energy, max_energy,
                                                        energy_bin_width):
    '''
    helper function for loading all response functions right away, so only have to do this once especially if looping
    over XFCT data
    :param dir_response_fn: str
    :param filename_1: str, "responsefunction"
    :param filename_2: str, "keV5x5.npy"
    :param min_energy: float
    :param max_energy: float or int
    :param energy_bin_width: float
    :return: giant object array
    '''
    MASTER_RESPONSE = np.zeros(max_energy/energy_bin_width, dtype=object)
    for j in np.arange(min_energy, max_energy, step=energy_bin_width):
        energy, edep = np.load(dir_response_fn + filename_1 + str(j + 0.2) + filename_2)
        MASTER_RESPONSE[int(j * 5)] = edep
    return MASTER_RESPONSE

def readXMLhistogram(filename, max_energy, bin_width):
    '''
    Reads TOPAS XML file and outputs the histogram bin number and entry in a ntuple

    :return: energy bin in keV, histogram
    '''
    def missing_elements(L, start, end):
        return sorted(set(xrange(start, end + 1)).difference(L))

    length = int(int(max_energy) / bin_width)
    energy1 = np.array([])
    energydep1 = np.array([])

    e = xml.etree.ElementTree.parse(filename)

    for btype in e.findall('./histogram1d/data1d/bin1d'):
        energy1 = np.append(energy1, float(btype.get('binNum')))
        energydep1 = np.append(energydep1, float(btype.get('entries')))

    missing_indices = missing_elements(energy1, 0, length - 1)
    for i in missing_indices:
        energy1 = np.insert(energy1, i, i)
        energydep1 = np.insert(energydep1, i, 0)
    energy = energy1 * bin_width

    energydep = energydep1

    return energy, energydep

def stripEnergySpectrum(spectrum, MASTER_RESPONSE, min_energy, max_energy, width_energy, dirfile_feae,
                        number_of_primaries):
    '''
    Corrects a energy-deposited spectrum using the stripping method and returns the spectrum of x-rays that were
    incident on the CdTe crystal
    :return: corrected x-ray spectrum of len(spectrum)
    '''

    # load the full energy absorption efficiency
    topas_energy, topas_abseff = np.load(dirfile_feae)

    max_index = max_energy/width_energy
    n_true = np.zeros(int(max_index))

    for i in np.arange(min_energy, max_energy, step=width_energy)[::-1]:
        i_index = int(i/width_energy)
        n_det = np.take(spectrum, i_index)
        j = i + width_energy
        j_index = int(j/width_energy)

        while j_index < max_index:
            edep = np.take(MASTER_RESPONSE, j_index)
            energy = np.arange(0, j, step=width_energy)

            edep_norm = edep / float(number_of_primaries)
            energy = np.around(energy, decimals=1)
            i = np.around(i, decimals=1)
            resp = np.take(edep_norm, np.where(np.equal(energy, i))[0][0])

            n_det = n_det - (resp * np.take(n_true, j_index))
            j_index = j_index + 1
            j = j + width_energy

        topas_energy = np.around(topas_energy, 1)
        i_width = np.around(i + width_energy, 1)
        index_feae = np.where(np.equal(topas_energy, i_width))[0][0]
        if n_det < 0:
            n_true[i_index] = 0
        else:
            n_true[i_index] = 100 * n_det / np.take(topas_abseff,  index_feae)

    return n_true

if __name__ == "__main__":
    main()