import numpy as np

from onezone.constants import asym_to_anum, anum_to_asym


def generate_nomoto(filepath = './nomoto_latex_table.dat'):
    """
    Takes the Z=0 portion of Table 2 of Nomoto et. al. 2006
    and sums over all elemental isotopes to combine a table
    that gives yields for all elements from H to Bi and includes
    a total yield and metal yield value
    """

    raw_nomoto_yields = np.genfromtxt('nomoto_latex_table.dat',
                                  delimiter = '&', names = True,
                                  dtype = "|U2,f8,f8,f8,f8,f8,f8,f8")

    nomoto_e          = np.unique(raw_nomoto_yields['element'])
    masses            = np.array([float(y) for y in raw_nomoto_yields.dtype.names[1:]])
    nomoto_yields     = np.zeros(
                           (np.size(nomoto_e),7)
                         )

    j = 0
    for i in np.arange( np.size(nomoto_e)):

        for j in np.arange(np.size(raw_nomoto_yields['element'])):

            if raw_nomoto_yields['element'][j] == nomoto_e[i]:
                nomoto_yields[i] +=  np.array(list(raw_nomoto_yields[j])[1:])


    # get atomic numbers that exist
    nomoto_anum    = np.array([asym_to_anum[str(x).strip()] for x in nomoto_e])
    all_anum       = np.arange(1,84,1) # H to Bi

    final_yields = np.zeros( (np.size(all_anum)+2, 7))

    for i in np.arange(np.size(all_anum)):
        if all_anum[i] in nomoto_anum:
            final_yields[i + 2] = nomoto_yields[nomoto_anum == all_anum[i]]

    final_yields[0] = np.sum(final_yields, axis = 0)
    final_yields[1] = np.sum(final_yields[4:], axis = 0)


    return final_yields, masses


def generate_heger(filepath = './heger_yields.dat'):
    """
    Heger & Woosley 2002 - Table 3

    """
# anum element 65 70 75 80 85 90 95 100 105 110 115 120 125 130
    raw_heger_yields = np.genfromtxt('heger_latex_table.dat',
                                     delimiter = '&', names = True,
                                     dtype = "i2,|U4,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8,f8")

    heger_e      = [str(x).strip() for x in np.unique(raw_heger_yields['element'])]
    masses       = np.array([float(y) for y in raw_heger_yields.dtype.names[2:]])
    heger_yields = np.zeros( (np.size(heger_e),np.size(masses)))

    j = 0
    for i in np.arange(np.size(heger_e)):
        for j in np.arange(np.size( raw_heger_yields['element'])):
            if str(raw_heger_yields['element'][j]).strip() == heger_e[i]:
                heger_yields[i] += np.array(list(raw_heger_yields[j])[2:])

    # get atomic numbers
    heger_anum = np.array([asym_to_anum[str(x).strip()] for x in heger_e])
    all_anum   = np.arange(1,84,1) # H to Bi

    final_yields = np.zeros( (np.size(all_anum)+2,np.size(masses)))

    for i in np.arange(np.size(all_anum)):
        if all_anum[i] in heger_anum:
            final_yields[i+2] = heger_yields[ heger_anum == all_anum[i] ]

    final_yields[0] = np.sum(final_yields, axis=0)       # all elements
    final_yields[1] = np.sum(final_yields[4:], axis = 0) # all metals

    return final_yields, masses

def generate_table(outname = './popIII_yields.dat'):

    nomoto_yields, nomoto_masses = generate_nomoto()

    heger_yields,  heger_masses = generate_heger()

    # combine into a single table

    final_table  = np.concatenate( [nomoto_yields, heger_yields], axis=1)
    masses_final = np.concatenate( [nomoto_masses, heger_masses], axis=0)

    header = "# M Z m_tot m_metal " + ' '.join( [anum_to_asym[x] for x in np.arange(1,84,1)]) + '\n'

    output = open(outname, 'w')

    output.write(header)

    metallicity = 0.0
    for j in np.arange(np.size(masses_final)):

        output.write("%.3f %.3f "%(masses_final[j], metallicity))

        for i in np.arange(np.size(final_table[:,0])):
            output.write("%8.8E "%( final_table[i,j]))

        output.write("\n")

    output.close()

    return final_table, masses_final


if __name__ == '__main__':

    generate_table()
