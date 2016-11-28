import numpy as np
import datarail.utils.plate_fcts as pltfct

def _fingerprint_position(plate_dims):

    # check for proper plate dimension
    assert (1.*plate_dims[1]/plate_dims[0])==1.5, 'Unexpected plate dimension'
    assert np.log2(plate_dims[0])%1==0, 'Unexpected plate dimension'

    if plate_dims[0]==16:
        # 384-well plate --> 4 positions on each long edge, 2 on each side
        well_groups = [ ['A%02i' % i for i in range(6*j+1,6*j+7)] for j in range(0,4)] +\
                      [ ['P%02i' % i for i in range(6*j+1,6*j+7)] for j in range(0,4)] +\
                      [ ['%s01' % chr(i+66+j*6) for i in range(0,6)] for j in range(0,2)] +\
                      [ ['%s24' % chr(i+66+j*6) for i in range(0,6)] for j in range(0,2)]
    elif plate_dims[0]==8:
        # 96-well plate --> 2 positions on each long edge, 1 on each side
        well_groups = [ ['A%02i' % i for i in range(6*j+1,6*j+7)] for j in range(0,2)] +\
                      [ ['H%02i' % i for i in range(6*j+1,6*j+7)] for j in range(0,2)] +\
                      [ ['%s01' % chr(i+66) for i in range(0,6)] ] +\
                      [ ['%s12' % chr(i+66) for i in range(0,6)] ]

    return well_groups


def _check_sum(words):
    chksum = [True]*len(words[0])
    for w in words:
        chksum = [chksum[i] != (w[i]=='1') for i in range(0,len(w))]
    return ''.join([chr(48+i) for i in chksum])



def encode_fingerprint(fingerprint, plate_dims=[16,24]):
    ''' Take a string and convert it into a 6-bit code and output
    well names corresponding to the '1' values (positive controls)
    The string is converted to upper case and maps the ascii characters 32 (space) until 95 (underscore)'''


    dec_val = [ord(d)-32 for d in fingerprint]
    assert all([d<64 for d in dec_val]),  'Non-valid character in the fingerprint: use uppercase'
    assert all([d>=0 for d in dec_val]),  'Non-valid character in the fingerprint: use alphanumeric'
    bin_val = ['{0:06b}'.format(d) for d in dec_val]

    wells = _fingerprint_position(plate_dims)


    # save the last one for check sum
    assert len(bin_val)<len(wells), 'Too long fingerprint'
    if len(bin_val)<(len(wells)-1):
        # padding with spaces
        bin_val += ['000000']*(len(wells)-1-len(bin_val))
    # adding the check sum as the last word
    bin_val += [_check_sum(bin_val)]

    pos_wells = [ [wells[i][j] for j in range(0,len(bin_val[i])) if bin_val[i][j]=='1'] \
                  for i in range(0,len(bin_val))]
    #pos_wells2 = [sum(pos_wells[i],[]) for i in range(len(pos_wells))]
    pos_wells2 = [i for j in pos_wells for i in j]
    # return pos_wells2
    return pos_wells2



def decode_fingerprint(pos_wells, plate_dims=[16,24]):
    ''' decode_fingerprint takes a list of positive control wells corresponding to a
    6-bit code and output the fingerprint as a string '''

    wells = _fingerprint_position(plate_dims)

    # convert a list of list in a single list
    if type(pos_wells[1]) is list:
        pos_wells = sum(pos_wells, [])

    # convert to 'A01' from 'A1' if needed
    pos_wells = [pltfct.well_as_2digit(w) for w in pos_wells]

    bin_val = [ ''.join([chr(48+(w in pos_wells)) for w in wells[i]]) for i in range(0,len(wells)) ]

    chksum = _check_sum(bin_val[:-1])

    chksum_valid = all([ bin_val[-1][i] == chksum[i] for i in range(0,len(chksum)) ])
    if not chksum_valid:
        print 'Error with check sum'

    fingerprint = ''.join([ chr(int(b,2)+32) for b in bin_val[:-1] ])

    return fingerprint, chksum_valid
