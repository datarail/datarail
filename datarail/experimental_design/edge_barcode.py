import numpy as np

def _barcode_position(plate_dims):

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



def encode_barcode(barcode, plate_dims=[16,24]):
    ''' encode_barcode takes a string and convert it into a 6-bit code and output
    well names corresponding to the '1' values (positive controls)
    The string is converted to upper case and maps the ascii characters 32 (space) until 95 (underscore)'''


    dec_val = [ord(d)-32 for d in barcode]
    assert all([d<64 for d in dec_val]),  'Non-valid character in the barcode: use uppercase'
    assert all([d>=0 for d in dec_val]),  'Non-valid character in the barcode: use alphanumeric'
    bin_val = ['{0:06b}'.format(d) for d in dec_val]

    wells = _barcode_position(plate_dims)


    # save the last one for check sum
    assert len(bin_val)<len(wells), 'Too long barcode'
    if len(bin_val)<(len(wells)-1):
        # padding with spaces
        bin_val += ['000000']*(len(wells)-1-len(bin_val))
    # adding the check sum as the last word
    bin_val += [_check_sum(bin_val)]

    pos_wells = [ [wells[i][j] for j in range(0,len(bin_val[i])) if bin_val[i][j]=='1'] \
                  for i in range(0,len(bin_val))]

    return (pos_wells, bin_val)




def decode_barcode(pos_wells, plate_dims=[16,24]):
    ''' decode_barcode takes a list of positive control wells corresponding to a
    6-bit code and output the barcode as a string '''

    wells = _barcode_position(plate_dims)

    # convert a list of list in a single list
    if type(pos_wells[1]) is list:
        pos_wells = sum(pos_wells, [])

    bin_val = [ ''.join([chr(48+(w in pos_wells)) for w in wells[i]]) for i in range(0,len(wells)) ]

    chksum = _check_sum(bin_val[:-1])

    assert all([ bin_val[-1][i] == chksum[i] for i in range(0,len(chksum)) ]),\
        'Error with check sum'

    barcode = ''.join([ chr(int(b,2)+32) for b in bin_val[:-1] ])

    return barcode
