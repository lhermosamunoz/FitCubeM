import numpy as np

def DefRegionLines(l,l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14,sig_inicial,data_cor):
    newx1 = l[np.where(l>l1)[0][0]:np.where(l<l2)[0][-1]+1]             # SII2
    newy1 = data_cor[np.where(l>l1)[0][0]:np.where(l<l2)[0][-1]+1]
    newx2 = l[np.where(l>l3)[0][0]:np.where(l<l4)[0][-1]+1]             # SII1
    newy2 = data_cor[np.where(l>l3)[0][0]:np.where(l<l4)[0][-1]+1]
    newx3 = l[np.where(l>l5)[0][0]:np.where(l<l6)[0][-1]+1]             # NII2
    newy3 = data_cor[np.where(l>l5)[0][0]:np.where(l<l6)[0][-1]+1]
    newx4 = l[np.where(l>l7)[0][0]:np.where(l<l8)[0][-1]+1]             # Halpha
    newy4 = data_cor[np.where(l>l7)[0][0]:np.where(l<l8)[0][-1]+1]
    newx5 = l[np.where(l>l9)[0][0]:np.where(l<l10)[0][-1]+1]            # NII1
    newy5 = data_cor[np.where(l>l9)[0][0]:np.where(l<l10)[0][-1]+1]
    newx6 = l[np.where(l>l11)[0][0]:np.where(l<l12)[0][-1]+1]           # OI1
    newy6 = data_cor[np.where(l>l11)[0][0]:np.where(l<l12)[0][-1]+1]
    newx7 = l[np.where(l>l13)[0][0]:np.where(l<l14)[0][-1]+1]           # OI2
    newy7 = data_cor[np.where(l>l13)[0][0]:np.where(l<l14)[0][-1]+1]


    # Initial guesses of the fitting parameters
    sig0 = sig_inicial                  # SII2
    sig20 = 6.
    sig30 = 3.
    mu0  = newx1[np.argmax(newy1)]
    amp0 = max(newy1)
    amp20 = max(newy1)/2.
    amp30 = max(newy1)/1.5
    sig1 = sig_inicial                  # SII1
    sig21 = 6.
    sig31 = 3.
    mu1 = newx2[np.argmax(newy2)]
    amp1 = max(newy2)
    amp21 = max(newy2)/2.
    amp31 = max(newy2)/1.5
    sig2 = sig_inicial                  # NII2
    sig22 = 6.
    sig32 = 3.
    mu2 = newx3[np.argmax(newy3)]
    amp2 = max(newy3)
    amp22 = max(newy3)/2.
    amp32 = max(newy3)/1.5
    sig3 = sig_inicial                  # Halpha
    sig23 = 6.
    sig33 = 3.
    mu3 = newx4[np.argmax(newy4)]
    amp3 = max(newy4)
    amp23 = max(newy4)/2.
    amp33 = max(newy4)/1.5
    sig4 = sig_inicial                  # NII1
    sig24 = 6.
    sig34 = 3.
    mu4 = newx5[np.argmax(newy5)]
    amp4 = max(newy5)
    amp24 = max(newy5)/2.
    amp34 = max(newy5)/1.5
    sig5 = sig_inicial                  # OI1
    sig25 = 6.
    sig35 = 3.
    mu5  = newx6[np.argmax(newy6)]
    amp5 = max(newy6)
    amp25 = max(newy6)/2.
    amp35 = max(newy6)/1.5
    sig6 = sig_inicial                  # OI2
    sig26 = 6.
    sig36 = 3.
    mu6  = newx7[np.argmax(newy7)]
    amp6 = max(newy7)
    amp26 = max(newy7)/2.
    amp36 = max(newy7)/1.5

    return sig0,sig20,sig30,mu0,amp0,amp20,amp30,sig1,sig21,sig31,mu1,amp1,amp21,amp31,sig2,sig22,sig32,mu2,amp2,amp22,amp32,sig3,sig23,sig33,mu3,amp3,amp23,amp33,sig4,sig24,sig34,mu4,amp4,amp24,amp34,sig5,sig25,sig35,mu5,amp5,amp25,amp35,sig6,sig26,sig36,mu6,amp6,amp26,amp36

