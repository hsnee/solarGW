def show_res(starttime):
    import numpy as np
    from matplotlib.backends.backend_pdf import PdfPages
    import matplotlib.pyplot as plt
    fname = str(starttime)+'.pdf'
    h0 = np.array(np.loadtxt('h0'+str(starttime)+'.txt',dtype='f8'))
    p  = np.array(np.loadtxt('p'+str(starttime)+'.txt',dtype='f8'))
    p = np.exp(p-np.max(p))
    pathtointersect = '../../../intersect.txt'
    timearray = np.array(np.loadtxt(pathtointersect,dtype='f8'))
    StartTimes = timearray[:,0]
    EndTimes   = timearray[:,1]
    for i in range(len(StartTimes)):
        if starttime == StartTimes[i]:
            idx=i
        else:
            print 'Error: starttime unrecognised', starttime
    duration = EndTimes[idx] - StartTimes[idx]
    with PdfPages(fname) as pdf:
        fig1 = plt.figure()
        plt.plot(h0,p)
        plt.plot(h0,p,'+')
        
        plt.title('Probability Distribution for '+str(duration/60)+' minutes')
        plt.xlabel(r'$h_0$')
        plt.ylabel('Normalised probability')
        pdf.savefig(fig1)
        plt.close()
