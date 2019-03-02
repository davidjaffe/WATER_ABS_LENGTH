#!/usr/bin/env python
'''
re-implement model of water absorption length in dyb-doc-992-v3
units mm, MeV, nm
20190302
'''
import math
import sys
import random

import datetime

import numpy

import matplotlib.pyplot as plt

class abslen_model():
    def __init__(self):

        self.debug = 1

        self.dataDir = 'data_water_abslen/'
        self.figDir  = 'figures_water_abslen/'
        
        self.wlRange = [200., 800. ] # range of wavelengths

        self.modelPoints = None

        self.createdModels = {}

        self.filePrefixes = []
        for i in range(1,14):
            self.filePrefixes.append('abs'+ str(i))

        return
    def read_all(self):
        '''
        read all the data files
        '''
        all = {}
        for fpre in self.filePrefixes:
            L,authors = self.read_ald(fpre)
            all[fpre] = [L,authors]
        if self.debug>1: print all
        return all
    def read_ald(self,fpre):
        '''
        return lists of lists of wavelength and absorption coeff 1/cm from data files for water
        also return author names as string
        '''
        fn = fpre + '.vec'
        filename = self.dataDir + fn 
        f = open(filename,'r')
        if self.debug>1: print 'read_ald fpre',fpre
        wl,al = [],[]
        authors = None
        for line in f:
            if authors is None:
                last = line.find('"')
                if last<0: last = len(line)
                authors = line[1:last]
            if '*' not in line and len(line)>1:
                vars = line.split()
                if self.debug>1:
                    print 'read_ald line',line[:-1],'len(line)',len(line)
                    print 'read_ald vars',vars
                wl.append(float(vars[0]))
                al.append(float(vars[1]))

        f.close()
        return [wl,al], authors
    def work(self):
        '''
        get all data in dict
        dict[name] = [ [ [wl], [abslen] ], authors ]
        '''
        allData = self.read_all()
        for title in allData:
            L, authors = allData[title]
            wl,al = L
            self.drawIt(wl,al,'Wavelength(nm)','Absorption(1/cm)',title + " " + authors)
        return
    def makeModelPoints(self,drawModel=True):
        '''
        Return via global self.modelPoints a list of wavelengths and absorption.
        Use data points to make model of water absorption based on dyb-doc-992-v3.
        The model uses data from the following references in the specified wavelength ranges
        Quickenden and Irvin (200-320 nm), Pope and Fry  (380-727.5 nm) and Sullivan  (728 - 790 nm).
        The Sullivan absorption(790nm) is duplicated in the model at 800nm. ( absorption(800nm) = Sullivan absorption(790nm) )

        Smooth the Quickenden and Irvin data to remove fluctuations for 250-320nm

        The smoothed data points are written out to a file.
        '''

        QandI = 'abs4'
        PandF = 'abs2'
        Sull  = 'abs9'
        WLrange = { QandI: [200.,320.], PandF: [380.,727.5], Sull: [728.,790.] }
        WLend = 800.
        self.refsUsed = titles = [QandI, PandF, Sull]
        
        allData = self.read_all()
        
        Wavelength,Absorption = [],[]
        for title in titles:
            L, authors = allData[title]
            print 'abslen_model.makeModelPoints use',title,authors
            if title==QandI:
                for k in [1,2,3,4]:
                    LL = self.smooth(L[1],4)
                    LL[-1] = L[1][-1]
                    LL[0:40]= L[1][0:40]
                    L[1]=LL
                print 'abs_model.makeModelPoints smooth',title
            deltaWL,bestAL = 1.e20,None
            for wl,al in zip(L[0],L[1]):
                if WLrange[title][0] <= wl <= WLrange[title][1]:
                    Wavelength.append(wl)
                    Absorption.append(al)
                if title==Sull: # model adds 790nm point of Sullivan as 800nm
                    if abs(wl-WLend)<deltaWL :
                        bestAL = al
                        deltaWL = abs(wl-WLend)
            if bestAL is not None:
                Wavelength.append(WLend)
                Absorption.append(bestAL)
        # sort in ascending wavelength
        a = [x for _,x in sorted(zip(Wavelength,Absorption))]
        smooth = False
        if smooth:
            ns = 3
            w = self.smooth(sorted(Wavelength),ns)
            print 'abslen_model.makeModelPoints smooth with ns',ns
        else:
            w = sorted(Wavelength)
        
        if drawModel: self.drawIt(w,a,'Wavelength(nm)','Absorption(1/cm)','Model points')
        
        self.modelPoints = [w,a]

        outfn = self.dataDir + 'abslen_modelpoints.txt'
        f = open(outfn,'w')
        import datetime
        f.write('* ' + str(datetime.datetime.now())+'\n')
        f.write('* abslen_model.py model points\n')
        s = '* '
        for ref in self.refsUsed:
            authors = allData[ref][1]
            s += authors + ' '
        s += '\n'
        f.write(s)
        f.write('* wavelength(nm) absorption(1/cm)\n')
        for wl,al in zip(self.modelPoints[0],self.modelPoints[1]):
            f.write('  ' + str(wl) + '    ' + str(al) + '\n')
        f.close()
        print 'abslen_model.py Wrote out modelPoints to',outfn
    
        
        return
    def createModel(self,A200=0.001):
        '''
        absorption = a
        a = amin + ai, ai = ai(200nm)/(WL/200nm)^4
        amin = interpolated absorption from modelPoints
        ai(200nm) = input parameter of absorption in 1/cm

        results are returned in global self.createdModels
        '''

        # if model exists, then don't recreate it
        if A200 in self.createdModels: return

        # get data points if needed
        if self.modelPoints is None: self.makeModelPoints(drawModel=False)
            
        WL,AL = self.modelPoints
        wlModel = [float(q) for q in range(int(self.wlRange[0]), int(self.wlRange[1]), 1)]
        alModel = []

        lenWL = len(WL)
        for wl in wlModel:
            j = min(range(len(WL)), key=lambda i: abs(WL[i]-wl))
            if WL[j]<wl:
                k = min(j + 1,lenWL)
            elif WL[j]>wl:
                k = max(0,j-1)
            else:
                k = j
            if k==j:
                amin = AL[k]
            else:
                amin = ((AL[k]-AL[j])*wl + (AL[j]*WL[k]-AL[k]*WL[j]))/(WL[k]-WL[j])
            ai = A200/pow(wl/200.,4)
            alModel.append(amin + ai)
        self.createdModels[A200] = [ wlModel, alModel ]
        if self.debug>1: self.drawIt(wlModel,alModel,'Wavelength(nm)','Absorption(1/cm)','Model for ai(200nm)='+str(A200))
        return
    def smooth(self,y, box_pts):
        box = numpy.ones(box_pts)/box_pts
        y_smooth = numpy.convolve(y, box, mode='same')
        return y_smooth

    def drawExample(self):
        '''
        draw data points and model lines
        '''
        A200 = [0., 0.001, 0.010, 0.05]
        lstyle=[':','-','-.','--']
        cstyle=['r','g','c','m']
        for a200 in A200:
            self.createModel(A200=a200)

        plt.clf()
        fig,ax = plt.subplots()

        

        ax.set_title('Water absorption data and model curves')
        ax.set_xlabel('Wavelength(nm)')
        ax.set_ylabel('Absorption(1/cm)')
        ax.set_yscale('log')
        ax.set_xlim(200.,800.)
        ax.set_ylim(1.e-5, 1.e-1)
        special = False
        if special:
            ax.set_xlim(200.,400.)
            ax.set_ylim(1.e-4,1.e-3)
        
        WL,AL = self.modelPoints
        X,Y = numpy.array(WL),numpy.array(AL)
        ax.plot(X,Y,'o',color='y',label='Model points')

        allData = self.read_all()
        for i,title in enumerate(self.refsUsed):
            L, authors = allData[title]
            X,Y = numpy.array(L[0]),numpy.array(L[1])
            ax.plot(X,Y,'x',color=cstyle[i],label='Ref: '+authors)
            print 'abslen_model.drawExample title',title,'color',cstyle[i]
        

        for i,a200 in enumerate(A200):
            WL,AL = self.createdModels[a200]
            X,Y = numpy.array(WL),numpy.array(AL)
            minAL = min(AL)
            L = 1./minAL/100.
            minWL = WL[AL.index(minAL)]
            s = '({0:.1f}m, {1}nm) '.format(L,minWL)
            ax.plot(X,Y,linestyle=lstyle[i],color=cstyle[i],label=s+'a(200nm)='+str(a200))

        ax.legend(loc='upper center',prop={'size':8})
        ax.grid()
        show = False
        if show:
            plt.show()
        else:
            figpdf = self.figDir + 'example.pdf'
            plt.savefig(figpdf)
            print 'abslen_model.drawExample wrote',figpdf
        return
    def drawIt(self,x,y,xtitle,ytitle,title,figDir=None,ylog=True,xlims=[200.,800.],ylims=[1.e-5,1.e-1]):
        '''
        draw graph defined by x,y

        '''
        plt.clf()
        plt.grid()
        plt.title(title)
        figpdf = 'FIG_'+title.replace(' ','_') + '.pdf'

        X = numpy.array(x)
        Y = numpy.array(y)
        plt.plot(X,Y,'o-')
        plt.xlabel(xtitle)
        plt.ylabel(ytitle)
        if ylog : plt.yscale('log')
        plt.xlim(xlims)
        plt.ylim(ylims)

        if figDir is not None:
            figpdf = figDir + figpdf
            plt.savefig(figpdf)
            print 'abslen_model.drawIt wrote',figpdf
        else:
            plt.show()
        return
        
if __name__ == '__main__' :
    if len(sys.argv)>1:
        if 'help' in sys.argv[1].lower():
            print 'abslen_model.py '
            sys.exit()

    alm = abslen_model()
    #alm.read_all()
    #alm.work()
    #alm.makeModelPoints()
    #alm.createModel(A200=0.)
    alm.drawExample()
    sys.exit('abslen_model over')

