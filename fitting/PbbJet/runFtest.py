import ROOT as r,sys,math,os
from ROOT import TFile, TTree, TChain, gPad, gDirectory
from multiprocessing import Process
from optparse import OptionParser
from operator import add
import math
import sys
import time
import array


def exec_me(command, dryRun=False):
    print command
    if not dryRun:
        os.system(command)

        
if __name__ == "__main__":
    parser = OptionParser()
    parser.add_option('-m','--mass'   ,action='store',type='int',dest='mass'   ,default=125, help='mass')
    parser.add_option('--nr1' ,action='store',type='int',dest='NR1'   ,default=1, help='order of rho polynomial for model 1')
    parser.add_option('--np1' ,action='store',type='int',dest='NP1'   ,default=1, help='order of pt polynomial for model 1')
    parser.add_option('--nr2' ,action='store',type='int',dest='NR2'   ,default=2, help='order of rho polynomial for model 2')
    parser.add_option('--np2' ,action='store',type='int',dest='NP2'   ,default=1, help='order of pt polynomial for model 2')
    parser.add_option('--scale',dest='scale', default=1,type='float',help='scale factor to scale MC (assuming only using a fraction of the data)')
    parser.add_option('-l','--lumi'   ,action='store',type='float',dest='lumi'   ,default=36.4, help='lumi')
    parser.add_option('-i','--ifile', dest='ifile', default = 'hist_1DZbb.root',help='file with histogram inputs', metavar='ifile')
    parser.add_option('-t','--toys'   ,action='store',type='int',dest='toys'   ,default=200, help='number of toys')
    
    parser.add_option('-n','--n' ,action='store',type='int',dest='n'   ,default=5*20, help='number of bins')
    parser.add_option('--just-plot', action='store_true', dest='justPlot', default=False, help='just plot')
    parser.add_option('--pseudo', action='store_true', dest='pseudo', default=False, help='run on asimov dataset')
    parser.add_option('--blind', action='store_true', dest='blind', default=False, help='run on blinded dataset')
    parser.add_option('--dry-run',dest="dryRun",default=False,action='store_true',
                  help="Just print out commands to run")    


    (options,args) = parser.parse_args()
    if options.pseudo:
        cardsDir1 = 'cards_mc_r%ip%i'%(options.NR1,options.NP1)
        cardsDir2 = 'cards_mc_r%ip%i'%(options.NR2,options.NP2)
    else:        
        cardsDir1 = 'cards_r%ip%i'%(options.NR1,options.NP1)
        cardsDir2 = 'cards_r%ip%i'%(options.NR2,options.NP2)

    pseudoString = ''
    if options.pseudo:
        pseudoString = '--pseudo'

    blindString = ''
    if options.blind:
        blindString = '--blind'
    
    fillString = ''
    if options.box=='CA15':
        fillString = '--fillCA15'
    if not options.pseudo:
        exec_me('python writeMuonCRDatacard.py -i ./ -o %s/ -c %s --mass %s %s'%(sigDir1, cut, options.mass, fillString),options.dryRun)
        exec_me('python writeMuonCRDatacard.py -i ./ -o %s/ -c %s --mass %s %s'%(sigDir2, cut, options.mass, fillString),options.dryRun)
    
    
    if not options.pseudo:
        exec_me('combineCards.py cat1=%s/card_rhalphabet_cat1.txt cat2=%s/card_rhalphabet_cat2.txt cat3=%s/card_rhalphabet_cat3.txt cat4=%s/card_rhalphabet_cat4.txt cat5=%s/card_rhalphabet_cat5.txt cat6=%s/card_rhalphabet_cat6.txt %s/datacard_muonCR.txt > %s/card_rhalphabet_muonCR_r%ip%i.txt'%(sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,options.NR1, options.NP1),options.dryRun)
        exec_me('combineCards.py cat1=%s/card_rhalphabet_cat1.txt cat2=%s/card_rhalphabet_cat2.txt cat3=%s/card_rhalphabet_cat3.txt cat4=%s/card_rhalphabet_cat4.txt cat5=%s/card_rhalphabet_cat5.txt cat6=%s/card_rhalphabet_cat6.txt %s/datacard_muonCR.txt > %s/card_rhalphabet_muonCR_r%ip%i.txt'%(sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,options.NR2, options.NP2),options.dryRun)
        #exec_me("text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125  --PO verbose --PO 'map=.*/*hqq125:r[1,0,20]' --PO 'map=.*/zqq:r_z[1,0,20]' %s/card_rhalphabet_muonCR_r%ip%i.txt -o %s/card_rhalphabet_muonCR_floatZ_r%ip%i.root"%(sigDir1, options.NR1, options.NP1, sigDir1, options.NR1, options.NP1))
        #exec_me("text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125  --PO verbose --PO 'map=.*/*hqq125:r[1,0,20]' --PO 'map=.*/zqq:r_z[1,0,20]' %s/card_rhalphabet_muonCR_r%ip%i.txt -o %s/card_rhalphabet_muonCR_floatZ_r%ip%i.root"%(sigDir2, options.NR2, options.NP2, sigDir2, options.NR2, options.NP2))
    else:
        exec_me('combineCards.py cat1=%s/card_rhalphabet_cat1.txt cat2=%s/card_rhalphabet_cat2.txt cat3=%s/card_rhalphabet_cat3.txt cat4=%s/card_rhalphabet_cat4.txt cat5=%s/card_rhalphabet_cat5.txt cat6=%s/card_rhalphabet_cat6.txt > %s/card_rhalphabet_r%ip%i.txt'%(sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,sigDir1,options.NR1, options.NP1),options.dryRun)
        exec_me('combineCards.py cat1=%s/card_rhalphabet_cat1.txt cat2=%s/card_rhalphabet_cat2.txt cat3=%s/card_rhalphabet_cat3.txt cat4=%s/card_rhalphabet_cat4.txt cat5=%s/card_rhalphabet_cat5.txt cat6=%s/card_rhalphabet_cat6.txt > %s/card_rhalphabet_r%ip%i.txt'%(sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,sigDir2,options.NR2, options.NP2),options.dryRun)
        #exec_me("text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125  --PO verbose --PO 'map=.*/*hqq125:r[1,0,20]' --PO 'map=.*/zqq:r_z[1,0,20]' %s/card_rhalphabet_r%ip%i.txt -o %s/card_rhalphabet_r%ip%i.txt"%(sigDir1, options.NR1, options.NP1, sigDir1, options.NR1, options.NP1))
        #exec_me("text2workspace.py -P HiggsAnalysis.CombinedLimit.PhysicsModel:multiSignalModel -m 125  --PO verbose --PO 'map=.*/*hqq125:r[1,0,20]' --PO 'map=.*/zqq:r_z[1,0,20]' %s/card_rhalphabet_r%ip%i.txt -o %s/card_rhalphabet_r%ip%i.txt"%(sigDir2, options.NR2, options.NP2, sigDir2, options.NR2, options.NP2))

    p1 = int((options.NR1+1)*(options.NP1+1)) # free paramaters
    p2 = int((options.NR2+1)*(options.NP2+1)) # free parameters
    

    dataString = ''
    if not options.pseudo:
        dataString = '--data'

    if options.pseudo:
        exec_me('python limit.py -M FTest --datacard %s/card_rhalphabet_r%ip%i.txt --datacard-alt %s/card_rhalphabet_r%ip%i.txt -o %s -n %i --p1 %i --p2 %i -t %i --lumi %f %s -r %f --freezeNuisances tqqeffSF,tqqnormSF'%(sigDir1, options.NR1, options.NP1, sigDir2, options.NR2, options.NP2, ftestDir, options.n, p1, p2, options.toys, options.lumi, dataString, options.r),options.dryRun)
    else:
        exec_me('python limit.py -M FTest --datacard %s/card_rhalphabet_muonCR_r%ip%i.txt --datacard-alt %s/card_rhalphabet_muonCR_r%ip%i.txt -o %s -n %i --p1 %i --p2 %i -t %i --lumi %f %s -r %f'%(sigDir1, options.NR1, options.NP1, sigDir2, options.NR2, options.NP2, ftestDir_muonCR, options.n, p1, p2, options.toys, options.lumi, dataString, options.r),options.dryRun)
     

     
