lumi = 19700.
import os
import math
from array import array
import optparse
import ROOT
from ROOT import *
import scipy
from CutOnTree import *

# TW = "(1.0*(1-0.2*0.19361022101362746)*2.71828^(-0.0013*(1-0.2*0.19771016591081025)*0.5*(MCantitoppt+MCtoppt)))" # central value with real N/a inserted
TW = "1.0"

FileName = "/home/dfehling/work/compare/data_0708.root"

tFileName = ["ttjets7_0708", "ttjets10_0708"]
txs = [245*0.074,245*0.014]
tn = [3082812,1249111]
tFilePrefix = "/home/dfehling/work/compare/"
# Vars
var = "jet1mass"
varb = [40,100,500]

# var2 = "((jet1tau3/jet1tau2)*(jet2topTagged))"
var2 = "(jet1tau3/jet1tau2)"
varb2 = [20, 0, 1]

# Preselection Cuts:
# PreSel = "jet1tau2/jet1tau1>0.1&&jet1minMass>50&&jet1nSubj>2"
PreSel = "jet1tau2/jet1tau1>0.1&&jet1minMass>50&&jet1nSubj>2&&jet2topTagged"

# plots:
tPlot = TH2F("tPlotM", "", varb[0],varb[1],varb[2],varb2[0],varb2[1],varb2[2])
dPlot = TH2F("dPlotM", "", varb[0],varb[1],varb[2],varb2[0],varb2[1],varb2[2])

# Fill Plots:

write2dplot(FileName, 1.0, dPlot, var, var2, PreSel, "1.0")
for i in range(len(tFileName)):
	print "xsec * lumi / Nevents = ",txs[i]*lumi/tn[i]
	write2dplot(tFilePrefix+tFileName[i]+'.root', txs[i]*lumi/tn[i], tPlot, var, var2, "("+PreSel+")", TW)

# Now that the plots are filled, we can fit on them.
# First we subtract the non-non-top from the field: ttbar and single top are removed:

dPlot.Add(tPlot,-1.0)
dPlot.Draw("colz")
c1.SaveAs("test.png")

# Start Fitting:
#First arrange the plots in a fittable way: separate by mass bins and compute P/F
# Do this for Muons:
x = []
y = []
exl = []
eyl = []
exh = []
eyh = []
hx = []
hy = []
ehx = []
ehy = []

bins = [[100,120],[120,140],[250,270],[270,400]]
# bins = [[100,120],[120,140],[140,195],[195,250],[250,270],[270,400]]
for b in bins:
	passed = 0
	failed = 0
	for i in range(dPlot.GetNbinsX()): # this is slightly crude, but it works well and doens't have to be run too many times. It loops trhough the bins I define and fills them bin-by-bin (on the histogram) with pass/fail events
		for j in range(dPlot.GetNbinsY()):
			if dPlot.GetXaxis().GetBinCenter(i) < b[1] and dPlot.GetXaxis().GetBinCenter(i) > b[0]:
				if dPlot.GetYaxis().GetBinCenter(j) > 0.55:
					failed = failed + dPlot.GetBinContent(i,j)
				else:
					passed = passed + dPlot.GetBinContent(i,j)
	# If we have low statistics subtracting the ttbar and single top can leave us with "gaps" or just negative events in the least populated bins. This is of course bad, so we just set them to zero. If you want to change the binning it should be such that this is minimized,but this his here ot protect you if need be.
	if passed < 0:
		passed = 0
	if failed < 0:
		failed = 0
	if passed == 0 or failed == 0:
		continue


	x.append((float((b[0]+b[1])/2)-170.))
	exl.append(float((b[1]-b[0])/2))
	exh.append(float((b[1]-b[0])/2))
	y.append(passed/(failed))
	mep = math.sqrt(passed)
	mef = math.sqrt(failed)
	err = (passed/(failed))*math.sqrt((mep/passed)+(mef/(passed))**2)
	# print "\n\nbin = ",b
	# print "x = ",x
	# print "exl = ",exl
	# print "exh = ",exh
	# print "y = ",y
	# print "mep = ",mep
	# print "passed = ", passed
	# print "mef = ",mef
	# print "failed = ", failed
	# print "err = ",err
	eyh.append(err)
	if (passed/failed) - err > 0.:
		eyl.append(err)
	else:
		eyl.append(passed/failed)
G = TGraphAsymmErrors(len(x), scipy.array(x), scipy.array(y), scipy.array(exl), scipy.array(exh), scipy.array(eyl), scipy.array(eyh))

# We now have our three lG Graphs and will not need to use any of the previous files.
#Mu:
funclin = TF1("mfitting_function_linear", "[0]+ [1]*x",-70,230)
funclin.SetParameter(0, 0.5) 
funclin.SetParameter(1, 0.5)
G.Fit(funclin, "EMRN") # This is the fitting step BE SURE TO INCLUDE "E" AS AN OPTION. "E" has the fit take the error into account, thus weighing low statistics bin less than high statistic bins.
# G.Fit(funclin, "E")
# G.Fit(funclin, "EMRNQ")
funclin.SetLineColor(kBlue)

mfitter = TVirtualFitter.GetFitter() # careful if you change the order of things, TVirtualFitter will remember the last fitter used.
mcov = mfitter.GetCovarianceMatrixElement(0,1)

funclinup = TF1("mfitting_function_linear_up", "[0]+ [1]*x + sqrt((x*x*[3]*[3])+(x*2*[4])+([2]*[2]))",-70,230)
funclinup.SetParameter(0, funclin.GetParameter(0))
funclinup.SetParameter(1, funclin.GetParameter(1))
funclinup.SetParameter(2, funclin.GetParErrors()[0])
funclinup.SetParameter(3, funclin.GetParErrors()[1])
funclinup.SetParameter(4, mcov)

funclindn = TF1("mfitting_function_linear_dn", "[0]+ [1]*x - sqrt((x*x*[3]*[3])+(x*2*[4])+([2]*[2]))",-70,230)
funclindn.SetParameter(0, funclin.GetParameter(0))
funclindn.SetParameter(1, funclin.GetParameter(1))
funclindn.SetParameter(2, funclin.GetParErrors()[0])
funclindn.SetParameter(3, funclin.GetParErrors()[1])
funclindn.SetParameter(4, mcov)

funclinup.SetLineColor(kBlue)
funclindn.SetLineColor(kBlue)
funclinup.SetLineStyle(2)
funclindn.SetLineStyle(2)


# Now some beautification:
for G in [G]:
	G.SetMarkerStyle(21)
	G.GetXaxis().SetTitle("#Delta(jet - top)_{mass} [GeV/c^{2}]")
	G.GetYaxis().SetTitle("N_{passed}/N_{failed}")
	G.GetYaxis().SetTitleSize(0.05)
	G.GetYaxis().SetTitleOffset(0.9)
	G.GetXaxis().SetTitleOffset(0.87)
	G.GetXaxis().SetTitleSize(0.04)

leg = TLegend(0.2,0.7,0.45,0.89)
leg.SetFillColor(0)
leg.SetLineColor(0)
leg.AddEntry(funclin, "linear fit", "L")
leg.AddEntry(funclinup, "error in fit", "L")
leg.AddEntry(G, "data", "PL")

# bleg = TLegend(0.2,0.6,0.55,0.89)
# bleg.SetFillColor(0)
# bleg.SetLineColor(0)
# bleg.AddEntry(efunclin, "linear fit (e channel)", "L")
# bleg.AddEntry(funclin, "linear fit (#mu channel)", "L")
# bleg.AddEntry(bfunclin, "linear fir (both channels)", "L")
# bleg.AddEntry(bfunclinup, "error in fit", "L")
# bleg.AddEntry(bG, "data (both channels)", "PL")

C1 = TCanvas("c1", "c1", 1200, 600)
C1.cd()
G.Draw("AP")
# G.GetYaxis().SetRangeUser(-0.,0.5)
funclin.Draw("same")
# bfunclin.Draw("same")
funclinup.Draw("same")
funclindn.Draw("same")
leg.Draw()

# C2 = TCanvas("c2", "c2", 600,600)
# C2.cd()
# bG.Draw("AP")
# bG.GetYaxis().SetRangeUser(-0.,0.6)
# efunclin.Draw("same")
# funclin.Draw("same")
# bfunclin.Draw("same")
# bfunclinup.Draw("same")
# bfunclindn.Draw("same")
# bleg.Draw()

funclin.Print()
funclinup.Print()
funclindn.Print()

print "first = ",funclin.GetParameter(0),";"
print "firstErr = ",funclin.GetParErrors()[0],";"
print "second = ",funclin.GetParameter(1),";"
print "secondErr = ",funclin.GetParErrors()[1],";"

# funclinup.SetParameter(0, funclin.GetParameter(0))
# funclinup.SetParameter(1, funclin.GetParameter(1))
# funclinup.SetParameter(2, funclin.GetParErrors()[0])
# funclinup.SetParameter(3, funclin.GetParErrors()[1])
