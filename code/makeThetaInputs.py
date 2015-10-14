from ROOT import *
import sys
gROOT.ProcessLine('.L makeThetaInputs_newMistag.C++')

sys.path.append("/home/dfehling/work/compare")

ptBins=[["_pt400-500",400,500],["_pt500-600",500,600],["_pt600-700",600,700],["_pt700-800",700,800],["_pt800-inf",800,99999]]
print ptBins

for ptBin in ptBins:
	print "Working on ptBin",ptBin[0]
	makeThetaInputs(ptBin[0], ptBin[1], ptBin[2])
# makeThetaInputs("test");
