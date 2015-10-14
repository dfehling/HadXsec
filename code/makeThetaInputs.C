#include <stdlib.h>
#include <sys/stat.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <TTree.h>
#include <TH1F.h>
#include <TCanvas.h>
#include <TRandom3.h>
#include <TFile.h>
#include <THStack.h>
#include <TLegend.h>
#include <TLatex.h>
#include <TF1.h>

std::vector<std::string> split(const std::string &s, char delim);

void runOverTree(TTree * tree, std::vector<std::vector<double> > nSubCuts, std::vector<std::vector<TH1F*> > mistag, 
                               std::vector<std::vector<std::vector<std::vector<TH1F*> > > > histos, int sample,
                               double ptLow, double ptHigh, double systCoefficient=0.0);

void scaleHistos(std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos, double systCoefficient=0.0);
void combineHistos(std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos);
void plotHistos(TH1F *dataH, TH1F *ttbarH, TH1F *qcdH, string tags_label, string nsub_label, string output, string syst_label="",
                const TString& ptLabel="");


void makeThetaInputs(const TString& ptLabel="", double ptLow=0, double ptHigh=99999) {
// string var=" ", int tags=3, int nsub=4, string input_label = "0708_mistag", string label="", 
    // string add_selection ="1", string minMassLabel="lowMPMnSubj", int nBins=0, float min=0, float max=0) {

    TString test = "test";
    test + ptLabel;
    std::vector<double> ptCuts = {ptLow,ptHigh};

    string var = "taggedMass";
    string folder_png, folder_root, file_out;

    string prefix_path = "/home/dfehling/work/compare/";
    // // system("mkdir -p '$var'");
    string full_path = prefix_path + var;
    // std::cout << full_path << std::endl;
    mkdir(full_path.c_str(),0777);
    folder_png = prefix_path + var + "/PNG/";
    folder_root = prefix_path + var + "/ROOT/";
    mkdir(folder_png.c_str(),0777);
    mkdir(folder_root.c_str(),0777);

    // gROOT->Reset();
    // gStyle->SetOptStat(0);
    // gStyle->SetOptTitle(0);

    std::vector<string> btagLabels = {"0", "1", "2", "all"};
    // std::vector<string> nsubLabels = {"00-04", "04-07", "07-10", "00-07", "00-10"};
    std::vector<string> nsubLabels = {"000-055"};
    // std::vector<string> nsubLabels = {"000-040"};//,"00-10"};
    std::vector<string>  etaLabels = {"", "_10eta24"};

    std::vector<string> systLabel = {""};
    std::vector<double> systCoef  = {0.0};
    // std::vector<string> systLabel = {"", "__scale__plus", "__scale__minus"};
    // std::vector<double> systCoef  = {0.0, 1.0, -1.0};

    string sample_label[3][2] = { {"__DATA",    "__QCD"},
                                  {"__TTBAR",   "__ttbar7mis"},
                                  {"__ttbar10", "__ttbar10mis"} };


    // const int test = nsubLabels.size();

    std::vector<int>    btagCuts;
    // std::vector<string> nsubCuts_s[5];
    std::vector<std::vector<string> > nsubCuts_s;
    // std::vector<double> nsubCuts[5];
    std::vector<std::vector<double> > nsubCuts;

    for(unsigned int i=0;i<btagLabels.size();i++) {
        if (btagLabels[i] == "all") {
            btagCuts.push_back(-1);
            continue;
        }
        else{
            btagCuts.push_back(std::strtol((btagLabels[i]).c_str(),NULL,0));
        }
    }

    for(unsigned int i=0;i<nsubLabels.size();i++) {
        // std::cout << i << std::endl;
        nsubCuts_s.push_back(split(nsubLabels[i],'-'));
        // std::cout << nsubCuts_s[i][0] << std::endl;
        // nsubCuts[i].push_back(std::strtod((nsubCuts_s[i][0]).c_str(),NULL)/10.0);
        // nsubCuts[i].push_back(std::strtod((nsubCuts_s[i][1]).c_str(),NULL)/10.0);
        std::vector<double> tempVector;

        tempVector.push_back(std::strtod((nsubCuts_s[i][0]).c_str(),NULL)/100.0);
        tempVector.push_back(std::strtod((nsubCuts_s[i][1]).c_str(),NULL)/100.0);
        nsubCuts.push_back(tempVector);
    }

    int nBins=10;
    int min=140;
    int max=250;

    string input_label = "0708";

    TFile data(Form("%sdata_%s.root", prefix_path.c_str() ,input_label.c_str()));
    TFile ttbar7(Form("%sttjets7_%s.root", prefix_path.c_str() ,input_label.c_str()));
    TFile ttbar10(Form("%sttjets10_%s.root", prefix_path.c_str() ,input_label.c_str()));

    TString outname = Form("%s%s",folder_root.c_str(), var.c_str()) + ptLabel + Form("_%s.root", nsubLabels[0].c_str());
    TFile *outfile = new TFile(outname, "recreate");
    // TFile *outfile = new TFile(outname, "update");
    // TFile *outfile = new TFile(Form("%s%s%s_%s.root",folder_root.c_str(), var.c_str(), nsubLabels[0].c_str()),"recreate");
    outfile->cd();
    // TFile mistagFile("mistag_MPM_0915.root");
    // TFile mistagFile("mistag_MPM_noNsubjReq_0916.root");
    // TFile mistagFile("mistag_invertedtau32.root");
    // TFile mistagFile("mistag_MPM_noNsubjReq_divideTotal_0917.root");
    // TFile mistagFile("mistag_MPM_divideTotal_0917.root");

    TTree * trees[3] = { (TTree*)data.Get("treeVars"),
                         (TTree*)ttbar7.Get("treeVars"),
                         (TTree*)ttbar10.Get("treeVars") };

    // We will have an array of histograms set up like this: histos[btagBin][nSubjettinessBin][SelectionBin][SampleBin]
    // TH1F *histos[btagLabels.size()][nsubLabels.size()][3][2];
    std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos[3];
    // std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos_up;
    // std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos_down;
    // TH1F *mistagHistos[btagLabels.size()][nsubLabels.size()];

    std::vector <std::vector <TH1F*> > mistagHistos;

    for(unsigned int sys=0; sys<systLabel.size(); sys++) {
        for(unsigned int sample=0;sample<3;sample++) {
            std::vector<std::vector<std::vector<std::vector<TH1F*> > > > tempHisto4;
            for(unsigned int btags=0;btags<btagLabels.size();btags++) {
                std::vector<std::vector<std::vector<TH1F*> > > tempHisto3;
                // std::vector<TH1F*> tempMistagHisto;
                for(unsigned int nsub=0;nsub<nsubLabels.size();nsub++) {
                    std::vector<std::vector<TH1F*> > tempHisto2;
                    // TString mistagname = Form("MR_%sbtag_%stau", btagLabels[btags].c_str(), nsubLabels[nsub].c_str());
                    // if(sample==0) {
                    //     tempMistagHisto.push_back((TH1F*)mistagFile.Get(mistagname)->Clone());
                    //     std::cout << mistagname << std::endl;
                    // }
                    // std::cout << mistagname << std::endl;
                    // if (sample==0) mistagHistos[btags][nsub] = (TH1F*)mistagFile.Get(mistagname)->Clone();
                    for(unsigned int eta=0;eta<2;eta++) {
                        std::vector<TH1F*> tempHisto;
                        for(unsigned int sel=0;sel<2;sel++) {
                            TString histname = Form("%s_%sbtag_%stau%s%s%s", var.c_str(),btagLabels[btags].c_str(), nsubLabels[nsub].c_str(), etaLabels[eta].c_str(), sample_label[sample][sel].c_str(),systLabel[sys].c_str());
                            TH1F * tempHisto_ = new TH1F(histname,histname,nBins,min,max);
                            // std::cout << histname << std::endl;
                            tempHisto_->Sumw2();
                            tempHisto.push_back(tempHisto_);
                        }
                        tempHisto2.push_back(tempHisto);
                    }
                    tempHisto3.push_back(tempHisto2);
                }
                tempHisto4.push_back(tempHisto3);
                // if(sample==0) mistagHistos.push_back(tempMistagHisto);
            }
            histos[sys].push_back(tempHisto4);
            // histos_up.push_back(tempHisto4);
            // histos_down.push_back(tempHisto4);   
        }
    }

    for(unsigned int sys=0; sys<systLabel.size();sys++) {
        for(unsigned int sample=0;sample<3;sample++) {
            runOverTree(trees[sample],nsubCuts,mistagHistos,histos[sys][sample],sample,ptLow,ptHigh,systCoef[sys]);
        }
    scaleHistos(histos[sys],systCoef[sys]);
    combineHistos(histos[sys]);
    }


    // outfile->cd();
    // for(unsigned int sample=0;sample<2;sample++) {

    data.Close();
    ttbar7.Close();
    ttbar10.Close();
    
    for(unsigned int sys=0; sys<systLabel.size();sys++) {
        for(unsigned int btags=0;btags<btagLabels.size();btags++) {
        // for(unsigned int btags=1;btags<3;btags++) {
            for(unsigned int nsub=0;nsub<nsubLabels.size();nsub++) {
                TH1F * data  = histos[sys][0][btags][nsub][0][0];
                TH1F * ttbar = histos[sys][1][btags][nsub][0][0];
                TH1F * qcd   = histos[sys][0][btags][nsub][0][1];
                // histos[0][btags][nsub][0][0]->Write();
                // histos[0][btags][nsub][0][1]->Write();
                // histos[1][btags][nsub][0][0]->Write();
                plotHistos(data,ttbar,qcd,btagLabels[btags],nsubLabels[nsub],folder_png,systLabel[sys],ptLabel);
                // outfile->cd();
                if(sys==0) data->Write();
                ttbar->Write();
                qcd->Write();
                // histos[0][btags][nsub][1][0]->Write();
                // histos[0][btags][nsub][1][1]->Write();
                // histos[1][btags][nsub][1][0]->Write();
            }
        }
    }


    // for(unsigned int sample=0;sample<3;sample++) {
    //     for(unsigned int btags=0;btags<btagLabels.size();btags++) {
    //         for(unsigned int nsub=0;nsub<nsubLabels.size();nsub++) {
    //             for(unsigned int eta=0;eta<2;eta++) {
    //                 for(unsigned int sel=0;sel<2;sel++) {
    //                     histos[sample][btags][nsub][eta][sel]->Write();
    //                 }
    //             }
    //         }
    //     }
    // }

    outfile->Close();
}

void runOverTree(TTree * tree, std::vector<std::vector<double> > nSubCuts, std::vector<std::vector<TH1F*> > mistag, 
                               std::vector<std::vector<std::vector<std::vector<TH1F*> > > > histos, int sample, 
                               double ptLow, double ptHigh, double systCoefficient) {

    int total=0;
    float totalWt=0;

    int nentries = (int) tree->GetEntries();

    float taggedPt, taggedEta, taggedMass, taggedTau32, jet1pt, jet2pt, jet1eta, jet2eta, jet1mass, jet2mass, jet1tau32, jet2tau32, x;
    float mistagWt[2] = {1.0, -1.0};
    int index, mistagBin, jet1bTagged, jet2bTagged, jet1topTagged, jet2topTagged;
    float jet1minMass;
    int jet1nSubj;

    TRandom3 *random_mt = new TRandom3(170); // need to make this 0 when actually running

    std::cout << nentries << std::endl;

    tree->SetBranchAddress("index", &index);
    tree->SetBranchAddress("jet1pt", &jet1pt);
    tree->SetBranchAddress("jet2pt", &jet2pt);
    tree->SetBranchAddress("jet1eta", &jet1eta);
    tree->SetBranchAddress("jet2eta", &jet2eta);
    tree->SetBranchAddress("jet1mass", &jet1mass);
    tree->SetBranchAddress("jet2mass", &jet2mass);
    tree->SetBranchAddress("jet1topTagged", &jet1topTagged);
    tree->SetBranchAddress("jet2topTagged", &jet2topTagged);
    tree->SetBranchAddress("jet1bTagged", &jet1bTagged);                // in the latest ntuples these have changed to _sub
    tree->SetBranchAddress("jet2bTagged", &jet2bTagged);
    tree->SetBranchAddress("jet1tau32", &jet1tau32);
    tree->SetBranchAddress("jet2tau32", &jet2tau32);
    tree->SetBranchAddress("jet1minMass", &jet1minMass);
    tree->SetBranchAddress("jet1nSubj", &jet1nSubj);

    
    std::cout << ptLow << " " << ptHigh << std::endl;
    for(int i=0; i<nentries; i++) {
        tree->GetEntry(i);

        if(jet1pt<ptLow || jet1pt >= ptHigh) continue;

        //This is to ensure that we pick randomly and are not biased. We also are only taking signal and background with at least 1 toptag at this point
        float x = random_mt->Rndm();
        taggedMass=0;
        if(x < 0.5) {  // SHOULD I RANDOMLY PICK FOR SIGNAL AND BACKGROUND SEPARATELY?
            taggedPt = jet1pt;
            taggedEta = jet1eta;
            taggedTau32 = jet1tau32;                                // maybe use number of tau tags?
            // if(jet1topTagged) taggedMass = jet1mass;
            taggedMass = jet1mass;
            // else continue;
        }
        else {
            taggedPt = jet2pt;
            taggedEta = jet2eta;
            taggedTau32 = jet2tau32;
            // if(jet2topTagged) taggedMass = jet2mass;
            taggedMass = jet2mass;
            // else continue;
        }
        // if(taggedMass==0) continue;

        int bTags=0;
        if(jet1bTagged) bTags++;
        if(jet2bTagged) bTags++;

        int tTags=0;
        if(jet1topTagged) tTags++;
        if(jet2topTagged) tTags++;

        int etaBin=-1;
        if(abs(taggedEta)<1.0) etaBin=0;
        else if(abs(taggedEta)>=1.0&&abs(taggedEta)<2.4) etaBin=1;
        else continue;

        bool passPre = jet1minMass>50.0 && jet1nSubj>2;
        // bool passOther = jet1minMass>50.0 && jet1nSubj>2 && jet2topTagged;

        // jet 2 toptagged and jet1tau32 > 0.55
          //  1  p0           1.88637e-01   1.43763e-02   1.70921e-05   7.90597e-05
          //  2  p1           4.96022e-04   1.87321e-04   1.87321e-04  -3.37276e-01
        // mistagWt[1]=.188637+(.000496022*(jet1mass-170.));

        // jet1tau32 > 0.55
           // 1  p0           1.27072e-01   1.07958e-02   1.29764e-05   7.56836e-05
           // 2  p1           3.17809e-04   1.33020e-04   1.33020e-04  -2.71077e-01

        // mistagWt[1]=.127072+(.000317809*(jet1mass-170.));

        // jet2 toptagged in preselection and jet1tau32 > 0.55
       // 1  p0           1.22354e-01   2.71476e-02  -4.11305e-04   2.18165e-04
       // 2  p1           2.93623e-04   2.82489e-04   2.82489e-04  -5.99508e-02

        double first, firstErr, second, secondErr;
        // first = .122354;
        // firstErr = .0271476;
        // second = .000293623;
        // secondErr = .000282489;

        first =  0.122353729882 ;
        firstErr =  0.027147579072 ;
        second =  0.000293623479694 ;
        secondErr =  0.000282489486328 ;


           // 1  p0           1.08503e-02   7.80048e-03   2.23706e-04   3.62545e-03
           // 2  p1           6.01435e-05   1.08938e-04   1.08938e-04  -1.39310e+00



        // first     = .0108503;
        // firstErr  = .00780048;
        // second    = .0000601435;
        // secondErr = .000108938;


        first+= systCoefficient*firstErr;
        second+=systCoefficient*secondErr;

        // mistagWt[1]=.122354+(.000293623*(jet1mass-170.));
        mistagWt[1]=first+(second*(jet1mass-170.));
        // std::cout << mistagWt[1] << std::endl;

        for(unsigned int nSub=0;nSub<nSubCuts.size();nSub++) {

            float allBmistagWt = mistagWt[1];
 
            // bool passTaggedTau = taggedTau32>nSubCuts[nSub][0]&&taggedTau32<=nSubCuts[nSub][1] && passOther;
            bool passTaggedTau = jet1tau32>nSubCuts[nSub][0]&&jet1tau32<=nSubCuts[nSub][1];
            bool passBothTau = (jet1tau32>nSubCuts[nSub][0]&&jet1tau32<=nSubCuts[nSub][1])&&(jet2tau32>nSubCuts[nSub][0]&&jet2tau32<=nSubCuts[nSub][1]);

            if(passPre && jet2topTagged) {
                //Don't weight as this is our signal region
                if(passTaggedTau) {
                    histos[bTags][nSub][etaBin][0]->Fill(taggedMass,mistagWt[0]);
                    histos[3][nSub][etaBin][0]->Fill(taggedMass,mistagWt[0]);
                }
                //Now check our Tau32 cut to fill the qcd
                // if(!passTaggedTau) {
                // else if(jet2topTagged){
                else {
                    histos[bTags][nSub][etaBin][1]->Fill(taggedMass,mistagWt[1]);
                    histos[3][nSub][etaBin][1]->Fill(taggedMass,mistagWt[1]);
                }
            }
        }
    }
}

void scaleHistos(std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos, double systCoefficient) {

    // sample | btag | nsub | eta | sel

    float topTagScale[2] = {1.173, 0.705};
    float nsubScale[2]   = {0.990, 0.845};
    float btagScale[2]   = {0.915, 1.52};
    float topTagScaleErr[2] = {0.085, 0.138};
    float nsubScaleErr[2]   = {0.079, 0.181};
    float btagScaleErr[2]   = {0.074, 0.232};

    // for(int i=0; i<2;i++) {
    //     topTagScale[i]+=systCoefficient*topTagScaleErr[i];
    //     nsubScale[i]+=systCoefficient*nsubScaleErr[i];
    //     btagScale[i]+=systCoefficient*btagScaleErr[i];
    // }

    float topXsec = 245;
    float lumi    = 19700;

    float sampleScale[3] = {1, 0.074/3082812., 0.014/1249111.};
    // float sampleScale[3] = {1, 1, 1};

    float scales[3];
    for(unsigned int sample=1;sample<histos.size();sample++) {
        for(unsigned int btags=0;btags<histos[0].size()-1;btags++) {
            for(unsigned int nsub=0;nsub<histos[0][0].size();nsub++) {
                for(unsigned int eta=0;eta<histos[0][0][0].size();eta++) {
                    for(unsigned int sel=0;sel<histos[0][0][0][0].size();sel++) {
                    scales[0] = sampleScale[sample]*topTagScale[eta]*pow(btagScale[eta],btags)*nsubScale[eta];
                    scales[1] = sampleScale[sample]*lumi*topXsec*topTagScale[eta]*pow(btagScale[eta],btags)*nsubScale[eta];
                    scales[2] = sampleScale[sample]*lumi*topXsec*topTagScale[eta]*nsubScale[eta];
                        if(sample==0&&sel==0) continue;     // Don't scale Data
                        // float scale = sampleScale[sample]*lumi*topXsec*topTagScale[eta]*pow(btagScale[eta],btags)*nsubScale[eta];
                        // float allBscale = sampleScale[sample]*lumi*topXsec*topTagScale[eta]*nsubScale[eta];
                        // float qcdscale = sampleScale[sampe]*topTagScale[eta]*pow(btagScale[eta],btags)*nsubScale[eta];
                        // histos[sample][btags][nsub][eta][0]->Scale(scales[1]);
                        // histos[sample][btags][nsub][eta][1]->Scale(scales[2]);
                        histos[sample][btags][nsub][eta][sel]->Scale(scales[sel+1]);
                        if(btags==0&&sel==0) {
                            histos[sample][3][nsub][eta][0]->Scale(scales[2]);
                            histos[sample][3][nsub][eta][1]->Scale(scales[2]);
                        }
                    }
                }
            }
        }
    }
}

void combineHistos(std::vector<std::vector<std::vector<std::vector<std::vector<TH1F*> > > > > histos) {
    //Combine the eta-separated histograms
    for(unsigned int sample=0;sample<histos.size();sample++) {
        for(unsigned int btags=0;btags<histos[0].size();btags++) {
            for(unsigned int nsub=0;nsub<histos[0][0].size();nsub++) {
                for(unsigned int sel=0;sel<histos[0][0][0][0].size();sel++) {
                    histos[sample][btags][nsub][0][sel]->Add(histos[sample][btags][nsub][1][sel]);
                }
            }
        }
    }
    //Combine the ttbar
    for(unsigned int btags=0;btags<histos[0].size();btags++) {
        for(unsigned int nsub=0;nsub<histos[0][0].size();nsub++) {
            //Subtract mistagged ttbar contribution from ttbar MC
            histos[1][btags][nsub][0][0]->Add(histos[1][btags][nsub][0][1],-1.0);
            histos[2][btags][nsub][0][0]->Add(histos[2][btags][nsub][0][1],-1.0);
            //Add this to the QCD estimate
            histos[0][btags][nsub][0][1]->Add(histos[1][btags][nsub][0][1]);
            histos[0][btags][nsub][0][1]->Add(histos[2][btags][nsub][0][1]);
            //Combine the ttbar MC
            histos[1][btags][nsub][0][0]->Add(histos[2][btags][nsub][0][0]);
            // std::cout << histos[1][btags][nsub][0][0]->Integral() << std::endl;
        }
    }
}

void plotHistos(TH1F *dataH, TH1F *ttbarH, TH1F *qcdH, string tags_label, string nsub_label, string output, string syst_label, const TString& ptLabel) {
    dataH->SetLineColor(kBlack);
    dataH->SetMarkerStyle(20);
    ttbarH->SetLineColor(kRed);
    ttbarH->SetFillColor(kRed);
    qcdH->SetLineColor(kYellow);
    qcdH->SetFillColor(kYellow);

    std::cout << "plotting " << syst_label << std::endl;

    THStack *stack = new THStack("stack", "stack");
    stack->Add(qcdH);
    stack->Add(ttbarH);

    cout << "ttbar    : " << ttbarH->Integral() << endl;
    cout << "QCD      : " << qcdH->Integral() << endl;
    cout << "DATA     : " << dataH->Integral() << endl;


    TCanvas *c1 = new TCanvas("c1", "c1",600,500);
    // c1->Range(0,0,1,1);
    // c1->Draw();
    // c1->cd();

    TPad *c1_1 = new TPad("c1_1", "newpad",0.01,0.01,0.99,0.32);
    c1_1->Draw();
    TPad *c1_2 = new TPad("c1_2", "newpad",0.01,0.33,0.99,0.99);
    c1_2->Draw(); 

    c1->cd();
    c1_2->cd();
    c1_2->SetTopMargin(0.1);
    c1_2->SetBottomMargin(0.01);
    c1_2->SetRightMargin(0.05);
    c1_2->SetLeftMargin(0.1);
    c1_2->SetFillStyle(0);
    // c1_2->SetLogy();

    double maximum = max(dataH->GetMaximum(),(ttbarH->GetMaximum()+qcdH->GetMaximum()));
    // std::cout << maximum << std::endl;

    dataH->SetMaximum( 1.2* maximum );
    dataH->SetMinimum( 0.0 );
    stack->SetMaximum( 1.2* maximum );

    dataH->GetYaxis()->SetTitle("Events");
    dataH->GetYaxis()->SetTitleOffset(0.7);

    // ROOT::SetOwnership(stack,kFALSE);

    // dataH->Draw("axis");
    // stack->SetTitle(Form("%sbtags_%stau",tags_label.c_str(), nsub_label.c_str()));
    stack->SetTitle("");
    stack->Draw("hist");
    dataH->Draw("same");
    // dataH->Draw("same");
    // stack->Draw("hist same");


    TLegend *leg = new TLegend(0.75,0.6,0.9,0.88);
    // leg->AddEntry(dataH, Form("%sbtags_%stau",tags_label.c_str(), nsub_label.c_str()), "");
    leg->AddEntry(dataH, "Data", "p");
    leg->AddEntry(qcdH, "QCD", "f");
    leg->AddEntry(ttbarH, "t#bar{t}", "f");
    leg->SetFillColor(0);
    leg->Draw("same");

    TLatex *cmsLabel = new TLatex();
    cmsLabel->SetNDC();
    cmsLabel->DrawLatex(0.1,0.9, "CMS Preliminary, #sqrt{s} = 8 TeV, 19.7 fb^{-1}");
    // cmsLabel->DrawLatex(0.1,0.9, "CMS Preliminary, #sqrt{s} = 8 TeV, 7.051 fb^{-1}");
    cmsLabel->DrawLatex(0.7, 0.9, Form("%sbtags  tau #in%s",tags_label.c_str(), nsub_label.c_str()));

    // gPad->RedrawAxis();

    TH1F *totalH = new TH1F();
    totalH = (TH1F *) qcdH->Clone("total");
    totalH->Add(ttbarH);
    // totalH->Add(ttbar7H);
    // // cout << totalH->Integral() << endl;
    // // cout << dataH->Integral() << endl;

    TH1F *ratioH = new TH1F();
    ratioH->Sumw2();
    ratioH = (TH1F*) dataH->Clone("ratio");
    ratioH->Divide(totalH);

    c1_1->cd();
    c1_1->SetTopMargin(0.01);
    c1_1->SetBottomMargin(0.3);
    c1_1->SetRightMargin(0.05);
    c1_1->SetLeftMargin(0.1);
    c1_1->SetFillStyle(0);

    ratioH->GetYaxis()->SetRangeUser(0.,2.);
    ratioH->GetYaxis()->SetTitle("Data / BG Ratio");
    ratioH->GetYaxis()->SetTitleOffset(0.4);
    ratioH->GetYaxis()->SetTitleSize(0.11);
    ratioH->GetXaxis()->SetLabelSize(0.11);
    ratioH->GetXaxis()->SetTitleSize(0.11);
    ratioH->GetXaxis()->SetTitle( "taggedMass");
    ratioH->Draw("E");

    TF1 *line = new TF1("line", "1", 140, 250);
    line->SetLineColor(kBlack);

    c1->Draw();
    line->Draw("same");

    TString label = Form("%staggedMass", output.c_str()) + ptLabel + Form("_%sbtag_%stau%s.png",tags_label.c_str(),nsub_label.c_str(),syst_label.c_str());

    std::cout << label << std::endl;
    // c1->SaveAs("./test.png");
    c1->Print(label);


    // delete stack;
    delete c1;



}


std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}
