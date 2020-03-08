#include "TCutG.h"
#include "TCut.h"
#include "TMath.h"

 Double_t xray_peak(Double_t *x,Double_t *par)
 {
   return par[3]+par[4]*x[0]+par[2]*TMath::Exp(-((x[0]-par[0])*(x[0]-par[0]))/(2*par[1]*par[1]));
 } 

Double_t peak_ran(Double_t *x,Double_t *par)
{
  return par[2]*TMath::Exp(-((x[0]-par[0])*(x[0]-par[0]))/(2*par[1]*par[1]));
} 

Double_t linear(Double_t *x,Double_t *par)
{
  return par[0]+par[1]*x[0];
} 






void fitting_prog()
{  
  printf("OPTQ 2 Calibration\n\n");
  
  
  char fname[100];




  Double_t par[9];
  Double_t real_energy_of_xray[3]={0};

    
     //Optimization
     sprintf(fname,"for_fit/for_calib_190627_max_he_optq.root");
     TFile *fin=new TFile(fname);
     TH1F *clbr_max_he=(TH1F*)fin->Get("clbr_max_he");
     TCanvas *c2 = new TCanvas("c2");
     gPad->SetLogy(1);
     clbr_max_he->SetLineWidth(2);
     clbr_max_he->Rebin(1);
     int first_peak_lb = 3900000;
     int first_peak_ub = 4200000;
     int fourth_lb = 4700000;
     int fourth_ub = 4900000;
     int six_lb = 5100000;
     int six_ub = 5300000;
     real_energy_of_xray[0]=302;
     real_energy_of_xray[1]=356;
     real_energy_of_xray[2]=383;
    

    /*
     //Before Max Peak
     sprintf(fname,"for_fit/for_calib_190627_max_he_before_pq.root");
     TFile *fin=new TFile(fname);
     TH1F *clbr_max_he=(TH1F*)fin->Get("clbr_max_he");
     TCanvas *c2 = new TCanvas("c2");
     gPad->SetLogy(1);
     clbr_max_he->SetLineWidth(2);
     clbr_max_he->Rebin(1);
     int first_peak_lb = 1200000;
     int first_peak_ub = 1400000;
     int fourth_lb = 1400000;
     int fourth_ub = 1600000;
     int six_lb = 1700000;
     int six_ub = 2100000;
     real_energy_of_xray[0]=302;
     real_energy_of_xray[1]=356;
     real_energy_of_xray[2]=383;
    */


    /*
     //After the peak
     sprintf(fname,"for_fit/for_calib_190627_max_he_after_pq.root");
     TFile *fin=new TFile(fname);
     TH1F *clbr_max_he=(TH1F*)fin->Get("clbr_max_he");
     TCanvas *c2 = new TCanvas("c2");
     gPad->SetLogy(1);
     clbr_max_he->SetLineWidth(2);
     clbr_max_he->Rebin(1);
     int first_peak_lb = 1800000;
     int first_peak_ub = 2000000;
     int fourth_lb = 2200000;
     int fourth_ub = 2500000;
     int six_lb = 2500000;
     int six_ub = 2700000;
     real_energy_of_xray[0]=302;
     real_energy_of_xray[1]=356;
     real_energy_of_xray[2]=383;
    */

  /*
  //Max_Peak
     sprintf(fname,"for_fit/for_calib_190627_max_he.root");
     TFile *fin=new TFile(fname);
     TH1F *clbr_max_he=(TH1F*)fin->Get("clbr_max_he");
    TCanvas *c2 = new TCanvas("c2");
    gPad->SetLogy(1);
    clbr_max_he->SetLineWidth(2);
    clbr_max_he->Rebin(1);
  int first_peak_lb = 1280;
  int first_peak_ub = 1440;
  int fourth_lb = 6340;
  int fourth_ub = 6740;
  int six_lb = 6920;
  int six_ub = 7260;
  real_energy_of_xray[0]=81;
  real_energy_of_xray[1]=356;
  real_energy_of_xray[2]=383;
   */

    /*
  //Q_factor (Total integral)
    sprintf(fname,"for_fit/for_calib_190627_max_he_q.root");
    TFile *fin=new TFile(fname);
    TH1F *clbr_max_he=(TH1F*)fin->Get("clbr_max_he");
    TCanvas *c2 = new TCanvas("c2");
    gPad->SetLogy(1);
    clbr_max_he->SetLineWidth(2);
    clbr_max_he->Rebin(1);
    int first_peak_lb = 3300000;
    int first_peak_ub = 3500000;
    int fourth_lb = 4100000;
    int fourth_ub = 4399000;
    int six_lb = 4500000;
    int six_ub = 4700000;
    real_energy_of_xray[0]=302;
    real_energy_of_xray[1]=356;
    real_energy_of_xray[2]=383;
    */
    
  TF1 *Gaus1    =new TF1("gaus1","gaus",first_peak_lb,first_peak_ub);
  clbr_max_he->Fit(Gaus1,"R+");
    Gaus1->GetParameters(&par[0]);
  cout << "Constant: " << par[0] << endl;
  cout << "Mean: " << par[1] << endl;
  cout << "Error: " << par[2] << endl;
  TF1 *fitfun01=new TF1("fitfun01",xray_peak,first_peak_lb,first_peak_ub,5);
  fitfun01->SetLineColor(2);
  fitfun01->SetParameter(0,par[1]);
  fitfun01->SetParameter(1,par[2]);
  fitfun01->SetParameter(2,par[0]);
  fitfun01->SetParLimits(3,-1,1);
  fitfun01->SetParLimits(4,-1,1);  
  clbr_max_he->Fit("fitfun01","","",first_peak_lb,first_peak_ub);  


    TF1 *Gaus2    =new TF1("gaus2","gaus",fourth_lb,fourth_ub);
    clbr_max_he->Fit(Gaus2,"R+");
    Gaus2->GetParameters(&par[3]);
    cout << "Constant: " << par[3] << endl;
    cout << "Mean: " << par[4] << endl;
    cout << "Error: " << par[5] << endl;
  TF1 *fitfun02=new TF1("fitfun02",xray_peak,fourth_lb,fourth_ub,5); 
  fitfun02->SetLineColor(2);
  fitfun02->SetParameter(0,par[4]);
  fitfun02->SetParameter(1,par[5]);
  fitfun02->SetParameter(2,par[3]);
  fitfun02->SetParLimits(3,-1,1);
  fitfun02->SetParLimits(4,-1,1);  
  clbr_max_he->Fit("fitfun02","","",fourth_lb,fourth_ub);  


   TF1 *Gaus3    =new TF1("gaus3","gaus",six_lb,six_ub);
    clbr_max_he->Fit(Gaus3,"R+");
    Gaus3->GetParameters(&par[6]);
    cout << "Constant: " << par[6] << endl;
    cout << "Mean: " << par[7] << endl;
    cout << "Error: " << par[8] << endl;
  TF1 *fitfun03=new TF1("fitfun03",xray_peak,six_lb,six_ub,5); 
  fitfun03->SetLineColor(2);
  fitfun03->SetParameter(0,par[7]);  // mean
  fitfun03->SetParameter(1,par[8]);  // sigma
  fitfun03->SetParameter(2,par[6]);  // constant
  fitfun03->SetParLimits(3,-1,1);
  fitfun03->SetParLimits(4,-1,1);  
  clbr_max_he->Fit("fitfun03","","",six_lb,six_ub); 


  clbr_max_he->Draw("");
  fitfun01->Draw("same");
  fitfun02->Draw("same");
  fitfun03->Draw("same");
  //

  //GET FIT PARAMETER
  Double_t mean[3]={0};
  Double_t sigma[3]={0};
  
  mean[0]=fitfun01->GetParameter(0);  
  sigma[0]=fitfun01->GetParError(0);
  
  mean[1]=fitfun02->GetParameter(0);  
  sigma[1]=fitfun02->GetParError(0);

  mean[2]=fitfun03->GetParameter(0);  
  sigma[2]=fitfun03->GetParError(0);


  for(int ii =0;ii<3;ii++)
    {
      printf("%f  %f   %f \n",real_energy_of_xray[ii],mean[ii],sigma[ii]); 
    } 
 
  TF1 *linefitfun=new TF1("linefitfun",linear,-2E+06,32E+06,2); 
  linefitfun->SetLineColor(4);
  linefitfun->SetParameter(0,1.19047e-01);  
  linefitfun->SetParameter(1,4.62746e-04);

  TCanvas *c3 = new TCanvas("c3");
  gStyle->SetOptFit(1);
  TGraphErrors *cali = new TGraphErrors(3,mean,real_energy_of_xray,sigma);
  cali->SetMarkerColor(2);
  cali->SetMarkerStyle(21);
  cali->Fit("linefitfun","","",-2E+06,32E+06);  
  cali->Draw("ap"); 


  printf("const %.8e\t slope  %.8e\n",linefitfun->GetParameter(0),linefitfun->GetParameter(1));
    
}
