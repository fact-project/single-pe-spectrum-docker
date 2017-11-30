#include <iomanip>
#include <iostream>

#include <TF1.h>
#include <TH2.h>
#include <TProfile2D.h>
#include <TFile.h>
#include <TFitResult.h>
#include <TStyle.h>

#include "MLog.h"
#include "MLogManip.h"
#include "MStatusArray.h"
#include "MStatusDisplay.h"
#include "MHCamera.h"
#include "MGeomCamFACT.h"
#include "MParameters.h"
#include "MArrayI.h"
#include "MRawRunHeader.h"
#include "PixelMap.h"

using namespace std;

// --------------------------------------------------------------------------

// Fit function for a single pe spectrum
Double_t fcn_g(Double_t *xx, Double_t *par)
{
    const Double_t ampl  = par[0];
    const Double_t gain  = par[1];
    const Double_t sigma = par[2]*gain;
    const Double_t cross = par[3];
    const Double_t shift = par[4];
    const Double_t noise = par[5]<0 ? sigma : par[5];
    const Double_t expo  = par[6];

    const Double_t P = cross*TMath::Exp(-cross);

    Double_t y = 0;
    for (int N=1; N<14; N++)
    {
        const Double_t muN  = N*gain + shift;
        const Double_t sigN = TMath::Sqrt(N*sigma*sigma + noise*noise);

        const Double_t p =
            TMath::Power(N*P, N-1)/TMath::Power(TMath::Factorial(N-1), expo);

        y += TMath::Gaus(xx[0], muN, sigN) * p / sigN;
    }

    const Double_t sig1 = TMath::Sqrt(sigma*sigma + noise*noise);
    return ampl*sig1*y;
}

// Calculate the crosstalk from the function parameters
Double_t xtalk(TF1 &f)
{
    Double_t cross = f.GetParameter(3);
    Double_t expo  = f.GetParameter(6);

    const Double_t P = cross*TMath::Exp(-cross);

    Double_t y = 0;
    for (int N=2; N<14; N++)
        y +=
            TMath::Power(N*P, N-1)/TMath::Power(TMath::Factorial(N-1), expo);

    return y / (y + 1);
}

// calculate the integral in units per millisecond
double integral(TF1 &func, TH1 &hist)
{
    const Double_t sigma = func.GetParameter(2)*func.GetParameter(1);
    const Double_t cross = func.GetParameter(3);
    const Double_t noise = func.GetParameter(5)<0 ? sigma : func.GetParameter(5);
    const Double_t expo  = func.GetParameter(6);

    const Double_t P = cross*TMath::Exp(-cross);

    Double_t sum = 0;
    for (int N=1; N<14; N++)
        sum +=
            TMath::Power(N*P, N-1)/TMath::Power(TMath::Factorial(N-1), expo);

    const Double_t scale = hist.GetBinWidth(1);
    const Double_t sig1  = TMath::Sqrt(sigma*sigma + noise*noise);
    const Double_t integ = func.GetParameter(0)*sum*sig1*sqrt(TMath::TwoPi())/scale;

    return integ/1e-9/1e6;
}


// --------------------------------------------------------------------------

// Print function parameters
void PrintFunc(TF1 &f, float integration_window=30)
{
    //cout << "--------------------" << endl;
    cout << "Ampl:       " << setw(8) << f.GetParameter(0) << " +/- " << f.GetParError(0) << endl;
    cout << "Gain:       " << setw(8) << f.GetParameter(1) << " +/- " << f.GetParError(1) << endl;
    cout << "Rel.sigma:  " << setw(8) << f.GetParameter(2) << " +/- " << f.GetParError(2) << endl;
    cout << "Baseline:   " << setw(8) << f.GetParameter(4)/integration_window << " +/- " << f.GetParError(4)/integration_window << endl;
    cout << "Crosstalk:  " << setw(8) << f.GetParameter(3) << " +/- " << f.GetParError(3) << endl;
    cout << "Pcrosstalk: " << setw(8) << xtalk(f) << endl;
    if (f.GetParameter(5)>=0)
        cout << "Noise:      " << setw(8) << f.GetParameter(5)/sqrt(integration_window) << " +/- " << f.GetParError(5)/sqrt(integration_window) << endl;
    cout << "Expo:       " << setw(8) << f.GetParameter(6) << " +/- " << f.GetParError(6) << endl;
    //cout << "--------------------" << endl;

}

// --------------------------------------------------------------------------

// The parameters for the function are the filenam from the gain extraction
// and the output filename
int fit_spectra(const char *filename = "20130222_018_018.root",
                     const char *outfile  = "20130222_018_018-fit.root",
                     bool fixednoise=true)
{
    // Read the mapping between bias channels and hardware channels
    PixelMap pmap;
    if (!pmap.Read("FACTmapV5a.txt"))
    {
        cout << "FACTmapV5a.txt not found." << endl;
        return 1;
    }

    // It is more correct to fit the integral, but this is super
    // slow, especially for 1440 pixel. To get that, one would have
    // to analytically integrate and fit this function.
    // (Currently the integral is switched off for the 1440 individual
    // spectra in all cases)
    bool fast = false; // Switch off using integral

    // Values which should be read from the file but not available atm
    Int_t integration_window = 30;

    // Map for which pixel shall be plotted and which not
    TArrayC usePixel(1440);
    memset(usePixel.GetArray(), 1, 1440);

    // List of Pixel that should be ignored in camera view
    usePixel[424]  = 0;
    usePixel[923]  = 0;
    usePixel[1208] = 0;
    usePixel[583]  = 0;
    usePixel[830]  = 0;
    usePixel[1399] = 0;
    usePixel[113]  = 0;
    usePixel[115]  = 0;
    usePixel[354]  = 0;
    usePixel[423]  = 0;
    usePixel[1195] = 0;
    usePixel[1393] = 0;

    cout << setprecision(3);

    // ======================================================
    // Read data and histograms from file

    TFile file(filename);
    if (file.IsZombie())
    {
        gLog << err << "Opening file '" << filename << "' failed." << endl;
        return 1;
    }

    MStatusArray arr;
    if (arr.Read()<=0)
    {
        gLog << err << "Reading of MStatusArray from '" << filename << "' failed." << endl;
        return 2;
    }

    TH2 *hsignal = (TH2*)arr.FindObjectInCanvas("Signal", "TH2F", "MHSingles");
    if (!hsignal)
    {
        gLog << err << "Histogram Signal not found in '" << filename << "'." << endl;
        return 3;
    }
    TH2 *htime = (TH2*)arr.FindObjectInCanvas("Time", "TH2F", "MHSingles");
    if (!htime)
    {
        gLog << err << "Histogram Time not found in '" << filename << "'." << endl;
        return 4;
    }
    TProfile2D *hpulse = (TProfile2D*)arr.FindObjectInCanvas("Pulse", "TProfile2D", "MHSingles");
    if (!hpulse)
    {
        gLog << err << "Histogram Pulse not found in '" << filename << "'." << endl;
        return 5;
    }
    TH2F *hbase = (TH2F*)arr.FindObjectInCanvas("Baseline", "TH2F", "MHBaseline");
    if (!hbase)
    {
        gLog << err << "Histogram Baseline not found in '" << filename << "'." << endl;
        return 6;
    }

    MRawRunHeader header;
    if (header.Read()<=0)
    {
        gLog << err << "MRawRunheader not found in '" << filename << "'." << endl;
        return 7;
    }

    MParameterI par("NumEvents");
    if (par.Read()<=0)
    {
        gLog << err << "NumEvents not found in '" << filename << "'." << endl;
        return 8;
    }

    MArrayI ext;
    if (ext.Read("ExtractionRange")<=0)
    {
        gLog << err << "ExtractionRange not found in '" << filename << "'." << endl;
        return 9;

    }

    // ======================================================

    MStatusDisplay *d = new MStatusDisplay;

    // Camera geometry for displays
    MGeomCamFACT fact;

    // ------------------ Spectrum Fit ---------------
    // Instantiate the display histograms
    MHCamera cRate(fact);
    MHCamera cGain(fact);
    MHCamera cRelSigma(fact);
    MHCamera cCrosstalk(fact);
    MHCamera cBaseline(fact);
    MHCamera cNoise(fact);
    MHCamera cChi2(fact);
    MHCamera cNormGain(fact);
    MHCamera cFitProb(fact);
    MHCamera cCrosstalkP(fact);
    MHCamera cCoeffR(fact);

    // Set name and title for the histograms
    cRate.SetNameTitle      ("Rate",      "Dark count rate");
    cGain.SetNameTitle      ("Gain",      "Gain distribution");
    cRelSigma.SetNameTitle  ("RelSigma",  "Rel. Sigma");
    cCrosstalk.SetNameTitle ("Crosstalk", "Crosstalk probability");
    cBaseline.SetNameTitle  ("Baseline",  "Baseline per sample");
    cNoise.SetNameTitle     ("Noise",     "Noise per sample");
    cChi2.SetNameTitle      ("Chi2",      "\\Chi^2");
    cNormGain.SetNameTitle  ("NormGain",  "Normalized gain");
    cFitProb.SetNameTitle   ("FitProb",   "Root's fit probability");
    cCrosstalkP.SetNameTitle("Pxtalk",    "Crosstalk coeff. P");
    cCoeffR.SetNameTitle    ("CoeffR",    "Coefficient R");

    // Instantiate 1D histograms for the distributions
    // including TM channels
    TH1F hRate1      ("Rate1",      "Dark count rate",       150,  0,    15);
    TH1F hGain1      ("Gain1",      "Gain distribution",     100,  0,   400);
    TH1F hRelSigma1  ("RelSigma1",  "Rel. Sigma",            160,  0,  0.40);
    TH1F hCrosstalk1 ("Crosstalk1", "Crosstalk probability",  90,  0,  0.30);
    TH1F hBaseline1  ("Baseline1",  "Baseline per sample",    75, -7.5, 7.5);
    TH1F hNoise1     ("Noise1",     "Noise per sample",       60,  0,    30);
    TH1F hChiSq1     ("ChiSq1",     "\\Chi^2",               200,  0,     4);
    TH1F hNormGain1  ("NormGain1",  "Normalized gain",        51,  0.5, 1.5);
    TH1F hFitProb1   ("FitProb1",   "FitProb distribution",  100,  0,     1);
    TH1F hCrosstalkP1("Pxtalk1",    "Crosstalk coeff.",       90,  0,   0.3);
    TH1F hCoeffR1    ("CoeffR1",    "Coefficient R",          90, -1,     2);

    // excluding TM channels
    TH1F hRate2      ("Rate2",      "Dark count rate",       150,  0,    15);
    TH1F hGain2      ("Gain2",      "Gain distribution",     100,  0,   400);
    TH1F hRelSigma2  ("RelSigma2",  "Rel. Sigma",            160,  0,  0.40);
    TH1F hCrosstalk2 ("Crosstalk2", "Crosstalk probability",  90,  0,  0.30);
    TH1F hBaseline2  ("Baseline2",  "Baseline per sample",    75, -7.5, 7.5);
    TH1F hNoise2     ("Noise2",     "Noise per sample",       60,  0,    30);
    TH1F hChiSq2     ("ChiSq2",     "\\Chi^2",               200,  0,     4);
    TH1F hNormGain2  ("NormGain2",  "Normalized gain",        51,  0.5, 1.5);
    TH1F hFitProb2   ("FitProb2",   "FitProb distribution",  100,  0,     1);
    TH1F hCrosstalkP2("Pxtalk2",    "Crosstalk coeff.",       90,  0,   0.3);
    TH1F hCoeffR2    ("CoeffR2",    "Coefficient R",          90, -1,     2);

    // Histigram for the sum of all spectrums
    TH1F hSum("Sum1", "Signal spectrum of all pixels",
              hsignal->GetNbinsY(), hsignal->GetYaxis()->GetXmin(),  hsignal->GetYaxis()->GetXmax());
    hSum.SetXTitle("pulse integral [mV sample]");
    hSum.SetYTitle("Counts");
    hSum.SetStats(false);
    hSum.Sumw2();

    // Histogram for the sum of all pixels excluding the ones with faulty fits
    TH1F hSumClear1("SumC1", "Signal spectrum of all pixels (incl TM)",
                    hsignal->GetNbinsY(), hsignal->GetYaxis()->GetXmin(),  hsignal->GetYaxis()->GetXmax());
    hSumClear1.SetXTitle("pulse integral [mV sample]");
    hSumClear1.SetYTitle("Counts");
    hSumClear1.SetStats(false);
    hSumClear1.SetLineColor(kBlue);
    hSumClear1.Sumw2();

    TH1F hSumClear2("SumC2", "Signal spectrum of all pixels (excp TM)",
                    hsignal->GetNbinsY(), hsignal->GetYaxis()->GetXmin(),  hsignal->GetYaxis()->GetXmax());
    hSumClear2.SetXTitle("pulse integral [mV sample]");
    hSumClear2.SetYTitle("Counts");
    hSumClear2.SetStats(false);
    hSumClear2.SetLineColor(kBlue);
    hSumClear2.Sumw2();

    // Arrival time spectrum of the extracted pulses
    TH1F hTime("Time", "Arrival time spectrum", htime->GetNbinsY(), htime->GetYaxis()->GetXmin(), htime->GetYaxis()->GetXmax());
    hTime.SetXTitle("pulse arrival time [sample]");
    hTime.SetYTitle("Counts");
    hTime.SetStats(false);

    // average pulse shape of the extracted pulses
    TH1F hPulse("Puls", "Average pulse", hpulse->GetNbinsY(), hpulse->GetYaxis()->GetXmin(), hpulse->GetYaxis()->GetXmax());
    hPulse.SetXTitle("pulse arrival time [sample]");
    hPulse.SetYTitle("Counts");
    hPulse.SetStats(false);

    // Spektrum for the normalized individual spectra
    TH1F hSumScale1("SumScale1", "Signal spectrum of all pixels (incl TM)",  120, 0.05, 12.05);
    hSumScale1.SetXTitle("pulse integral [pe]");
    hSumScale1.SetYTitle("Rate");
    hSumScale1.SetStats(false);
    hSumScale1.Sumw2();

    TH1F hSumScale2("SumScale2", "Signal spectrum of all pixels (excl TM)",  120, 0.05, 12.05);
    hSumScale2.SetXTitle("pulse integral [pe]");
    hSumScale2.SetYTitle("Rate");
    hSumScale2.SetStats(false);
    hSumScale2.Sumw2();

    // define fit function for Amplitudespectrum
    TF1 func("spektrum", fcn_g, 0, hSum.GetXaxis()->GetXmax(), 7);
    func.SetNpx(2000);
    func.SetParNames("Maximum", "Gain", "Sigma", "XtalkProb", "Offset", "Noise", "Expo");
    func.SetLineColor(kRed);

    //--------------------- fill sum spectrum --------------------------------

    d->SetStatusLine("Calculating sum spectrum", 0);

    // Begin of Loop over Pixels
    for (Int_t pixel = 0; pixel < hsignal->GetNbinsX(); pixel++)
    {
        //jump to next pixel if the current is marked as faulty
        if (usePixel[pixel]==0)
            continue;

        TH1D *hist = hsignal->ProjectionY("proj", pixel+1, pixel+1);
        hSum.Add(hist);
        delete hist;
    }

    //----------------- get starting values -------------------------------

    hSum.GetXaxis()->SetRangeUser(150, hSum.GetXaxis()->GetXmax());

    const Int_t    maxbin   = hSum.GetMaximumBin();
    const Double_t maxpos   = hSum.GetBinCenter(maxbin);
    const Double_t binwidth = hSum.GetBinWidth(maxbin);
    const Double_t ampl     = hSum.GetBinContent(maxbin);

    double fwhmSum = 0;

    //Calculate full width half Maximum
    for (int i=1; i<maxbin; i++)
        if (hSum.GetBinContent(maxbin-i)+hSum.GetBinContent(maxbin+i)<ampl)
        {
            fwhmSum = 2*(i-0.5)*binwidth;
            break;
        }

    if (fwhmSum==0)
    {
        gLog << warn << "Could not determine start value for sigma." << endl;
    }

    Double_t sigma_est = fwhmSum/2.3548; // FWHM = 2*sqrt(2*ln(2))*sigma

    Double_t fitmin = maxpos-3*sigma_est; // was 3
    Double_t fitmax = hSum.GetXaxis()->GetXmax();

    // ------------------- fit --------------------------------

    //Fit and draw spectrum
    func.SetParLimits(0, 0, 2*ampl);       // Amplitude
    func.SetParLimits(1, 0, 2*maxpos);     // Gain
    func.SetParLimits(2, 0, 1);            // Sigma
    func.SetParLimits(3, 0, 1);            // Crosstalk
    if (!fixednoise)
        func.SetParLimits(5, 0, 150);      // Noise
    func.SetParLimits(6, 0, 2);            // Expo


    func.SetParameter(0, ampl);                         // Amplitude
    func.SetParameter(1, maxpos);                       // Gain
    func.SetParameter(2, 0.1);                          // Sigma
    func.SetParameter(3, 0.16);                         // Crosstalk
    func.SetParameter(4, 0*integration_window);         // Baseline
    if (fixednoise)
        func.FixParameter(5, -1);                       // Noise
    else
        func.SetParameter(5, 0.1*maxpos);               // Noise

    func.SetParameter(6, 0.95);                            // Expo

    func.SetRange(fitmin, fitmax);
    hSum.Fit(&func, fast?"N0QSR":"IN0QSR");

    Double_t res_par[7];
    func.GetParameters(res_par);

    //func.FixParameter(6, func.GetParameter(6));                          // Expo

    // ------------------ display result -------------------------------

    cout << "--------------------" << endl;
    cout << "AmplEst:    " << ampl      << endl;
    cout << "GainEst:    " << maxpos    << endl;
    cout << "SigmaEst:   " << sigma_est << endl;
    PrintFunc(func, integration_window);
    cout << "--------------------" << endl;

    gROOT->SetSelectedPad(0);
    TCanvas &c11 = d->AddTab("SumHist");
    c11.cd();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();
    hSum.GetXaxis()->SetRange();
    hSum.DrawCopy("hist");
    func.DrawCopy("same");

    // ===================================================================
    //  Gain Calculation for each Pixel
    // ===================================================================

    // counter for number of processed pixel
    int count_ok = 0;

    // Begin of Loop over Pixels
    for (Int_t pixel=0; pixel<hsignal->GetNbinsX(); pixel++)
    {
        // User information
        d->SetStatusLine(Form("Fitting pixel #%d", pixel), 0);
        d->SetProgressBarPosition((pixel+1.)/hsignal->GetNbinsX(), 1);

        // Skip pixels known to be faulty
        if (usePixel[pixel]==0)
            continue;

        //Projectipon of a certain Pixel out of Ampl.Spectrum
        TH1D *hist = hsignal->ProjectionY("proj", pixel+1, pixel+1);
        hist->SetDirectory(0);

        if (hist->GetEntries()<100)
        {
            gLog << warn << pixel << " ...histogram empty." << endl;
            usePixel[pixel] = 0;
            delete hist;
            continue;
        }

        //Rebin Projection
        hist->Rebin(2);

        // Fit range
        hist->GetXaxis()->SetRangeUser(150, hist->GetXaxis()->GetXmax());

        // Determine start values
        const Int_t    maxBin = hist->GetMaximumBin();
        const Double_t maxPos = hist->GetBinCenter(maxBin);

        const Double_t gain    = res_par[1];
        const Double_t GainRMS = res_par[2];

        const double fit_min = maxPos-GainRMS*gain*2.5;
        const double fit_max = fitmax;//maxPos+gain*(maxOrder-0.5);

        TArrayD cpy_par(7, res_par);

        cpy_par[0] = hist->GetBinContent(maxBin);
        cpy_par[1] = maxPos-res_par[4];  // correct position for avg baseline

        func.SetParameters(cpy_par.GetArray());
        func.SetParLimits(0, 0, 2*cpy_par[0]);
        func.SetParLimits(1, 0, 2*cpy_par[1]);

        // For individual spectra, the average fit yields 1 anyway
        //func.SetParameter(6, 0); // Expo

        // ----------- Fit Pixels spectrum ---------------

        const TFitResultPtr rc = hist->Fit(&func, /*fast?*/"LLN0QS"/*:"LLIN0QS"*/, "", fit_min, fit_max);

        // ----------- Calculate quality parameter ---------------

        Int_t b1 = hist->GetXaxis()->FindFixBin(fit_min);
        Int_t b2 = hist->GetXaxis()->FindFixBin(fit_max);

        Double_t chi2 = 0;
        Int_t    cnt  = 0;
        for (int i=b1; i<=b2; i++)
        {
            if (hist->GetBinContent(i)<1.5 || func.Eval(hist->GetBinCenter(i))<1.5)
                continue;

            const Double_t y = func.Integral(hist->GetBinLowEdge(i), hist->GetBinLowEdge(i+1));
            const Double_t v = hist->GetBinContent(i)*hist->GetBinWidth(i);

            const Double_t chi = (v-y)/v;

            chi2 += chi*chi;
            cnt ++;
        }

        chi2 = cnt==0 ? 0 : sqrt(chi2/cnt);

        // ----------------- Fit result --------------------

        const double fit_prob   = rc->Prob();

        const float fRate      = integral(func, *hist)/(ext[pixel]*0.5);
        const float fGain      = func.GetParameter(1);
        const float fGainRMS   = func.GetParameter(2);
        const float fCrosstalkP= func.GetParameter(3);
        const float fCrosstlk  = xtalk(func);
        const float fOffset    = func.GetParameter(4);
        const float fNoise     = func.GetParameter(5)<0 ? fGainRMS*fGain/sqrt(integration_window) : func.GetParameter(5)/sqrt(integration_window);
        const float fCoeffR    = func.GetParameter(6);

        // Fill histograms with result values
        cRate.SetBinContent(      pixel+1, fRate);
        cGain.SetBinContent(      pixel+1, fGain);
        cRelSigma.SetBinContent(  pixel+1, fGainRMS);
        cCrosstalk.SetBinContent( pixel+1, fCrosstlk);
        cBaseline.SetBinContent(  pixel+1, fOffset/integration_window);
        cNoise.SetBinContent(     pixel+1, fNoise);
        cChi2.SetBinContent(      pixel+1, chi2);
        cNormGain.SetBinContent(  pixel+1, fGain/gain);
        cFitProb.SetBinContent(   pixel+1, fit_prob);
        cCrosstalkP.SetBinContent(pixel+1, fCrosstalkP);
        cCoeffR.SetBinContent(    pixel+1, fCoeffR);

        // ======================================================

        // Try to determine faulty fits

        bool ok = int(rc)==0;

        // mark pixels suspicious with failed fit
        if (!ok)
            gLog << warn <<  pixel << " ...fit failed!" << endl;

        // mark pixels suspicious with negative GainRMS
        if (fabs(fGain/gain-1)>0.3)
        {
            gLog << warn <<  pixel << " ...gain deviates more than 30% from sum-gain." << endl;
            ok = 0;
        }

        if (fabs(fOffset/integration_window)>3)
        {
            gLog << warn <<  pixel << " ...baseline deviates." << endl;
            ok = 0;
        }

        // cancel out pixel where the fit was not succsessfull
        usePixel[pixel] = ok;

        // Plot pixel 0 and 5 (TM) and all faulty fits
        if (pixel==0 || pixel==5 || !ok)
        {
            TCanvas &c = d->AddTab(Form("Pix%d", pixel));
            c.cd();
            gPad->SetLogy();
            gPad->SetGridx();
            gPad->SetGridy();

            hist->SetStats(false);
            hist->SetXTitle("Extracted signal");
            hist->SetYTitle("Counts");
            hist->SetName(Form("Pix%d", pixel));
            hist->GetXaxis()->SetRange();
            hist->DrawCopy("hist")->SetDirectory(0);
            func.DrawCopy("SAME")->SetRange(fit_min, fit_max);

            cout << "--------------------" << endl;
            cout << "Pixel:      "   << pixel  << endl;
            cout << "fit prob:   "   << fit_prob  << endl;
            cout << "AmplEst:    "   << cpy_par[0] << endl;
            cout << "GainEst:    "   << cpy_par[1] << endl;
            PrintFunc(func, integration_window);
            cout << "--------------------" << endl;
        }

        if (!ok)
        {
            delete hist;
            continue;
        }

        // Fill Parameters into histograms
        hRate1.Fill(      fRate);
        hGain1.Fill(      fGain);
        hRelSigma1.Fill(  fGainRMS);
        hCrosstalk1.Fill( fCrosstlk);
        hBaseline1.Fill(  fOffset/integration_window);
        hNoise1.Fill(     fNoise);
        hChiSq1.Fill(     chi2);
        hNormGain1.Fill(  fGain/gain);
        hFitProb1.Fill(   fit_prob);
        hCrosstalkP1.Fill(fCrosstalkP);
        hCoeffR1.Fill(    fCoeffR);

        if (!pmap.index(pixel).isTM())
        {
            hRate2.Fill(      fRate);
            hGain2.Fill(      fGain);
            hRelSigma2.Fill(  fGainRMS);
            hCrosstalk2.Fill( fCrosstlk);
            hBaseline2.Fill(  fOffset/integration_window);
            hNoise2.Fill(     fNoise);
            hChiSq2.Fill(     chi2);
            hNormGain2.Fill(  fGain/gain);
            hFitProb2.Fill(   fit_prob);
            hCrosstalkP2.Fill(fCrosstalkP);
            hCoeffR2.Fill(    fCoeffR);
        }

        // Fill sum spectrum
        for (int b=1; b<=hist->GetNbinsX(); b++)
            hSumScale1.Fill((hist->GetBinCenter(b)-fOffset)/fGain, hist->GetBinContent(b));

        if (!pmap.index(pixel).isTM())
            for (int b=1; b<=hist->GetNbinsX(); b++)
                hSumScale2.Fill((hist->GetBinCenter(b)-fOffset)/fGain, hist->GetBinContent(b));

        delete hist;

        // Because of the rebinning...
        hist = hsignal->ProjectionY("proj", pixel+1, pixel+1);
        hSumClear1.Add(hist);
        if (!pmap.index(pixel).isTM())
            hSumClear2.Add(hist);
        delete hist;

        hist = htime->ProjectionY("proj", pixel+1, pixel+1);
        hTime.Add(hist);
        delete hist;

        hist = hpulse->ProjectionY("proj", pixel+1, pixel+1);
        hPulse.Add(hist);
        delete hist;

        count_ok++;
    }

    //------------------fill histograms-----------------------
    // Display only pixels used and with valid fits

    cRate.SetUsed(usePixel);
    cGain.SetUsed(usePixel);
    cRelSigma.SetUsed(usePixel);
    cCrosstalk.SetUsed(usePixel);
    cBaseline.SetUsed(usePixel);
    cNoise.SetUsed(usePixel);
    cChi2.SetUsed(usePixel);
    cNormGain.SetUsed(usePixel);
    cFitProb.SetUsed(usePixel);
    cCrosstalkP.SetUsed(usePixel);
    cCoeffR.SetUsed(usePixel);

    // --------------------------------------------------------
    // Display data

    TCanvas *canv = &d->AddTab("Cams1");
    canv->Divide(3,2);

    canv->cd(1);
    cRate.DrawCopy();

    canv->cd(2);
    cGain.DrawCopy();

    canv->cd(3);
    cBaseline.DrawCopy();

    canv->cd(4);
    cRelSigma.DrawCopy();

    canv->cd(5);
    cCrosstalk.DrawCopy();

    canv->cd(6);
    cNoise.DrawCopy();


    canv = &d->AddTab("Cams2");
    canv->Divide(3,2);

    canv->cd(1);
    cFitProb.DrawCopy();

    canv->cd(2);
    cChi2.DrawCopy();

    canv->cd(4);
    cCoeffR.DrawCopy();

    canv->cd(5);
    cCrosstalkP.DrawCopy();

    // --------------------------------------------------------

    gStyle->SetOptFit(1);

    canv = &d->AddTab("Hists1");
    canv->Divide(3,2);

    TH1 *hh = 0;

    canv->cd(1);
    hh = hRate1.DrawCopy();
    hh = hRate2.DrawCopy("same");

    canv->cd(2);
    hh = hGain1.DrawCopy();
    hh = hGain2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(3);
    hh = hBaseline1.DrawCopy();
    hh = hBaseline2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(4);
    hh = hRelSigma1.DrawCopy();
    hh = hRelSigma2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(5);
    hh = hCrosstalk1.DrawCopy();
    hh = hCrosstalk2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(6);
    hh = hNoise1.DrawCopy();
    hh = hNoise2.DrawCopy("same");
    hh->Fit("gaus");

    // --------------------------------------------------------

    canv = &d->AddTab("Hists2");
    canv->Divide(3,2);

    canv->cd(1);
    gPad->SetLogy();
    hh = hFitProb1.DrawCopy();
    hh = hFitProb2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(2);
    hChiSq1.DrawCopy();
    hChiSq2.DrawCopy("same");

    canv->cd(4);
    hh = hCoeffR1.DrawCopy();
    hh = hCoeffR2.DrawCopy("same");
    hh->Fit("gaus");

    canv->cd(5);
    hh = hCrosstalkP1.DrawCopy();
    hh = hCrosstalkP2.DrawCopy("same");
    hh->Fit("gaus");

    // --------------------------------------------------------

    canv = &d->AddTab("NormGain");
    canv->Divide(2,1);

    canv->cd(1);
    cNormGain.SetMinimum(0.8);
    cNormGain.SetMaximum(1.2);
    cNormGain.DrawCopy();

    canv->cd(2);
    gPad->SetLogy();
    hh = hNormGain1.DrawCopy();
    hh = hNormGain2.DrawCopy("same");
    hh->Fit("gaus");

    //------------------ Draw gain corrected sum specetrum -------------------
    gROOT->SetSelectedPad(0);
    c11.cd();
    hSumClear1.DrawCopy("hist same");

    //-------------------- fit gain corrected sum spectrum -------------------

    gROOT->SetSelectedPad(0);
    TCanvas &c11b = d->AddTab("CleanHist1");
    c11b.cd();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    const Int_t    maxbin1 = hSumClear1.GetMaximumBin();
    const Double_t ampl1   = hSumClear1.GetBinContent(maxbin1);

    func.SetParameters(res_par);
    func.SetParLimits(0, 0, 2*ampl1);
    func.SetParameter(0, ampl1);
    func.ReleaseParameter(6);

    hSumClear1.Fit(&func, fast?"LN0QSR":"LIN0QSR");

    hSumClear1.DrawCopy();
    func.DrawCopy("same");

    cout << "--------------------" << endl;
    PrintFunc(func, integration_window);
    cout << "--------------------" << endl;

    //-------------------- fit gain corrected sum spectrum -------------------

    gROOT->SetSelectedPad(0);
    TCanvas &c11c = d->AddTab("CleanHist2");
    c11c.cd();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    const Int_t    maxbin1b = hSumClear2.GetMaximumBin();
    const Double_t ampl1b   = hSumClear2.GetBinContent(maxbin1b);

    func.SetParameters(res_par);
    func.SetParLimits(0, 0, 2*ampl1b);
    func.SetParameter(0, ampl1b);
    func.ReleaseParameter(6);

    hSumClear2.Fit(&func, fast?"LN0QSR":"LIN0QSR");

    hSumClear2.DrawCopy();
    func.DrawCopy("same");

    cout << "--------------------" << endl;
    PrintFunc(func, integration_window);
    cout << "--------------------" << endl;

    //-------------------- fit gain corrected sum spectrum -------------------

    gROOT->SetSelectedPad(0);
    TCanvas &c12 = d->AddTab("GainHist1");
    c12.cd();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    const Int_t    maxbin2 = hSumScale1.GetMaximumBin();
    const Double_t ampl2   = hSumScale1.GetBinContent(maxbin2);

    //Set fit parameters
    Double_t par2[7] =
    {
        ampl2, 1, 0.1, res_par[3], 0, res_par[5]<0 ? -1 : res_par[5]/res_par[1], res_par[6]
    };

    func.SetParameters(par2);
    func.SetParLimits(0, 0, 2*ampl2);
    func.FixParameter(1, 1);
    func.FixParameter(4, 0);

    func.SetRange(0.62, 9);
    hSumScale1.Fit(&func, fast?"LN0QSR":"LIN0QSR");

    hSumScale1.DrawCopy();
    func.DrawCopy("same");

    cout << "--------------------" << endl;
    PrintFunc(func, integration_window);
    cout << "--------------------" << endl;

    //-------------------- fit gain corrected sum spectrum -------------------

    gROOT->SetSelectedPad(0);
    TCanvas &c12b = d->AddTab("GainHist2");
    c12b.cd();
    gPad->SetLogy();
    gPad->SetGridx();
    gPad->SetGridy();

    const Int_t    maxbin2b = hSumScale2.GetMaximumBin();
    const Double_t ampl2b   = hSumScale2.GetBinContent(maxbin2b);

    //Set fit parameters
    Double_t par2b[7] =
    {
        ampl2b, 1, 0.1, res_par[3], 0, res_par[5]<0 ? -1 : res_par[5]/res_par[1], res_par[6]
    };

    func.SetParameters(par2b);
    func.SetParLimits(0, 0, 2*ampl2b);
    func.FixParameter(1, 1);
    func.FixParameter(4, 0);

    func.SetRange(0.62, 9);
    hSumScale2.Fit(&func, fast?"LN0QSR":"LIN0QSR");

    hSumScale2.DrawCopy();
    func.DrawCopy("same");

    cout << "--------------------" << endl;
    PrintFunc(func, integration_window);
    cout << "--------------------" << endl;

    //--------fit gausses to peaks of gain corrected sum specetrum -----------

    d->AddTab("ArrTime");
    gPad->SetGrid();
    hTime.DrawCopy();

    // -----------------------------------------------------------------

    d->AddTab("Pulse");
    gPad->SetGrid();
    hPulse.DrawCopy();

    // ================================================================

    cout << "Saving results to '" << outfile << "'" << endl;
    d->SaveAs(outfile);
    cout << "..success!" << endl;

    TFile f(outfile, "UPDATE");
    par.Write();
    ext.Write("ExtractionRange");
    header.Write();

    return 0;
}
