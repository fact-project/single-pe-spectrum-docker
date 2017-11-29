#include <TStyle.h>
#include <TCanvas.h>
#include <TSystem.h>
#include <TF1.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TMath.h>
#include <TGraph.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>
#include <TFile.h>

#include <cstdio>
#include <stdio.h>
#include <stdint.h>

#include "MH.h"
#include "MArrayI.h"
#include "MLog.h"
#include "MLogManip.h"
#include "MDirIter.h"
#include "MFillH.h"
#include "MEvtLoop.h"
#include "MCamEvent.h"
#include "MHCamEvent.h"
#include "MGeomApply.h"
#include "MTaskList.h"
#include "MParList.h"
#include "MContinue.h"
#include "MBinning.h"
#include "MDrsCalibApply.h"
#include "MDrsCalibration.h"
#include "MRawFitsRead.h"
#include "MBadPixelsCam.h"
#include "MStatusDisplay.h"
#include "MTaskInteractive.h"
#include "MPedestalSubtractedEvt.h"
#include "MHCamera.h"
#include "MGeomCamFACT.h"
#include "MRawRunHeader.h"
#include "MPedestalCam.h"
#include "MPedestalPix.h"
#include "MParameters.h"

using namespace std;
using namespace TMath;

// Structure to store the extracted value and the extracted time
struct Single
{
    float    fSignal;
    float    fTime;
};

//Storage class to kKeep a list of the extracted single
class MSingles : public MParContainer, public MCamEvent
{
    Int_t fIntegrationWindow;

    vector<vector<Single>> fData;

public:
    MSingles(const char *name=NULL, const char *title=NULL) : fIntegrationWindow(30)
    {
        fName = name ? name : "MSingles";
        fName = title ? title : "Storeage for vector of singles";
    }

    void InitSize(const UInt_t i)
    {
        fData.resize(i);
    }

    vector<Single> &operator[](UInt_t i) { return fData[i]; }
    vector<Single> &GetVector(UInt_t i)  { return fData[i]; }

    UInt_t GetNumPixels() const { return fData.size(); }

    void SetIntegrationWindow(Int_t w) { fIntegrationWindow = w; }
    Int_t GetIntegrationWindow() const { return fIntegrationWindow; }

    Bool_t GetPixelContent(Double_t &, Int_t , const MGeomCam &, Int_t) const
    {
        return kTRUE;
    }
    void   DrawPixelContent(Int_t) const { }

    ClassDef(MSingles, 1)
};

// Histogram class to extract the baseline
class MHBaseline : public MH
{
    TH2F fBaseline;

    MPedestalCam *fPedestal;

    // The baseline is only extracted where also the signal is extracted
    // FIXME: Make sure this is consistent with MExtractSingles
    UShort_t fSkipStart;
    UShort_t fSkipEnd;

public:
    MHBaseline() : fPedestal(0), fSkipStart(20), fSkipEnd(10)
    {
        fName = "MHBaseline";

        // Setup the histogram
        fBaseline.SetName("Baseline");
        fBaseline.SetTitle("Median spectrum");
        fBaseline.SetXTitle("Pixel [idx]");
        fBaseline.SetYTitle("Median baseline [mV]");
        fBaseline.SetDirectory(NULL);
    }

    Bool_t ReInit(MParList *plist)
    {
        fPedestal = (MPedestalCam*)plist->FindCreateObj("MPedestalCam");
        if (!fPedestal)
            return kFALSE;

        const MRawRunHeader *header = (MRawRunHeader*)plist->FindObject("MRawRunHeader");
        if (!header)
        {
            *fLog << err << "MRawRunHeader not found... abort." << endl;
            return kFALSE;
        }

        // Bin width should be around 1 dac count which is about 0.5mV
        MBinning binsx, binsy;
        binsx.SetEdges(header->GetNumNormalPixels(), -0.5, header->GetNumNormalPixels()-0.5);
        binsy.SetEdges(100, -20.5, 29.5);

        // Setup binnning
        MH::SetBinning(fBaseline, binsx, binsy);

        return kTRUE;
    }

    // Fill the samples into the histogram
    Int_t Fill(const MParContainer *par, const Stat_t)
    {
        const MPedestalSubtractedEvt *evt = dynamic_cast<const MPedestalSubtractedEvt*>(par);

        const Int_t n = evt->GetNumSamples()-fSkipStart-fSkipEnd;

        // Loop over all pixels
        for (int pix=0; pix<evt->GetNumPixels(); pix++)
        {
            // Get samples for each pixel
            const Float_t *ptr = evt->GetSamples(pix);

            // Average two consecutive samples
            for (int i=0; i<n; i+=2)
            {
                const Double_t val = 0.5*ptr[i+fSkipStart]+0.5*ptr[i+1+fSkipStart];
                fBaseline.Fill(pix, val);
            }
        }

        return kTRUE;
    }

    // Extract the baseline value from the distrbutions
    Bool_t Finalize()
    {
        if (!fPedestal)
            return kTRUE;

        fPedestal->InitSize(fBaseline.GetNbinsX());
        fPedestal->SetNumEvents(GetNumExecutions());

        Int_t    cnt = 0;
        Double_t avg = 0;
        Double_t rms = 0;

        // Loop over all 'pixels'
        for (int x=0; x<fBaseline.GetNbinsX(); x++)
        {
            // Get the corresponding slice from the histogram
            TH1D *hist = fBaseline.ProjectionY("proj", x+1, x+1);

            // Get the maximum bin
            const Int_t bin  = hist->GetMaximumBin();

            // Calculate a parabola through this and the surrounding points
            // on logarithmic values (that's a gaussian)

            //
            // Quadratic interpolation
            //
            // calculate the parameters of a parabula such that
            //    y(i)  = a + b*x(i) +   c*x(i)^2
            //    y'(i) =     b      + 2*c*x(i)
            //
            //

            // -------------------------------------------------------------------------
            // a = y2;
            // b = (y3-y1)/2;
            // c = (y3+y1)/2 - y2;

            const Double_t v1 = hist->GetBinContent(bin-1);
            const Double_t v2 = hist->GetBinContent(bin);
            const Double_t v3 = hist->GetBinContent(bin+1);
            if (v1<=0 || v2<=0 || v3<=0)
                continue;

            const Double_t y1 = TMath::Log(v1);
            const Double_t y2 = TMath::Log(v2);
            const Double_t y3 = TMath::Log(v3);

            //Double_t a = y2;
            const Double_t b = (y3-y1)/2;
            const Double_t c = (y1+y3)/2 - y2;
            if (c>=0)
                continue;

            const Double_t w  = -1./(2*c);
            const Double_t dx =  b*w;

            if (dx<-1 || dx>1)
                continue;

            // y = exp( - (x-k)^2 / s^2 / 2 )
            //
            // -2*s^2 * log(y) = x^2 - 2*k*x + k^2
            //
            // a = (k/s0)^2/2
            // b = k/s^2
            // c = -1/(2s^2)      -->    s = sqrt(-1/2c)

            const Double_t binx = hist->GetBinCenter(bin);
            const Double_t binw = hist->GetBinWidth(bin);

            //const Double_t p = hist->GetBinCenter(hist->GetMaximumBin());
            const Double_t p = binx + dx*binw;

            avg += p;
            rms += p*p;

            cnt++;

            // Store baseline and sigma
            MPedestalPix &pix = (*fPedestal)[x];

            pix.SetPedestal(p);
            pix.SetPedestalRms(TMath::Sqrt(w)*binw);

            delete hist;
        }

        avg /= cnt;
        rms /= cnt;

        cout << "Baseline(" << cnt << "): " << avg << " +- " << sqrt(rms-avg*avg) << endl;

        return kTRUE;
    }

    // Draw histogram
    void Draw(Option_t *)
    {
        TVirtualPad *pad = gPad ? gPad : MakeDefCanvas(this);

        AppendPad("");

        pad->SetBorderMode(0);
        pad->SetFrameBorderMode(0);
        fBaseline.Draw("colz");
    }


    ClassDef(MHBaseline, 1);
};

// Histogram class for the signal and time distribution as
// well as the pulse shape
class MHSingles : public MH
{
    TH2F       fSignal;
    TH2F       fTime;
    TProfile2D fPulse;

    UInt_t fNumEntries;

    MSingles               *fSingles;   //!
    MPedestalSubtractedEvt *fEvent;     //!
    MBadPixelsCam          *fBadPix;    //!

public:
    MHSingles() : fNumEntries(0), fSingles(0), fEvent(0)
    {
        fName = "MHSingles";

        // Setup histograms
        fSignal.SetName("Signal");
        fSignal.SetTitle("pulse integral spectrum");
        fSignal.SetXTitle("pixel [SoftID]");
        fSignal.SetYTitle("time [a.u.]");
        fSignal.SetDirectory(NULL);

        fTime.SetName("Time");
        fTime.SetTitle("arival time spectrum");
        fTime.SetXTitle("pixel [SoftID]");
        fTime.SetYTitle("time [a.u.]");
        fTime.SetDirectory(NULL);

        fPulse.SetName("Pulse");
        fPulse.SetTitle("average pulse shape spectrum");
        fPulse.SetXTitle("pixel [SoftID]");
        fPulse.SetYTitle("time [a.u.]");
        fPulse.SetDirectory(NULL);
    }

    Bool_t SetupFill(const MParList *plist)
    {
        fSingles = (MSingles*)plist->FindObject("MSingles");
        if (!fSingles)
        {
            *fLog << err << "MSingles not found... abort." << endl;
            return kFALSE;
        }

        fEvent = (MPedestalSubtractedEvt*)plist->FindObject("MPedestalSubtractedEvt");
        if (!fEvent)
        {
            *fLog << err << "MPedestalSubtractedEvt not found... abort." << endl;
            return kFALSE;
        }

        fBadPix = (MBadPixelsCam*)plist->FindObject("MBadPixelsCam");
        if (!fBadPix)
        {
            *fLog << err << "MBadPixelsCam not found... abort." << endl;
            return kFALSE;
        }

        fNumEntries = 0;

        return kTRUE;
    }

    Bool_t ReInit(MParList *plist)
    {
        const MRawRunHeader *header = (MRawRunHeader*)plist->FindObject("MRawRunHeader");
        if (!header)
        {
            *fLog << err << "MRawRunHeader not found... abort." << endl;
            return kFALSE;
        }

        // Setup binning
        const Int_t w = fSingles->GetIntegrationWindow();

        MBinning binsx, binsy;
        binsx.SetEdges(fSingles->GetNumPixels(), -0.5, fSingles->GetNumPixels()-0.5);
        binsy.SetEdges(22*w, -10*w, 100*w);

        MH::SetBinning(fSignal, binsx, binsy);

        const UShort_t roi = header->GetNumSamples();

        // Correct for DRS timing!!
        MBinning binst(roi, -0.5, roi-0.5);
        MH::SetBinning(fTime, binsx, binst);

        MBinning binspy(2*roi, -roi-0.5, roi-0.5);
        MH::SetBinning(fPulse, binsx, binspy);

        return kTRUE;
    }

    // Fill singles into histograms
    Int_t Fill(const MParContainer *, const Stat_t)
    {
        // Get pointer to samples to fill samples
        const Float_t *ptr = fEvent->GetSamples(0);
        const Short_t  roi = fEvent->GetNumSamples();

        // Loop over all pixels
        for (unsigned int i=0; i<fSingles->GetNumPixels(); i++)
        {
            if ((*fBadPix)[i].IsUnsuitable())
                continue;

            // loop over all singles
            const vector<Single> &result = fSingles->GetVector(i);

            for (unsigned int j=0; j<result.size(); j++)
            {
                // Fill signal and time
                fSignal.Fill(i, result[j].fSignal);
                fTime.Fill  (i, result[j].fTime);

                if (!ptr)
                    continue;

                // Fill pulse shape
                const Short_t tm = floor(result[j].fTime);

                for (int s=0; s<roi; s++)
                    fPulse.Fill(i, s-tm, ptr[s]);
            }

             ptr+=roi;
        }

        fNumEntries++;

        return kTRUE;
    }

    // Getter for histograms
    const TH2 *GetSignal() const { return &fSignal; }
    const TH2 *GetTime() const   { return &fTime; }
    const TH2 *GetPulse() const  { return &fPulse; }

    UInt_t GetNumEntries() const { return fNumEntries; }

    void Draw(Option_t *)
    {
        TVirtualPad *pad = gPad ? gPad : MakeDefCanvas(this);

        AppendPad("");

        pad->Divide(2,2);

        pad->cd(1);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fSignal.Draw("colz");

        pad->cd(2);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fTime.Draw("colz");

        pad->cd(3);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fPulse.Draw("colz");
    }

    void DrawCopy(Option_t *)
    {
        TVirtualPad *pad = gPad ? gPad : MakeDefCanvas(this);

        AppendPad("");

        pad->Divide(2,2);

        pad->cd(1);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fSignal.DrawCopy("colz");

        pad->cd(2);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fTime.DrawCopy("colz");

        pad->cd(3);
        gPad->SetBorderMode(0);
        gPad->SetFrameBorderMode(0);
        fPulse.DrawCopy("colz");
    }

    ClassDef(MHSingles, 1)
};

// Task to extract the singles
class MExtractSingles : public MTask
{
    MSingles *fSingles;
    MPedestalCam *fPedestal;
    MPedestalSubtractedEvt *fEvent;

    // On time for each pixel in samples
    MArrayI fExtractionRange;

    // Number of samples for sliding average
    size_t fNumAverage;

    // The range in which the singles have to fit in
    // FIXME: Make sure this is consistent with MHBaseline
    UInt_t fStartSkip;
    UInt_t fEndSkip;

    UInt_t fIntegrationSize;
    UInt_t fMaxSearchWindow;

    Int_t    fMaxDist;
    Double_t fThreshold;

public:
    MExtractSingles() : fSingles(0), fPedestal(0), fEvent(0),
        fExtractionRange(1440),
        fNumAverage(10), fStartSkip(20), fEndSkip(10),
        fIntegrationSize(3*10), fMaxSearchWindow(30)
    {
    }

    void SetMaxDist(Int_t m) { fMaxDist = m; }
    void SetThreshold(Int_t th) { fThreshold = th; }

    UInt_t GetIntegrationSize() const { return fIntegrationSize; }

    const MArrayI &GetExtractionRange() const { return fExtractionRange; }


    Int_t PreProcess(MParList *plist)
    {
        fSingles = (MSingles*)plist->FindCreateObj("MSingles");
        if (!fSingles)
            return kFALSE;

        fEvent = (MPedestalSubtractedEvt*)plist->FindObject("MPedestalSubtractedEvt");
        if (!fEvent)
        {
            *fLog << err << "MPedestalSubtractedEvt not found... abort." << endl;
            return kFALSE;
        }

        fPedestal = (MPedestalCam*)plist->FindObject("MPedestalCam");
        if (!fPedestal)
        {
            *fLog << err << "MPedestalCam not found... abort." << endl;
            return kFALSE;
        }

        return kTRUE;
    }

    Int_t Process()
    {
        const UInt_t roi = fEvent->GetNumSamples();

        const size_t navg              = fNumAverage;
        const float  threshold         = fThreshold;
        const UInt_t start_skip        = fStartSkip;
        const UInt_t end_skip          = fEndSkip;
        const UInt_t integration_size  = fIntegrationSize;
        const UInt_t max_search_window = fMaxSearchWindow;
        const UInt_t max_dist          = fMaxDist;

        vector<float> val(roi-navg);//result of first sliding average
        for (int pix=0; pix<fEvent->GetNumPixels(); pix++)
        {
            // Total number of samples as 'on time'
            fExtractionRange[pix] += roi-start_skip-navg-end_skip-integration_size;

            // Clear singles for this pixel
            vector<Single> &result = fSingles->GetVector(pix);
            result.clear();

            // Get pointer to samples
            const Float_t *ptr = fEvent->GetSamples(pix);

            // Get Baseline
            const Double_t ped = (*fPedestal)[pix].GetPedestal();

            // Subtract baseline and do a sliding average over
            // the samples to remove noise (mainly the 200MHz noise)
            for (size_t i=0; i<roi-navg; i++)
            {
                val[i] = 0;
                for (size_t j=i; j<i+navg; j++)
                    val[i] += ptr[j];
                val[i] /= navg;
                val[i] -= ped;
            }

            // peak finding via threshold
            UInt_t imax = roi-navg-end_skip-integration_size;
            for (UInt_t i=start_skip; i<imax; i++)
            {
                //search for threshold crossings
                if (val[i+0]>threshold ||
                    val[i+4]>threshold ||

                    val[i+5]<threshold ||
                    val[i+9]<threshold)
                    continue;

                //search for maximum after threshold crossing
                UInt_t k_max = i+5;
                for (UInt_t k=i+6; k<i+max_search_window; k++)
                {
                    if (val[k] > val[k_max])
                        k_max = k;
                }

                // Check if the findings make sense
                if (k_max == i+5 || k_max == i + max_search_window-1)
                    continue;

                //search for half maximum before maximum
                UInt_t k_half_max = 0;
                for (UInt_t k=k_max; k>k_max-25; k--)
                {
                    if (val[k-1] < val[k_max]/2 &&
                        val[k] >= val[k_max]/2 )
                    {
                        k_half_max = k;
                        break;
                    }
                }
                // Check if the finding makes sense
                if (k_half_max+integration_size > roi-navg-end_skip)
                    continue;
                if (k_half_max == 0)
                    continue;
                if (k_max - k_half_max > max_dist)
                    continue;

                // FIXME: Evaluate arrival time more precisely!!!
                // FIXME: Get a better integration window

                // Compile "single"
                Single single;
                single.fSignal = 0;
                single.fTime   = k_half_max;

                // Crossing of the threshold is at 5
                // Due to the averaging whihc is used to dtermine the
                // pulse, the k_half_max is roughly where the rising
                // edge starts, by adding 2.5ns we get to half the leading edge.
                for (UInt_t j=k_half_max+5; j<k_half_max+integration_size+5; j++)
                    single.fSignal += ptr[j];

                single.fSignal -= integration_size*ped;

                // Add single to list
                result.push_back(single);

                // skip falling edge
                i += 10+integration_size;

                // Remove skipped samples from 'on time'
                fExtractionRange[pix] -= i>imax ? 10+integration_size-(i-imax) : 10+integration_size;
            }
        }
        return kTRUE;
    }
};

int extract_singles_mc(
                    const char  *rundir    = "/fact/raw/2012/07/23/",
                    const char  *filesuffix = "Events.fits.gz",
                    int          firstRunID =  6,
                    int          lastRunID  =  7,
                    const char  *drsfile    = "/fact/raw/2012/07/23/20120723_004.drs.fits.gz",
                    const char  *outfile    = ".",
                    Int_t        max_dist   = 14,
                    Double_t     threshold  =  5
                   )
{
    // ======================================================

    MDirIter iter;

    gLog << "Output file:" << outfile << endl;
    if (!gSystem->AccessPathName(outfile)){
        gLog << "Output file exists allready" << endl;
        return 0;
      }

    // add input files to directory iterator
    for (Int_t fileNr = firstRunID; fileNr <= lastRunID; fileNr++)
    {
        TString filename = rundir;
        filename += Form("%08d", fileNr);
        filename += filesuffix;
        iter.AddDirectory(gSystem->DirName(rundir), gSystem->BaseName(filename+".fz"));
        iter.AddDirectory(gSystem->DirName(filename), gSystem->BaseName(filename+".gz"));
        iter.AddDirectory(gSystem->DirName(filename), gSystem->BaseName(filename));
    }
    // ======================================================

    // true:  Display correctly mapped pixels in the camera displays
    //        but the value-vs-index plot is in software/spiral indices
    // false: Display pixels in hardware/linear indices,
    //        but the order is the camera display is distorted
    bool usemap = true;

    // map file to use (get that from La Palma!)
    const char *map = usemap ? "FACTmapV5a.txt" : NULL;

    // ------------------------------------------------------

    Long_t max0 = 0;
    Long_t max1 = 0;

    // ======================================================

    if (map && gSystem->AccessPathName(map, kFileExists))
    {
        gLog << "ERROR - Cannot access mapping file '" << map << "'" << endl;
        return 1;
    }
    if (gSystem->AccessPathName(drsfile, kFileExists))
    {
        gLog << "ERROR - Cannot access drsfile file '" << drsfile << "'" << endl;
        return 2;
    }

    // --------------------------------------------------------------------------------

    gLog.Separator("Singles");
    gLog << "Extract singles with '" << drsfile << "'" << endl;
    gLog << endl;

    // ------------------------------------------------------

    MDrsCalibration drscalib300;
    if (!drscalib300.ReadFits(drsfile))
        return 3;

    // ------------------------------------------------------

    iter.Sort();
    iter.Print();

    TString title;
    title =  iter.Next();
    title += "; ";
    title += max1;
    iter.Reset();
    gLog << "Title: " << title << endl;

    // ======================================================

    MStatusDisplay *d = new MStatusDisplay;

    MBadPixelsCam badpixels;
    badpixels.InitSize(1440);
    badpixels[ 424].SetUnsuitable(MBadPixelsPix::kUnsuitable);
    badpixels[ 583].SetUnsuitable(MBadPixelsPix::kUnsuitable);
    badpixels[ 830].SetUnsuitable(MBadPixelsPix::kUnsuitable);
    badpixels[ 923].SetUnsuitable(MBadPixelsPix::kUnsuitable);
    badpixels[1208].SetUnsuitable(MBadPixelsPix::kUnsuitable);
    badpixels[1399].SetUnsuitable(MBadPixelsPix::kUnsuitable);

    //  Twin pixel
    //     113
    //     115
    //     354
    //     423
    //    1195
    //    1393

    MH::SetPalette();

    MContinue cont("(MRawEvtHeader.GetTriggerID&0xff00)!=0x100", "SelectCal");

    MGeomApply apply;

    MDrsCalibApply drsapply;

    MRawFitsRead read;
    read.LoadMap(map);
    read.AddFiles(iter);

    // ======================================================

    gLog << endl;
    gLog.Separator("Calculating baseline");

    MTaskList tlist0;

    MRawRunHeader header;
    MPedestalCam  pedcam;

    MParList plist;
    plist.AddToList(&tlist0);
    plist.AddToList(&drscalib300);
    plist.AddToList(&header);
    plist.AddToList(&pedcam);

    // ------------------ Setup the tasks ---------------

    MFillH fill0("MHBaseline", "MPedestalSubtractedEvt");

    drsapply.SetRemoveSpikes(4);

    // ------------------ Setup eventloop and run analysis ---------------

    tlist0.AddToList(&read);
    tlist0.AddToList(&apply);
    tlist0.AddToList(&drsapply);
    tlist0.AddToList(&fill0);

    // ------------------ Setup and run eventloop ---------------

    MEvtLoop loop0(title);
    loop0.SetDisplay(d);
    loop0.SetParList(&plist);

    if (!loop0.Eventloop(max0))
        return 4;

    if (!loop0.GetDisplay())
        return 5;

    // ----------------------------------------------------------------

    MGeomCamFACT fact;
    MHCamera hped(fact);
    hped.SetCamContent(pedcam);
    hped.SetCamError(pedcam, 4);
    hped.SetAllUsed();
    MHCamEvent display;
    display.SetHist(hped);

    d->AddTab("Baseline");
    display.Clone()->Draw();

    // ======================================================

    gLog << endl;
    gLog.Separator("Extracting singles");

    MTaskList tlist1;

    MSingles singles;

    plist.Replace(&tlist1);
    plist.AddToList(&badpixels);
    plist.AddToList(&singles);

    // ------------------ Setup the tasks ---------------

    MExtractSingles extract;
    extract.SetMaxDist(max_dist);
    extract.SetThreshold(threshold);

    MFillH fill("MHSingles");

    // ------------------ Setup eventloop and run analysis ---------------

    tlist1.AddToList(&read);
    tlist1.AddToList(&apply);
    tlist1.AddToList(&drsapply);
    tlist1.AddToList(&extract);
    tlist1.AddToList(&fill);

    // ------------------ Setup and run eventloop ---------------

    MEvtLoop loop1(title);
    loop1.SetDisplay(d);
    loop1.SetParList(&plist);

    if (!loop1.Eventloop(max1))
        return 6;

    if (!loop1.GetDisplay())
        return 7;

    if (fill.GetNumExecutions()==0)
        return 8;

    // =============================================================

    gStyle->SetOptFit(kTRUE);

    MHSingles* h = (MHSingles*) plist.FindObject("MHSingles");
    if (!h)
        return 9;

    const TH2 *htime   = h->GetTime();
    const TH2 *hpulse  = h->GetPulse();

    d->AddTab("Time");
    TH1 *ht = htime->ProjectionY();
    ht->SetYTitle("Counts [a.u.]");
    ht->Scale(1./1440);
    ht->Draw();

    d->AddTab("Pulse");
    TH1 *hp = hpulse->ProjectionY();
    hp->SetYTitle("Counts [a.u.]");
    hp->Scale(1./1440);
    hp->Draw();

    d->SaveAs(outfile);

    TFile f(outfile, "UPDATE");

    MParameterI par("NumEvents");
    par.SetVal(fill.GetNumExecutions());
    par.Write();

    MParameterI win("IntegrationWindow");
    win.SetVal(extract.GetIntegrationSize());
    win.Write();

    extract.GetExtractionRange().Write("ExtractionRange");

    if (firstRunID==lastRunID)
        header.Write();

    return 0;
}
