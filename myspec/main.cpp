#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <thread>
#include <unistd.h>
#include <Accelerate/Accelerate.h>
#include "writeOutImage.hpp"
#include "lancos.hpp"
#include "waveAu.hpp"
#include "utils.hpp"
using namespace std;
typedef long long ll;
namespace fs=std::__fs::filesystem;
rgb jetMap[256];
mutex mtx;

void processAud(ll workerID,ll *job,waveAu* masterSrc, ll *totBinned, ll *masterBin, ll bins, ll binL, ll binR, double *frame,double *hammingWindow,FFTSetup fftSetup,double scale,bool output,string outputPath, ll resizeX, ll resizeY, double lfreq,double rfreq, ll dftWindowSample,ll interpolationSample,ll frameIncSample,ll coefsNum,ll chunkSample,ll strideSample) { //transform the audio
    ll log2n=log2(interpolationSample);
    ll *myBin;
    if (bins != -1) {
        myBin=new ll[bins];
        for (ll i=0;i<bins;i++) myBin[i]=0;
    }
    //do cut
    ll numWhat=0;
    
    for (ll startSample=0;;startSample+=strideSample) {
        zeroD:;
        if (startSample+chunkSample>masterSrc->datasize) break;
        //ll endSample=startSample+chunkSample;
        //from startSample..<endSample
        waveAu wavSrc=waveAu(masterSrc->is8,masterSrc->samplerate,chunkSample,masterSrc->channels,false); //inherit from master
        //copy from masterSrc
        wavSrc.rawData16=&masterSrc->rawData16[startSample];
        wavSrc.hasData=true;
        
        vector<double*>spectr;
        for (ll i=dftWindowSample/2;1;i+=frameIncSample) { //if dftWindowSample is odd then take the lower bound!
            double *coefs=new double[coefsNum];
            ll leftEdg=i-dftWindowSample/2;
            ll rightEdg=leftEdg+dftWindowSample-1;
            if (rightEdg>=wavSrc.datasize) {
                delete[] coefs;
                break;
            }
            
            bool isZeroTrap=true;
            for (ll j=leftEdg;j<=rightEdg&&isZeroTrap;j++) if (wavSrc.rawData16[j]!=0) isZeroTrap=false;
            if (isZeroTrap) {
                startSample+=strideSample;
                for (ll i=0;i<spectr.size();i++) delete[] spectr[i];
                delete[] coefs;
                goto zeroD;
            }
            
            
            //leftEdg...rightEdg (inclusive)
            for (ll j=leftEdg;j<=rightEdg;j++) {
                frame[j-leftEdg]=wavSrc.rawData16[j]/32768.0;
                //hamming window
                frame[j-leftEdg]*=hammingWindow[j-leftEdg];
            }
            for (ll j=dftWindowSample;j<interpolationSample;j++) frame[j]=0;
            //execute fourier transform
                    
            float *input=new float[interpolationSample];
            float *output=new float[interpolationSample];
            for (ll i=0;i<interpolationSample;i++) input[i]=frame[i];
            COMPLEX_SPLIT A;
            A.realp=new float[interpolationSample/2];
            A.imagp=new float[interpolationSample/2];
            
            vDSP_ctoz((COMPLEX*)input,2,&A,1,interpolationSample/2);
            vDSP_fft_zrip(fftSetup,&A,1,log2n,FFT_FORWARD);
            
            coefs[0]=A.realp[0]*A.realp[0];
            for (ll j=1;j<coefsNum;j++) {
                coefs[j]=A.realp[j]*A.realp[j]+A.imagp[j]*A.imagp[j];
                coefs[j]/=4;
            }
            
            for (ll j=0;j<coefsNum;j++) {
                coefs[j]=10*log10(coefs[j]/scale);
            }
            spectr.push_back(coefs);
            
            delete[] input;
            delete[] output;
            delete[] A.realp;
            delete[] A.imagp;
        }
        
        double **src=new double*[spectr.size()];
        for (ll i=0;i<spectr.size();i++) src[i]=spectr[i];
        
        if (bins != -1) {
            for (ll i=0;i<spectr.size();i++) {
                for (ll j=0;j<coefsNum;j++) {
                    if (src[i][j]<binL||src[i][j]>binR) {
                        cout<<"Out of bin "<<i<<" "<<j<<" "<<src[i][j]<<endl;
                    } else {
                        ll daBin=(src[i][j]-binL)/(binR-binL)*bins;
                        daBin=max(daBin,0ll);
                        daBin=min(daBin,bins-1);
                        myBin[daBin]++;
                    }
                }
            }
            mtx.lock();
            for (ll i=0;i<bins;i++) {
                masterBin[i]+=myBin[i];
                *totBinned+=myBin[i];
            }
            mtx.unlock();
        }
        
        if (output) {
            double **dst=new double*[resizeX];
            
            for (ll i=0;i<resizeX;i++) {
                dst[i]=new double[resizeY];
            }
//            lancos3Resample(src,dst,spectr.size(),coefsNum,resizeX,resizeY);
            for (ll i=0;i<resizeX;i++) {
                for (ll j=0;j<resizeY;j++) {
                    dst[i][j]=src[(int)((double)i/resizeX*spectr.size())][int((double)j/resizeY*coefsNum)];
//                    dst[i][j]=src[i][j];
                }
            }
            
            for (ll i=0;i<spectr.size();i++) delete[] spectr[i];
            delete[] src;
            
            png_bytepp pngdt=new png_bytep[resizeX];
            for (ll i=0;i<resizeX;i++) pngdt[i]=new png_byte[3*resizeY];
            double minRng=lfreq,maxRng=rfreq;
            for (ll i=0;i<resizeX;i++) {
                png_bytep curRow=pngdt[i];
                for (ll j=0;j<resizeY;j++) {
                    png_bytep dbt=&curRow[j*3];
                    double thePix=(dst[i][j]-minRng)/(maxRng-minRng);
                    ll fin=thePix*256;
                    fin=max(fin,0ll);
                    fin=min(fin,255ll);
                    rgb curPixel=jetMap[fin];
                    dbt[0]=curPixel.r;
                    dbt[1]=curPixel.g;
                    dbt[2]=curPixel.b;
                }
            }
            numWhat++;
            writeOutImage(pngdt,outputPath+"_"+to_string(numWhat)+".png",resizeX,resizeY);
            for (ll i=0;i<resizeX;i++) delete[] dst[i];
            delete[] dst;
            for (ll i=0;i<resizeX;i++) delete[] pngdt[i];
            delete[] pngdt;
        } else {
            for (ll i=0;i<spectr.size();i++) delete[] spectr[i];
            delete[] src;
        }
        
        //IMPORTANT: Erase this
        break;
    }
    if (bins != -1) {
        delete[] myBin;
    }
    delete[] masterSrc->rawData16;
    delete masterSrc;
//    cout<<job[workerID]<<" execution complete"<<endl;
    job[workerID]=-1;
    //important: deallocate everything
}

bool runJobsBatch(vector<analObj> toAnal, ll thrs, ll binL, ll binR, ll bins, double destroyPercentileL, double destroyPercentileR, double dftWindow, double frameInc, double chunkLen, double strideLen, ll resizeX, ll resizeY) {
    ll totBinned=0;
    
    bool hasInit=false;
    ll dftWindowSample,interpolationSample,frameIncSample,coefsNum,chunkSample,strideSample;
    dftWindowSample=interpolationSample=frameIncSample=coefsNum=chunkSample=strideSample=-1; //if you get -1 you know you fucked up
    FFTSetup fourierTransform[thrs];
    double *frame[thrs];
    double *hammingWindow;
    double scale;
    scale=-1;
    thread thr[thrs];
    ll* job=new ll[thrs];
    ll* daBins=new ll[bins];
    for (ll i=0;i<bins;i++) daBins[i]=0;
    for (ll i=0;i<thrs;i++) job[i]=-1;
    
    cout<<"Running Preprocess..."<<endl;
    for (ll times=0;times<toAnal.size();times++) {
//        cout<<"Preprocess Job "<<times<<" "<<toAnal[times].fpth<<endl;
        waveAu* daWav=new waveAu(toAnal[times].inPath);
        if (!daWav->hasData) {
            cout<<"Read error at "<<toAnal[times].inPath<<", skipping..."<<endl;
            continue;
        }
        daWav->toMono();
        daWav->bitrate48t16();
        if (!hasInit) {
            cout<<"Initializing based on file "<<toAnal[times].inPath<<endl;
            if (daWav->is8) {
                cout<<"No 8-bit!"<<endl;
                return false;
            }
            dftWindowSample=dftWindow*daWav->samplerate;
            interpolationSample=next2Pow(4*dftWindowSample); //interpolate!
            frameIncSample=frameInc*daWav->samplerate;
            chunkSample=chunkLen*daWav->samplerate;
            strideSample=strideLen*daWav->samplerate;
            cout<<"System configuration\nExecution threads:"<<thrs<<"\n\nData configuration\nChunk samples:"<<chunkSample<<"\nStride samples:"<<strideSample<<"\nOutput image x:"<<resizeX<<"\nOutput image y:"<<resizeY<<"\nBins:"<<binL<<"..."<<binR<<" ("<<bins<<" bins)\nBin Percentiles(L R):"<<destroyPercentileL<<" "<<destroyPercentileR<<"\n\nFeature configuration\nDiscrete Fourier Transform samples:"<<dftWindowSample<<"\nInterpolated (zero padding) Discrete Fourier Transform window:"<<interpolationSample<<"\nFrame increment:"<<frameIncSample<<"\n";
            
            //initialize everything
            for (ll i=0;i<thrs;i++) frame[i]=new double[interpolationSample];
            hammingWindow=new double[dftWindowSample];
            uint32_t log2n=log2(interpolationSample);
            for (ll i=0;i<thrs;i++) fourierTransform[i]=vDSP_create_fftsetup(log2n,FFT_RADIX2);
            coefsNum=floor(interpolationSample/2.0);
            
            double hammingCorrectCoefficient=0;
            for (ll i=0;i<dftWindowSample;i++) {
                hammingWindow[i]=0.53836-(1-0.53836)*cos(2.0*M_PI*i/(double)(dftWindowSample-1));
                hammingCorrectCoefficient+=hammingWindow[i]*hammingWindow[i];
            }
            scale=0.5*daWav->samplerate*hammingCorrectCoefficient;
            
            hasInit=true;
            cout<<"\nInitialization complete. Executing..."<<endl;
        }
        
        //find one
        ll freeWorker=-1;
        while (true) {
            for (ll i=0;i<thrs;i++) {
                if (job[i]==-1) {
                    freeWorker=i;
                    break;
                }
            }
            if (freeWorker!=-1) break;
            usleep(1e4);
        }
        job[freeWorker]=times;
        if (thr[freeWorker].joinable()) thr[freeWorker].join();

        thr[freeWorker]=thread(processAud,freeWorker,job,daWav, &totBinned, daBins, bins, binL, binR, frame[freeWorker],hammingWindow,fourierTransform[freeWorker],scale,false,"", -1, -1, -1, -1, dftWindowSample, interpolationSample, frameIncSample, coefsNum, chunkSample, strideSample);
    }
    //merge with the threads
    for (ll i=0;i<thrs;i++) if (thr[i].joinable()) thr[i].join();
    
    if (totBinned==0) {
        cout<<"Fatal error: no elements in bin!"<<endl;
        if (hasInit) {
            for (ll i=0;i<thrs;i++) {
                delete[] frame[i];
                vDSP_destroy_fftsetup(fourierTransform[i]);
            }
            delete[] hammingWindow;
        }

        delete[] daBins;
        delete[] job;

        return false;
    }
    //determine L and R
    ll tryBin=totBinned*destroyPercentileL/100.0;
    ll bound=0;
    ll theTot=0;
    for (;;bound++) {
        if (theTot+daBins[bound]>tryBin) break;
        theTot+=daBins[bound];
    }
    double lfreq=bound*(binR-binL)/(double)bins+binL+(tryBin-theTot)/(double)daBins[bound]*(binR-binL)/bins; //linear interpolate left frequency
    
    theTot=0;
    bound=bins-1;
    tryBin=totBinned*destroyPercentileR/100.0;;
    for (;;bound--) {
        if (theTot+daBins[bound]>tryBin) break;
        theTot+=daBins[bound];
    }
    double rfreq=bound*(binR-binL)/(double)bins+binL-(tryBin-theTot)/(double)daBins[bound]*(binR-binL)/bins; //linear interpolate left frequency
    cout<<"Frequency Range "<<lfreq<<" "<<rfreq<<endl;
    
    cout<<"Outputting images"<<endl;
    for (ll i=0;i<thrs;i++) job[i]=-1;
    for (ll times=0;times<toAnal.size();times++) {
//        if (toAnal[times].fpth.parent_path().filename()!="Malesad   ") continue;
//        cout<<"Output Job "<<times<<" "<<toAnal[times].fpth<<endl;
        waveAu* daWav=new waveAu(toAnal[times].inPath);
        if (!daWav->hasData) {
            cout<<"Read error at "<<toAnal[times].inPath<<", skipping..."<<endl;
            continue;
        }
        daWav->toMono();
        daWav->bitrate48t16();
        
        //find one
        ll freeWorker=-1;
        while (true) {
            for (ll i=0;i<thrs;i++) {
                if (job[i]==-1) {
                    freeWorker=i;
                    break;
                }
            }
            if (freeWorker!=-1) break;
            usleep(1e4);
        }
        job[freeWorker]=times;
        if (thr[freeWorker].joinable()) thr[freeWorker].join();
        
        thr[freeWorker]=thread(processAud,freeWorker,job,daWav, nullptr, nullptr, -1, -1, -1, frame[freeWorker],hammingWindow,fourierTransform[freeWorker],scale,true,toAnal[times].outPath.string(), resizeX, resizeY, lfreq, rfreq, dftWindowSample, interpolationSample, frameIncSample, coefsNum, chunkSample, strideSample);
    }
    
    for (ll i=0;i<thrs;i++) if (thr[i].joinable()) thr[i].join();
    
    if (hasInit) {
        for (ll i=0;i<thrs;i++) {
            delete[] frame[i];
            vDSP_destroy_fftsetup(fourierTransform[i]);
        }
        delete[] hammingWindow;
    }
    
    delete[] daBins;
    delete[] job;
    
    return true;
}

int main() {
    chrono::steady_clock::time_point grs=chrono::steady_clock::now();
    
    ifstream in("/Users/legitmichel777/Developer/Orate/Code/Spectrogram Generation/jet.csv");
    for (ll i=0;i<256;i++) {
        char take;
        double tmp;
        in>>tmp;
        jetMap[i].r=tmp*255;
        in>>take;
        in>>tmp;
        jetMap[i].g=tmp*255;
        in>>take;
        in>>tmp;
        jetMap[i].b=tmp*255;
    }
    
//    for (double dftWindow=0.001; dftWindow<0.1;dftWindow+=0.0001) {
//        vector<analObj>toAnal;
//        fs::path toDir="/Users/legitmichel777/Developer/Orate/Datasets/windowsize-eval/"+to_string(dftWindow)+".wav";
//        toAnal.push_back((analObj){"/Users/legitmichel777/Developer/Orate/Datasets/emoDB/08a04Tb.wav", toDir.string()});
//        runJobsBatch(toAnal, 8, -230, 0, 10000, 1, 1, dftWindow, dftWindow/4, 1.0, 0.1, 225, 225);
//    }
    
    vector<double>dftWindows;
    dftWindows.push_back(0.0071);
    dftWindows.push_back(0.01533898305);
    dftWindows.push_back(0.0254);
    
    struct csvEntry {
        string image1;
        string image2;
        string image3;
        string emotion;
    };
    
//    {
//        fs::recursive_directory_iterator itr("/Users/legitmichel777/Developer/Orate/Datasets/emoDB-sorted");
//        vector<analObj>toAnal;
//        for (const fs::directory_entry& files:itr) {
//            if (files.path().filename().extension()==".wav") {
//                string s=files.path().filename().string();
//                fs::path outPath="/Volumes/smerge/"+files.path().parent_path().filename().string();
//                if (!fs::exists(outPath)) fs::create_directory(outPath);
//                outPath/=files.path().filename();
//                toAnal.push_back((analObj){files.path(),outPath});
//            }
//        }
//    }
    
    fs::recursive_directory_iterator itr("/Users/legitmichel777/Developer/Orate/Datasets/emoDB-sorted");

    fs::path toDir="/Users/legitmichel777/Developer/Orate/Datasets/temp-gen";
    if (fs::exists(toDir)) {
        fs::remove_all(toDir);
    }
    fs::create_directory(toDir);
    vector<analObj>toAnal;
    for (const fs::directory_entry& files:itr) {
        if (files.path().filename().extension()==".wav") {
            string s=files.path().filename().string();
            fs::path outPath=toDir/files.path().parent_path().filename();
            if (!fs::exists(outPath)) fs::create_directory(outPath);
            outPath/=files.path().filename();
            toAnal.push_back((analObj){files.path(),outPath});
        }
    }
    cout<<"Found "<<toAnal.size()<<" files for analysis."<<endl;

    double dftWindow=0.01533898305;
    runJobsBatch(toAnal, 8, -230, 0, 10000, 1, 1, dftWindow, dftWindow/4, 1.0, 0.1, 225, 225);
    
    cout<<chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - grs).count()<<endl;
}

