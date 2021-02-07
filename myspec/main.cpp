#include <iostream>
#include <fstream>
#include <png.h>
#include <vector>
#include <filesystem>
#include <thread>
#include <unistd.h>
#include <Accelerate/Accelerate.h>
using namespace std;
typedef long long ll;
namespace fs=std::__fs::filesystem;
class waveAu {
public:
    int samplerate,datasize;
    short channels;
    short *rawData16;
    char *rawData8;
    bool is8;
    bool hasData;
    waveAu(bool cis8,int csamplerate,int cdatasize,int cchannels,bool initWithData=true) {
        is8=cis8;
        samplerate=csamplerate;
        datasize=cdatasize;
        channels=cchannels;
        if (!initWithData) return;
        if (is8) {
            rawData8=new char[datasize];
        } else {
            rawData16=new short[datasize];
        }
        hasData=true;
    }
    waveAu(string file) {
        hasData=false;
        ifstream in(file);
        if (!in.good()) {cerr<<"File not found error!"<<endl;return;}
        
        char verify[4]; //will hold RIFF, WAVE, and fmt to verify that this is indeed a wave file.
        in.read(verify,4);
        
        if (strncmp(verify,"RIFF",4)) {cerr<<"File format error! (RIFF)"<<endl;return;}
        int fsize;
        //next up: the size
        in.read((char*)&fsize,4); // evil bit hack lol
        
        //MARK: Verifications
        in.read(verify,4);
        if (strncmp(verify,"WAVE",4)) {cerr<<"File format error! (VERFWAVE)"<<endl;return;}
//        for (ll i=0;i<610;i++) {
//            in.read(verify,1);
//        }
        in.read(verify,4);
        
        if (strncmp(verify,"fmt ",4)) {cerr<<"File format error! (VERFfmt)"<<endl;return;}
        
        int formatChunkSize;
        in.read((char*)&formatChunkSize,4); // evil bit hack lol
        if (formatChunkSize!=16) {cerr<<"Format not supported :/"<<endl;return;}
        
        short formatTag,blockAlign,bitspsample; //formatTag: way data is stored, channels: number of channels (1 for mono, 2 for stereo), block alignment: bits/sample/8*channels  bits per sample: 8 bit or 16 bit sound
        int avgbytespsec;
        //number of samples per second and how many bytes per second
        in.read((char*)&formatTag,2);in.read((char*)&channels,2);in.read((char*)&samplerate,4);in.read((char*)&avgbytespsec,4);in.read((char*)&blockAlign,2);in.read((char*)&bitspsample,2);
        if (channels!=1&&channels!=2) {cerr<<"Channel count not supported :/"<<endl;return;}
//        for (ll i=0;i<4052;i++) {
//            in.read(verify,1);
//        }
        in.read(verify,4);
        if (strncmp(verify,"data",4)) {cerr<<"File format error! (VERFDATA)"<<endl;return;}
        in.read((char*)&datasize,4);
        
        if (bitspsample==8) {
            is8=true;
            rawData8=new char[datasize];
            in.read((char*)rawData8,datasize);
        } else {
            is8=false;
            rawData16=new short[datasize/2];
            in.read((char*)rawData16,datasize);
            datasize=datasize/2;
        }
        hasData=true;
    }
    void destroy() {
        if (!hasData) return;
        if (is8) delete[] rawData8;
        else delete[] rawData16;
        hasData=false;
    }
    void toMono() {
        if (hasData&&channels==2) {
            if (is8) {
                char *monoDt=new char[datasize/2];
                for (ll i=0;i<datasize/2;i++) monoDt[i]=((int)rawData8[2*i]+rawData8[2*i+1])/2;
                delete[] rawData8;
                rawData8=monoDt;
            } else {
                short *monoDt=new short[datasize/2];
                for (ll i=0;i<datasize/2;i++) monoDt[i]=((int)rawData16[2*i]+rawData16[2*i+1])/2;
                delete[] rawData16;
                rawData16=monoDt;
            }
            channels=1;
            datasize=datasize/2;
        }
    }
    void bitrate48t16() {
        ll finDtSz=datasize/3;
        short *newDt=new short[finDtSz];
        for (ll i=0;i<finDtSz;i++) {
            newDt[i]=rawData16[i*3];
        }
        delete[] rawData16;
        rawData16=newDt;
        datasize=finDtSz;
    }
    double *toDouble() {
        double* rturn=new double[datasize];
        for (ll i=0;i<datasize;i++) {
            if (is8) rturn[i]=rawData8[i];
            else rturn[i]=rawData16[i];
        }
        return rturn;
    }
};

double const dftWindow=0.01533898305;
double const frameInc=0.003813559322;
double const chunkLen=1.0;
double const strideLen=0.1;
ll const resizeX=256;
ll const resizeY=256;
struct rgb {
    ll r,g,b;
};
rgb jetMap[256];
bool writeOutImg(png_bytepp pngdt,string file,ll height,ll width) {
    //this is mostly copied and formatted from someone else's code lmao
    png_structp png_ptr=png_create_write_struct(PNG_LIBPNG_VER_STRING,NULL,NULL,NULL);
    if (png_ptr==NULL) { return 0; }
    png_infop info_ptr=png_create_info_struct(png_ptr);
    if (info_ptr==NULL) {png_destroy_write_struct(&png_ptr,(png_infopp)NULL); return 0;};
    
    if (setjmp(png_jmpbuf(png_ptr))) {png_destroy_write_struct(&png_ptr,&info_ptr); return 0;}
    
    FILE *ofs=fopen(file.c_str(),"wb");
    if (ofs==NULL) {
        return 0;
    }
    png_init_io(png_ptr,ofs);
    
    png_set_IHDR(png_ptr, info_ptr,width,height,8,PNG_COLOR_TYPE_RGB,PNG_INTERLACE_NONE,PNG_COMPRESSION_TYPE_DEFAULT,PNG_FILTER_TYPE_DEFAULT);
    
    png_write_info(png_ptr,info_ptr);
    
    png_write_image(png_ptr,pngdt);
    
    png_write_end(png_ptr,NULL);
    
    png_destroy_write_struct(&png_ptr,&info_ptr);
    
    fclose(ofs);
    return true;
}

//Image interpolation algorithm
double sinc(double x) {
    x*=M_PI;
    if (x<0.01f && x>-0.01f) return 1.0+x*x*(-1.0/6.0+x*x*1.0/120.0);
    return sin(x)/x;
}
double clip(double t) {
    const double eps = 0.0000125;
    if (fabs(t)<eps) return 0.0;
    return t;
}
double lancos(double t) {
    if (t<0.0f) t=-t;

    if (t<3.0f) return clip(sinc(t)*sinc(t/3.0));
    else return (0.0f);
}
//Interpolation in the x direction:
static inline float lancos3_resample_x(double** arr,ll height,ll width,ll x,ll y,double xscale) { //recoded from programmersought.com/article/9978929119/
    double scaleReach;
    if (xscale>1.0) scaleReach=3.0; //if it's downscaling, then at least use 3 blocks for the kernel!
    else scaleReach=3.0/xscale;

    double xcenter=(double)x/xscale; //corresponding xcenter in the source image

    if (y<0) y=0; //clip within bounds
    if (y>=height) y=height-1;
    if (xscale>1.0f) xscale=1.0f;
    
    double rturn=0;
    double coef_sum=0.0;
    for (ll i=floor(xcenter-scaleReach);i<=ceil(xcenter+scaleReach);i++) { //move Lancos filter across
        x=i;
        if (i<0) x=0;
        if (i>=width) x=width-1; //clip within bounds!
        double pix=arr[y][x];
        double coef=lancos((xcenter-i)*xscale);
        rturn+=pix*coef;
        coef_sum+=coef;
    }
    rturn/=coef_sum;
    return rturn;
}
void lancos3Resample(double** src,double** dst,ll src_rows,ll src_cols,ll dst_rows,ll dst_cols) {
    double xratio = dst_cols/(double)src_cols;
    double yratio = dst_rows/(double)src_rows;
    
    double scaleReach;
    double scale=0.0;
    if (yratio>1.0) {
        scaleReach=3.0;
        scale=1.0;
    } else {
        scaleReach=3.0/yratio;
        scale=yratio;
    }
    
    for (ll i=0;i<dst_rows;i++) {
        for (ll j=0;j<dst_cols;j++) {
            double summ=0;
            double coef_sum=0.0;
            
            double c = (double)i/yratio;
            // Interpolate in the x direction first, then interpolate in the y direction.
            for (ll k=floor(c-scaleReach);k<=ceil(c+scaleReach);k++) {
                double pix = lancos3_resample_x(src, src_rows,src_cols,j,k,xratio);
                double coef = lancos((c - k)*scale);
                coef_sum += coef;
                pix *= coef;
                summ += pix;
            }
            double val=summ/coef_sum;
            dst[i][j]=val;
        }
    }
}
inline ll next2Pow(ll x) {
    //get the next power of 2
    if ((x&(-x))==x) return x;
    for (ll i=1;i<=32;(i<<=1)) { //bit haccs
        x|=(x>>i);
    }
    return x+1;
}
fs::path analPth="/Users/legitmichel777/Desktop/Orate/Datasets/Retired/sortedRAVDESS";
fs::path toDir="/Users/legitmichel777/Desktop/Orate/Datasets/superPlan-RAVDESSdd";
ll thrs=8;
ll binL=-230,binR=0;
ll bins=10000;
ll dftWindowSample,interpolationSample,frameIncSample,coefsNum,chunkSample,strideSample;
ll totBinned;
double destroyPercentileL=1;
double destroyPercentileR=1;
mutex mtx;
void processAud(ll workerID,ll *job,waveAu* masterSrc,ll *masterBin,double *frame,double *hammingWindow,FFTSetup fftSetup,double *outCoefs,double scale,bool output,string outputPath,double lfreq,double rfreq) { //transform the audio
    ll log2n=log2(interpolationSample);
    ll *myBin=new ll[bins];
    for (ll i=0;i<bins;i++) myBin[i]=0;
    //do cut
    ll resNum=0;
    ll numWhat=0;
    for (ll startSample=0;;startSample+=strideSample) {
        chrono::steady_clock::time_point grs=chrono::steady_clock::now();
        zeroD:;
        if (startSample+chunkSample>masterSrc->datasize) break;
        //ll endSample=startSample+chunkSample;
        //from startSample..<endSample
        waveAu wavSrc=waveAu(masterSrc->is8,masterSrc->samplerate,chunkSample,masterSrc->channels,false); //inherit from master
        //copy from masterSrc
        wavSrc.rawData16=&masterSrc->rawData16[startSample];
        wavSrc.hasData=true;
        
        //signal process
        vector<double*>spectr;
        for (ll i=dftWindowSample/2;1;i+=frameIncSample) { //if dftWindowSample is odd then take the lower bound!
            double *coefs=new double[coefsNum];
            ll leftEdg=i-dftWindowSample/2;
            ll rightEdg=(frameIncSample%2?i+dftWindowSample/2:i+dftWindowSample/2+1);
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
//        cout<<chrono::duration_cast<chrono::milliseconds>(chrono::steady_clock::now() - grs).count()<<" ms"<<endl;
        
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
            writeOutImg(pngdt,outputPath+"_"+to_string(numWhat)+".png",resizeX,resizeY);
            for (ll i=0;i<resizeX;i++) delete[] dst[i];
            delete[] dst;
            for (ll i=0;i<resizeX;i++) delete[] pngdt[i];
            delete[] pngdt;
        } else {
            for (ll i=0;i<spectr.size();i++) delete[] spectr[i];
            delete[] src;
        }
    }
    mtx.lock();
    for (ll i=0;i<bins;i++) {
        masterBin[i]+=myBin[i];
        totBinned+=myBin[i];
    }
    mtx.unlock();
    delete[] myBin;
    delete[] masterSrc->rawData16;
    delete masterSrc;
    cout<<job[workerID]<<" execution complete"<<endl;
    job[workerID]=-1;
    //important: deallocate everything
}
int main() {
    ifstream in("/Users/legitmichel777/Desktop/Orate/Code/jet.csv");
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
    
    if (fs::exists(toDir)) {
        fs::remove_all(toDir);
    }
    fs::create_directory(toDir);
    
    bool hasInit=false;
    fs::recursive_directory_iterator itr(analPth);
    dftWindowSample=interpolationSample=frameIncSample=coefsNum=chunkSample=strideSample=-1; //if you get -1 you know you fucked up
    FFTSetup fourierTransform[thrs];
    double *frame[thrs];
    double *outCoefs[thrs];
    double *hammingWindow;
    double scale;
    scale=-1;
    struct analObj {
        fs::path fpth,emo;
    };
    vector<analObj>toAnal;
    for (const fs::directory_entry& files:itr) {
        if (files.path().filename().extension()==".wav") {
            string s=files.path().filename().string();
            if (s.size()>4) {
                if (s.substr(0,4)=="conv") {
                    toAnal.push_back((analObj){files.path(),files.path().parent_path().filename()});
                }
            }
        }
    }
    cout<<"Found "<<toAnal.size()<<" files for analysis."<<endl;
    thread thr[thrs];
    ll* job=new ll[thrs];
    ll *daBins=new ll[bins];
    for (ll i=0;i<bins;i++) daBins[i]=0;
    for (ll i=0;i<thrs;i++) job[i]=-1;
    cout<<"Running Fourier Transform..."<<endl;
    for (ll times=0;times<toAnal.size();times++) {
        cout<<"Job "<<times<<" "<<toAnal[times].fpth<<endl;
        waveAu* daWav=new waveAu(toAnal[times].fpth);
        if (!daWav->hasData) {
            cout<<"Read error at "<<toAnal[times].fpth<<", skipping..."<<endl;
            continue;
        }
        daWav->toMono();
        daWav->bitrate48t16();
        if (!hasInit) {
            cout<<"Initializing based on file "<<toAnal[times].fpth<<endl;
            if (daWav->is8) {
                cout<<"No 8-bit!"<<endl;
                return 0;
            }
            dftWindowSample=dftWindow*daWav->samplerate;
            interpolationSample=next2Pow(4*dftWindowSample); //interpolate!
            frameIncSample=frameInc*daWav->samplerate;
            chunkSample=chunkLen*daWav->samplerate;
            strideSample=strideLen*daWav->samplerate;
            cout<<"System configuration\nExecution threads:"<<thrs<<"\n\nData configuration\nChunk samples:"<<chunkSample<<"\nStride samples:"<<strideSample<<"\nOutput image x:"<<resizeX<<"\nOutput image y:"<<resizeY<<"\nBins:"<<binL<<"..."<<binR<<" ("<<bins<<" bins)\nBin Percentiles(L R):"<<destroyPercentileL<<" "<<destroyPercentileR<<"\n\nFeature configuration\nDiscrete Fourier Transform samples:"<<dftWindowSample<<"\nInterpolated (zero padding) Discrete Fourier Transform window:"<<interpolationSample<<"\nFrame increment:"<<frameIncSample<<"\n";
            
            //initialize everything
            for (ll i=0;i<thrs;i++) frame[i]=new double[interpolationSample];
            for (ll i=0;i<thrs;i++) outCoefs[i]=new double[interpolationSample];
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
        
        thr[freeWorker]=thread(processAud,freeWorker,job,daWav,daBins,frame[freeWorker],hammingWindow,fourierTransform[freeWorker],outCoefs[freeWorker],scale,false,"",0,0);
    }
    //merge with the threads
    for (ll i=0;i<thrs;i++) if (thr[i].joinable()) thr[i].join();
    
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
        cout<<"Job "<<times<<" "<<toAnal[times].fpth<<endl;
        waveAu* daWav=new waveAu(toAnal[times].fpth);
        if (!daWav->hasData) {
            cout<<"Read error at "<<toAnal[times].fpth<<", skipping..."<<endl;
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
        
        fs::path goTo=toDir/toAnal[times].emo;
        if (!fs::exists(goTo)) fs::create_directory(goTo);
        goTo/=toAnal[times].fpth.filename();
        
        thr[freeWorker]=thread(processAud,freeWorker,job,daWav,daBins,frame[freeWorker],hammingWindow,fourierTransform[freeWorker],outCoefs[freeWorker],scale,true,goTo.string(),lfreq,rfreq);
    }
    
    for (ll i=0;i<thrs;i++) if (thr[i].joinable()) thr[i].join();
    
    for (ll i=0;i<thrs;i++) {
        delete[] frame[i];
        delete[] outCoefs[i];
        vDSP_destroy_fftsetup(fourierTransform[i]);
    }
    delete[] hammingWindow;
    
    delete[] daBins;
    delete[] job;
}

