//
//  waveAu.hpp
//  waveAu
//
//  Created by LegitMichel777 on 2021/8/8.
//

#ifndef waveAu_h
#define waveAu_h

#include <iostream>
using namespace std;
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
        
        if (strncmp(verify,"RIFF",4)) {
            cerr<<"File format error! (RIFF)"<<endl;
            in.close();
            return;
        }
        int fsize;
        //next up: the size
        in.read((char*)&fsize,4); // evil bit hack lol
        
        //MARK: Verifications
        in.read(verify,4);
        if (strncmp(verify,"WAVE",4)) {
            cerr<<"File format error! (VERFWAVE)"<<endl;
            in.close();
            return;
        }
        
        in.read(verify,4);
        
        if (strncmp(verify,"fmt ",4)) {
            cerr<<"File format error! (VERFfmt)"<<endl;
            in.close();
            return;
        }
        
        int formatChunkSize;
        in.read((char*)&formatChunkSize,4); // evil bit hack lol
        if (formatChunkSize!=16) {
            cerr<<"Format not supported :/"<<endl;
            in.close();
            return;
        }
        
        short formatTag,blockAlign,bitspsample; //formatTag: way data is stored, channels: number of channels (1 for mono, 2 for stereo), block alignment: bits/sample/8*channels  bits per sample: 8 bit or 16 bit sound
        int avgbytespsec;
        //number of samples per second and how many bytes per second
        in.read((char*)&formatTag,2);in.read((char*)&channels,2);in.read((char*)&samplerate,4);in.read((char*)&avgbytespsec,4);in.read((char*)&blockAlign,2);in.read((char*)&bitspsample,2);
        if (channels!=1&&channels!=2) {
            cerr<<"Channel count not supported :/"<<endl;
            in.close();
            return;
        }
//        for (ll i=0;i<4052;i++) {
//            in.read(verify,1);
//        }
        in.read(verify,4);
        if (strncmp(verify,"data",4)) {
            cout<<"File format error! (VERFDATA) Expected data, got "<<verify[0]<<verify[1]<<verify[2]<<verify[3]<<endl;
            in.close();
            return;
        }
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
        in.close();
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

#endif /* waveAu_h */
