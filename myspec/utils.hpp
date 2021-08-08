//
//  utils.hpp
//  utils
//
//  Created by LegitMichel777 on 2021/8/8.
//

#ifndef utils_h
#define utils_h

#include <filesystem>
namespace fs=std::__fs::filesystem;

struct rgb {
    ll r,g,b;
};

inline ll next2Pow(ll x) {
    //get the next power of 2
    if ((x&(-x))==x) return x;
    for (ll i=1;i<=32;(i<<=1)) { //bit haccs
        x|=(x>>i);
    }
    return x+1;
}

struct analObj {
    fs::path inPath,outPath;
};

void replaceAll(string &str, const string& from, const string& to) {
    size_t start_pos = 0;
    while((start_pos = str.find(from, start_pos)) != std::string::npos) {
        str.replace(start_pos, from.length(), to);
        start_pos += to.length(); // Handles case where 'to' is a substring of 'from'
    }
}

#endif /* utils_h */
