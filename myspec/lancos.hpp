//
//  lanczos.hpp
//  lanczos
//
//  Created by LegitMichel777 on 2021/8/8.
//

#ifndef lanczos_h
#define lanczos_h

typedef long long ll;

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

#endif /* lanczos_h */
