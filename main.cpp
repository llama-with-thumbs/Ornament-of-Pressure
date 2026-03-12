#include <SDL2/SDL.h>
#include <GL/glew.h>
#include <cmath>
#include <vector>
#include <unordered_map>
#include <string>

static const float PI    = 3.14159265358979323846f;
static const float SQRT3 = 1.7320508075688772f;
static const int   FFT_N = 1024;
static const int   NUM_ORNAMENTS = 6;
static const int   THUMB_SIZE = 64;      // thumbnail pixel size
static const int   SIDEBAR_W  = 86;      // sidebar width in pixels
static const int   THUMB_PAD  = 8;       // padding between thumbnails

// ============================================================
// Color layers
// ============================================================
#define NUM_LAYERS 7
enum Layer {
    L_STAR_BOUNDARY = 0,  // gold
    L_SPOKES,             // warm white
    L_GRID_MAIN,          // teal
    L_GRID_FINE,          // soft blue
    L_MODULES,            // purple
    L_STARS8,             // crimson
    L_ROTATING,           // deep magenta
};

struct LayerColor { float r, g, b; };
static LayerColor baseColors[NUM_LAYERS] = {
    { 1.00f, 0.82f, 0.28f },  // gold
    { 1.00f, 0.95f, 0.85f },  // warm white
    { 0.20f, 0.88f, 0.95f },  // teal
    { 0.40f, 0.65f, 1.00f },  // bright blue
    { 0.75f, 0.38f, 1.00f },  // purple
    { 1.00f, 0.25f, 0.35f },  // crimson
    { 0.95f, 0.32f, 0.78f },  // magenta
};

// ============================================================
// FFT
// ============================================================
struct Complex { float re, im; };
static void fft(Complex* x, int n) {
    for (int i=1,j=0;i<n;i++){int b=n>>1;for(;j&b;b>>=1)j^=b;j^=b;
        if(i<j){Complex t=x[i];x[i]=x[j];x[j]=t;}}
    for (int len=2;len<=n;len<<=1){
        float ang=-2.f*PI/len; Complex wl={cosf(ang),sinf(ang)};
        for(int i=0;i<n;i+=len){Complex w={1,0};
            for(int j=0;j<len/2;j++){
                Complex u=x[i+j];
                Complex v={w.re*x[i+j+len/2].re-w.im*x[i+j+len/2].im,
                           w.re*x[i+j+len/2].im+w.im*x[i+j+len/2].re};
                x[i+j]={u.re+v.re,u.im+v.im};
                x[i+j+len/2]={u.re-v.re,u.im-v.im};
                Complex nw={w.re*wl.re-w.im*wl.im,w.re*wl.im+w.im*wl.re};w=nw;}}}
}

// ============================================================
// Audio
// ============================================================
struct AudioState { float ring[FFT_N*4]; int wp; SDL_AudioDeviceID dev; };
static AudioState g_au={};
static void audioCB(void*,Uint8*s,int len){
    float*f=(float*)s;int c=len/(int)sizeof(float);
    for(int i=0;i<c;i++){g_au.ring[g_au.wp%(FFT_N*4)]=f[i];g_au.wp++;}}

// ============================================================
// Geometry
// ============================================================
struct Vec2 { float x, y; };

static void hankinPoly(const Vec2*v,int nv,float ca,std::vector<Vec2>&s){
    float cx=0,cy=0;
    for(int i=0;i<nv;i++){cx+=v[i].x;cy+=v[i].y;}cx/=nv;cy/=nv;
    struct Ray{float ox,oy,dx,dy;int e;};
    Ray rays[48];int rc=0;
    for(int e=0;e<nv;e++){
        Vec2 a=v[e],b=v[(e+1)%nv];
        float mx=(a.x+b.x)*.5f,my=(a.y+b.y)*.5f;
        float ea=atan2f(b.y-a.y,b.x-a.x),perp=ea+PI*.5f;
        float tx=mx+cosf(perp)*.01f,ty=my+sinf(perp)*.01f;
        float inw=((tx-cx)*(tx-cx)+(ty-cy)*(ty-cy)<(mx-cx)*(mx-cx)+(my-cy)*(my-cy))?perp:perp+PI;
        rays[rc++]={mx,my,cosf(inw+ca),sinf(inw+ca),e};
        rays[rc++]={mx,my,cosf(inw-ca),sinf(inw-ca),e};}
    for(int i=0;i<rc;i++){
        float bt=1e9f,bx=0,by=0;bool ok=false;
        for(int j=0;j<rc;j++){if(rays[j].e==rays[i].e)continue;
            float dn=rays[i].dx*rays[j].dy-rays[i].dy*rays[j].dx;
            if(fabsf(dn)<1e-10f)continue;
            float ex=rays[j].ox-rays[i].ox,ey=rays[j].oy-rays[i].oy;
            float t=(ex*rays[j].dy-ey*rays[j].dx)/dn;
            float u=(ex*rays[i].dy-ey*rays[i].dx)/dn;
            if(t>.001f&&u>.001f&&t<bt){bt=t;bx=rays[i].ox+t*rays[i].dx;by=rays[i].oy+t*rays[i].dy;ok=true;}}
        if(ok){s.push_back({rays[i].ox,rays[i].oy});s.push_back({bx,by});}}
}

static void sqHankin(float sc,float cx,float cy,float maxR,float caDeg,float rot,std::vector<Vec2>&s){
    float ca=caDeg*PI/180.f;
    int n=(int)ceilf(maxR*2.5f/sc)+1;float h=n*sc*.5f;
    float cR=cosf(rot),sR=sinf(rot);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        float lx=-h+i*sc,ly=-h+j*sc;
        float off[4][2]={{0,0},{sc,0},{sc,sc},{0,sc}};
        Vec2 q[4];for(int k=0;k<4;k++){float px=lx+off[k][0],py=ly+off[k][1];
            q[k]={px*cR-py*sR+cx,px*sR+py*cR+cy};}
        float qx=(q[0].x+q[2].x)*.5f,qy=(q[0].y+q[2].y)*.5f;
        if(sqrtf((qx-cx)*(qx-cx)+(qy-cy)*(qy-cy))>maxR*1.2f)continue;
        hankinPoly(q,4,ca,s);}
}

// ============================================================
// Ornament generator — outputs to per-layer segment vectors
// ============================================================
// --- Helper: add a star polygon {n/k} as segments ---
static void starPoly(float cx,float cy,float r,int n,int skip,float rot,std::vector<Vec2>&s){
    Vec2 pts[64];
    for(int i=0;i<n;i++){float a=2*PI*i/n+rot;pts[i]={cx+r*cosf(a),cy+r*sinf(a)};}
    for(int i=0;i<n;i++){s.push_back(pts[i]);s.push_back(pts[(i+skip)%n]);}}

// --- Helper: add a regular polygon as segments ---
static void regPoly(float cx,float cy,float r,int n,float rot,std::vector<Vec2>&s){
    for(int i=0;i<n;i++){
        float a1=2*PI*i/n+rot,a2=2*PI*(i+1)/n+rot;
        s.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
        s.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}}

// --- Helper: polar star boundary r(θ) = R0 + R1*cos(nθ) ---
static void polarStar(float cx,float cy,float R0,float R1,int folds,int res,std::vector<Vec2>&s){
    for(int i=0;i<res;i++){
        float a1=2*PI*i/res,a2=2*PI*(i+1)/res;
        float r1=R0+R1*cosf(folds*a1),r2=R0+R1*cosf(folds*a2);
        s.push_back({cx+r1*cosf(a1),cy+r1*sinf(a1)});
        s.push_back({cx+r2*cosf(a2),cy+r2*sinf(a2)});}}

// --- Helper: Hankin on triangle grid centered ---
static void triHankin(float sc,float cx,float cy,float maxR,float caDeg,std::vector<Vec2>&s){
    float ca=caDeg*PI/180.f;
    float e1x=sc,e2x=sc*.5f,e2y=sc*SQRT3*.5f;
    int ni=(int)ceilf(maxR*2.5f/sc)+2,nj=(int)ceilf(maxR*2.5f/(sc*SQRT3*.5f))+2;
    float ox=cx-(ni*e1x+nj*e2x)*.5f,oy=cy-nj*e2y*.5f;
    auto vtx=[&](int i,int j)->Vec2{return{i*e1x+j*e2x+ox,j*e2y+oy};};
    for(int i=-1;i<=ni;i++)for(int j=-1;j<=nj;j++){
        Vec2 v00=vtx(i,j),v10=vtx(i+1,j),v01=vtx(i,j+1),v11=vtx(i+1,j+1);
        float tx=(v00.x+v10.x+v01.x)/3,ty=(v00.y+v10.y+v01.y)/3;
        if(sqrtf((tx-cx)*(tx-cx)+(ty-cy)*(ty-cy))<maxR){Vec2 t[3]={v00,v10,v01};hankinPoly(t,3,ca,s);}
        tx=(v10.x+v11.x+v01.x)/3;ty=(v10.y+v11.y+v01.y)/3;
        if(sqrtf((tx-cx)*(tx-cx)+(ty-cy)*(ty-cy))<maxR){Vec2 t[3]={v10,v11,v01};hankinPoly(t,3,ca,s);}
    }
}

// --- Helper: Hankin on hexagonal cells ---
static void hexHankin(float sc,float cx,float cy,float maxR,float caDeg,std::vector<Vec2>&s){
    float ca=caDeg*PI/180.f;
    float colStep=sc*SQRT3,rowStep=sc*1.5f;
    int cols=(int)ceilf(maxR*2/colStep)+3,rows=(int)ceilf(maxR*2/rowStep)+3;
    for(int r=-rows;r<=rows;r++)for(int c=-cols;c<=cols;c++){
        float hx=c*colStep+(r%2!=0?colStep*.5f:0)+cx;
        float hy=r*rowStep+cy;
        if(sqrtf((hx-cx)*(hx-cx)+(hy-cy)*(hy-cy))>maxR*1.2f)continue;
        Vec2 hex[6];
        for(int i=0;i<6;i++){float a=PI*i/3-PI/6;hex[i]={hx+sc*cosf(a),hy+sc*sinf(a)};}
        hankinPoly(hex,6,ca,s);
    }
}

// ============================================================
// Ornament 0: 12-fold Mandala (original)
// ============================================================
static void generateOrnament_0(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float R1=R0*starDepth;
    float a=R0*0.16f;
    float modDist=R0*.78f, modW=R0*.24f, modH=R0*.40f, cham=modW*.35f;

    auto inside=[&](float x,float y)->bool{
        float dx=x-cx,dy=y-cy,r=sqrtf(dx*dx+dy*dy),th=atan2f(dy,dx);
        if(r<R0+R1*cosf(12*th)+a*.4f)return true;
        for(int k=0;k<4;k++){float rot=PI*.5f*k,cA=cosf(rot),sA=sinf(rot);
            float lx=dx*cA+dy*sA,ly=-dx*sA+dy*cA;
            if(lx>modDist-modH-a*.35f&&lx<modDist+modH+a*.35f&&fabsf(ly)<modW+a*.35f)return true;}
        return false;};

    auto clip=[&](std::vector<Vec2>&raw,std::vector<Vec2>&out){
        for(size_t i=0;i+1<raw.size();i+=2)
            if(inside(raw[i].x,raw[i].y)&&inside(raw[i+1].x,raw[i+1].y)){
                out.push_back(raw[i]);out.push_back(raw[i+1]);}};

    // ==========================================================
    // L_STAR_BOUNDARY: multi-ring 12-fold star boundaries
    // ==========================================================
    {
        auto&L=layers[L_STAR_BOUNDARY];
        int SN=480;
        // Primary 12-fold star
        polarStar(cx,cy,R0,R1,12,SN,L);
        // Secondary 12-fold ring
        polarStar(cx,cy,R0*.58f,R1*.58f,12,SN,L);
        // Tertiary 12-fold ring
        polarStar(cx,cy,R0*.35f,R1*.35f,12,SN,L);
        // Inner circle
        regPoly(cx,cy,R0*.2f,SN,0,L);
        // Outer decorative 24-fold ripple
        polarStar(cx,cy,R0*.78f,R1*.15f,24,SN,L);
        // 12-pointed star polygon {12/5} at primary radius
        starPoly(cx,cy,R0*.92f,12,5,0,L);
        // Inner {12/5}
        starPoly(cx,cy,R0*.48f,12,5,PI/12,L);
        // Tiny {12/5} at center
        starPoly(cx,cy,R0*.15f,12,5,0,L);
    }

    // ==========================================================
    // L_SPOKES: 12 + 24 radial axes
    // ==========================================================
    {
        auto&L=layers[L_SPOKES];
        // 12 primary spokes (as 6 diameters)
        for(int k=0;k<6;k++){
            float th=2*PI*k/12;
            float rA=R0+R1*cosf(12*th),rB=R0+R1*cosf(12*(th+PI));
            L.push_back({cx+rA*cosf(th),cy+rA*sinf(th)});
            L.push_back({cx+rB*cosf(th+PI),cy+rB*sinf(th+PI)});}
        // 12 secondary spokes (between primary, shorter)
        for(int k=0;k<12;k++){
            float th=2*PI*k/12+PI/12;
            float rr=R0*.58f+R1*.58f*cosf(12*th);
            L.push_back({cx+R0*.2f*cosf(th),cy+R0*.2f*sinf(th)});
            L.push_back({cx+rr*cosf(th),cy+rr*sinf(th)});}
        // Connecting arcs between spoke tips and star boundary
        // (radial tick marks at each spoke intersection with secondary ring)
        for(int k=0;k<12;k++){
            float th=2*PI*k/12;
            float r1=R0*.52f,r2=R0*.62f;
            float da=PI/48;
            L.push_back({cx+r1*cosf(th-da),cy+r1*sinf(th-da)});
            L.push_back({cx+r2*cosf(th),cy+r2*sinf(th)});
            L.push_back({cx+r2*cosf(th),cy+r2*sinf(th)});
            L.push_back({cx+r1*cosf(th+da),cy+r1*sinf(th+da)});}
    }

    // ==========================================================
    // L_GRID_MAIN: Hankin on square + triangle + hex grids
    // ==========================================================
    {
        std::vector<Vec2> raw;raw.reserve(80000);
        // Square grid at 0° and 45°
        sqHankin(a,      cx,cy,R0*1.25f,caDeg,0,raw);
        sqHankin(a,      cx,cy,R0*1.25f,caDeg,PI/4,raw);
        // Additional square grid at 22.5° (creates 16-fold intersections)
        sqHankin(a*1.2f, cx,cy,R0*1.0f, caDeg*.95f,PI/8,raw);
        // Hexagonal Hankin layer (6-fold detail)
        hexHankin(a*.9f, cx,cy,R0*.85f, caDeg*.85f,raw);
        clip(raw,layers[L_GRID_MAIN]);
    }

    // ==========================================================
    // L_GRID_FINE: multi-scale fine detail
    // ==========================================================
    {
        std::vector<Vec2> raw;raw.reserve(40000);
        // Fine square grid in center
        sqHankin(a*.4f,  cx,cy,R0*.38f,caDeg*.9f,0,raw);
        sqHankin(a*.4f,  cx,cy,R0*.38f,caDeg*.9f,PI/4,raw);
        // Fine triangle grid (6-fold rosettes)
        triHankin(a*.5f, cx,cy,R0*.4f, caDeg*.85f,raw);
        // Micro detail at very center
        sqHankin(a*.22f, cx,cy,R0*.18f,caDeg*.8f,0,raw);
        sqHankin(a*.22f, cx,cy,R0*.18f,caDeg*.8f,PI/4,raw);
        clip(raw,layers[L_GRID_FINE]);
    }

    // ==========================================================
    // L_ROTATING: audio-reactive rotating layers
    // ==========================================================
    {
        std::vector<Vec2> raw;raw.reserve(30000);
        sqHankin(a*1.4f, cx,cy,R0*.9f, caDeg,layerRot,raw);
        sqHankin(a*1.0f, cx,cy,R0*.6f, caDeg*.9f,layerRot*1.5f,raw);
        // Rotating hex layer
        hexHankin(a*1.2f,cx,cy,R0*.7f, caDeg*.8f,raw);
        clip(raw,layers[L_ROTATING]);
    }

    // ==========================================================
    // L_MODULES: module borders + internal structure
    // ==========================================================
    {
        auto&L=layers[L_MODULES];
        Vec2 baseP[8]={
            {modDist-modH+cham,-modW},{modDist+modH-cham,-modW},
            {modDist+modH,-modW+cham},{modDist+modH,modW-cham},
            {modDist+modH-cham,modW},{modDist-modH+cham,modW},
            {modDist-modH,modW-cham},{modDist-modH,-modW+cham}};

        for(int k=0;k<4;k++){
            float rot=PI*.5f*k,cA=cosf(rot),sA=sinf(rot);
            float sCx=modDist*cA+cx,sCy=modDist*sA+cy;

            // Rotate helper
            auto rotV=[&](float x,float y)->Vec2{return{x*cA-y*sA+cx,x*sA+y*cA+cy};};

            // Outer module border
            Vec2 rP[8];
            for(int i=0;i<8;i++) rP[i]=rotV(baseP[i].x,baseP[i].y);
            for(int i=0;i<8;i++){L.push_back(rP[i]);L.push_back(rP[(i+1)%8]);}

            // Inner border (75%)
            Vec2 inner[8];
            for(int i=0;i<8;i++){
                float px=(baseP[i].x-modDist)*.75f+modDist,py=baseP[i].y*.75f;
                inner[i]=rotV(px,py);}
            for(int i=0;i<8;i++){L.push_back(inner[i]);L.push_back(inner[(i+1)%8]);}

            // Innermost border (45%)
            Vec2 inner2[8];
            for(int i=0;i<8;i++){
                float px=(baseP[i].x-modDist)*.45f+modDist,py=baseP[i].y*.45f;
                inner2[i]=rotV(px,py);}
            for(int i=0;i<8;i++){L.push_back(inner2[i]);L.push_back(inner2[(i+1)%8]);}

            // Cross-hatching: connect outer vertices to star vertices
            float r2=modW*.7f;
            Vec2 s8[8];
            for(int i=0;i<8;i++){float sa=2*PI*i/8+rot;
                s8[i]={sCx+r2*cosf(sa),sCy+r2*sinf(sa)};}
            for(int i=0;i<8;i++){L.push_back(rP[i]);L.push_back(s8[i%8]);}
            // Also connect inner border to star
            for(int i=0;i<8;i++){L.push_back(inner[i]);L.push_back(s8[i%8]);}

            // Diamond lattice inside module
            for(int i=0;i<8;i++){
                L.push_back(rP[i]);L.push_back(inner[(i+1)%8]);
                L.push_back(inner[i]);L.push_back(rP[(i+1)%8]);}

            // Small circles at module tips
            float tipX=modDist+modH*.85f,tipY=0;
            Vec2 tp=rotV(tipX,tipY);
            regPoly(tp.x,tp.y,modW*.18f,8,rot,L);
            // Circle at inner end
            float tipX2=modDist-modH*.75f,tipY2=0;
            Vec2 tp2=rotV(tipX2,tipY2);
            regPoly(tp2.x,tp2.y,modW*.15f,8,rot,L);
        }

        // Diagonal inter-module connectors (between adjacent arms)
        for(int k=0;k<4;k++){
            float th1=PI*.5f*k+PI/4; // midway between arms
            float rd=R0*.65f;
            float dth=PI/14;
            Vec2 a1={cx+rd*cosf(th1-dth),cy+rd*sinf(th1-dth)};
            Vec2 a2={cx+rd*cosf(th1+dth),cy+rd*sinf(th1+dth)};
            Vec2 b1={cx+rd*1.15f*cosf(th1),cy+rd*1.15f*sinf(th1)};
            Vec2 b2={cx+rd*.85f*cosf(th1),cy+rd*.85f*sinf(th1)};
            L.push_back(a1);L.push_back(b1);
            L.push_back(b1);L.push_back(a2);
            L.push_back(a2);L.push_back(b2);
            L.push_back(b2);L.push_back(a1);
        }
    }

    // ==========================================================
    // L_STARS8: {8/3} + {6/2} + {10/3} star polygons
    // ==========================================================
    {
        auto&L=layers[L_STARS8];
        for(int k=0;k<4;k++){
            float rot=PI*.5f*k,cA=cosf(rot),sA=sinf(rot);
            float sCx=modDist*cA+cx,sCy=modDist*sA+cy;
            float r2=modW*.72f;

            // Primary {8/3}
            starPoly(sCx,sCy,r2,8,3,rot,L);
            // Inner {8/3} rotated
            starPoly(sCx,sCy,r2*.42f,8,3,rot+PI/8,L);
            // Tiny {8/2} (octagram) at center of module
            starPoly(sCx,sCy,r2*.2f,8,2,rot,L);
            // Circumscribed octagon
            regPoly(sCx,sCy,r2*.88f,8,rot+PI/8,L);
        }

        // {12/5} stars at diagonal positions (between modules)
        for(int k=0;k<4;k++){
            float th=PI*.5f*k+PI/4;
            float rd=R0*.62f;
            float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
            starPoly(sx,sy,a*.55f,12,5,th,L);
            regPoly(sx,sy,a*.65f,12,th+PI/12,L);
        }

        // Ring of small {6/2} stars around center
        for(int k=0;k<12;k++){
            float th=2*PI*k/12;
            float rd=R0*.38f;
            float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
            starPoly(sx,sy,a*.32f,6,2,th,L);
        }

        // {10/3} stars at secondary ring intersections
        for(int k=0;k<12;k++){
            float th=2*PI*k/12+PI/12;
            float rd=R0*.55f;
            float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
            starPoly(sx,sy,a*.25f,10,3,th,L);
        }
    }
}

// ============================================================
// Ornament 1: Hexagonal Lattice (6-fold)
// ============================================================
static void generateOrnament_1(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float a=R0*0.14f;
    (void)starDepth;

    // L_STAR_BOUNDARY: concentric hexagonal rings
    {
        auto&L=layers[L_STAR_BOUNDARY];
        for(int ring=1;ring<=5;ring++){
            float r=R0*ring*.18f;
            regPoly(cx,cy,r,6,PI/6,L);
        }
        polarStar(cx,cy,R0*.92f,R0*.08f,6,360,L);
    }
    // L_SPOKES: 6 radial spokes
    {
        auto&L=layers[L_SPOKES];
        for(int k=0;k<6;k++){
            float th=PI*k/3+PI/6;
            L.push_back({cx,cy});
            L.push_back({cx+R0*.92f*cosf(th),cy+R0*.92f*sinf(th)});
        }
    }
    // L_GRID_MAIN: hexagonal Hankin tiling
    {
        std::vector<Vec2> raw;raw.reserve(60000);
        hexHankin(a,cx,cy,R0*1.1f,caDeg,raw);
        layers[L_GRID_MAIN]=raw;
    }
    // L_GRID_FINE: triangular Hankin overlay
    {
        std::vector<Vec2> raw;raw.reserve(40000);
        triHankin(a*.7f,cx,cy,R0*.7f,caDeg*.85f,raw);
        layers[L_GRID_FINE]=raw;
    }
    // L_MODULES: {6/2} stars at hex vertices
    {
        auto&L=layers[L_MODULES];
        for(int ring=1;ring<=3;ring++){
            float rd=R0*ring*.25f;
            for(int k=0;k<6;k++){
                float th=PI*k/3+PI/6;
                float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
                starPoly(sx,sy,a*.6f,6,2,th,L);
                regPoly(sx,sy,a*.75f,6,th+PI/6,L);
            }
        }
    }
    // L_STARS8: central rosette
    {
        auto&L=layers[L_STARS8];
        starPoly(cx,cy,R0*.35f,12,5,0,L);
        starPoly(cx,cy,R0*.2f,6,2,PI/6,L);
        regPoly(cx,cy,R0*.4f,12,PI/12,L);
    }
    // L_ROTATING: rotating hex Hankin
    {
        std::vector<Vec2> raw;raw.reserve(20000);
        hexHankin(a*1.5f,cx,cy,R0*.8f,caDeg*.9f,raw);
        // manual rotate
        float c=cosf(layerRot),s=sinf(layerRot);
        for(auto&v:raw){float dx=v.x-cx,dy=v.y-cy;v.x=dx*c-dy*s+cx;v.y=dx*s+dy*c+cy;}
        layers[L_ROTATING]=raw;
    }
}

// ============================================================
// Ornament 2: Square Kufic / Rectilinear
// ============================================================
static void generateOrnament_2(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float a=R0*0.12f;
    (void)starDepth;

    // L_STAR_BOUNDARY: square boundary rings
    {
        auto&L=layers[L_STAR_BOUNDARY];
        for(int ring=1;ring<=4;ring++){
            float r=R0*ring*.22f;
            regPoly(cx,cy,r,4,PI/4,L);
            regPoly(cx,cy,r*.92f,4,0,L); // rotated 45°
        }
    }
    // L_SPOKES: diagonal + axis lines
    {
        auto&L=layers[L_SPOKES];
        for(int k=0;k<4;k++){
            float th=PI*k/4;
            L.push_back({cx+R0*.15f*cosf(th),cy+R0*.15f*sinf(th)});
            L.push_back({cx+R0*.88f*cosf(th),cy+R0*.88f*sinf(th)});
        }
    }
    // L_GRID_MAIN: dense square Hankin at 0°
    {
        std::vector<Vec2> raw;raw.reserve(60000);
        sqHankin(a,cx,cy,R0*1.1f,caDeg,0,raw);
        layers[L_GRID_MAIN]=raw;
    }
    // L_GRID_FINE: overlaid 45° square Hankin
    {
        std::vector<Vec2> raw;raw.reserve(60000);
        sqHankin(a,cx,cy,R0*1.1f,caDeg,PI/4,raw);
        layers[L_GRID_FINE]=raw;
    }
    // L_MODULES: {8/3} stars at grid intersections
    {
        auto&L=layers[L_MODULES];
        for(int i=-3;i<=3;i++)for(int j=-3;j<=3;j++){
            float sx=cx+i*a*2.5f,sy=cy+j*a*2.5f;
            float d=sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy));
            if(d>R0*.9f)continue;
            starPoly(sx,sy,a*.7f,8,3,0,L);
        }
    }
    // L_STARS8: larger {8/2} at corners
    {
        auto&L=layers[L_STARS8];
        float pos[]={-1,-1, 1,-1, 1,1, -1,1};
        for(int k=0;k<4;k++){
            float sx=cx+pos[k*2]*R0*.55f,sy=cy+pos[k*2+1]*R0*.55f;
            starPoly(sx,sy,a*1.2f,8,2,PI/8,L);
            regPoly(sx,sy,a*1.4f,8,0,L);
        }
        starPoly(cx,cy,R0*.25f,8,3,0,L);
    }
    // L_ROTATING
    {
        std::vector<Vec2> raw;raw.reserve(30000);
        sqHankin(a*1.8f,cx,cy,R0*.85f,caDeg*.9f,layerRot,raw);
        layers[L_ROTATING]=raw;
    }
}

// ============================================================
// Ornament 3: 8-fold Star (Octagonal)
// ============================================================
static void generateOrnament_3(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float a=R0*0.15f;
    float R1=R0*starDepth;

    // L_STAR_BOUNDARY: 8-fold polar star
    {
        auto&L=layers[L_STAR_BOUNDARY];
        polarStar(cx,cy,R0,R1,8,480,L);
        polarStar(cx,cy,R0*.6f,R1*.5f,8,480,L);
        regPoly(cx,cy,R0*.25f,8,PI/8,L);
        starPoly(cx,cy,R0*.85f,8,3,0,L);
    }
    // L_SPOKES: 8 radial lines
    {
        auto&L=layers[L_SPOKES];
        for(int k=0;k<8;k++){
            float th=2*PI*k/8;
            float rr=R0+R1*cosf(8*th);
            L.push_back({cx+R0*.2f*cosf(th),cy+R0*.2f*sinf(th)});
            L.push_back({cx+rr*cosf(th),cy+rr*sinf(th)});
        }
    }
    // L_GRID_MAIN: sq Hankin 0° and 45°
    {
        std::vector<Vec2> raw;raw.reserve(60000);
        sqHankin(a,cx,cy,R0*1.2f,caDeg,0,raw);
        sqHankin(a,cx,cy,R0*1.2f,caDeg,PI/4,raw);
        layers[L_GRID_MAIN]=raw;
    }
    // L_GRID_FINE: fine 22.5° grid
    {
        std::vector<Vec2> raw;raw.reserve(30000);
        sqHankin(a*.6f,cx,cy,R0*.5f,caDeg*.9f,PI/8,raw);
        layers[L_GRID_FINE]=raw;
    }
    // L_MODULES: octagonal frames
    {
        auto&L=layers[L_MODULES];
        for(int ring=1;ring<=3;ring++){
            float rd=R0*ring*.28f;
            for(int k=0;k<8;k++){
                float th=2*PI*k/8;
                float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
                regPoly(sx,sy,a*.5f,8,th,L);
            }
        }
    }
    // L_STARS8: {8/3} constellation
    {
        auto&L=layers[L_STARS8];
        starPoly(cx,cy,R0*.45f,8,3,0,L);
        for(int k=0;k<8;k++){
            float th=2*PI*k/8+PI/8;
            float sx=cx+R0*.65f*cosf(th),sy=cy+R0*.65f*sinf(th);
            starPoly(sx,sy,a*.55f,8,3,th,L);
        }
    }
    // L_ROTATING
    {
        std::vector<Vec2> raw;raw.reserve(20000);
        sqHankin(a*1.3f,cx,cy,R0*.75f,caDeg*.85f,layerRot,raw);
        layers[L_ROTATING]=raw;
    }
}

// ============================================================
// Ornament 4: 10-fold Decagonal
// ============================================================
static void generateOrnament_4(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float a=R0*0.14f;
    float R1=R0*starDepth;

    // L_STAR_BOUNDARY: 10-fold star
    {
        auto&L=layers[L_STAR_BOUNDARY];
        polarStar(cx,cy,R0,R1,10,500,L);
        polarStar(cx,cy,R0*.55f,R1*.45f,10,500,L);
        starPoly(cx,cy,R0*.88f,10,3,0,L);
        starPoly(cx,cy,R0*.5f,10,4,PI/10,L);
    }
    // L_SPOKES: 10 radial
    {
        auto&L=layers[L_SPOKES];
        for(int k=0;k<10;k++){
            float th=2*PI*k/10;
            float rr=R0+R1*cosf(10*th);
            L.push_back({cx+R0*.18f*cosf(th),cy+R0*.18f*sinf(th)});
            L.push_back({cx+rr*cosf(th),cy+rr*sinf(th)});
        }
        for(int k=0;k<10;k++){
            float th=2*PI*k/10+PI/10;
            L.push_back({cx+R0*.3f*cosf(th),cy+R0*.3f*sinf(th)});
            L.push_back({cx+R0*.55f*cosf(th),cy+R0*.55f*sinf(th)});
        }
    }
    // L_GRID_MAIN: 5 overlapping square grids at 36° increments
    {
        std::vector<Vec2> raw;raw.reserve(80000);
        for(int g=0;g<5;g++)
            sqHankin(a,cx,cy,R0*1.1f,caDeg,g*PI/5,raw);
        layers[L_GRID_MAIN]=raw;
    }
    // L_GRID_FINE: fine penrose-like detail
    {
        std::vector<Vec2> raw;raw.reserve(30000);
        sqHankin(a*.5f,cx,cy,R0*.4f,caDeg*.9f,PI/10,raw);
        triHankin(a*.6f,cx,cy,R0*.45f,caDeg*.8f,raw);
        layers[L_GRID_FINE]=raw;
    }
    // L_MODULES: decagonal frames
    {
        auto&L=layers[L_MODULES];
        for(int k=0;k<10;k++){
            float th=2*PI*k/10;
            float sx=cx+R0*.7f*cosf(th),sy=cy+R0*.7f*sinf(th);
            regPoly(sx,sy,a*.5f,10,th,L);
            regPoly(sx,sy,a*.35f,5,th+PI/10,L);
        }
    }
    // L_STARS8: {10/3} and {10/4} stars
    {
        auto&L=layers[L_STARS8];
        starPoly(cx,cy,R0*.35f,10,3,0,L);
        for(int k=0;k<10;k++){
            float th=2*PI*k/10+PI/10;
            float sx=cx+R0*.55f*cosf(th),sy=cy+R0*.55f*sinf(th);
            starPoly(sx,sy,a*.4f,10,4,th,L);
        }
        for(int k=0;k<5;k++){
            float th=2*PI*k/5;
            float sx=cx+R0*.4f*cosf(th),sy=cy+R0*.4f*sinf(th);
            starPoly(sx,sy,a*.3f,5,2,th,L);
        }
    }
    // L_ROTATING
    {
        std::vector<Vec2> raw;raw.reserve(20000);
        sqHankin(a*1.5f,cx,cy,R0*.7f,caDeg*.85f,layerRot,raw);
        layers[L_ROTATING]=raw;
    }
}

// ============================================================
// Ornament 5: Muqarnas Cascade (radial nested)
// ============================================================
static void generateOrnament_5(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float a=R0*0.13f;
    float R1=R0*starDepth;

    // L_STAR_BOUNDARY: many concentric 12-fold rings
    {
        auto&L=layers[L_STAR_BOUNDARY];
        for(int ring=2;ring<=8;ring++){
            float r=R0*ring*.12f;
            float depth=R1*ring*.1f;
            polarStar(cx,cy,r,depth,12,480,L);
        }
    }
    // L_SPOKES: 12 + 12 interleaved
    {
        auto&L=layers[L_SPOKES];
        for(int k=0;k<12;k++){
            float th=2*PI*k/12;
            L.push_back({cx+R0*.1f*cosf(th),cy+R0*.1f*sinf(th)});
            L.push_back({cx+R0*.96f*cosf(th),cy+R0*.96f*sinf(th)});
        }
        for(int k=0;k<12;k++){
            float th=2*PI*k/12+PI/12;
            L.push_back({cx+R0*.25f*cosf(th),cy+R0*.25f*sinf(th)});
            L.push_back({cx+R0*.75f*cosf(th),cy+R0*.75f*sinf(th)});
        }
    }
    // L_GRID_MAIN: hex Hankin
    {
        std::vector<Vec2> raw;raw.reserve(60000);
        hexHankin(a,cx,cy,R0*1.1f,caDeg,raw);
        hexHankin(a*1.5f,cx,cy,R0*1.1f,caDeg*.95f,raw);
        layers[L_GRID_MAIN]=raw;
    }
    // L_GRID_FINE: tri + sq overlay
    {
        std::vector<Vec2> raw;raw.reserve(40000);
        triHankin(a*.6f,cx,cy,R0*.6f,caDeg*.85f,raw);
        sqHankin(a*.5f,cx,cy,R0*.5f,caDeg*.9f,0,raw);
        layers[L_GRID_FINE]=raw;
    }
    // L_MODULES: nested {12/5} and {6/2}
    {
        auto&L=layers[L_MODULES];
        for(int ring=1;ring<=4;ring++){
            float rd=R0*ring*.2f;
            starPoly(cx,cy,rd,12,5,ring*PI/24,L);
            for(int k=0;k<12;k++){
                float th=2*PI*k/12+ring*PI/12;
                float sx=cx+rd*cosf(th),sy=cy+rd*sinf(th);
                starPoly(sx,sy,a*.35f,6,2,th,L);
            }
        }
    }
    // L_STARS8: {8/3} ring
    {
        auto&L=layers[L_STARS8];
        for(int k=0;k<12;k++){
            float th=2*PI*k/12;
            float sx=cx+R0*.65f*cosf(th),sy=cy+R0*.65f*sinf(th);
            starPoly(sx,sy,a*.5f,8,3,th,L);
        }
        starPoly(cx,cy,R0*.2f,12,5,0,L);
    }
    // L_ROTATING
    {
        std::vector<Vec2> raw;raw.reserve(20000);
        hexHankin(a*1.3f,cx,cy,R0*.75f,caDeg*.8f,raw);
        float c=cosf(layerRot),s=sinf(layerRot);
        for(auto&v:raw){float dx=v.x-cx,dy=v.y-cy;v.x=dx*c-dy*s+cx;v.y=dx*s+dy*c+cy;}
        layers[L_ROTATING]=raw;
    }
}

// ============================================================
// Ornament dispatcher
// ============================================================
typedef void(*OrnamentFunc)(float,float,float,float,float,float,std::vector<Vec2>[NUM_LAYERS]);
static OrnamentFunc ornamentFuncs[NUM_ORNAMENTS] = {
    generateOrnament_0, generateOrnament_1, generateOrnament_2,
    generateOrnament_3, generateOrnament_4, generateOrnament_5
};
static const char* ornamentNames[NUM_ORNAMENTS] = {
    "12-fold Mandala", "Hexagonal", "Square Kufic",
    "Octagonal", "Decagonal", "Muqarnas"
};

// ============================================================
// Generate simplified thumbnail geometry (few segments, no Hankin)
// ============================================================
static void generateThumb(int idx, float sz, std::vector<Vec2>& out) {
    float cx=sz*.5f, cy=sz*.5f, r=sz*.38f;
    switch(idx) {
    case 0: // 12-fold mandala: 12-pointed star + circle
        for(int i=0;i<12;i++){
            float a1=2*PI*i/12,a2=2*PI*((i+5)%12)/12;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<24;i++){
            float a1=2*PI*i/24,a2=2*PI*(i+1)/24;float rr=r*.55f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 1: // Hexagonal: hex + inner {6/2}
        for(int i=0;i<6;i++){
            float a1=PI*i/3+PI/6,a2=PI*(i+1)/3+PI/6;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<6;i++){
            float a1=PI*i/3+PI/6,a2=PI*((i+2)%6)/3+PI/6;float rr=r*.6f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 2: // Square: overlapping squares
        for(int q=0;q<2;q++){
            float rot=q*PI/4;float rr=r*(q==0?1.f:.72f);
            for(int i=0;i<4;i++){
                float a1=PI*i/2+PI/4+rot,a2=PI*(i+1)/2+PI/4+rot;
                out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
                out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}}
        // {8/3} in center
        for(int i=0;i<8;i++){
            float a1=2*PI*i/8,a2=2*PI*((i+3)%8)/8;float rr=r*.4f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 3: // Octagonal: octagon + {8/3}
        for(int i=0;i<8;i++){
            float a1=2*PI*i/8,a2=2*PI*(i+1)/8;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<8;i++){
            float a1=2*PI*i/8,a2=2*PI*((i+3)%8)/8;float rr=r*.7f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 4: // Decagonal: decagon + {10/3}
        for(int i=0;i<10;i++){
            float a1=2*PI*i/10,a2=2*PI*(i+1)/10;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<10;i++){
            float a1=2*PI*i/10,a2=2*PI*((i+3)%10)/10;float rr=r*.65f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 5: // Muqarnas: concentric circles + 12 spokes
        for(int ring=1;ring<=3;ring++){
            float rr=r*ring*.32f;
            for(int i=0;i<24;i++){
                float a1=2*PI*i/24,a2=2*PI*(i+1)/24;
                out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
                out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}}
        for(int k=0;k<12;k++){
            float th=2*PI*k/12;
            out.push_back({cx+r*.15f*cosf(th),cy+r*.15f*sinf(th)});
            out.push_back({cx+r*.95f*cosf(th),cy+r*.95f*sinf(th)});}
        break;
    }
}

// ============================================================
// Chain builder — closed contours only
// ============================================================
struct ChainBuilder {
    std::vector<float> vx,vy;
    std::vector<std::vector<int>> adj;
    std::unordered_map<int64_t,std::vector<int>> grid;
    float eps;

    void clear(float e){vx.clear();vy.clear();adj.clear();grid.clear();eps=e;}
    int64_t ck(float x,float y){
        int gx=(int)floorf(x/(eps*2)),gy=(int)floorf(y/(eps*2));
        return((int64_t)(gx+500000)<<21)^(int64_t)(gy+500000);}
    int findOrAdd(float x,float y){
        float e2=eps*eps;int gx=(int)floorf(x/(eps*2)),gy=(int)floorf(y/(eps*2));
        for(int dx=-1;dx<=1;dx++)for(int dy=-1;dy<=1;dy++){
            int64_t k=((int64_t)(gx+dx+500000)<<21)^(int64_t)(gy+dy+500000);
            auto it=grid.find(k);if(it!=grid.end())
                for(int idx:it->second){float ddx=vx[idx]-x,ddy=vy[idx]-y;
                    if(ddx*ddx+ddy*ddy<e2)return idx;}}
        int idx=(int)vx.size();vx.push_back(x);vy.push_back(y);
        adj.push_back({});grid[ck(x,y)].push_back(idx);return idx;}
    void addSeg(Vec2 a,Vec2 b){
        int ai=findOrAdd(a.x,a.y),bi=findOrAdd(b.x,b.y);
        if(ai==bi)return;adj[ai].push_back(bi);adj[bi].push_back(ai);}

    // Iteratively remove all degree-1 vertices (dangling ends)
    void prune(){
        std::vector<int> q;
        for(int v=0;v<(int)adj.size();v++)
            if((int)adj[v].size()==1) q.push_back(v);
        while(!q.empty()){
            int v=q.back();q.pop_back();
            if(adj[v].size()!=1)continue;
            int other=adj[v][0]; adj[v].clear();
            auto&oa=adj[other];
            for(int i=0;i<(int)oa.size();i++)
                if(oa[i]==v){oa.erase(oa.begin()+i);break;}
            if(oa.size()==1)q.push_back(other);
        }
    }

    // Trace only closed loops (GL_LINE_LOOP)
    void build(std::vector<Vec2>&oV,std::vector<GLint>&oS,std::vector<GLsizei>&oC){
        int nv=(int)vx.size();
        std::vector<std::vector<bool>>used(nv);
        for(int i=0;i<nv;i++)used[i].assign(adj[i].size(),false);

        auto mark=[&](int f,int ei){int to=adj[f][ei];used[f][ei]=true;
            for(int i=0;i<(int)adj[to].size();i++)
                if(adj[to][i]==f&&!used[to][i]){used[to][i]=true;break;}};
        auto unused=[&](int v)->int{
            for(int i=0;i<(int)adj[v].size();i++)if(!used[v][i])return i;return-1;};

        // At junctions, pick the straightest continuation
        auto best=[&](int cur,int prev)->int{
            if(prev<0)return unused(cur);
            float inA=atan2f(vy[cur]-vy[prev],vx[cur]-vx[prev]);
            int b=-1;float bs=-999;
            for(int i=0;i<(int)adj[cur].size();i++){if(used[cur][i])continue;
                int nx=adj[cur][i];
                float sc=cosf(atan2f(vy[nx]-vy[cur],vx[nx]-vx[cur])-inA);
                if(sc>bs){bs=sc;b=i;}}return b;};

        // Trace closed contour from start vertex
        auto trace=[&](int start,int fe){
            GLint si=(GLint)oV.size();
            oV.push_back({vx[start],vy[start]});
            int prev=-1,cur=start,ce=fe;
            bool closed=false;
            while(ce>=0){
                int nx=adj[cur][ce];
                mark(cur,ce);
                if(nx==start){closed=true;break;}  // loop closed — don't duplicate start
                oV.push_back({vx[nx],vy[nx]});
                prev=cur;cur=nx;
                ce=best(cur,prev);}
            GLsizei cnt=(GLsizei)((int)oV.size()-si);
            if(closed&&cnt>=3){
                oS.push_back(si);oC.push_back(cnt);
            }else{oV.resize(si);}  // discard open paths
        };

        // Trace from all vertices with unused edges
        for(int v=0;v<nv;v++){
            int e=unused(v);
            while(e>=0){trace(v,e);e=unused(v);}
        }
    }
};

// (morph is now a cross-fade, no geometry distortion needed)

// ============================================================
// Main
// ============================================================
int main(int argc,char*argv[]){
    (void)argc;(void)argv;
    SDL_Init(SDL_INIT_VIDEO|SDL_INIT_AUDIO);
    SDL_GL_SetAttribute(SDL_GL_DOUBLEBUFFER,1);
    SDL_GL_SetAttribute(SDL_GL_STENCIL_SIZE,8);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MAJOR_VERSION,2);
    SDL_GL_SetAttribute(SDL_GL_CONTEXT_MINOR_VERSION,1);

    SDL_Window*win=SDL_CreateWindow("Ornament of Pressure",
        SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,1024,1024,
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext gl=SDL_GL_CreateContext(win);
    glewInit();SDL_GL_SetSwapInterval(1);

    GLuint vbo;glGenBuffers(1,&vbo);

    SDL_AudioSpec want={},have={};
    want.freq=44100;want.format=AUDIO_F32SYS;want.channels=1;want.samples=1024;want.callback=audioCB;
    g_au.dev=SDL_OpenAudioDevice(NULL,1,&want,&have,SDL_AUDIO_ALLOW_FORMAT_CHANGE);
    if(g_au.dev>0)SDL_PauseAudioDevice(g_au.dev,0);

    float sBass=.5f,sMid=.5f,sHigh=0,layerRot=0;
    float sensitivity=1.0f;  // audio sensitivity: 0.1 .. 10.0, controlled by scroll/arrow keys
    Complex fftBuf[FFT_N];float hann[FFT_N];
    for(int i=0;i<FFT_N;i++)hann[i]=.5f*(1-cosf(2*PI*i/(FFT_N-1)));

    // Per-layer rendering data
    struct LayerRender {
        std::vector<Vec2> verts;
        std::vector<GLint> starts;
        std::vector<GLsizei> counts;
    };

    std::vector<Vec2> layerSegs[NUM_LAYERS];
    LayerRender lr[NUM_LAYERS];
    ChainBuilder cb;

    // Ornament selection and morphing state
    int curOrnament = 0;
    int targetOrnament = -1;  // -1 = no transition
    float morphT = 0;         // 0..1: 0=fully old, 1=fully new
    static const float MORPH_SPEED = 0.018f;  // per frame (~3.3s transition)

    // Second set of layer data for cross-fade
    std::vector<Vec2> layerSegs2[NUM_LAYERS];
    LayerRender lr2[NUM_LAYERS];

    // Pre-generate thumbnail geometry
    std::vector<Vec2> thumbGeom[NUM_ORNAMENTS];
    for (int i = 0; i < NUM_ORNAMENTS; i++)
        generateThumb(i, (float)THUMB_SIZE, thumbGeom[i]);

    bool run=true;
    while(run){
        SDL_Event ev;
        while(SDL_PollEvent(&ev)){
            if(ev.type==SDL_QUIT)run=false;
            if(ev.type==SDL_KEYDOWN&&ev.key.keysym.sym==SDLK_ESCAPE)run=false;

            // Click on sidebar thumbnails
            if(ev.type==SDL_MOUSEBUTTONDOWN && ev.button.button==SDL_BUTTON_LEFT){
                int mx=ev.button.x, my=ev.button.y;
                // Get window size (logical pixels)
                int winW,winH;SDL_GetWindowSize(win,&winW,&winH);
                // Get DPI scale
                int drawW,drawH;SDL_GL_GetDrawableSize(win,&drawW,&drawH);
                float dpiScale=(float)drawW/(float)winW;
                // Sidebar is on the left, SIDEBAR_W logical pixels
                if(mx < SIDEBAR_W){
                    for(int i=0;i<NUM_ORNAMENTS;i++){
                        int ty = THUMB_PAD + i*(THUMB_SIZE+THUMB_PAD) + (SIDEBAR_W-THUMB_SIZE)/2;
                        int tx = (SIDEBAR_W-THUMB_SIZE)/2;
                        if(mx>=tx && mx<tx+THUMB_SIZE && my>=ty && my<ty+THUMB_SIZE){
                            if(i!=curOrnament && targetOrnament<0){
                                targetOrnament=i;
                                morphT=0;
                            }
                            break;
                        }
                    }
                }
            }

            // Number keys 1-6 to switch
            if(ev.type==SDL_KEYDOWN && ev.key.keysym.sym>=SDLK_1 && ev.key.keysym.sym<=SDLK_6){
                int idx=ev.key.keysym.sym-SDLK_1;
                if(idx!=curOrnament && targetOrnament<0){
                    targetOrnament=idx;
                    morphT=0;
                }
            }

            // Sensitivity: scroll wheel or up/down arrows
            if(ev.type==SDL_MOUSEWHEEL){
                sensitivity*=(ev.wheel.y>0)?1.15f:0.87f;
                if(sensitivity<0.1f)sensitivity=0.1f;
                if(sensitivity>10.f)sensitivity=10.f;
            }
            if(ev.type==SDL_KEYDOWN&&ev.key.keysym.sym==SDLK_UP){
                sensitivity*=1.2f;if(sensitivity>10.f)sensitivity=10.f;}
            if(ev.type==SDL_KEYDOWN&&ev.key.keysym.sym==SDLK_DOWN){
                sensitivity*=0.83f;if(sensitivity<0.1f)sensitivity=0.1f;}
        }

        // Advance morph
        if(targetOrnament>=0){
            morphT+=MORPH_SPEED;
            if(morphT>=1.0f){
                curOrnament=targetOrnament;
                targetOrnament=-1;
                morphT=0;
            }
        }

        if(g_au.dev>0){
            SDL_LockAudioDevice(g_au.dev);int pos=g_au.wp;SDL_UnlockAudioDevice(g_au.dev);
            for(int i=0;i<FFT_N;i++){int idx=(pos-FFT_N+i+FFT_N*4)%(FFT_N*4);
                fftBuf[i]={g_au.ring[idx]*hann[i],0};}
            fft(fftBuf,FFT_N);
            float bass=0,mid=0,high=0;
            for(int k=1;k<=8;k++)bass+=fftBuf[k].re*fftBuf[k].re+fftBuf[k].im*fftBuf[k].im;
            for(int k=9;k<=80;k++)mid+=fftBuf[k].re*fftBuf[k].re+fftBuf[k].im*fftBuf[k].im;
            for(int k=81;k<=300;k++)high+=fftBuf[k].re*fftBuf[k].re+fftBuf[k].im*fftBuf[k].im;
            float s=sensitivity;
            bass=fminf(sqrtf(bass/8)*3*s,1);mid=fminf(sqrtf(mid/72)*5*s,1);high=fminf(sqrtf(high/220)*6*s,1);
            sBass=sBass*.80f+bass*.20f;sMid=sMid*.82f+mid*.18f;sHigh=sHigh*.85f+high*.15f;
        }

        float R0=250+sBass*120, depth=.12f+sBass*.18f, angle=25+sMid*40;
        layerRot+=sHigh*.008f+.001f;

        int wW,wH;SDL_GL_GetDrawableSize(win,&wW,&wH);
        int winLW,winLH;SDL_GetWindowSize(win,&winLW,&winLH);
        float dpiScale=(float)wW/(float)winLW;
        float sidebarPx=SIDEBAR_W*dpiScale;

        // Generate per-layer segments for current ornament
        for(int i=0;i<NUM_LAYERS;i++) layerSegs[i].clear();
        ornamentFuncs[curOrnament]((float)wW,(float)wH,R0,depth,angle,layerRot,layerSegs);

        // Build chains per layer (current/old ornament)
        for(int i=0;i<NUM_LAYERS;i++){
            lr[i].verts.clear();lr[i].starts.clear();lr[i].counts.clear();
            if(layerSegs[i].empty())continue;
            cb.clear(0.5f);
            for(size_t j=0;j+1<layerSegs[i].size();j+=2)
                cb.addSeg(layerSegs[i][j],layerSegs[i][j+1]);
            cb.prune();
            lr[i].verts.reserve(layerSegs[i].size());
            cb.build(lr[i].verts,lr[i].starts,lr[i].counts);
        }

        // If morphing, also generate the target ornament
        bool morphing=(targetOrnament>=0);
        if(morphing){
            for(int i=0;i<NUM_LAYERS;i++) layerSegs2[i].clear();
            ornamentFuncs[targetOrnament]((float)wW,(float)wH,R0,depth,angle,layerRot,layerSegs2);
            for(int i=0;i<NUM_LAYERS;i++){
                lr2[i].verts.clear();lr2[i].starts.clear();lr2[i].counts.clear();
                if(layerSegs2[i].empty())continue;
                cb.clear(0.5f);
                for(size_t j=0;j+1<layerSegs2[i].size();j+=2)
                    cb.addSeg(layerSegs2[i][j],layerSegs2[i][j+1]);
                cb.prune();
                lr2[i].verts.reserve(layerSegs2[i].size());
                cb.build(lr2[i].verts,lr2[i].starts,lr2[i].counts);
            }
        }

        // --- Render ---
        glViewport(0,0,wW,wH);
        glClearColor(0.02f,0.02f,0.04f,1);
        glClearStencil(0);
        glClear(GL_COLOR_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

        // ---- Draw sidebar background ----
        glMatrixMode(GL_PROJECTION);glLoadIdentity();glOrtho(0,wW,wH,0,-1,1);
        glMatrixMode(GL_MODELVIEW);glLoadIdentity();
        glDisable(GL_LINE_SMOOTH);
        glBegin(GL_QUADS);
        glColor3f(0.06f,0.06f,0.08f);
        glVertex2f(0,0);glVertex2f(sidebarPx,0);
        glVertex2f(sidebarPx,(float)wH);glVertex2f(0,(float)wH);
        glEnd();
        // Sidebar divider line
        glBegin(GL_LINES);
        glColor3f(0.15f,0.15f,0.18f);
        glVertex2f(sidebarPx,0);glVertex2f(sidebarPx,(float)wH);
        glEnd();

        // ---- Draw thumbnails ----
        glEnable(GL_LINE_SMOOTH);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA,GL_ONE_MINUS_SRC_ALPHA);

        for(int ti=0;ti<NUM_ORNAMENTS;ti++){
            float tx=(SIDEBAR_W-THUMB_SIZE)/2*dpiScale;
            float ty=(THUMB_PAD + ti*(THUMB_SIZE+THUMB_PAD) + (SIDEBAR_W-THUMB_SIZE)/2)*dpiScale;
            float tsz=THUMB_SIZE*dpiScale;

            // Highlight box for selected ornament
            bool selected=(ti==curOrnament);
            bool isTarget=(ti==targetOrnament);
            if(selected||isTarget){
                glBegin(GL_LINE_LOOP);
                float br=selected?0.7f:0.4f;
                glColor3f(0.9f*br,0.7f*br,0.2f*br);
                glVertex2f(tx-2,ty-2);glVertex2f(tx+tsz+2,ty-2);
                glVertex2f(tx+tsz+2,ty+tsz+2);glVertex2f(tx-2,ty+tsz+2);
                glEnd();
            }

            // Dark thumbnail background
            glBegin(GL_QUADS);
            glColor3f(0.03f,0.03f,0.05f);
            glVertex2f(tx,ty);glVertex2f(tx+tsz,ty);
            glVertex2f(tx+tsz,ty+tsz);glVertex2f(tx,ty+tsz);
            glEnd();

            // Draw thumbnail geometry scaled to thumbnail viewport
            if(!thumbGeom[ti].empty()){
                // Create scaled copy
                std::vector<Vec2> scaled=thumbGeom[ti];
                for(auto&v:scaled){v.x=v.x*dpiScale+tx;v.y=v.y*dpiScale+ty;}

                glBindBuffer(GL_ARRAY_BUFFER,vbo);
                glBufferData(GL_ARRAY_BUFFER,scaled.size()*sizeof(Vec2),scaled.data(),GL_STREAM_DRAW);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(2,GL_FLOAT,0,0);
                float br=selected?1.0f:(isTarget?0.7f:0.45f);
                glColor3f(0.9f*br,0.75f*br,0.25f*br);
                glLineWidth(1.5f);
                glDrawArrays(GL_LINES,0,(GLsizei)scaled.size());
                glDisableClientState(GL_VERTEX_ARRAY);
            }
        }

        // ---- Draw sensitivity bar below thumbnails ----
        {
            float barX = 8*dpiScale;
            float barY = (THUMB_PAD + NUM_ORNAMENTS*(THUMB_SIZE+THUMB_PAD) + (SIDEBAR_W-THUMB_SIZE)/2 + 10)*dpiScale;
            float barW = (SIDEBAR_W - 16)*dpiScale;
            float barH = 6*dpiScale;
            // Background
            glBegin(GL_QUADS);
            glColor3f(0.12f,0.12f,0.15f);
            glVertex2f(barX,barY);glVertex2f(barX+barW,barY);
            glVertex2f(barX+barW,barY+barH);glVertex2f(barX,barY+barH);
            glEnd();
            // Fill: log scale, sensitivity 0.1..10 maps to 0..1
            float fill=fminf(fmaxf((logf(sensitivity)-logf(0.1f))/(logf(10.f)-logf(0.1f)),0.f),1.f);
            glBegin(GL_QUADS);
            glColor3f(1.0f*fill, 0.8f*fill, 0.2f);
            glVertex2f(barX,barY);glVertex2f(barX+barW*fill,barY);
            glVertex2f(barX+barW*fill,barY+barH);glVertex2f(barX,barY+barH);
            glEnd();
            // Border
            glBegin(GL_LINE_LOOP);
            glColor3f(0.3f,0.3f,0.35f);
            glVertex2f(barX,barY);glVertex2f(barX+barW,barY);
            glVertex2f(barX+barW,barY+barH);glVertex2f(barX,barY+barH);
            glEnd();
            // Small triangle icon above bar (speaker-like icon)
            float iconY=barY-12*dpiScale;
            float iconCx=barX+barW*0.5f;
            float s3=3*dpiScale;
            // Sound wave arcs to hint "sensitivity"
            glColor3f(0.5f,0.5f,0.55f);
            glLineWidth(1.0f);
            for(int arc=0;arc<3;arc++){
                float r=(6+arc*5)*dpiScale;
                float alpha0=fill*0.8f+0.2f;
                if(arc>(int)(fill*3))alpha0*=0.3f;
                glColor4f(0.6f,0.55f,0.3f,alpha0);
                glBegin(GL_LINE_STRIP);
                for(int a=-30;a<=30;a+=5){
                    float rad=a*PI/180.f;
                    glVertex2f(iconCx+r*sinf(rad), iconY-r*cosf(rad)+8*dpiScale);
                }
                glEnd();
            }
        }

        // ---- Draw main ornament (cross-fade) ----
        float pulse[NUM_LAYERS] = {
            0.90f + sBass * 0.10f,
            0.80f + sMid  * 0.20f,
            0.75f + sBass * 0.25f,
            0.65f + sHigh * 0.35f,
            0.80f + sMid  * 0.20f,
            0.85f + sBass * 0.15f,
            0.55f + sHigh * 0.45f,
        };
        float lw[NUM_LAYERS] = { 2.2f, 1.8f, 1.3f, 1.0f, 2.0f, 2.0f, 1.1f };
        // Fill alpha per layer (brighter fills)
        float fillAlpha[NUM_LAYERS] = { 0.35f, 0.22f, 0.28f, 0.20f, 0.30f, 0.38f, 0.18f };

        // Smoothstep blend factor
        float blend = 0;
        if(morphing){
            float t=morphT;
            blend=t*t*(3-2*t);
        }
        float alphaOld = morphing ? (1.0f - blend) : 1.0f;
        float alphaNew = blend;

        // Lambda: draw one ornament's layers (fills + lines)
        auto drawOrnament = [&](LayerRender lrArr[], float alpha) {
            // Pass 1: filled polygons for every other contour (stencil even-odd)
            for(int i=0;i<NUM_LAYERS;i++){
                if(lrArr[i].verts.empty())continue;
                int nc=(int)lrArr[i].starts.size();
                if(nc==0)continue;

                glBindBuffer(GL_ARRAY_BUFFER,vbo);
                glBufferData(GL_ARRAY_BUFFER,lrArr[i].verts.size()*sizeof(Vec2),
                             lrArr[i].verts.data(),GL_STREAM_DRAW);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(2,GL_FLOAT,0,0);

                float p=pulse[i]*alpha;
                float fa=fillAlpha[i]*alpha;

                // Fill every other contour using stencil for correct non-convex fill
                for(int c=0;c<nc;c+=2){
                    GLint start=lrArr[i].starts[c];
                    GLsizei cnt=lrArr[i].counts[c];
                    if(cnt<3)continue;

                    // Step 1: write stencil with triangle fan (invert bits)
                    glEnable(GL_STENCIL_TEST);
                    glClear(GL_STENCIL_BUFFER_BIT);
                    glColorMask(GL_FALSE,GL_FALSE,GL_FALSE,GL_FALSE);
                    glStencilFunc(GL_ALWAYS,0,0xFF);
                    glStencilOp(GL_KEEP,GL_KEEP,GL_INVERT);
                    glDrawArrays(GL_TRIANGLE_FAN,start,cnt);

                    // Step 2: draw filled quad where stencil!=0
                    glColorMask(GL_TRUE,GL_TRUE,GL_TRUE,GL_TRUE);
                    glStencilFunc(GL_NOTEQUAL,0,0xFF);
                    glStencilOp(GL_KEEP,GL_KEEP,GL_KEEP);
                    glColor4f(baseColors[i].r*p*0.5f, baseColors[i].g*p*0.5f,
                              baseColors[i].b*p*0.5f, fa);
                    // Draw a screen-filling quad (the stencil clips it)
                    glBindBuffer(GL_ARRAY_BUFFER,0);
                    glDisableClientState(GL_VERTEX_ARRAY);
                    glBegin(GL_QUADS);
                    glVertex2f(0,0);glVertex2f((float)wW,0);
                    glVertex2f((float)wW,(float)wH);glVertex2f(0,(float)wH);
                    glEnd();

                    glDisable(GL_STENCIL_TEST);

                    // Re-bind for next contour
                    glBindBuffer(GL_ARRAY_BUFFER,vbo);
                    glEnableClientState(GL_VERTEX_ARRAY);
                    glVertexPointer(2,GL_FLOAT,0,0);
                }
                glDisableClientState(GL_VERTEX_ARRAY);
            }

            // Pass 2: line outlines for ALL contours
            for(int i=0;i<NUM_LAYERS;i++){
                if(lrArr[i].verts.empty())continue;
                glBindBuffer(GL_ARRAY_BUFFER,vbo);
                glBufferData(GL_ARRAY_BUFFER,lrArr[i].verts.size()*sizeof(Vec2),
                             lrArr[i].verts.data(),GL_STREAM_DRAW);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(2,GL_FLOAT,0,0);

                float p=pulse[i]*alpha;
                glColor4f(baseColors[i].r*p, baseColors[i].g*p, baseColors[i].b*p, alpha);
                glLineWidth(lw[i]);
                glMultiDrawArrays(GL_LINE_LOOP,lrArr[i].starts.data(),
                                  lrArr[i].counts.data(),(GLsizei)lrArr[i].starts.size());
                glDisableClientState(GL_VERTEX_ARRAY);
            }
        };

        // Draw current ornament
        drawOrnament(lr, alphaOld);

        // Draw target ornament fading in
        if(morphing)
            drawOrnament(lr2, alphaNew);

        glBindBuffer(GL_ARRAY_BUFFER,0);

        SDL_GL_SwapWindow(win);
    }

    if(g_au.dev>0)SDL_CloseAudioDevice(g_au.dev);
    SDL_GL_DeleteContext(gl);SDL_DestroyWindow(win);SDL_Quit();return 0;
}
