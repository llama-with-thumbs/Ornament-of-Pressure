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
static const int   THUMB_SIZE = 56;      // thumbnail pixel size
static const int   SIDEBAR_W  = 130;     // sidebar width in pixels
static const int   THUMB_PAD  = 24;      // padding between thumbnails (room for name text)

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
    { 1.00f, 1.00f, 1.00f },  // white
    { 0.85f, 0.85f, 0.85f },  // light gray
    { 0.95f, 0.95f, 0.95f },  // near white
    { 0.75f, 0.75f, 0.75f },  // medium gray
    { 0.90f, 0.90f, 0.90f },  // off white
    { 1.00f, 1.00f, 1.00f },  // white
    { 0.65f, 0.65f, 0.65f },  // dim gray
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
// Helpers for book-authentic constructions
// ============================================================

// Tiling helper: hex grid — calls fn(cx,cy) for each hex center within radius
template<typename F>
static void hexGrid(float sc,float cx,float cy,float maxR,F fn){
    float colStep=sc*SQRT3,rowStep=sc*1.5f;
    int cols=(int)ceilf(maxR*2/colStep)+3,rows=(int)ceilf(maxR*2/rowStep)+3;
    for(int r=-rows;r<=rows;r++)for(int c=-cols;c<=cols;c++){
        float hx=c*colStep+(r%2!=0?colStep*.5f:0)+cx;
        float hy=r*rowStep+cy;
        if(sqrtf((hx-cx)*(hx-cx)+(hy-cy)*(hy-cy))>maxR*1.2f)continue;
        fn(hx,hy);
    }
}

// Tiling helper: square grid — calls fn(cx,cy) for each cell center
template<typename F>
static void sqGrid(float sc,float cx,float cy,float maxR,float rot,F fn){
    int n=(int)ceilf(maxR*2.5f/sc)+1;float h=n*sc*.5f;
    float cR=cosf(rot),sR=sinf(rot);
    for(int i=0;i<n;i++)for(int j=0;j<n;j++){
        float lx=-h+(i+.5f)*sc,ly=-h+(j+.5f)*sc;
        float px=lx*cR-ly*sR+cx,py=lx*sR+ly*cR+cy;
        if(sqrtf((px-cx)*(px-cx)+(py-cy)*(py-cy))>maxR*1.2f)continue;
        fn(px,py,lx,ly,cR,sR);
    }
}

// Double triangle (Seal of Solomon) at hex midpoints — THE star-and-hexagon pattern (p.2-3)
static void starHexPattern(float sc,float cx,float cy,float maxR,std::vector<Vec2>&s){
    hexGrid(sc,cx,cy,maxR,[&](float hx,float hy){
        // Draw two overlapping equilateral triangles (Seal of Solomon)
        for(int t=0;t<2;t++){
            float rot0=t*PI/6;  // 0° and 30°
            for(int i=0;i<3;i++){
                float a1=2*PI*i/3+rot0-PI/6, a2=2*PI*(i+1)/3+rot0-PI/6;
                float r=sc*SQRT3/3;
                s.push_back({hx+r*cosf(a1),hy+r*sinf(a1)});
                s.push_back({hx+r*cosf(a2),hy+r*sinf(a2)});
            }
        }
    });
}

// Breath of the Compassionate pattern (p.8-9) — stars and crosses from square grid
static void breathPattern(float sc,float cx,float cy,float maxR,float breathe,std::vector<Vec2>&s){
    // breathe: 0=standard, positive=stars expand, negative=contract
    int n=(int)ceilf(maxR*2.5f/sc)+2;float h=n*sc*.5f;
    // For each square cell, draw the khaātam (double square = 8-pointed star)
    for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sc,sy=cy+j*sc;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>maxR*1.2f)continue;
        float r=sc*.5f*(1.0f+breathe*.15f);
        // Two overlapping squares at 45°
        for(int q=0;q<2;q++){
            float rot=q*PI/4;
            for(int k=0;k<4;k++){
                float a1=PI*k/2+PI/4+rot, a2=PI*(k+1)/2+PI/4+rot;
                s.push_back({sx+r*cosf(a1),sy+r*sinf(a1)});
                s.push_back({sx+r*cosf(a2),sy+r*sinf(a2)});
            }
        }
    }
}

// Rosette construction (p.10-11): octagonal rosette with petals
static void rosette8(float cx,float cy,float r,float rot,std::vector<Vec2>&s){
    // Central 8-pointed star
    float starR=r*.42f;
    for(int i=0;i<8;i++){
        float a1=2*PI*i/8+rot,a2=2*PI*((i+3)%8)/8+rot;
        s.push_back({cx+starR*cosf(a1),cy+starR*sinf(a1)});
        s.push_back({cx+starR*cosf(a2),cy+starR*sinf(a2)});
    }
    // Petals: kite shapes between star points and outer octagon
    for(int i=0;i<8;i++){
        float a=2*PI*i/8+rot;
        float aL=a-PI/8, aR=a+PI/8;
        // Outer point of petal
        Vec2 outer={cx+r*cosf(a),cy+r*sinf(a)};
        // Inner point (star tip)
        Vec2 inner={cx+starR*.7f*cosf(a),cy+starR*.7f*sinf(a)};
        // Side points
        float sideR=r*.68f;
        Vec2 left={cx+sideR*cosf(aL),cy+sideR*sinf(aL)};
        Vec2 right={cx+sideR*cosf(aR),cy+sideR*sinf(aR)};
        s.push_back(inner);s.push_back(left);
        s.push_back(left);s.push_back(outer);
        s.push_back(outer);s.push_back(right);
        s.push_back(right);s.push_back(inner);
    }
    // Outer octagon
    regPoly(cx,cy,r,8,rot+PI/8,s);
}

// Dodecagon-based 12-fold star (p.16-19): stars at midpoints of dodecagon edges
static void dodecaStar(float cx,float cy,float r,float rot,std::vector<Vec2>&s){
    // Regular dodecagon
    Vec2 dv[12];
    for(int i=0;i<12;i++){
        float a=2*PI*i/12+rot;
        dv[i]={cx+r*cosf(a),cy+r*sinf(a)};
    }
    // Dodecagon outline
    for(int i=0;i<12;i++){s.push_back(dv[i]);s.push_back(dv[(i+1)%12]);}
    // 12-pointed star from midpoints: two overlapping hexagons
    Vec2 mid[12];
    for(int i=0;i<12;i++){
        mid[i]={(dv[i].x+dv[(i+1)%12].x)*.5f,(dv[i].y+dv[(i+1)%12].y)*.5f};
    }
    // Connect midpoints to form the 12-fold star
    for(int i=0;i<12;i++){
        s.push_back(mid[i]);s.push_back(mid[(i+5)%12]);
    }
    // Inner 12-pointed star
    for(int i=0;i<12;i++){
        s.push_back(mid[i]);s.push_back(mid[(i+4)%12]);
    }
}

// Decagon-based 10-fold: Umm al-Girih (p.34-35)
static void ummAlGirih(float cx,float cy,float r,float rot,std::vector<Vec2>&s){
    // Regular decagon
    Vec2 dv[10];
    for(int i=0;i<10;i++){
        float a=2*PI*i/10+rot;
        dv[i]={cx+r*cosf(a),cy+r*sinf(a)};
    }
    // Decagon outline
    for(int i=0;i<10;i++){s.push_back(dv[i]);s.push_back(dv[(i+1)%10]);}
    // Stars from midpoints of decagon edges (the girih construction)
    Vec2 mid[10];
    for(int i=0;i<10;i++){
        mid[i]={(dv[i].x+dv[(i+1)%10].x)*.5f,(dv[i].y+dv[(i+1)%10].y)*.5f};
    }
    // 10-pointed star: connect midpoints skipping 3
    for(int i=0;i<10;i++){
        s.push_back(mid[i]);s.push_back(mid[(i+3)%10]);
    }
    // Inner pentagons at corners
    for(int i=0;i<10;i+=2){
        float a=2*PI*i/10+rot;
        float pr=r*.28f;
        regPoly(cx+r*.78f*cosf(a),cy+r*.78f*sinf(a),pr,5,a,s);
    }
}

// Khaātam (8-fold seal) for zillīj (p.26-29)
static void khaatam(float cx,float cy,float r,float rot,std::vector<Vec2>&s){
    // Octagram: connect every 3rd vertex of octagon ({8/3})
    for(int i=0;i<8;i++){
        float a1=2*PI*i/8+rot, a2=2*PI*((i+3)%8)/8+rot;
        s.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
        s.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});
    }
    // Inner octagon connecting the star intersections
    float ir=r*cosf(PI/8)/cosf(PI/8-PI*3/8)*cosf(PI/8); // intersection radius
    ir=r*0.541f; // simplified: cos(pi/8)*tan(pi/8) ratio
    regPoly(cx,cy,ir,8,rot+PI/8,s);
    // Outer octagon
    regPoly(cx,cy,r,8,rot+PI/8,s);
}

// ============================================================
// Ornament 0: Star & Hexagon — Seal of Solomon (p.2-5)
// ============================================================
static void generateOrnament_0(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.16f;
    // p.2-5: Star-and-hexagon (Seal of Solomon) on hex grid
    // L_STAR_BOUNDARY: the hex grid star-and-hexagon pattern
    starHexPattern(sc,cx,cy,R0*1.1f,layers[L_STAR_BOUNDARY]);
    // L_SPOKES: 6 radial axes
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<6;k++){float th=PI*k/3;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*cosf(th),cy+R0*sinf(th)});}}
    // L_GRID_MAIN: Hankin on hex grid (the substructure)
    hexHankin(sc,cx,cy,R0*1.1f,caDeg,layers[L_GRID_MAIN]);
    // L_GRID_FINE: finer triangle Hankin overlay
    triHankin(sc*.6f,cx,cy,R0*.7f,caDeg*.85f,layers[L_GRID_FINE]);
    // L_MODULES: concentric hexagonal rings (p.3 construction circles)
    { auto&L=layers[L_MODULES];
      for(int ring=1;ring<=4;ring++) regPoly(cx,cy,R0*ring*.24f,6,PI/6,L);
      // Seal of Solomon at center, larger
      starPoly(cx,cy,R0*.35f,6,2,0,L);
      starPoly(cx,cy,R0*.35f,6,2,PI/6,L); }
    // L_STARS8: {6/2} stars at hex vertices (p.3 bottom)
    { auto&L=layers[L_STARS8];
      hexGrid(sc*2,cx,cy,R0*.85f,[&](float hx,float hy){
        starPoly(hx,hy,sc*.45f,6,2,0,L);});
      starPoly(cx,cy,R0*.6f,12,5,0,L); }
    // L_ROTATING: rotating hex Hankin layer
    { std::vector<Vec2> raw;
      hexHankin(sc*1.5f,cx,cy,R0*.85f,caDeg*.9f,raw);
      float c=cosf(layerRot),s=sinf(layerRot);
      for(auto&v:raw){float dx=v.x-cx,dy=v.y-cy;v.x=dx*c-dy*s+cx;v.y=dx*s+dy*c+cy;}
      layers[L_ROTATING]=raw; }
}

// ============================================================
// Ornament 1: Breath of the Compassionate (p.8-9)
// ============================================================
static void generateOrnament_1(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.14f;
    // p.8-9: Stars and crosses from square grid (double square = khaātam)
    // L_STAR_BOUNDARY: the breath pattern — expanding/contracting double squares
    breathPattern(sc,cx,cy,R0*1.1f,starDepth,layers[L_STAR_BOUNDARY]);
    // L_SPOKES: 4+4 radial axes
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<8;k++){float th=PI*k/4;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*.95f*cosf(th),cy+R0*.95f*sinf(th)});}}
    // L_GRID_MAIN: Hankin on square grid at 0°
    sqHankin(sc,cx,cy,R0*1.1f,caDeg,0,layers[L_GRID_MAIN]);
    // L_GRID_FINE: Hankin on square grid at 45°
    sqHankin(sc,cx,cy,R0*1.1f,caDeg,PI/4,layers[L_GRID_FINE]);
    // L_MODULES: khaātam at key intersections (p.9 bottom)
    { auto&L=layers[L_MODULES];
      int n=(int)ceilf(R0/(sc*3))+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sc*3,sy=cy+j*sc*3;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0)continue;
        khaatam(sx,sy,sc*.9f,0,L);}}
    // L_STARS8: crosses between stars (p.9)
    { auto&L=layers[L_STARS8];
      int n=(int)ceilf(R0/(sc*3))+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+(i+.5f)*sc*3,sy=cy+(j+.5f)*sc*3;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0)continue;
        // Cross shape
        float cr=sc*.6f;
        for(int k=0;k<4;k++){float th=PI*k/2;
          L.push_back({sx,sy});
          L.push_back({sx+cr*cosf(th),sy+cr*sinf(th)});}
        regPoly(sx,sy,sc*.45f,4,PI/4,L);}}
    // L_ROTATING: rotating square Hankin
    { std::vector<Vec2> raw;
      sqHankin(sc*1.6f,cx,cy,R0*.85f,caDeg*.9f,layerRot,raw);
      layers[L_ROTATING]=raw; }
}

// ============================================================
// Ornament 2: Eight-Fold Rosettes (p.10-11)
// ============================================================
static void generateOrnament_2(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.18f;
    (void)starDepth;
    // p.10-11: Octagonal rosettes with petals on square grid
    // L_STAR_BOUNDARY: rosettes at grid points
    { auto&L=layers[L_STAR_BOUNDARY];
      int n=(int)ceilf(R0/(sc*2.5f))+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sc*2.5f,sy=cy+j*sc*2.5f;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0*1.1f)continue;
        rosette8(sx,sy,sc*1.1f,0,L);}}
    // L_SPOKES: 8 radial
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<8;k++){float th=PI*k/4;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*.95f*cosf(th),cy+R0*.95f*sinf(th)});}}
    // L_GRID_MAIN: Hankin 0° + 45°
    { std::vector<Vec2> raw;
      sqHankin(sc,cx,cy,R0*1.1f,caDeg,0,raw);
      sqHankin(sc,cx,cy,R0*1.1f,caDeg,PI/4,raw);
      layers[L_GRID_MAIN]=raw; }
    // L_GRID_FINE: fine detail at 22.5°
    sqHankin(sc*.5f,cx,cy,R0*.5f,caDeg*.9f,PI/8,layers[L_GRID_FINE]);
    // L_MODULES: octagonal frames at half-grid offsets
    { auto&L=layers[L_MODULES];
      int n=(int)ceilf(R0/(sc*2.5f))+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+(i+.5f)*sc*2.5f,sy=cy+(j+.5f)*sc*2.5f;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0)continue;
        regPoly(sx,sy,sc*.7f,8,PI/8,L);
        starPoly(sx,sy,sc*.5f,8,3,0,L);}}
    // L_STARS8: small {8/2} at rosette centers
    { auto&L=layers[L_STARS8];
      int n=(int)ceilf(R0/(sc*2.5f))+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sc*2.5f,sy=cy+j*sc*2.5f;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0*1.05f)continue;
        starPoly(sx,sy,sc*.35f,8,2,PI/8,L);}}
    // L_ROTATING
    { std::vector<Vec2> raw;
      sqHankin(sc*1.8f,cx,cy,R0*.8f,caDeg*.85f,layerRot,raw);
      layers[L_ROTATING]=raw; }
}

// ============================================================
// Ornament 3: Six of One / Three Times Four — 12-fold (p.16-19)
// ============================================================
static void generateOrnament_3(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.16f;
    (void)starDepth;
    // p.16-19: 12-fold stars from dodecagons + hexagons + squares
    // L_STAR_BOUNDARY: dodecagonal stars at hex grid points
    { auto&L=layers[L_STAR_BOUNDARY];
      hexGrid(sc*3.5f,cx,cy,R0*1.1f,[&](float hx,float hy){
        dodecaStar(hx,hy,sc*1.6f,0,L);});}
    // L_SPOKES: 12 radial axes
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<12;k++){float th=2*PI*k/12;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*.95f*cosf(th),cy+R0*.95f*sinf(th)});}}
    // L_GRID_MAIN: hex Hankin
    hexHankin(sc*1.2f,cx,cy,R0*1.1f,caDeg,layers[L_GRID_MAIN]);
    // L_GRID_FINE: triangle Hankin overlay
    triHankin(sc*.7f,cx,cy,R0*.7f,caDeg*.85f,layers[L_GRID_FINE]);
    // L_MODULES: hexagonal frames between dodecagons
    { auto&L=layers[L_MODULES];
      hexGrid(sc*3.5f,cx,cy,R0,[&](float hx,float hy){
        regPoly(hx,hy,sc*1.75f,6,PI/6,L);
        regPoly(hx,hy,sc*1.3f,6,0,L);
        // Small squares between dodecagons
        for(int k=0;k<6;k++){float th=PI*k/3+PI/6;
          float sx=hx+sc*2.4f*cosf(th),sy=hy+sc*2.4f*sinf(th);
          regPoly(sx,sy,sc*.4f,4,PI/4,L);}});
      // Central {12/5}
      starPoly(cx,cy,R0*.4f,12,5,0,L); }
    // L_STARS8: {6/2} between dodecagons
    { auto&L=layers[L_STARS8];
      hexGrid(sc*3.5f,cx,cy,R0*.9f,[&](float hx,float hy){
        for(int k=0;k<6;k++){float th=PI*k/3+PI/6;
          float sx=hx+sc*2.0f*cosf(th),sy=hy+sc*2.0f*sinf(th);
          starPoly(sx,sy,sc*.5f,6,2,th,L);}});
      starPoly(cx,cy,R0*.25f,12,5,PI/12,L); }
    // L_ROTATING
    { std::vector<Vec2> raw;
      hexHankin(sc*2.0f,cx,cy,R0*.85f,caDeg*.8f,raw);
      float c=cosf(layerRot),s=sinf(layerRot);
      for(auto&v:raw){float dx=v.x-cx,dy=v.y-cy;v.x=dx*c-dy*s+cx;v.y=dx*s+dy*c+cy;}
      layers[L_ROTATING]=raw; }
}

// ============================================================
// Ornament 4: Umm al-Girih — 10-fold (p.34-37)
// ============================================================
static void generateOrnament_4(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.15f;
    float R1=R0*starDepth;
    // p.34-37: "Mother of patterns" — decagons edge-to-edge with bowtie hexagons
    // L_STAR_BOUNDARY: decagonal stars (Umm al-Girih construction)
    { auto&L=layers[L_STAR_BOUNDARY];
      // Central large decagon star
      ummAlGirih(cx,cy,R0*.65f,0,L);
      // Ring of smaller decagons
      for(int k=0;k<10;k++){float th=2*PI*k/10;
        float d=R0*.78f;
        ummAlGirih(cx+d*cosf(th),cy+d*sinf(th),R0*.25f,th,L);}
      polarStar(cx,cy,R0,R1,10,500,L); }
    // L_SPOKES: 10 radial
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<10;k++){float th=2*PI*k/10;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*cosf(th),cy+R0*sinf(th)});}
      for(int k=0;k<10;k++){float th=2*PI*k/10+PI/10;
        L.push_back({cx+R0*.3f*cosf(th),cy+R0*.3f*sinf(th)});
        L.push_back({cx+R0*.6f*cosf(th),cy+R0*.6f*sinf(th)});}}
    // L_GRID_MAIN: 5 overlapping square grids at 36° (quasi-periodic)
    { std::vector<Vec2> raw;raw.reserve(80000);
      for(int g=0;g<5;g++) sqHankin(sc,cx,cy,R0*1.1f,caDeg,g*PI/5,raw);
      layers[L_GRID_MAIN]=raw; }
    // L_GRID_FINE: fine detail
    { std::vector<Vec2> raw;
      sqHankin(sc*.5f,cx,cy,R0*.45f,caDeg*.9f,PI/10,raw);
      triHankin(sc*.6f,cx,cy,R0*.5f,caDeg*.8f,raw);
      layers[L_GRID_FINE]=raw; }
    // L_MODULES: pentagons and bowties (p.34 shapes)
    { auto&L=layers[L_MODULES];
      for(int k=0;k<10;k++){float th=2*PI*k/10;
        float d=R0*.52f;
        regPoly(cx+d*cosf(th),cy+d*sinf(th),sc*.6f,5,th,L);}
      for(int k=0;k<10;k++){float th=2*PI*k/10+PI/10;
        float d=R0*.42f;
        regPoly(cx+d*cosf(th),cy+d*sinf(th),sc*.5f,10,th,L);}
      // {10/3} at center
      starPoly(cx,cy,R0*.3f,10,3,0,L); }
    // L_STARS8: {10/4} pentagrammaton stars (p.36-37)
    { auto&L=layers[L_STARS8];
      starPoly(cx,cy,R0*.55f,10,4,0,L);
      for(int k=0;k<5;k++){float th=2*PI*k/5;
        float d=R0*.65f;
        starPoly(cx+d*cosf(th),cy+d*sinf(th),sc*.55f,5,2,th,L);}
      for(int k=0;k<10;k++){float th=2*PI*k/10+PI/10;
        float d=R0*.48f;
        starPoly(cx+d*cosf(th),cy+d*sinf(th),sc*.35f,10,3,th,L);}}
    // L_ROTATING
    { std::vector<Vec2> raw;
      sqHankin(sc*1.5f,cx,cy,R0*.75f,caDeg*.85f,layerRot,raw);
      layers[L_ROTATING]=raw; }
}

// ============================================================
// Ornament 5: Zillīj — 8-fold khaātam compositions (p.26-29)
// ============================================================
static void generateOrnament_5(float w, float h, float R0, float starDepth,
    float caDeg, float layerRot, std::vector<Vec2> layers[NUM_LAYERS])
{
    float cx=w*.5f, cy=h*.5f;
    float sc=R0*0.15f;
    (void)starDepth;
    // p.26-29: Zillīj 8-fold with khaātam (octagram seal) and saft (hexagon)
    // L_STAR_BOUNDARY: grid of khaātam + saft pattern
    { auto&L=layers[L_STAR_BOUNDARY];
      float sp=sc*2.4f; // spacing for khaātam grid
      int n=(int)ceilf(R0*1.1f/sp)+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sp,sy=cy+j*sp;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0*1.1f)continue;
        khaatam(sx,sy,sc*.95f,0,L);}
      // Saft (elongated hexagons) between khaātam
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+(i+.5f)*sp,sy=cy+j*sp;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0)continue;
        // Vertical saft
        float sh=sc*.55f,sw=sc*.25f;
        L.push_back({sx-sw,sy-sh});L.push_back({sx+sw,sy-sh});
        L.push_back({sx+sw,sy-sh});L.push_back({sx+sw,sy+sh});
        L.push_back({sx+sw,sy+sh});L.push_back({sx-sw,sy+sh});
        L.push_back({sx-sw,sy+sh});L.push_back({sx-sw,sy-sh});
        // Horizontal saft
        float sx2=cx+i*sp,sy2=cy+(j+.5f)*sp;
        if(sqrtf((sx2-cx)*(sx2-cx)+(sy2-cy)*(sy2-cy))>R0)continue;
        L.push_back({sx2-sh,sy2-sw});L.push_back({sx2+sh,sy2-sw});
        L.push_back({sx2+sh,sy2-sw});L.push_back({sx2+sh,sy2+sw});
        L.push_back({sx2+sh,sy2+sw});L.push_back({sx2-sh,sy2+sw});
        L.push_back({sx2-sh,sy2+sw});L.push_back({sx2-sh,sy2-sw});}}
    // L_SPOKES: 8 radial
    { auto&L=layers[L_SPOKES];
      for(int k=0;k<8;k++){float th=PI*k/4;
        L.push_back({cx+sc*cosf(th),cy+sc*sinf(th)});
        L.push_back({cx+R0*.95f*cosf(th),cy+R0*.95f*sinf(th)});}}
    // L_GRID_MAIN: Hankin 0°
    sqHankin(sc,cx,cy,R0*1.1f,caDeg,0,layers[L_GRID_MAIN]);
    // L_GRID_FINE: Hankin 45°
    sqHankin(sc,cx,cy,R0*1.1f,caDeg,PI/4,layers[L_GRID_FINE]);
    // L_MODULES: larger framing octagons (p.27 compositions)
    { auto&L=layers[L_MODULES];
      float sp=sc*2.4f;
      int n=(int)ceilf(R0/sp)+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+i*sp,sy=cy+j*sp;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0*1.05f)continue;
        regPoly(sx,sy,sc*1.15f,8,PI/8,L);}}
    // L_STARS8: {8/2} stars between khaātam
    { auto&L=layers[L_STARS8];
      float sp=sc*2.4f;
      int n=(int)ceilf(R0/sp)+1;
      for(int i=-n;i<=n;i++)for(int j=-n;j<=n;j++){
        float sx=cx+(i+.5f)*sp,sy=cy+(j+.5f)*sp;
        if(sqrtf((sx-cx)*(sx-cx)+(sy-cy)*(sy-cy))>R0)continue;
        starPoly(sx,sy,sc*.55f,8,2,PI/8,L);
        regPoly(sx,sy,sc*.65f,4,0,L);}}
    // L_ROTATING
    { std::vector<Vec2> raw;
      sqHankin(sc*1.8f,cx,cy,R0*.85f,caDeg*.9f,layerRot,raw);
      layers[L_ROTATING]=raw; }
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
    "Star & Hexagon", "Breath of Compassionate", "Eight-Fold Rosettes",
    "Six of One (12-fold)", "Umm al-Girih (10-fold)", "Zillij (8-fold)"
};

// ============================================================
// Thumbnail pictograms matching book patterns
// ============================================================
static void generateThumb(int idx, float sz, std::vector<Vec2>& out) {
    float cx=sz*.5f, cy=sz*.5f, r=sz*.38f;
    switch(idx) {
    case 0: // Star & Hexagon: hex + double triangle (Seal of Solomon)
        for(int i=0;i<6;i++){float a1=PI*i/3+PI/6,a2=PI*(i+1)/3+PI/6;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int t=0;t<2;t++){float rot=t*PI/6;
            for(int i=0;i<3;i++){float a1=2*PI*i/3+rot-PI/6,a2=2*PI*(i+1)/3+rot-PI/6;float rr=r*.65f;
                out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
                out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}}
        break;
    case 1: // Breath: two overlapping squares (khaātam)
        for(int q=0;q<2;q++){float rot=q*PI/4;
            for(int i=0;i<4;i++){float a1=PI*i/2+PI/4+rot,a2=PI*(i+1)/2+PI/4+rot;
                out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
                out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}}
        // Small cross
        for(int k=0;k<4;k++){float th=PI*k/2;float rr=r*.3f;
            out.push_back({cx,cy});out.push_back({cx+rr*cosf(th),cy+rr*sinf(th)});}
        break;
    case 2: // Rosette: octagon + petals
        for(int i=0;i<8;i++){float a1=2*PI*i/8,a2=2*PI*(i+1)/8;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<8;i++){float a=2*PI*i/8;float rr=r*.5f;
            float aL=a-PI/8,aR=a+PI/8;
            out.push_back({cx+rr*.5f*cosf(a),cy+rr*.5f*sinf(a)});
            out.push_back({cx+r*.75f*cosf(aL),cy+r*.75f*sinf(aL)});
            out.push_back({cx+r*.75f*cosf(aL),cy+r*.75f*sinf(aL)});
            out.push_back({cx+r*cosf(a),cy+r*sinf(a)});
            out.push_back({cx+r*cosf(a),cy+r*sinf(a)});
            out.push_back({cx+r*.75f*cosf(aR),cy+r*.75f*sinf(aR)});
            out.push_back({cx+r*.75f*cosf(aR),cy+r*.75f*sinf(aR)});
            out.push_back({cx+rr*.5f*cosf(a),cy+rr*.5f*sinf(a)});}
        break;
    case 3: // 12-fold: dodecagon + inner star
        for(int i=0;i<12;i++){float a1=2*PI*i/12,a2=2*PI*(i+1)/12;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<12;i++){float a1=2*PI*i/12,a2=2*PI*((i+5)%12)/12;float rr=r*.65f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 4: // Umm al-Girih: decagon + {10/3} star
        for(int i=0;i<10;i++){float a1=2*PI*i/10,a2=2*PI*(i+1)/10;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<10;i++){float a1=2*PI*i/10,a2=2*PI*((i+3)%10)/10;float rr=r*.6f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        for(int i=0;i<5;i++){float a1=2*PI*i/5+PI/10,a2=2*PI*(i+1)/5+PI/10;float rr=r*.3f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
        break;
    case 5: // Zillīj: khaātam ({8/3} + octagon)
        for(int i=0;i<8;i++){float a1=2*PI*i/8,a2=2*PI*((i+3)%8)/8;
            out.push_back({cx+r*cosf(a1),cy+r*sinf(a1)});
            out.push_back({cx+r*cosf(a2),cy+r*sinf(a2)});}
        for(int i=0;i<8;i++){float a1=2*PI*i/8+PI/8,a2=2*PI*(i+1)/8+PI/8;float rr=r*.75f;
            out.push_back({cx+rr*cosf(a1),cy+rr*sinf(a1)});
            out.push_back({cx+rr*cosf(a2),cy+rr*sinf(a2)});}
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

// Per-layer rendering data (at file scope for morphLayers)
struct LayerRender {
    std::vector<Vec2> verts;
    std::vector<GLint> starts;
    std::vector<GLsizei> counts;
};

// ============================================================
// Minimal 3x5 stroke font — renders text as GL_LINES
// Each glyph: pairs of (x1,y1,x2,y2) normalized to 0..3 x 0..5
// ============================================================
static void drawStrokeChar(float x, float y, float s, char ch) {
    // Simple segment data for common chars (uppercase + digits + &-/)
    struct Seg { float x1,y1,x2,y2; };
    static const Seg empty[] = {{0,0,0,0}};
    #define G(name) static const Seg name
    G(A)[]={{0,5,0,1},{0,1,1,0},{1,0,2,0},{2,0,3,1},{3,1,3,5},{0,3,3,3}};
    G(B)[]={{0,0,0,5},{0,0,2,0},{2,0,3,1},{3,1,2,2.5f},{2,2.5f,0,2.5f},{2,2.5f,3,3.5f},{3,3.5f,2,5},{2,5,0,5}};
    G(C)[]={{3,0,1,0},{1,0,0,1},{0,1,0,4},{0,4,1,5},{1,5,3,5}};
    G(D)[]={{0,0,0,5},{0,0,2,0},{2,0,3,1},{3,1,3,4},{3,4,2,5},{2,5,0,5}};
    G(E)[]={{3,0,0,0},{0,0,0,5},{0,5,3,5},{0,2.5f,2,2.5f}};
    G(F)[]={{3,0,0,0},{0,0,0,5},{0,2.5f,2,2.5f}};
    G(GG)[]={{3,1,2,0},{2,0,1,0},{1,0,0,1},{0,1,0,4},{0,4,1,5},{1,5,3,5},{3,5,3,3},{3,3,2,3}};
    G(H)[]={{0,0,0,5},{3,0,3,5},{0,2.5f,3,2.5f}};
    G(I)[]={{1,0,2,0},{1.5f,0,1.5f,5},{1,5,2,5}};
    G(K)[]={{0,0,0,5},{3,0,0,2.5f},{0,2.5f,3,5}};
    G(L)[]={{0,0,0,5},{0,5,3,5}};
    G(M)[]={{0,5,0,0},{0,0,1.5f,2},{1.5f,2,3,0},{3,0,3,5}};
    G(N)[]={{0,5,0,0},{0,0,3,5},{3,5,3,0}};
    G(O)[]={{1,0,2,0},{2,0,3,1},{3,1,3,4},{3,4,2,5},{2,5,1,5},{1,5,0,4},{0,4,0,1},{0,1,1,0}};
    G(P)[]={{0,0,0,5},{0,0,2,0},{2,0,3,1},{3,1,3,2},{3,2,2,2.5f},{2,2.5f,0,2.5f}};
    G(R)[]={{0,0,0,5},{0,0,2,0},{2,0,3,1},{3,1,2,2.5f},{2,2.5f,0,2.5f},{1.5f,2.5f,3,5}};
    G(S)[]={{3,0,1,0},{1,0,0,1},{0,1,0,2},{0,2,1,2.5f},{1,2.5f,2,2.5f},{2,2.5f,3,3.5f},{3,3.5f,3,4},{3,4,2,5},{2,5,0,5}};
    G(T)[]={{0,0,3,0},{1.5f,0,1.5f,5}};
    G(U)[]={{0,0,0,4},{0,4,1,5},{1,5,2,5},{2,5,3,4},{3,4,3,0}};
    G(X)[]={{0,0,3,5},{3,0,0,5}};
    G(Z)[]={{0,0,3,0},{3,0,0,5},{0,5,3,5}};
    G(n0)[]={{1,0,2,0},{2,0,3,1},{3,1,3,4},{3,4,2,5},{2,5,1,5},{1,5,0,4},{0,4,0,1},{0,1,1,0}};
    G(n1)[]={{1,1,1.5f,0},{1.5f,0,1.5f,5},{0.5f,5,2.5f,5}};
    G(n2)[]={{0,1,1,0},{1,0,2,0},{2,0,3,1},{3,1,3,2},{3,2,0,5},{0,5,3,5}};
    G(n5)[]={{3,0,0,0},{0,0,0,2.5f},{0,2.5f,2,2.5f},{2,2.5f,3,3.5f},{3,3.5f,2,5},{2,5,0,5}};
    G(n6)[]={{2,0,1,0},{1,0,0,1},{0,1,0,4},{0,4,1,5},{1,5,2,5},{2,5,3,4},{3,4,2,2.5f},{2,2.5f,0,2.5f}};
    G(n8)[]={{1,0,2,0},{2,0,3,1},{3,1,2,2.5f},{2,2.5f,1,2.5f},{1,2.5f,0,1},{0,1,1,0},{1,2.5f,0,3.5f},{0,3.5f,1,5},{1,5,2,5},{2,5,3,3.5f},{3,3.5f,2,2.5f}};
    G(AMP)[]={{2,1,1,0},{1,0,0,1},{0,1,1,2.5f},{1,2.5f,0,4},{0,4,1,5},{1,5,3,3},{3,3,3,4},{3,4,3,5}};
    G(DASH)[]={{0,2.5f,3,2.5f}};
    G(SPC)[]={{0,0,0,0}};
    #undef G

    const Seg* segs=empty; int ns=0;
    #define MAP(c,arr) case c: segs=arr;ns=sizeof(arr)/sizeof(Seg);break;
    switch(ch>='a'&&ch<='z'?ch-32:ch){
        MAP('A',A) MAP('B',B) MAP('C',C) MAP('D',D) MAP('E',E) MAP('F',F)
        MAP('G',GG) MAP('H',H) MAP('I',I) MAP('K',K) MAP('L',L) MAP('M',M)
        MAP('N',N) MAP('O',O) MAP('P',P) MAP('R',R) MAP('S',S) MAP('T',T)
        MAP('U',U) MAP('X',X) MAP('Z',Z)
        MAP('0',n0) MAP('1',n1) MAP('2',n2) MAP('5',n5) MAP('6',n6) MAP('8',n8)
        MAP('&',AMP) MAP('-',DASH) MAP(' ',SPC)
        default: break;
    }
    #undef MAP
    for(int i=0;i<ns;i++){
        glVertex2f(x+segs[i].x1*s, y+segs[i].y1*s);
        glVertex2f(x+segs[i].x2*s, y+segs[i].y2*s);
    }
}

static void drawStrokeText(float x, float y, float scale, const char* text, float r, float g, float b, float a) {
    glColor4f(r,g,b,a);
    glLineWidth(1.0f);
    glBegin(GL_LINES);
    float cx=x;
    for(const char*p=text;*p;p++){
        drawStrokeChar(cx,y,scale,*p);
        cx+=4*scale; // char width + spacing
    }
    glEnd();
}

// Resample a closed chain to exactly N evenly-spaced vertices
static std::vector<Vec2> resampleChain(const Vec2* pts, int n, int targetN) {
    if(n<2||targetN<2) return std::vector<Vec2>(targetN,pts[0]);
    // Compute cumulative arc lengths
    std::vector<float> cum(n+1,0);
    for(int i=1;i<=n;i++){
        float dx=pts[i%n].x-pts[(i-1)%n].x, dy=pts[i%n].y-pts[(i-1)%n].y;
        cum[i]=cum[i-1]+sqrtf(dx*dx+dy*dy);
    }
    float totalLen=cum[n];
    if(totalLen<1e-6f) return std::vector<Vec2>(targetN,pts[0]);
    std::vector<Vec2> out(targetN);
    int seg=0;
    for(int i=0;i<targetN;i++){
        float t=(float)i/(float)targetN*totalLen;
        while(seg<n-1 && cum[seg+1]<t) seg++;
        float segLen=cum[seg+1]-cum[seg];
        float frac=(segLen>1e-8f)?(t-cum[seg])/segLen:0;
        out[i]={pts[seg%n].x+(pts[(seg+1)%n].x-pts[seg%n].x)*frac,
                pts[seg%n].y+(pts[(seg+1)%n].y-pts[seg%n].y)*frac};
    }
    return out;
}

// Interpolate two LayerRender arrays into a merged result
static void morphLayers(LayerRender src[], LayerRender dst[], LayerRender out[],
                        int numLayers, float t, float cx, float cy) {
    float s = t*t*(3-2*t); // smoothstep
    for(int i=0;i<numLayers;i++){
        out[i].verts.clear(); out[i].starts.clear(); out[i].counts.clear();
        int nSrc=(int)src[i].starts.size(), nDst=(int)dst[i].starts.size();
        int nMax=nSrc>nDst?nSrc:nDst;
        for(int c=0;c<nMax;c++){
            bool hasSrc=(c<nSrc), hasDst=(c<nDst);
            if(hasSrc && hasDst){
                // Both exist — resample to same count and lerp
                int cntS=src[i].counts[c], cntD=dst[i].counts[c];
                int target=cntS>cntD?cntS:cntD;
                if(target<3) target=3;
                auto rS=resampleChain(&src[i].verts[src[i].starts[c]],cntS,target);
                auto rD=resampleChain(&dst[i].verts[dst[i].starts[c]],cntD,target);
                GLint start=(GLint)out[i].verts.size();
                for(int v=0;v<target;v++){
                    out[i].verts.push_back({rS[v].x*(1-s)+rD[v].x*s,
                                            rS[v].y*(1-s)+rD[v].y*s});
                }
                out[i].starts.push_back(start);
                out[i].counts.push_back(target);
            } else if(hasSrc) {
                // Source only — shrink toward center
                GLint start=(GLint)out[i].verts.size();
                int cnt=src[i].counts[c];
                const Vec2* p=&src[i].verts[src[i].starts[c]];
                for(int v=0;v<cnt;v++){
                    out[i].verts.push_back({p[v].x*(1-s)+cx*s, p[v].y*(1-s)+cy*s});
                }
                out[i].starts.push_back(start);
                out[i].counts.push_back(cnt);
            } else {
                // Dest only — grow from center
                GLint start=(GLint)out[i].verts.size();
                int cnt=dst[i].counts[c];
                const Vec2* p=&dst[i].verts[dst[i].starts[c]];
                for(int v=0;v<cnt;v++){
                    out[i].verts.push_back({cx*(1-s)+p[v].x*s, cy*(1-s)+p[v].y*s});
                }
                out[i].starts.push_back(start);
                out[i].counts.push_back(cnt);
            }
        }
    }
}

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
        SDL_WINDOWPOS_CENTERED,SDL_WINDOWPOS_CENTERED,700,700,
        SDL_WINDOW_OPENGL|SDL_WINDOW_RESIZABLE|SDL_WINDOW_ALLOW_HIGHDPI);
    SDL_GLContext gl=SDL_GL_CreateContext(win);
    glewInit();SDL_GL_SetSwapInterval(1);

    GLuint vbo;glGenBuffers(1,&vbo);

    SDL_AudioSpec want={},have={};
    want.freq=44100;want.format=AUDIO_F32SYS;want.channels=1;want.samples=1024;want.callback=audioCB;
    g_au.dev=SDL_OpenAudioDevice(NULL,1,&want,&have,SDL_AUDIO_ALLOW_FORMAT_CHANGE);
    if(g_au.dev>0)SDL_PauseAudioDevice(g_au.dev,0);

    float sBass=.5f,sMid=.5f,sHigh=0,layerRot=0;
    float sensitivity=4.0f;  // audio sensitivity: 0.1 .. 10.0, controlled by scroll/arrow keys
    Complex fftBuf[FFT_N];float hann[FFT_N];
    for(int i=0;i<FFT_N;i++)hann[i]=.5f*(1-cosf(2*PI*i/(FFT_N-1)));

    std::vector<Vec2> layerSegs[NUM_LAYERS];
    LayerRender lr[NUM_LAYERS];
    ChainBuilder cb;

    // Ornament selection and morphing state
    int curOrnament = 0;
    int targetOrnament = -1;  // -1 = no transition
    float morphT = 0;         // 0..1: 0=fully old, 1=fully new
    static const float MORPH_SPEED = 0.018f;  // per frame (~3.3s transition)

    // Second set of layer data for morph target
    std::vector<Vec2> layerSegs2[NUM_LAYERS];
    LayerRender lr2[NUM_LAYERS];
    LayerRender lrMorph[NUM_LAYERS]; // interpolated result

    // Reveal animation — ornament paints outward from center
    float revealRadius = 0;           // current reveal radius in pixels
    static const float REVEAL_SPEED = 8.0f; // pixels per frame
    bool revealing = true;            // true while reveal is growing

    // Pre-generate thumbnail geometry
    std::vector<Vec2> thumbGeom[NUM_ORNAMENTS];
    for (int i = 0; i < NUM_ORNAMENTS; i++)
        generateThumb(i, (float)THUMB_SIZE, thumbGeom[i]);

    // Set initial window title
    { std::string t = std::string("Ornament of Pressure — ") + ornamentNames[0];
      SDL_SetWindowTitle(win, t.c_str()); }

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
                revealRadius=0; revealing=true; // restart reveal for new ornament
                // Update window title with pattern name
                std::string title = std::string("Ornament of Pressure — ") + ornamentNames[curOrnament];
                SDL_SetWindowTitle(win, title.c_str());
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

        float R0=306; // fixed size — audio drives geometry, not scale
        float depth=.08f+sBass*.25f;    // bass → star depth / breath amount
        float angle=20+sMid*50;          // mids → contact angle (pattern openness)
        layerRot+=sHigh*.012f+.001f;     // highs → rotation speed

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
            // Draw pattern name below thumbnail
            {
                static const char* shortNames[NUM_ORNAMENTS] = {
                    "Star & Hex", "Breath", "8-Fold", "Six of One", "Umm al-Girih", "Zillij"
                };
                float nameScale = 3.0f * dpiScale;
                const char* nm = shortNames[ti];
                int nchars = (int)strlen(nm);
                // Each char advances 4*scale; last char is 3*scale wide
                float tw = nchars > 0 ? (nchars - 1) * 4.0f * nameScale + 3.0f * nameScale : 0;
                float nx = tx + tsz * 0.5f - tw * 0.5f;
                float ny = ty + tsz + 3 * dpiScale;
                float tbr = selected ? 0.9f : (isTarget ? 0.6f : 0.35f);
                drawStrokeText(nx, ny, nameScale, nm, tbr, tbr * 0.85f, tbr * 0.3f, 1.0f);
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
        float lw[NUM_LAYERS] = { 2.0f, 1.5f, 1.2f, 0.9f, 1.8f, 1.8f, 1.0f };
        // Fill alpha per layer — heavy white fills
        float fillAlpha[NUM_LAYERS] = { 0.55f, 0.40f, 0.48f, 0.35f, 0.50f, 0.55f, 0.30f };

        // Advance reveal radius
        float maxReveal = sqrtf((float)(wW*wW+wH*wH)); // diagonal
        if(revealing){
            revealRadius += REVEAL_SPEED;
            if(revealRadius >= maxReveal) { revealRadius = maxReveal; revealing = false; }
        }

        float rcx = sidebarPx + (wW - sidebarPx) * 0.5f;
        float rcy = wH * 0.5f;

        // If morphing, interpolate geometry; otherwise just draw current
        LayerRender* drawLR = lr;
        if(morphing){
            float ocx = sidebarPx + (wW - sidebarPx) * 0.5f;
            float ocy = wH * 0.5f;
            morphLayers(lr, lr2, lrMorph, NUM_LAYERS, morphT, ocx, ocy);
            drawLR = lrMorph;
        }

        // Helper: compute centroid distance from reveal center for a contour
        auto contourDist = [&](const Vec2* verts, GLint start, GLsizei cnt) -> float {
            float cx=0,cy=0;
            for(int v=0;v<cnt;v++){cx+=verts[start+v].x;cy+=verts[start+v].y;}
            cx/=cnt; cy/=cnt;
            float dx=cx-rcx, dy=cy-rcy;
            return sqrtf(dx*dx+dy*dy);
        };

        // Reveal fade band width (pixels)
        float fadeBand = 40.0f;

        // Draw ornament (fills + lines)
        {
            // Pass 1: filled polygons for every other contour (stencil even-odd)
            for(int i=0;i<NUM_LAYERS;i++){
                if(drawLR[i].verts.empty())continue;
                int nc=(int)drawLR[i].starts.size();
                if(nc==0)continue;

                glBindBuffer(GL_ARRAY_BUFFER,vbo);
                glBufferData(GL_ARRAY_BUFFER,drawLR[i].verts.size()*sizeof(Vec2),
                             drawLR[i].verts.data(),GL_STREAM_DRAW);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(2,GL_FLOAT,0,0);

                float p=pulse[i];
                float fa=fillAlpha[i];

                // Fill every other contour using stencil for correct non-convex fill
                for(int c=0;c<nc;c+=2){
                    GLint start=drawLR[i].starts[c];
                    GLsizei cnt=drawLR[i].counts[c];
                    if(cnt<3)continue;

                    // Reveal: skip contours beyond reveal radius
                    float dist=contourDist(drawLR[i].verts.data(),start,cnt);
                    if(dist > revealRadius) continue;
                    float revealAlpha = (revealRadius-dist < fadeBand) ?
                        (revealRadius-dist)/fadeBand : 1.0f;

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
                    glColor4f(baseColors[i].r*p*0.8f, baseColors[i].g*p*0.8f,
                              baseColors[i].b*p*0.8f, fa*revealAlpha);
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

            // Pass 2: line outlines for ALL contours (with reveal)
            for(int i=0;i<NUM_LAYERS;i++){
                if(drawLR[i].verts.empty())continue;
                int nc=(int)drawLR[i].starts.size();

                glBindBuffer(GL_ARRAY_BUFFER,vbo);
                glBufferData(GL_ARRAY_BUFFER,drawLR[i].verts.size()*sizeof(Vec2),
                             drawLR[i].verts.data(),GL_STREAM_DRAW);
                glEnableClientState(GL_VERTEX_ARRAY);
                glVertexPointer(2,GL_FLOAT,0,0);

                float p=pulse[i];
                glLineWidth(lw[i]);
                // Draw each contour individually with reveal alpha
                for(int c=0;c<nc;c++){
                    GLint start=drawLR[i].starts[c];
                    GLsizei cnt=drawLR[i].counts[c];
                    float dist=contourDist(drawLR[i].verts.data(),start,cnt);
                    if(dist > revealRadius) continue;
                    float revealAlpha = (revealRadius-dist < fadeBand) ?
                        (revealRadius-dist)/fadeBand : 1.0f;
                    glColor4f(baseColors[i].r*p, baseColors[i].g*p, baseColors[i].b*p, revealAlpha);
                    glDrawArrays(GL_LINE_LOOP,start,cnt);
                }
                glDisableClientState(GL_VERTEX_ARRAY);
            }
        }

        glBindBuffer(GL_ARRAY_BUFFER,0);

        SDL_GL_SwapWindow(win);
    }

    if(g_au.dev>0)SDL_CloseAudioDevice(g_au.dev);
    SDL_GL_DeleteContext(gl);SDL_DestroyWindow(win);SDL_Quit();return 0;
}
