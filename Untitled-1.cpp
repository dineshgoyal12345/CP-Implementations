#pragma GCC optimize("Ofast")
#pragma GCC target("avx,avx2,fma")
#include<bits/stdc++.h>
#include<ext/pb_ds/assoc_container.hpp>
#include<ext/pb_ds/tree_policy.hpp>
#define pi 3.141592653589793238
#define int long long
#define ll long long
#define ld long double
#define fi first
#define se second

using namespace __gnu_pbds;
using namespace std;
template <typename T>
using ordered_set = tree<T, null_type, less<T>, rb_tree_tag, tree_order_statistics_node_update>;
mt19937 rnd(chrono::high_resolution_clock::now().time_since_epoch().count());


long long powm(long long a, long long b,long long mod) {
long long res = 1;
while (b > 0) {
    if (b & 1)
        res = res * a %mod;
    a = a * a %mod;
    b >>= 1;
}

return res;
}

ll gcd(ll a, ll b)
{
if (b == 0)
    return a;
return gcd(b, a % b);      
}


namespace geometry{
    #ifndef M_PI
        #define M_PI acosl(-1)
    #endif
 
    typedef double T;
    struct pt {
        T x,y;
        pt operator+(pt p) {return {x+p.x, y+p.y};}
        pt operator-(pt p) {return {x-p.x, y-p.y};}
        pt operator*(T d) {return {x*d, y*d};}
        pt operator/(T d) {return {x/d, y/d};} // only for floating-point
    };
 
    bool operator==(pt a, pt b) {return a.x == b.x && a.y == b.y;}
    bool operator!=(pt a, pt b) {return !(a == b);}
 
    T sq(pt p) {return p.x*p.x + p.y*p.y;}
    double abs(pt p) {return sqrt(sq(p));}
 
 
    ostream& operator<<(ostream& os, pt p) {
        #ifdef DBG_LOCAL
            return os << "("<< p.x << "," << p.y << ")";
        #else
            return os<<p.x<<" "<<p.y;
        #endif
    }
 
    istream& operator>>(istream& os, pt& p) {
        return os>>p.x>>p.y;
    }
 
    template <typename T> int sgn(T x) {
        return (T(0) < x) - (x < T(0));
    }
 
    pt translate(pt v, pt p) {return p+v;}
 
    pt scale(pt c, double factor, pt p) {
        return c + (p-c)*factor;
    }
 
    pt rot(pt p, double a) {
        return {p.x*cos(a) - p.y*sin(a), p.x*sin(a) + p.y*cos(a)};
    }
 
    pt perp(pt p) {return {-p.y, p.x};}
 
    T dot(pt v, pt w) {return v.x*w.x + v.y*w.y;}
 
    T cross(pt v, pt w) {return v.x*w.y - v.y*w.x;}
 
    bool isPerp(pt v, pt w) {return dot(v,w) == 0;}
 
    pt linearTransfo(pt p, pt q, pt r, pt fp, pt fq) {
        pt pq = q-p, num{cross(pq, fq-fp), dot(pq, fq-fp)};
        return fp + pt{cross(r-p, num), dot(r-p, num)} / sq(pq);
    }
 
    double angle(pt v, pt w) {
        return acos(clamp(dot(v,w) / abs(v) / abs(w), (double)-1.0, (double)1.0));
    }
 
    T orient(pt a, pt b, pt c) {return cross(b-a,c-a);}
 
    bool inAngle(pt a, pt b, pt c, pt p) {
        assert(orient(a,b,c) != 0);
        if (orient(a,b,c) < 0) swap(b,c);
        return orient(a,b,p) >= 0 && orient(a,c,p) <= 0;
    }
 
 
    double orientedAngle(pt a, pt b, pt c) {
        if (orient(a,b,c) >= 0)
        return angle(b-a, c-a);
        else
        return 2*M_PI - angle(b-a, c-a);
    }
 
 
    bool isConvex(vector<pt> p) {
        bool hasPos=false, hasNeg=false;
        for (int i=0, n=p.size(); i<n; i++) {
        int o = orient(p[i], p[(i+1)%n], p[(i+2)%n]);
        if (o > 0) hasPos = true;
        if (o < 0) hasNeg = true;
        }
        return !(hasPos && hasNeg);
    }
    bool half(pt p, pt v = {0,1}) { // true if in blue half
        return cross(v,p) < 0 || (cross(v,p) == 0 && dot(v,p) <
0);
    }
 
    void polarSort(vector<pt> &v) {
        sort(v.begin(), v.end(), [](pt v, pt w) {
            return  make_tuple(half(v), 0, sq(v)) <
                   make_tuple(half(w), cross(v,w), sq(w));
        });
    }
 
    void polarSortAround(pt o, vector<pt> &v) {
        sort(v.begin(), v.end(), [o](pt v, pt w) {
            return  make_tuple(half(v-o), 0, sq(v-o)) <
              make_tuple(half(w-o), cross(v-o, w-o), sq(w-o));
        });
    }
 
 
    struct line {
        pt v; T c;
        // From direction vector v and offset c
        line(pt v, T c) : v(v), c(c) {}
        // From equation ax+by=c
        line(T a, T b, T c) : v({b,-a}), c(c) {}
        // From points P and Q
        line(pt p, pt q) : v(q-p), c(cross(v,p)) {}
        // Will be defined later:
        // - these work with T = int
        T side(pt p) {return cross(v,p)-c;}
        double dist(pt p) {return std::abs(side(p)) / abs(v);}
        double sqDist(pt p) {return side(p)*side(p) / (double)sq(v);}
        line perpThrough(pt p) {return {p, p + perp(v)};}
        bool cmpProj(pt p, pt q) {
            return dot(v,p) < dot(v,q);
        }
        line translate(pt t) {return {v, c + cross(v,t)};}
        // - these require T = double
        line shiftLeft(double dist) {return {v, c + dist*abs(v)};}
 
        pt proj(pt p) {return p - perp(v)*side(p)/sq(v);}
        pt refl(pt p) {return p - perp(v)*2*side(p)/sq(v);}
    };
 
    ostream& operator<<(ostream& os, line p) {
            return os<<-p.v.y<<" "<<p.v.x<<" "<<p.c<<endl;
    }
 
    bool inter(line l1, line l2, pt &out) {
        T d = cross(l1.v, l2.v);
        if (d == 0) return false;
        out = (l2.v*l1.c - l1.v*l2.c) / d; // requires floating-point coordinates
        return true;
    }
    line bisector(line l1, line l2, bool interior) {
        assert(cross(l1.v, l2.v) != 0); // l1 and l2 cannot be parallel!
        double sign = interior ? 1 : -1;
        return {l2.v/abs(l2.v) + l1.v/abs(l1.v) * sign,
                l2.c/abs(l2.v) + l1.c/abs(l1.v) * sign};
    }
    bool inDisk(pt a, pt b, pt p) {
        return dot(a-p, b-p) <= 0;
    }
    bool onSegment(pt a, pt b, pt p) {
        return orient(a,b,p) == 0 && inDisk(a,b,p);
    }
    bool properInter(pt a, pt b, pt c, pt d, pt &out) {
        double oa = orient(c,d,a),
        ob = orient(c,d,b),
        oc = orient(a,b,c),
        od = orient(a,b,d);
        // Proper intersection exists iff opposite signs
        if (oa*ob < 0 && oc*od < 0) {
            out = (a*ob - b*oa) / (ob-oa);
            return true;
        }
        return false;
    }
 
    struct cmpX {
        bool operator()(const pt& a, const pt& b) const {
            return make_pair(a.x, a.y) < make_pair(b.x, b.y);
        }
    };
    set<pt,cmpX> inters(pt a, pt b, pt c, pt d) {
        pt out;
        if (properInter(a,b,c,d,out)) return {out};
        set<pt,cmpX> s;
        if (onSegment(c,d,a)) s.insert(a);
        if (onSegment(c,d,b)) s.insert(b);
        if (onSegment(a,b,c)) s.insert(c);
        if (onSegment(a,b,d)) s.insert(d);
        return s;
    }
 
    double segPoint(pt a, pt b, pt p) {
        if (a != b) {
            line l(a,b);
            if (l.cmpProj(a,p) && l.cmpProj(p,b)) // if closest to projection
            return l.dist(p);
            // output distance to line
        }
        return min(abs(p-a), abs(p-b)); // otherwise distance to A or B
    }
    double segSeg(pt a, pt b, pt c, pt d) {
        pt dummy;
        if (properInter(a,b,c,d,dummy))
            return 0;
        return min({segPoint(a,b,c), segPoint(a,b,d),
                    segPoint(c,d,a), segPoint(c,d,b)});
    }
 
 
    double areaTriangle(pt a, pt b, pt c) {
        return std::abs(cross(b-a, c-a)) / 2.0;
    }
    double areaPolygon(const vector<pt>& p) {
        double area = 0.0;
        for (int i = 0, n = p.size(); i < n; i++) {
            area += cross(p[i], p[(i+1)%n]); // wrap back to 0 if i == n-1
        }
        return std::abs(area) / 2.0;
    }
    // true if P at least as high as A (blue part)
    bool above(pt a, pt p) {
        return p.y >= a.y;
    }
    // check if [PQ] crosses ray from A
    bool crossesRay(pt a, pt p, pt q) {
        return (above(a,q) - above(a,p)) * orient(a,p,q) > 0;
    }
    // if strict, returns false when A is on the boundary
    bool inPolygon(vector<pt> p, pt a, bool strict = true) {
        int numCrossings = 0;
        for (int i = 0, n = p.size(); i < n; i++) {
            if (onSegment(p[i], p[(i+1)%n], a))
                return !strict;
            numCrossings += crossesRay(a, p[i], p[(i+1)%n]);
        }
        return numCrossings & 1; // inside if odd number of crossings
    }
    // amplitude travelled around point A, from P to Q
    double angleTravelled(pt a, pt p, pt q) {
        double ampli = angle(p-a, q-a);
        if (orient(a,p,q) > 0) return ampli;
        else return -ampli;
    }
 
    int windingNumber(vector<pt> p, pt a) {
        double ampli = 0;
        for (int i = 0, n = p.size(); i < n; i++)
            ampli += angleTravelled(a, p[i], p[(i+1)%n]);
        return round(ampli / (2*M_PI));
    }
 
 
    //Circle
 
 
 
    pt circumCenter(pt a, pt b, pt c) {
        b = b-a, c = c-a; // consider coordinates relative to A
        assert(cross(b,c) != 0); // no circumcircle if A,B,C aligned
        return a + perp(b*sq(c) - c*sq(b))/cross(b,c)/2;
    }
 
 
    int circleLine(pt o, double r, line l, pair<pt,pt> &out) {
        double h2 = r*r - l.sqDist(o);
        if (h2 >= 0) { // the line touches the circle
            pt p = l.proj(o); // point P
            pt h = l.v*sqrt(h2)/abs(l.v); // vector parallel to l, of length h
            out = {p-h, p+h};
        }
        return 1 + sgn(h2);
    }
 
    int circleCircle(pt o1, double r1, pt o2, double r2, pair<pt,pt> &out) {
        pt d=o2-o1; double d2=sq(d);
        if (d2 == 0) {assert(r1 != r2); return 0;} // concentric circles
        double pd = (d2 + r1*r1 - r2*r2)/2; // = |O_1P| * d
        double h2 = r1*r1 - pd*pd/d2; // = hˆ2
        if (h2 >= 0) {
            pt p = o1 + d*pd/d2, h = perp(d)*sqrt(h2/d2);
            out = {p-h, p+h};
        }
        return 1 + sgn(h2);
    }
 
};
 
using namespace geometry;

double solve(pt a,pt b, pt c)
{
    line l(a,b);
    pt z=l.proj(c);
    if(dot(b-a,c-a)>=0)
        return l.dist(c);
    return abs(c-a);
}

signed main()
{

ios_base::sync_with_stdio(false);
cin.tie(NULL);
cout.tie(0);

#ifndef ONLINE_JUDGE
if(fopen("INPUT.txt","r"))
{
freopen("INPUT.txt","r",stdin);
freopen("OUTPUT.txt","w",stdout);
}
#endif

// freopen("raydist.in","r",stdin);
// freopen("raydist.out","w",stdout);
    
    
    pt a,b,c,d;
    cin>>a>>b>>c>>d;
    pt naya;
    line l1(a,b);
    line l2(c,d);

    bool f=inter(l1,l2,naya);

    if(dot(b-a,naya-a)<0||dot(d-c,naya-c)<0||f==0)
    {
        double mn=solve(a,b,c);
        mn=min(mn,solve(c,d,a));
        cout<<fixed<<setprecision(8)<<mn;
    }   
    else
        cout<<0; 

}