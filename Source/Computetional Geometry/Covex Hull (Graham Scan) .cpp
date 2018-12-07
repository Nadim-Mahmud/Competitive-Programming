
///    *** Convex Hull[Grahams Scan]  [nlog(n)]

#define ll long long
struct point{
    ll x,y;
}convex_points[MX],points[MX];;

/// global scope decraltion of min-left point of collcetion
point pivot;

///Distance calculation mathod
ll dist(point p1,point p2){
    return ((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
}

/**
*   https://www.geeksforgeeks.org/orientation-3-ordered-points/
*   calculating orientation based on slope (yi-yj)/(xi-xj)
*   by compering slope of p-q and q-r;
*   if p-q<q-r then counter-clockwise
*   @return 0 if colinear
*   @return 1 if clockwise (Right rotetion)
*   @return 2 if counter-clockwise (left rotetion)
*/
int orientation(point p,point q,point r){
    ll val = ((q.y-p.y)*(r.x-q.x) - (r.y-q.y)*(q.x-p.x));
    if(val==0) return 0;
    if(val>0) return 1;
    return 2;
}

/**
*   sorting by polor angle in counterclockwise order around point0.
*   If polor angle of two points is same, then put the nearest point first.
*/
bool cmp(point p1,point p2){
    ll o = orientation(pivot,p1,p2);
    if(o==0){
        return dist(pivot,p1) < dist(pivot,p2);
    }
    return (o==2);
}

/// returning previous value of top element
inline point nextToTop(stack<point>&st){
    point p,res;
    p = st.top();
    st.pop();
    res = st.top();
    st.push(p);
    return res;
}

int total;

/**
*   This function will calculate convexHull points
*   All arrays are in 0 based indexing
*   @param n total numbers of points
*/
bool convexHull(int n){
    ll i,pos=0,in=0,miny = points[0].y,minx = points[0].x;
    stack<point>st;

    /// Finding bottom-left most point
    for(i=0;i<n;i++){
        if((miny==points[i].y&&minx>points[i].x)||miny>points[i].y){
            minx = points[i].x;
            miny = points[i].y;
            pos = i;
        }
    }

    ///sorting element according to the criteria
    swap(points[0],points[pos]);
    pivot = points[0];
    sort(points+1,points+n,cmp);

    ///Now removing same angle point
    for(i=1;i<n;i++){
        while(i<n-1&&orientation(pivot,points[i],points[i+1])==0){
            i++;
        }
        points[++in] = points[i];
    }
    if(in<2) return 0;

    st.push(points[0]);
    st.push(points[1]);
    st.push(points[2]);
    for(i=3;i<=in;i++){
        ///only valid sequence is ant-clockwise
        while(orientation(nextToTop(st),st.top(),points[i])!=2){
            st.pop();
        }
        st.push(points[i]);
    }

    in = total = st.size();
    point tmp;
    /// storing convex points
    while(!st.empty()){
        tmp = st.top();
        st.pop();
        convex_points[--in] = tmp;
    }
    return 1;
}

