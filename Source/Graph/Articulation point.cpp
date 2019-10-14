
///         **** Articulation points ****
///   O(E+V)
///   intiate root with starting position of tree
///   by checking marked points we will get articulation points


#define MX 10005
vector<int>ed[MX];
int vs[MX],dt[MX],points[MX],low[MX],pr[MX],discovery_time,root;

void articulation_point(int n){
    int i,v,cn=0;
    vs[n] = 1;
    dt[n] = low[n] = ++discovery_time;
    for(i=0;i<ed[n].size();i++){
        v = ed[n][i];
        if(pr[n]==v) continue; // won parent visit
        // if backedge
        if(vs[v]) low[n] = min(low[n],dt[v]);
        else{
            pr[v] = n;
            articulation_point(v);
            // storing min discovery discovery_time over all chaild
            low[n] = min(low[n],low[v]);
            if(dt[n]<=low[v]&&n!=root){
                //then n is a articulation point
                // v is a root of a devided subtree
                points[n] = 1;
            }
            cn++;
        }
    }
    //if root has more than one child then it is also a articulation_point
    if(cn>1&&n==root) points[n]=1;
}

void free(){
    for(int i=0;i<MX;i++){
        ed[i].clear();
        vs[i]=low[i]=dt[i]=points[i]=0;
        pr[i]=-1;
    }
    discovery_time = 0;
}
