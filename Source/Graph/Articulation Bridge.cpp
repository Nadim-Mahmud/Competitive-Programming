
///         ****Articulation Bridge ****
///   O(E+V)
///   bridges are stored on artqb vector

#define MX 10005
vector<int>ed[MX], artqb[MX];
int vs[MX],dt[MX],low[MX],pr[MX],discovery_time,root;

void articulation_bridge(int n){
    int i,v;
    vs[n] = 1;
    dt[n] = low[n] = ++discovery_time;
    for(i=0;i<ed[n].size();i++){
        v = ed[n][i];
        if(pr[n]==v) continue; // won parent visit
        // if backedge
        if(vs[v]) low[n] = min(low[n],dt[v]);
        else{
            pr[v] = n;
            articulation_bridge(v);
            // storing min discovery discovery_time over all chaild
            low[n] = min(low[n],low[v]);
            if(dt[n]<low[v]){
                artqb[n].push_back(v);
                artqb[v].push_back(n);//for undirected.
                //u to v is articulation bridge here
            }
        }
    }
}

void free(){
    for(int i=0;i<MX;i++){
        ed[i].clear();
        artqb[i].clear();
        vs[i]=low[i]=dt[i]=0;
        mr[i]=stcol[i]=0;
        pr[i]=-1;
    }
    discovery_time = 0;
}
