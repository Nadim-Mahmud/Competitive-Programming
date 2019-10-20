/**
    Description : Find maximum matching. we can find bipartite graph by left or right array
                  lt = left, rt =right, ed = adjesency matrix
    Time Complexity : O(n^2*m)
*/


#define MX 55
vector<int>ed[MX];
int lt[MX],rt[MX];
bool vs[MX];

bool bpm(int u){
    int i,v;
    for(i=0;i<ed[u].size();i++){
        v = ed[u][i];
        if(vs[v]) continue;
        vs[v] = 1;
        if(rt[v]==-1||bpm(rt[v])){
            rt[v] = u;
            lt[u] = v;
            return 1;
        }
    }
    return 0;
}
// n = size of left set
int bipartiteMatching(int n){
    memset(lt,-1,sizeof lt);
    memset(rt,-1,sizeof rt);
    int mxBpm=0,i;
    for(i=1;i<=n;i++){
        memset(vs,0,sizeof vs);
        if(bpm(i)) mxBpm++;
    }
    return mxBpm;
}
