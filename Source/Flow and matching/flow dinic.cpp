/**
    Name : Flow Dinic
    Time Complexity : O(V^2*E)
    Use : if u to v have mutiple path dont add path multiple time.
    just add the path cost to capacity matrix
    * reverse edge must be push. if unidirectional graph then push reverse edge
    seting the reverse edge capacity zero
*/


typedef pair<int,int> pii;
vector<int>adj[MAX];
int cap[MAX][MAX];

int q[100000];
bool fl[MAX][MAX]; // u theke v te path ase kina

/**
    n = number of nodes
    s = starting node
    t = destination node
    returns maximum flow
*/

int dinic( int n,int s,int t ) {
    int prev[MAX], u, v, i, z, flow = 0, qh, qt, inc;
    while(1) {
        memset( prev, -1, sizeof( prev ) );
        qh = qt = 0;
        prev[s] = -2;
        q[qt++] = s;
        while( qt != qh && prev[t] == -1 ) {
            u = q[qh++];
            for(i = 0; i < adj[u].size(); i++) {
                v = adj[u][i];
                if( prev[v] == -1 && cap[u][v] ) {
                    prev[v] = u;
                    q[qt++] = v;
                }
            }
        }
        if(prev[t] == -1) break;
        for(z = 1; z <= n; z++) if( prev[z] !=- 1 && cap[z][t] ) {
            inc = cap[z][t];
            for( v = z, u = prev[v]; u >= 0; v = u, u=prev[v]) inc = min( inc, cap[u][v] );
            if( !inc ) continue;
            cap[z][t] -= inc;
            cap[t][z] += inc;
            for(v=z, u = prev[v]; u >= 0; v = u, u = prev[v]) {
                cap[u][v] -= inc;
                cap[v][u] += inc;
            }
            flow += inc;
        }
    }
    return flow;
}


        scanf("%d %d %d %d",&n,&s,&t,&c);
        CLR(cap);
        CLR(adj);
        CLR(fl);
        for(i=1;i<=c;i++){
            scanf("%d %d %d",&u,&v,&x);
            if(!fl[u][v]){
                adj[u].push_back(v);
                adj[v].push_back(u);
                fl[u][v] = fl[v][u];
            }
            cap[u][v] += x;
            cap[v][u] += x;
        }
        printf("Case %d: %d\n",cas,dinic(n,s,t));
