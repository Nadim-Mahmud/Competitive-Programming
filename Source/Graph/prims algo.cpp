/**
    Name : Minimum Spanning Tree Prims Algorithm

*/
#include <bits/stdc++.h>
using namespace std;
#define LL long long
#define EPS 0.00000001
#define PI 2*acos(0.0)
#define MOD 1000000007
#define ck(XX) cout<<XX<<endl
#define set(XX,POS) XX|(1<<POS)
#define reset(XX,POS) XX&(~(1<<POS))
#define check(XX,POS) (bool)(XX&(1<<POS))
#define toggle(XX,POS) (XX^(1<<POS))
#define Fin freopen("input.txt","r",stdin)
#define Fout freopen("output.txt","w",stdout)
#define valid(X,Y,R,C) X>=0 && X<R && Y>=0 && Y<C
#define MS(ARRAY,VALUE) memset(ARRAY,VALUE,sizeof(ARRAY))
#define RT printf("Run Time : %0.3lf seconds\n", clock()/(CLOCKS_PER_SEC*1.0))
struct Edge{
    int from;
    int to;
    int cost;
    Edge(){};
    Edge(int from_, int to_, int cost_){from=from_; to=to_; cost=cost_;}
};

bool operator < (Edge A, Edge B) {return A.cost > B.cost;}


int node, edge; //Total number of node and edge
int sow; //Sum Of Weight of MST
bool taken[100];
vector <Edge> adj[100]; //For storing given graph data
vector <Edge> tree[100]; //For creating new minimum spanning tree


void PMST(int source)
{
    sow = 0;
    MS(taken,0);

    priority_queue <Edge> Q;
    taken[source] = 1;
    for(int i=0; i<adj[source].size();i++) //Pushing all edge of source node
    {
        Q.push(Edge(adj[source][i].from, adj[source][i].to, adj[source][i].cost));
    }

    while(!Q.empty())
    {
        Edge u = Q.top(); //taking the minimum edge from taken nodes
        Q.pop();

        if(taken[u.to]) continue;

        taken[u.to] = 1;
        sow += u.cost;

        tree[u.from].push_back(u); //creating tree
        tree[u.to].push_back(Edge(u.to,u.from,u.cost)); //creating tree


        for(int i=0; i<adj[u.to].size();i++) //Pushing all edge of taken node
        {
            if(taken[adj[u.to][i].to]) continue;
            Q.push(Edge(adj[u.to][i].from, adj[u.to][i].to, adj[u.to][i].cost));
        }
    }
    return;
}


void Graph_Input(int x, int y, int cost)
{
    adj[x].push_back(Edge(x,y,cost));
    adj[y].push_back(Edge(y,x,cost));
    return;
}


void show_tree()
{
    for(int i=0; i< node; i++)
    {
        for(int j=0; j<tree[i].size(); j++)
        {
            printf("%d %d %d\n",tree[i][j].from, tree[i][j].to, tree[i][j].cost);
        }
    }
    printf("\n");
    return;
}


int main()
{
    int tc, cn=0;
    scanf("%d",&tc);
    while(tc--){
        for(int i=0; i<100; i++) {adj[i].clear(); tree[i].clear();}
        scanf("%d%d",&node, &edge);
        for(int i=0; i<edge; i++)
        {
            int x,y,z;
            scanf("%d%d%d",&x,&y,&z);
            Graph_Input(x,y,z);
        }
        PMST(1);

        printf("Case %d: %d\n", ++cn, sow);
        show_tree();
    }
    return 0;
}
