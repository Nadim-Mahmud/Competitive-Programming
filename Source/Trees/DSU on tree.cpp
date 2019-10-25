    ///1153 - Internet Bandwidth my code:
vector <int> ed[102];
int flw[102][102],cap[102][102],par[102],rcp[102];
char vis[101];
int s,t,d;
int dfs()
{
    queue <int> q;
    int v,f,i;
    q.push(s);
    MS(vis,0);
    vis[s]=1;
    rcp[s]=1E7;
    while(!q.empty()){
        f=q.front();
        q.pop();
        for(i=0;i<ed[f].size();i++){
            v=ed[f][i];
            if(!vis[v]&&(cap[f][v]-flw[f][v]>0)){
                vis[v]=1;
                par[v]=f;
                rcp[v]=min(cap[f][v]-flw[f][v],rcp[f]);
                if(v==d){
                    return 1;
                }
                q.push(v);
            }
        }
    }
    return 0;
}

int main()
{
    int i,j,a,b,ts,cn=0,n,m,w,cp,cfl,tf;
    //freopen("test.txt","r",stdin);
    scanf("%d",&ts);
    while(ts--){
        scanf("%d %d %d %d",&n,&s,&d,&m);
        for(i=0;i<=n;i++){
            ed[i].clear();
        }
        MS(cap,0)
        for(i=0;i<m;i++){
            scanf("%d %d %d",&a,&b,&w);
            if(w==0)continue;
            if(!cap[a][b]){
                ed[a].push_back(b);
                ed[b].push_back(a);
            }
            cap[a][b]+=w;
            cap[b][a]+=w;
        }
        tf=0;
        MS(flw,0);
        while(dfs()){
            cp=d;
            cfl=rcp[d];
            tf+=cfl;
            //P(tf)
            do{
                flw[par[cp]][cp]+=cfl;
                flw[cp][par[cp]]-=cfl;
                cp=par[cp];
            }while(cp!=s);
        }
        printf("Case %d: %d\n",++cn,tf);
    }
    return 0;
}
