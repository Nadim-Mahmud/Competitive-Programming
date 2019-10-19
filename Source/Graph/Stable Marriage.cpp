//// O(n*n)

int womens_priority[MX][MX],married[MX],couple[MX];
queue<int>mens_priority[MX];

void free(){
    for(int i=0;i<MX;i++){
        while(!mens_priority[i].empty()) mens_priority[i].pop();
        married[i] = couple[i] = 0;
        memset(womens_priority[i],0,sizeof womens_priority[i]);
    }
}

int main()
{
    int i,j,a,b,ts,cn=0,cas=0,n,m,x,y,sum=0,mn=INT_MAX,mx=0;
    //freopen("test.txt","r",stdin);
    cin>>ts;
    while(++cas<=ts){
        scanf("%d",&n);
        free();
        //womens
        for(i=1;i<=n;i++){
            for(j=1;j<=n;j++){
                scanf("%d",&x);
                womens_priority[i][x] = n-j;
            }
        }
        //mens
        for(i=1;i<=n;i++){
            for(j=1;j<=n;j++){
                scanf("%d",&x);
                mens_priority[n+i].push(x);
            }
        }
        cn=0;
        while(cn!=n){
            for(i=n+1;i<=2*n;i++){
                if(!married[i]){
                    x = mens_priority[i].front();
                    mens_priority[i].pop();
                    married[i] = 1;
                    if(couple[x]==0) {
                        couple[x] = i,cn++;
                    }
                    else{
                        if(womens_priority[x][couple[x]]<womens_priority[x][i]){
                            married[couple[x]] = 0;
                            couple[x] = i;
                        }
                        else married[i] = 0;
                    }
                }
            }
        }
        printf("Case %d:",cas);
        for(i=1;i<=n;i++){
            printf(" (%d %d)",i,couple[i]);
        }
        puts("");
    }
    return 0;
}
