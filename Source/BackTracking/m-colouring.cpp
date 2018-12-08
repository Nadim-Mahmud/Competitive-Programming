/**
*   this variation is to print maximum node with a spacific colour here 2
*/

///ncol is number of colours
vector<int>ed[MX];
int col[MX];
int n,m,tmp[MX],in,ncol=2;

bool safe(int pos,int cl){
    int x;
    for(int i=0;i<ed[pos].size();i++){
        /// here restricted adjecent colour is black||1
        x = ed[pos][i];
        if(col[x]==cl&&cl==1) return 0;
    }
    return 1;
}

void backtrack(int pos ,int cnt){
    if(pos>n){
        if(cnt>in){
            in = 0;
            for(int i=1;i<=n;i++){
                if(col[i]==1) tmp[++in] = i;
            }
        }
        return;
    }
    int i,cn=0;

    ///trying all colour in a spacific node

    for(i=1;i<=ncol;i++){

        if(i==1) cn = 1;
        else cn = 0;

        if(safe(pos,i)){
            col[pos] = i;
            backtrack(pos+1,cnt+cn);
            col[pos] = 0;
        }
    }
}

void free(){
    in = 0;
    for(int i=0;i<105;i++){
        col[i] = tmp[i] = 0;
        ed[i].clear();
    }
}

///calling with : backtrack(0,0);
