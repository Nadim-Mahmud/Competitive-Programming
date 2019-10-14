#include <bits/stdc++.h>
#define ll long long
#define P(X) cout<<"db "<<X<<endl;
#define P2(X,Y) cout<<"db2 "<<X<<" "<<Y<<endl;
#define rep(i,n) for(i=1;i<=n;i++)
#define FO freopen("t.txt","w",stdout);
#define MS(XX,YY) memset(XX,YY,sizeof(XX));
#define pii pair<int,int>
#define chk(n,i) (bool)(n&(1<<i))
#define on(n,i) (n|(1<<i))
#define off(n,i) n=n&(~(1<<i))
#define eps 10e-7
#define MX 1000005
using namespace std;

//ll ara[MX];
//ll ar1[MX],ar2[MX];
//char st[MX];
//char st1[MX],st2[MX];
int main()
{
    ll i,j,a,b,ts,cn=0,cas=0,n,m,x,y,sum=0,mn=INT_MAX,mx=0;
    //freopen("test.txt","r",stdin);
    cin>>ts;
    /*while(++cas<=ts){
        scanf("%lld",&);

        printf("Case %lld: \n",cas,);
    }*/
    return 0;
}


        ///******    Number Theory    ******///


/// Number of primes --->
/**
1 	10 	4
2 	100 	25
3 	1,000 	168
4 	10,000 	1,229
5 	100,000 	9,592
6 	1,000,000 	78,498
7 	10,000,000 	664,579
8 	100,000,000 	5,761,455
9 	1,000,000,000 	50,847,534
10 	10,000,000,000 	455,052,511
11 	100,000,000,000 	4,118,054,813
12 	1,000,000,000,000 	37,607,912,018
13 	10,000,000,000,000 	346,065,536,839
14 	100,000,000,000,000 	3,204,941,750,802
15 	1,000,000,000,000,000 	29,844,570,422,669
16 	10,000,000,000,000,000 	279,238,341,033,925
17 	100,000,000,000,000,000 	2,623,557,157,654,233
18 	1,000,000,000,000,000,000 	24,739,954,287,740,860
19 	10,000,000,000,000,000,000 	234,057,667,276,344,607
20 	100,000,000,000,000,000,000 	2,220,819,602,560,918,840
21 	1,000,000,000,000,000,000,000 	21,127,269,486,018,731,928
22 	10,000,000,000,000,000,000,000 	201,467,286,689,315,906,290
23 	100,000,000,000,000,000,000,000 	1,925,320,391,606,803,968,923
24 	1,000,000,000,000,000,000,000,000 	18,435,599,767,349,200,867,866
25 	10,000,000,000,000,000,000,000,000 	176,846,309,399,143,769,411,680
*/



/// *** ncr [n*r] [but it is dp]

ll cr[902][902];
ll ncr(int n,int r){
    if(n==r)return 1;
    if(r==1)return n;
    if(cr[n][r])return cr[n][r];
    return cr[n][r] = ncr(n-1,r)+ncr(n-1,r-1);
}

/// *** Is prime [sqrt(n)]

int isprime(int n){
    if(n==2)return 1;
    if(!(n%2)||n<2)return 0;
    int i,sq=sqrt(n)+2;
    for(i=3; i<sq; i+=2)if(!(n%i))return 0;
    return 1;
}

/// *** Normal Sieve [ n*sqrt(n) ]
///   [ mark can be check for is prime ]

#define SZ 100009
int pr[78600],in=0;
bool mr[SZ+3];

void sieve(){
    int i,j,sq,p;
    sq=sqrt(SZ)+2;
    mr[1]=1;
    for(i=2; i<sq; i++){
        if(!mr[i]){
            for(j=i*i;j<=SZ;j+=i){
                mr[j]=1;
            }
        }
    }
    for(i=2; i<SZ; i++){
        if(!mr[i]){
            pr[++in] = i;
        }
    }
}

/// *** Efficient Sieve [ n*log(n)~]
///   [mr array can't be use for prime check || Odd numbers only]

#define SZ 1000009
int pr[78600],in=0;
bool mr[SZ+3];

void sieve(){
    int i,j,sq,p;
    sq=sqrt(SZ)+2;
    for(i=3; i<sq; i+=2){
        if(!mr[i]){
            for(j=i*i; j<=SZ; j+=i<<1){
                mr[j]=1;
            }
        }
    }
    pr[++in] = 2;
    for(i=3; i<SZ; i+=2){
        if(!mr[i]){
            pr[++in] = i;
        }
    }
}

/// *** BitwiseSieve  [~ n log(n) may be ]
///  [It can be use for finding huge number of primes]

#define check(X) (mkr[X>>6]&(1<<((X&63)>>1)))
#define mark(X) mkr[X>>6]|=(1<<((X&63)>>1))
int mkr[10000900/64],SZ=10000020,pr[700000],in=0;

void bitwsiv()
{
    int i,j,rt=sqrt(SZ)+1;
    for(i=3; i<=rt; i+=2){
        if(!check(i)){
            for(j=i*i; j<=SZ; j+=i+i){
                mark(j);
            }
        }
    }
    pr[++in]=2;
    for(i=3; i<=SZ; i+=2){
        if(!check(i))pr[++in]=i;
    }
}

/// *** Segmented Sieve [ sqrt(up) + prime sieve ]

/** [l = lower limit, u = upper limit]
*   [first generate all prime upto sqrt(upper limit)]
*   [Checking prime
*   n = number into that segment]
*   if(!mark[n-l]) then it is prime
*/
bool mark[u-l];
void segsiv(ll l, ll u)
{
    ll i,j,lt;
    if(l==1) mark[0] = 1;
    for(i=1; i<=in&&pr[i]*pr[i]<=u; i++){
        lt = l/pr[i];
        lt *= pr[i];
        if(lt<l) lt += pr[i];
        if(pr[i]==lt) lt += pr[i];
        for(lt; lt<=u; lt+=pr[i]){
            mark[lt-l] = 1;
        }
    }
}


/// ***Number of divisors [Quberoot(n)]

/** There can be only two prime after qube root
*   so we factorize upto quberoot by the trial divison then
*   handle ramainig <=2 prime
*/

ll NOD(ll n){
    ll i,j,r,nod=1,cn;
    for(i=1;i<=in&&(pr[i]*pr[i]*pr[i])<=n;i++){
        cn = 1;
        while(n%pr[i]==0){
            n /= pr[i];
            cn++;
        }
        nod *= cn;
    }

    r = sqrtl(n);
    while((r+1)*(r+1)<=n) r++;

    if(MillerRabin(n,8)) nod *= 2;
    else if(r*r==n&&MillerRabin(r,8)) nod *= 3;
    else if(n!=1) nod *= 4;

    return nod;
}


/// ***Eulers totients Sieve [n*n]

bool mr[MX];
int pi[MX];

void phi(){
    int i,j;
    pi[1] = 1;
    for(i=2;i<MX;i++){
        if(!mr[i]){
            for(j=i;j<MX;j+=i){
                if(pi[j]==0) pi[j] = j;
                mr[j] = 1;
                pi[j] = pi[j]/i*(i-1);
            }
        }
    }
}


/// *** Modular Inverse [ log(n) ]

/**  [ large division = (upper%M)*modi(low)) ]
*    [ after modular division value should be moded ]
*    [ mainly it returns an iverse and multiplicable value ]
*    [ for loop == first loop throw all number(multiply) with simple mod -]
*    [- then modular inverse ]
*/

#define M 1000000007
#define pii pair<ll,ll>
pii extnuc(ll a,ll b)
{
    if(b==0)return pii(1,0);
    pii d=extnuc(b,a%b);
    return pii(d.second,d.first-d.second*(a/b));
}

ll modi(ll n) /// n must be moded value
{
    pii d=extnuc(n,M);
    return ((d.first%M)+M)%M;
}

///Now factorial & inverse factorial with MOD
void fact()
{
    ll i;
    fr[0]=fi[0]=1;
    for(i=1;i<2000007;i++){
        fr[i]=(fr[i-1]*i)%M;
        fi[i]=(fi[i-1]*modi(i))%M;
        //P(fr[i])
    }
}


/// ***Miller-Rabin Primality test [120]

/** usigned long long is not working sometimes
*   left shift == mult by 2 right shift == divide bye 2
*   multiplying two numbers (a*b)%c avoiding overflow
*/
ll mulmod(ll a, ll b, ll mod){
    ll x = 0,y = a % mod;
    while (b > 0){
        if (b&1) x = (x + y) % mod;
        y = (y<<1) % mod;
        b >>= 1;
    }
    return x % mod;
}

//Bigmod
ll modulo(ll n, ll r, ll mod){
    ll x = 1;
    while (r > 0){
        if (r&1) x = mulmod(x,n,mod);
        n = mulmod(n,n,mod);
        r >>= 1;
    }
    return x % mod;
}
/// higher value of "it" ensure higher percision [recomendation 7]

bool MillerRabin(ll p,int it)
{
    if (p < 2) return 0;
    if (p != 2 && p % 2==0) return 0;

    ll i,a,tmp,mod,s=p-1;

    while(s%2==0){
        s>>=1;
    }
    for(i=0;i<it;i++){
        a = rand()%(p-1)+1;
        tmp = s;
        mod = modulo(a, tmp, p);
        if(mod==1 || mod == p-1)continue;
        while (tmp!=p-1&&mod!=1&&mod!=p-1){
            mod = mulmod(mod, mod, p);
            tmp <<= 1;
        }
        if(mod!=p-1) return 0;
    }
    return 1;
}

/** this body of miller rabin giving correct ans upto 10^9
*   with lower time complexity
*/

for(i=1;i<=it;i++){
    a = rand() % (p - 1) + 1, tmp = s;
    mod = mod_pow(a, tmp, p);
    while(tmp != p - 1 && mod != 1 && mod != p - 1) {
        mod = mod_mul(mod, mod, p);
        tmp *= 2;
    }
    if (mod != p - 1 && tmp % 2 == 0) return false;
}


///       Gaussian Elemination [n^3]

void gauss(){
    for(i=0;i<n;i++){
        l = i;
        for(j=i+1;j<n;j++){
            if(abs(ara[j][i])>abs(ara[l][i])){
                l = j;
            }
        }
        if(i!=l){
            for(j=0;j<=n;j++){
                swap(ara[i][j], ara[l][j]);
            }
        }
        for(j=0;j<n;j++){
            if(j!=i){
                tmp = ara[j][i]/ara[i][i];
                for(k=i;k<=n;k++){
                    ara[j][k] -= tmp*ara[i][k];
                }
            }
        }
    }
    for(i=0;i<n;i++){
        ara[i][n] /= ara[i][i]; ///final result
    }
}


///           Catlan Number

///Time Complexity =~ onetime n for factorial
/// and in each query [log(n)]

#define ll long long
#define MX 1000005
#define MD 100000007

ll powr(ll n, ll p){
    if(p==0) return 1;
    if(p==1) return n;
    if(p&1LL) return (powr(n,p-1)*n)%MD;
    else{
        ll x = powr(n,p/2)%MD;
        return (x*x)%MD;
    }
}

ll inverse(ll n){
    return (powr(n,MD-2))%MD;
}

ll ft[MX];

void fact(){
    ll i;
    ft[0] = 1;
    for(i=1;i<MX;i++){
        ft[i] = (ft[i-1]*i)%MD;
    }
}

ll nCr(ll n,ll r){
    ll x = ft[n];
    ll y = inverse((ft[r]*ft[n-r])%MD)%MD;
    return (x*y)%MD;
}

ll catalan(ll n){
    ll  x = nCr(2*n,n);
    return (x*inverse((n+1)))%MD;
}


///   *** Matrix Exponentiation [2^3*log(n)]
///   by calling exponent this will calculate (ara^p)%md;
///   here n is the dimention of the matrix ara

#define ll long long
#define DIM 5

ll tmp[DIM][DIM],result[DIM][DIM];

void cpy(ll ar1[][DIM],ll ar2[][DIM],ll n){
    for(int i=1;i<=n;i++)
        for(int j=1;j<=n;j++)
            ar1[i][j] = ar2[i][j];
}

void mult(ll ar1[][DIM], ll ar2[][DIM], ll n, ll mod){
    int i,j,k;
    memset(tmp,0,sizeof tmp);
    for(i=1;i<=n;i++){
        for(j=1;j<=n;j++){
            for(k=1;k<=n;k++){
                tmp[i][j] = (tmp[i][j] + ar1[i][k]*ar2[k][j])%mod;
            }
        }
    }
    cpy(ar1,tmp,n);
}

void exponent(ll ara[][DIM], ll p,ll n, ll mod){
    int i,j;
    //initializing identity
    for(i=1;i<=n;i++)
        for(j=1;j<=n;j++)
            if(i==j) result[i][j] = 1;
            else result[i][j] = 0;

    while(p){
        if(p&1){
            mult(result,ara,n,mod);
        }
        mult(ara,ara,n,mod);
        p >>= 1;
    }
    cpy(ara,result,n);
}



             ///******    Graph    ******///


///  *** BFS [ total(node + edge) ]

vector <int> ed[MX];
int lev[MX],ms,par[MX];
bool vs[MX];
int bfs(int s){
    int i,j,d,f,v;
    queue <int> q;
    q.push(s);
    lev[s]=0;
    vs[s]=1;
    while(!q.empty()){
        f=q.front();
        q.pop();
        for(i=0;i<ed[f].size();i++){
            v=ed[f][i];
            if(!vs[v]){
                vs[v]=1;
                q.push(v);
                lev[v]=lev[f]+1;
            }
        }
    }
    return lev[destinition];
}

/// *** BFS in 2D Grid [edge + node]

int lev[MX][MX],m,n;
bool vs[MX][MX];
char st[MX][MX];
int dx[]={1,-1,0, 0, 1,1,};
int dy[]={0, 0,1,-1,-1, };

int bfs(int fx,int fy) /// starting position
{
    int i,v,x,y,md=0;
    queue <pii> q;
    pii pr;
    vs[fx][fy]=1;
    lev[fx][fy]=0;
    q.push(pii(fx,fy));
    while(!q.empty()){
        pr=q.front();
        fx=pr.first;
        fy=pr.second;
        q.pop();
        for(i=0;i<4;i++){
            x=fx+dx[i];
            y=fy+dy[i];
            if(x<0||x>=n||y<0||y>=m)continue; /// for zero based index
            if(!vs[x][y]&&st[x][y]!='#'){ /// # is blocked
                q.push(pii(x,y));
                vs[x][y]=1;
                lev[x][y]=lev[fx][fy]+1;
                if(st[x][y]=='d'){ /// d is destinition
                    md=max(md,lev[x][y]);
                }

            }
        }

    }
    return md;/// max distance of d

}

/// *** DFS [E+V]

vector<int>ed[MX];
bool vs[MX];
int lev[MX];

void dfs(int n)
{
    if(vs[n]) return;
    vs[n] = 1;
    int i,v;
    for(i=0;i<ed[n].size();i++){
        v = ed[n][i];
        dfs(v);
    }
}

///  *** IDDFS [v+E]

vector<int>ed[MX];
int vs[MX],pr[MX];

bool dfs(int f){
    stack<int>st;
    st.push(f);
    while(!st.empty()){
        f = st.top();
        st.pop();
        if(vs[f]) return 1;
        vs[f] = 1;
        for(int i=0;i<ed[f].size();i++){
            int v = ed[f][i];
            if(!vs[v]){
                st.push(v);
            }
        }
    }
    return 0;
}


///  *** Dijkstra [edge * log (node)]

struct node{
    int id,cost;
    node(){}
    node(int nid,int ncost)
    {
        id=nid;
        cost=ncost;
    }
    bool operator < (const node&x)const{
        return  cost>x.cost;
    }
};

vector <int> ed[MX],ec[MX];
int ds[MX];
void dxt(int s){
    priority_queue <node> q;
    q.push(node(s,0));
    ds[s]=0;
    node fn;
    int i,u,v;
    while(!q.empty()){
        fn=q.top();
        q.pop();
        u=fn.id;
        if(fn.cost!=ds[u])continue;
        for(i=0;i<ed[u].size();i++){
            v=ed[u][i];
            if(ds[v]>ds[u]+ec[u][i]){
                ds[v]=ds[u]+ec[u][i];
                q.push(node(v,ds[v]));
            }
        }
    }
}

/// *** Minimum Spanning Tree :: Kruskal [mlogm+m]

struct edge{
    int u,v,w;
    bool operator < (const edge& x)const {
        return w>x.w;
    }
}ab;

int ara[105],n;
vector<edge>vc;
/// vc.push_back(ab) make edge and push

int find(int n){
    if(ara[n] == n) return n;
    return ara[n] = find(ara[n]);
}

int mst()
{
    int sum=0,i,u,v,cn=0,mn=INT_MAX;
    for(i=0;i<=102;i++) ara[i] = i;
    sort(vc.begin(),vc.end());
    for(i=0;i<vc.size();i++){
        u = find(vc[i].u);
        v = find(vc[i].v);
        if(u!=v){
            ara[v] = u;
            sum+=vc[i].w;
            cn++;
            //mn = min(mn,vc[i].w);
            if(cn==n) break;
        }
    }
    return sum;
}


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




            ///     *** Data Structure ***    ///




///  *** Disjoint Set Union Find [n||1]

int parent(int n)
{
    if(rp[n]==n)return n;
    return rp[n]=parent(rp[n]);
}
void setUp(int a,int b){
    rp[parent(b)]=parent(a);
}


/// *** Squre Root Decomposition [Sqrt(n)*Query]

/** 1 based indexing
*   if index are 1 2 3 4 5 6 7 and block size is 3 then
*   1 2 is first block and 3 4 5 is in second block
*   this will reduce coding complexity
*/

int ara[MX],block[sqqrt(MX)],rt,in=0;

void update(){
    /// point update is simple
    /// segment update is like query
}

void creat(int n){
    in=0;
    int mn = INT_MAX,i;
    for(i=1;i<n;i++){
        mn = min(ara[i],mn);
        if(i%rt==0){
            block[++in] = mn;
            mn = INT_MAX;
        }
    }
    block[++in] = mn;
}

int query(int l,int r){
    int mn = ara[l];

    while(l%rt!=0&&l<=r&&l!=1){
        mn = min(ara[l],mn);
        l++;
    }
    while(l+rt<=r){
        l += rt;
        mn = min(mn,block[l/rt]);
    }
    while(l<=r){
        mn = min(ara[l],mn);
        l++;
    }
    return mn;
}

///   *** Merge Sort Tree [nlog(n) + query*log(n)*log(n)]


int ara[MX];
vector<int>seg[MX*4];
/// 1 based index

int bns(int n,int val){
    /// lower bound 1 2 2 2 3 4
    /// then returns 2
    int mid,i,j,low=0,high = seg[n].size()-1;
    while((high-low)>4){
        mid = (low+high)/2;
        if(seg[n][mid]<=val) low = mid;
        else high = mid - 1;
    }
    for(low;low<=high&&low<seg[n].size();low++){
        if(seg[n][low]>val) break;
    }
    return seg[n].size()-low; /// numbers greater than value
}

void mergee(int x,int y,int z){ /// merging 2 vector x and y to Z in sorted order
    int i,j,k,md,sz;
    sz = seg[x].size() + seg[y].size();
    for(i=0,j=0,k=0;k<sz;k++){
        if(i>=seg[x].size()) seg[z].push_back(seg[y][j++]);
        else if(j>=seg[y].size()) seg[z].push_back(seg[x][i++]);
        else if(seg[x][i]<seg[y][j]) seg[z].push_back(seg[x][i++]);
        else seg[z].push_back(seg[y][j++]);
    }
}

/** [low,high]  total range :: variable range
*   [qlow,qhigh] query range
*   pos = current position
*/
void creat(int low,int high,int pos){ /// creating merge sort tree
    if(low==high){
        seg[pos].push_back(ara[low]);
        return ;
    }
    int mid = (low+high)/2;
    creat(low,mid,pos*2);
    creat(mid+1,high,pos*2+1);
    mergee(pos*2,pos*2+1,pos);
    /// merge with stl
    /// merge(seg[pos*2].begin() , seg[pos*2].end(), seg[pos*2].begin(), seg[pos*2].end(),back_inserter(seg[pos]));
}

int query(int low,int high,int qlow,int qhigh,int pos,int val){
    if(qlow>qhigh) return 0;
    if(qlow>high||qhigh<low) return 0;
    if(qlow<=low&&qhigh>=high){
        return bns(pos,val);
    }
    int mid = (low + high)/2;
    return query(low,mid,qlow,qhigh,pos*2,val) + query(mid+1,high,qlow,qhigh,pos*2+1,val);
}

/// *** For Rnage orders statistics (find k'th number in sorted segment)

vector<pii>input;
vector<int>seg[MX*4];

void creat(int low,int high,int pos){ /// creating merge sort tree
    if(low==high){
        seg[pos].push_back(input[low-1].second); /// in is 0 based
        return ;
    }
    int mid = (low+high)/2;
    creat(low,mid,pos*2);
    creat(mid+1,high,pos*2+1);
    mergee(pos*2,pos*2+1,pos);
}

/** calculating total number in left range lower than the given index
*  if numbers are greater than equals to the searging value than look into left
*  searhing on right sub array and substracting left sub arrys given up values
*/

int query(int low,int high,int qlow,int qhigh,int pos,int val)
{
    if(low==high) return seg[pos][0];
    int mid = (low+high)>>1,left=pos<<1;

    int total = upper_bound(seg[left].begin(),seg[left].end(),qhigh) -
    lower_bound(seg[left].begin(),seg[left].end(),qlow);

    if(total>=val){
        return query(low,mid,qlow,qhigh,pos*2,val);
    }
    else{
        return query(mid+1,high,qlow,qhigh,pos*2+1,val-total);
    }
}

sort(input.begin(),input.end());




///   *** BIT [Log(n)] space [n]

/**   1 based index
*     which functions has inverse function that can be solve bye BIT
*     it works like consucative sums but in log(n)
*/

int n=SIZE of space;
void update(int idx,int val)//adding value val to idx index
{
    while(idx<=n){
        bitree[idx]+=val;
        idx+=idx&(-idx); // Last set of digit
    }
}
int query(int idx){// returns sum of 1 to idx index
    int sum=0;
    while(idx>0){
        sum+=bitree[idx];
        idx-=idx&(-idx);
    }
    return sum;
}

/// *** Segment Tree [log(total array size)*Query]

/** [ulow,uhigh] Query Range
*   [low,high] total range of root
*   [qlow,qhigh] Query Range
*   Currrent position = pos
*   0 based Index And Root is also 0
*/

int ara[MX],seg[4*MX],lazy[4*MX];

void creat(int low,int high,int pos)
{
    if(low==high){
        seg[pos] = ara[low]; // reached leaf and update
        return ;
    }
    int mid = (high+low)/2;
    creat(low,mid,pos*2+1);
    creat(mid+1,high,pos*2+2);
    seg[pos] += seg[pos*2+1] + seg[pos*2+2];
}

void update(int low,int high,int ulow,int uhigh,int val,int pos)
{
    if(low>high) return ;
    if(lazy[pos]!=0){ /// is not propagated yet
        seg[pos] += lazy[pos];
        if(low!=high){  ///if not leaf node
            lazy[pos*2+1] += lazy[pos];
            lazy[pos*2+2] += lazy[pos];
        }
        lazy[pos] = 0;
    }

    if(ulow>high||uhigh<low) return; ///No overlap
    if(ulow<=low&&uhigh>=high){ /// Total Overlap
        seg[pos] += val;
        if(low!=high){
            lazy[pos*2+1] += val;
            lazy[pos*2+2] += val;
        }
        return;
    }
    /// Partial overlap
    int mid = (high+low)/2;

    update(low,mid,ulow,uhigh,val,pos*2+1);
    update(mid+1,high,ulow,uhigh,val,pos*2+2);
    seg[pos] = seg[pos*2+1] + seg[pos*2+2]; /// Updating the intermediate node
}

int query(int low,int high,int qlow,int qhigh,int pos)
{
    if(low>high) return 0;
    if(lazy[pos]!=0){
        seg[pos] += lazy[pos];
        if(low!=high){
            lazy[pos*2+1] += lazy[pos];
            lazy[pos*2+2] += lazy[pos];
        }
        lazy[pos] = 0;
    }

    if(qlow>high||qhigh<low) return 0;

    if(qlow<=low&&qhigh>=high)
        return seg[pos];

    int mid = (high+low)/2;

    return query(low,mid,qlow,qhigh,pos*2+1) + query(mid+1,high,qlow,qhigh,pos*2+2);
}



///      *** Trie[(number of words)*(maximum lenght)]


/**
    This trie is for All Uppercase & Lowercase letters
*/

int mp(char ch){
    if(ch<95) return (int) (ch - 'A');
    return (int) (ch - 'a' + 26);
}

struct node{
    bool end;
    int next[55];
    int cnt; /// here cnt is counting how many element ends at this point
    void set(){
        cnt = 0;
        end = false;
        for(int i=0;i<=53;i++){
            next[i] = 0;
        }
    }
}ara[MX];

int ptr;

void insert(char *st){
    int i,in,x=0,y;
    for(i=0;st[i];i++){
        if(st[i]==' ') continue;
        in = mp(st[i]);
        /// If the chain is not exists
        /// allocating memory by ptr++ that means new array element
        if(!ara[x].next[in]){
            ///putting next arrays position at recent array's linked part
            ara[x].next[in] = ptr++;
            ara[ptr-1].set();
        }
        x = ara[x].next[in];
    }
    ///marking last element as last
    ara[x].end=1;
    ara[x].cnt++;
}


int search(char *st){
    int i,in,x=0,cn=0;
    for(i=0;st[i];i++){
        if(st[i]==' ')continue;
        cn++;
        in = mp(st[i]);
        ///if query element not exists then returning false
        if(!ara[x].next[in]) return 0;
        x = ara[x].next[in];
    }
    if(!cn) return 1;
    /// returning how many elements like "st"
    return ara[x].cnt;
}




///              String Mathcing


///    **** Knuth–Morris–Pratt(KMP) [n+k]

#define MX 1000005
///Largest possible suffix-prefix upto this index
int lps[MX];

/// prefix size is (j+1) as key i 0 based index

void lps_calc(char key[]){
    int i = 1,j = 0;
    while(key[i]){
        //cout<<"ff\n";
        if(key[i]==key[j]){
            lps[i] = ++j;
            i++;
        }
        else{
            if(j) j = lps[j-1];
            else lps[i++] = 0;
        }
    }
}
/**
*   @return -1 if no match found else return staring position of match
*/
int match(char txt[],char key[],int key_ln){
    int i=0,j=0,cn=0;
    while(txt[i]){
        if(txt[i]==key[j]){
            i++;
            j++;
            if(j==key_ln){
                return i - key_ln; //satring pos
                ///How many match
                //cn++;
                //j=lps[j-1];
            }
        }
        else{
            if(j) j = lps[j-1];
            else i++;
            /// j = 0 and pos are not equal then skipping this position
        }
    }
    return cn;
}

char txt[MX],key[MX];



///               COMPUTATIONAL GEOMETRY



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






///    ***Bigint Factorial ****

//package main;

import java.math.BigInteger;
import java.util.Scanner;

public class Main {

	public static Scanner sc;

	public static void main(String [] arg) {

		BigInteger [] fact = new BigInteger[200];

    	fact[0] = BigInteger.ONE;

		for(int i=1;i<=150;i++) {
			fact[i] = fact[i-1].multiply(new BigInteger(i + ""));
		}

        sc = new Scanner(System.in);
        int ts = sc.nextInt();
        int n,cas=0;

        while(++cas<=ts) {
        	n = sc.nextInt();
        	System.out.println(fact[n]);
        }
	}
}

