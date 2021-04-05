#include <bits/stdc++.h>

#define P(X) cout<<"db "<<X<<endl;
#define P2(X,Y) cout<<"d2 "<<X<<" "<<Y<<endl;
#define P3(X,Y,Z) cout<<"d3 "<<X<<" "<<Y<<" "<<Z<<endl;
#define SQ(x) ((x) * (x))

#define ll long long
#define pii pair<int,int>

#define bchk(n,i) (bool)(n&(1<<i))
#define bon(n,i) (n|(1<<i))
#define boff(n,i) n=n&(~(1<<i))

#define distance(a,b) (sq(a.x-b.x) + sq(a.y-b.y))
#define MAX3(a,b,c) max(a,max(b,c))
#define MS(XX,YY) memset(XX,YY,sizeof(XX));
#define FastIO ios_base::sync_with_stdio(0);cin.tie(nullptr);
#define eps 10e-9
#define MX 1000005

using namespace std;
int m,n;
int main()
{
    int i,j,test,cas=0;
    int a,b;
    freopen("test.txt","r",stdin);
    scanf("%d",&test);
    while(test--){
        scanf("%",&);

        printf("Case %d:\n",++cas,);
    }
    return 0;
}



/**

    DATA STRUCTURE! DATA STRUCTURE! DATA STRUCTURE!
    DATA STRUCTURE! DATA STRUCTURE! DATA STRUCTURE!
    DATA STRUCTURE! DATA STRUCTURE! DATA STRUCTURE!
    DATA STRUCTURE! DATA STRUCTURE! DATA STRUCTURE!

*/


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

/**
    Name : Sparse table(RMQ)
    Description : Find min/max
    Time Complexity : Build O(nlogn) Query O(1)
*/
#include <bits/stdc++.h>
using namespace std;
//0 Indexed
#define MX 10000
int spt[MX][22];
int n,ar[MX]={ 7, 2, 3, 0, 5, 10, 3, 12, 18 };
void buildST()
{
	for (int i = 0; i < n; i++) spt[i][0] = ar[i];

	for (int j = 1; (1 << j) <= n; j++) {
		for (int i = 0; (i + (1 << j) - 1) < n; i++) {
            spt[i][j] = min(spt[i + (1 << (j - 1))][j - 1] , spt[i][j - 1]);
		}
	}
}

int query(int l, int r)
{
    if(l>r) return INT_MAX;
	int j = (int)log2(r - l + 1);
	///j = 31 - __builtin_clz(r - l+1);
	return min (spt[l][j], spt[r - (1 << j) + 1][j]);
}

// Driver program
int main()
{

    n = 9;
	buildST();
	cout << query(4, 7) << endl;
	cout << query(7, 8) << endl;
	return 0;
}


///   *** BIT O(Log(n)) space O(n)

/**   1 based index
*     which functions has inverse function that can be solve bye BIT
*     it works like consucative sums but in log(n)
*/

int n=SIZE; //of space;
void update(int idx,int val)//adding value val to idx index
{
    while(idx<=n){
        bitree[idx]+=val;
        idx+=idx&(-idx); // Last set of digit
    }
}
int query(int idx){// returns sum of [1,idx] index
    int sum=0;
    while(idx>0){
        sum+=bitree[idx];
        idx-=idx&(-idx);
    }
    return sum;
}

/**
    Description : BIT range update range query
    Time Complexity : all log(n)
*/

/// Remember to use 1 based indexing
//const int MX = 100005;

ll query(ll *bit, int indx)
{
    ll sum = 0;
    while (indx) {
        sum += bit[indx];
        indx -= (indx & -indx);
    }
    return sum;
}

void update(ll *bit, int indx, ll x)
{
    while (indx < MX) {
        bit[indx] += x;
        indx += (indx & -indx);
    }
}
ll B1[MX],B2[MX];//set 0

void Rupdate(int l, int r, ll v){
    update(B1, l, v);
    update(B1, r+1, -v);
    update(B2, l, -((l-1)*v));
    update(B2, r+1, r*v);
}
ll Rquery1(int p){
    ll b1,b2;
    b1 = query(B1, p);
    b2 = query(B2, p);
    return b1 * p + b2;
}
ll Requery(int l,int r){
    return Rquery1(r)-Rquery1(l-1);
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
        seg[pos] += (high-low+1)*lazy[pos];
        if(low!=high){  ///if not leaf node
            lazy[pos*2+1] += lazy[pos];
            lazy[pos*2+2] += lazy[pos];
        }
        lazy[pos] = 0;
    }

    if(ulow>high||uhigh<low) return; ///No overlap
    if(ulow<=low&&uhigh>=high){ /// Total Overlap
        seg[pos] += (high-low+1)*val;
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
        seg[pos] += (high-low+1)*lazy[pos];
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


/**
	Description : Efficient and easy segment trees , Range [l, r)
	Time Complexity : O(logn)
	from: https://codeforces.com/blog/entry/18051
*/

const int N = 1e5;  // limit for array size
int n;  // array size
int tr[2 * N];
void build() {  // build the tree
  for (int i = n - 1; i > 0; --i) tr[i] = tr[i<<1] + tr[i<<1|1];
}

void modify(int p, int value) {  // set value at position p
  for (tr[p += n] = value; p > 1; p >>= 1) tr[p>>1] = tr[p] + tr[p^1];
}

int query(int l, int r) {  // sum on interval [l, r)
  int res = 0;
  for (l += n, r += n; l < r; l >>= 1, r >>= 1) {
    if (l&1) res += tr[l++];
    if (r&1) res += tr[--r];
  }
  return res;
}

int main() {
  //scanf("%d", &n);
  int ar[]={1,3,4,5,4,5,6,4,6,4,5,6,6,5};
  n=5;
  for (int i = 0; i < n; ++i) {
    tr[n+i]=ar[i];
  }
  build();
  printf("%d\n", query(0, 4));//that means [0,3]
  printf("%d\n", query(2, 3));
  modify(0, 5);
  printf("%d\n", query(0, 4));
  return 0;
}


/**
 * Name : Maximum subsegment sum
 * Description: Segment Tree with custom merge function. .
 * Usage: construct O(N), query O(lg(N)), update O(lg(N))
 * Source: https://github.com/dragonslayerx
 */

#include <iostream>
#include <cstdio>
using namespace std;

#define MAX 50100
#define INF -1000000000

struct node {
    int sum;
    int maxs, prefix, suffix;
    node(){
        sum = prefix = suffix = 0;
        maxs = INF;
    }

    node(int sum, int maxs, int prefix, int suffix) {
        setNode(sum, maxs, prefix, suffix);
    }

    void setNode(int sum, int maxs, int prefix, int suffix){
        this->sum =sum;
        this->maxs=maxs;
        this->prefix=prefix;
        this->suffix=suffix;
    }
};

int a[MAX];
node st[4*MAX];

node merge(node left, node right){
    node t;
    t.prefix = max(left.prefix, left.sum+right.prefix);
    t.suffix = max(right.suffix, right.sum+left.suffix);
    t.sum = left.sum+right.sum;
    t.maxs = left.maxs;
    t.maxs = max(t.maxs, right.maxs);
    t.maxs = max(t.maxs, left.suffix+right.prefix);
    return t;
}

node construct(int n, int ll, int rl){
    if (ll == rl) {
        st[n].setNode(a[ll], a[ll], a[ll], a[ll]);
    } else {
        node left = construct(2*n+1, ll, (ll+rl)/2);
        node right = construct(2*n+2, (ll+rl)/2+1, rl);
        st[n] = merge(left, right);
    }
    return st[n];
}

node query(int n, int ll, int rl, int x, int y){
    int mid = (ll+rl)/2;
    if (x==ll &&  y==rl) return st[n];
    else if (y <= mid) return query(2*n+1, ll, mid, x, y);
    else if (x > mid) return query(2*n+2, mid+1, rl, x, y);
    else {
        node left = query(2*n+1, ll, (ll+rl)/2, x, mid);
        node right = query(2*n+2, (ll+rl)/2+1, rl, mid+1, y);
        return merge(left, right);
    }
}

node update(int n, int ll, int rl, int p, int val){
    if (p < ll || p > rl) return st[n];
    if (p == ll &&  p == rl) {
        st[n].setNode(val, val, val, val);
        return st[n];
    } else {
        int mid = (ll+rl)/2;
        node left = update(2*n+1, ll, (ll+rl)/2, p, val);
        node right = update(2*n+2, (ll+rl)/2+1, rl, p, val);
        st[n] = merge(left, right);
    }
    return st[n];
}

int main()
{
    int n;
    scanf("%d", &n);
    for (int i = 0; i < n; i++) scanf("%d", a+i);
    construct(0, 0, n-1);
    int q;
    scanf("%d", &q);
    while (q--) {
        int x, y;
        scanf("%d%d", &x, &y);
        x--, y--;
        printf("%d\n", query(0, 0, n-1, x, y).maxs);
    }
}
© 2019 GitHub, Inc.

/**'
    Name : Disjoint set union find
    Description : Always check parent for size or anything
    Complexity : O(n) ~ O(log n) ~ O(1)
*/

#define MX 10000
int rp[MX],sz[MX];

int parent(int n){
    if(rp[n]==n)return n;
    return rp[n]=parent(rp[n]);
}

void setUp(int a,int b){
    a = parent(a);
    b = parent(b);
    if(a==b) return;
    if(sz[a]<sz[b]){
        rp[a] = rp[b];
        sz[b] += sz[a];
    }
    else{
        rp[b] = rp[a];
        sz[a] += sz[b];
    }
}

void init(){
    for(int i=0;i<=MX;i++)
        rp[i]=i,sz[i]=1;
}


/**
 * Name : Dijoint Set with undo
 * Description : DisjointSet (Makes a set of sets, merge sets, set membership, no. of sets, undo last operation,size of each component)
 * Time Complexity : parent O(lg(N)), setUp O(lg(N)), undo O(1),
 */

#define MX 10000
int rp[MX],sz[MX];
int compo;
int pts[MX*2],in=0;

int parent(int n){
    if(rp[n]==n)return n;
    return rp[n]=parent(rp[n]);
}

// additionally storing parent which is connected to another parents
void setUp(int a,int b){
    a = parent(a);
    b = parent(b);
    if(a==b){
        pts[++in]=-1;
        return;
    }
    if(sz[a]<sz[b]){
        rp[a] = rp[b];
        sz[b] += sz[a];
        pts[++in]=a;
    }
    else{
        rp[b] = rp[a];
        sz[a] += sz[b];
        pts[++in] = b;;
    }
    compo--;
}

void undo(){
    if(!in) return;
    int n = pts[in--];
    if(n!=-1) {
        sz[parent(rp[n])] -= sz[n];
        rp[n]=n;
        compo++;
    }
}

void init(int n){
    in=0;
    for(int i=0;i<=MX;i++){
        rp[i]=i;
        sz[i]=1;
    }
    compo=n;
}

/**
    Name : Trie with Dynamic memory
    Time complexity : o((number of words)*(maximum lenght))
*/


/// Trie form shafaetsplanet
struct node {
    bool endmark;
    node* next[26 + 1];
    node()
    {
        endmark = false;
        for (int i = 0; i < 26; i++)
            next[i] = NULL;
    }
} * root;
void insert(char* str, int len)
{
    node* curr = root;
    for (int i = 0; i < len; i++) {
        int id = str[i] - 'a';
        if (curr->next[id] == NULL)
            curr->next[id] = new node();
        curr = curr->next[id];
    }
    curr->endmark = true;
}
bool search(char* str, int len)
{
    node* curr = root;
    for (int i = 0; i < len; i++) {
        int id = str[i] - 'a';
        if (curr->next[id] == NULL)
            return false;
        curr = curr->next[id];
    }
    return curr->endmark;
}
void del(node* cur)
{
    for (int i = 0; i < 26; i++)
        if (cur->next[i])
            del(cur->next[i]);

    delete (cur);
}
int main()
{
    puts("ENTER NUMBER OF WORDS");
    root = new node();
    int num_word;
    cin >> num_word;
    for (int i = 1; i <= num_word; i++) {
        char str[50];
        scanf("%s", str);
        insert(str, strlen(str));
    }
    puts("ENTER NUMBER OF QUERY";);
    int query;
    cin >> query;
    for (int i = 1; i <= query; i++) {
        char str[50];
        scanf("%s", str);
        if (search(str, strlen(str)))
            puts("FOUND");
        else
            puts("NOT FOUND");
    }
    del(root);
    return 0;
}



///      *** Trie[(number of words)*(maximum lenght)]


/**
    Name : Trie with array
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



/**
    GRAPH! GRAPH! GRAPH! GRAPH! GRAPH! GRAPH!
    GRAPH! GRAPH! GRAPH! GRAPH! GRAPH! GRAPH!
    GRAPH! GRAPH! GRAPH! GRAPH! GRAPH! GRAPH!
    GRAPH! GRAPH! GRAPH! GRAPH! GRAPH! GRAPH!
*/

///  *** BFS [ total(node + edge) ]
/// level by level travarse

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

/**
    Description : dfs with coloring and visiting time
    Complexity : O(V+E)
*/

vector<vector<int>> adj; // graph represented as an adjacency list
int n; // number of vertices

vector<int> color;

vector<int> time_in, time_out;
int dfs_timer = 0;

void dfs(int v) {
    time_in[v] = dfs_timer++;
    color[v] = 1;
    for (int u : adj[v])
        if (color[u] == 0)
            dfs(u);
    color[v] = 2;
    time_out[v] = dfs_timer++;
}


///  *** Dijkstra O(edge * log (node))

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

/**
    Name : Warshal Algorithm

    he key idea to notice here is that, we go through all
    the nodes and for each node we try to make every path
    better by going through that. Hence, Floyd Warshall can add vertex online.
    n^2 loop can add this.

    Time Complexity : O(n^3)

*/


int mtx[102][102],n;//intialize with inf;
int next[102][102];//for finding path only
void wrsl()
{
    int i,j,k;
    for(i=1;i<=n;i++){//for finding path only
        for(j=1;j<=n;j++){
            next[i][j]=j;
        }
    }

    for(k=1;k<=n;k++){
        for(i=1;i<=n;i++){
            for(j=1;j<=n;j++){
                if(mtx[i][j]>mtx[i][k]+mtx[k][j]){
                    mtx[i][j]>mtx[i][k]+mtx[k][j];
                    next[i][j]=next[i][k];//for finding path only
                }
            }
        }
    }
}
//finding path using warshal, i to j
vector <int> path;
void findpath(int i,int j)
{
    path.clear();
    path.push_back(i);
    while(i!=j){
        i=next[i][j];
        path.push_back(i);
    }
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


/**
    Descripton :
    Source :https://blog.anudeep2011.com/heavy-light-decomposition/
*/


#include <cstdio>
#include <vector>
using namespace std;

#define root 0
#define N 10100
#define LN 14

vector <int> adj[N], costs[N], indexx[N];
int baseArray[N], ptr;
int chainNo, chainInd[N], chainHead[N], posInBase[N];
int depth[N], pa[LN][N], otherEnd[N], subsize[N];
int st[N*6], qt[N*6];

/*
 * make_tree:
 * Used to construct the segment tree. It uses the baseArray for construction
 */
void make_tree(int cur, int s, int e) {
	if(s == e-1) {
		st[cur] = baseArray[s];
		return;
	}
	int c1 = (cur<<1), c2 = c1 | 1, m = (s+e)>>1;
	make_tree(c1, s, m);
	make_tree(c2, m, e);
	st[cur] = st[c1] > st[c2] ? st[c1] : st[c2];
}

/*
 * update_tree:
 * Point update. Update a single element of the segment tree.
 */
void update_tree(int cur, int s, int e, int x, int val) {
	if(s > x || e <= x) return;
	if(s == x && s == e-1) {
		st[cur] = val;
		return;
	}
	int c1 = (cur<<1), c2 = c1 | 1, m = (s+e)>>1;
	update_tree(c1, s, m, x, val);
	update_tree(c2, m, e, x, val);
	st[cur] = st[c1] > st[c2] ? st[c1] : st[c2];
}

/*
 * query_tree:
 * Given S and E, it will return the maximum value in the range [S,E)
 */
void query_tree(int cur, int s, int e, int S, int E) {
	if(s >= E || e <= S) {
		qt[cur] = -1;
		return;
	}
	if(s >= S && e <= E) {
		qt[cur] = st[cur];
		return;
	}
	int c1 = (cur<<1), c2 = c1 | 1, m = (s+e)>>1;
	query_tree(c1, s, m, S, E);
	query_tree(c2, m, e, S, E);
	qt[cur] = qt[c1] > qt[c2] ? qt[c1] : qt[c2];
}

/*
 * query_up:
 * It takes two nodes u and v, condition is that v is an ancestor of u
 * We query the chain in which u is present till chain head, then move to next chain up
 * We do that way till u and v are in the same chain, we query for that part of chain and break
 */

int query_up(int u, int v) {
	if(u == v) return 0; // Trivial
	int uchain, vchain = chainInd[v], ans = -1;
	// uchain and vchain are chain numbers of u and v
	while(1) {
		uchain = chainInd[u];
		if(uchain == vchain) {
			// Both u and v are in the same chain, so we need to query from u to v, update answer and break.
			// We break because we came from u up till v, we are done
			if(u==v) break;
			query_tree(1, 0, ptr, posInBase[v]+1, posInBase[u]+1);
			// Above is call to segment tree query function
			if(qt[1] > ans) ans = qt[1]; // Update answer
			break;
		}
		query_tree(1, 0, ptr, posInBase[chainHead[uchain]], posInBase[u]+1);
		// Above is call to segment tree query function. We do from chainHead of u till u. That is the whole chain from
		// start till head. We then update the answer
		if(qt[1] > ans) ans = qt[1];
		u = chainHead[uchain]; // move u to u's chainHead
		u = pa[0][u]; //Then move to its parent, that means we changed chains
	}
	return ans;
}

/*
 * LCA:
 * Takes two nodes u, v and returns Lowest Common Ancestor of u, v
 */
int LCA(int u, int v) {
	if(depth[u] < depth[v]) swap(u,v);
	int diff = depth[u] - depth[v];
	for(int i=0; i<LN; i++) if( (diff>>i)&1 ) u = pa[i][u];
	if(u == v) return u;
	for(int i=LN-1; i>=0; i--) if(pa[i][u] != pa[i][v]) {
		u = pa[i][u];
		v = pa[i][v];
	}
	return pa[0][u];
}

void query(int u, int v) {
	/*
	 * We have a query from u to v, we break it into two queries, u to LCA(u,v) and LCA(u,v) to v
	 */
	int lca = LCA(u, v);
	int ans = query_up(u, lca); // One part of path
	int temp = query_up(v, lca); // another part of path
	if(temp > ans) ans = temp; // take the maximum of both paths
	printf("%d\n", ans);
}

/*
 * change:
 * We just need to find its position in segment tree and update it
 */
void change(int i, int val) {
	int u = otherEnd[i];
	update_tree(1, 0, ptr, posInBase[u], val);
}

/*
 * Actual HL-Decomposition part
 * Initially all entries of chainHead[] are set to -1.
 * So when ever a new chain is started, chain head is correctly assigned.
 * As we add a new node to chain, we will note its position in the baseArray.
 * In the first for loop we find the child node which has maximum sub-tree size.
 * The following if condition is failed for leaf nodes.
 * When the if condition passes, we expand the chain to special child.
 * In the second for loop we recursively call the function on all normal nodes.
 * chainNo++ ensures that we are creating a new chain for each normal child.
 */
void HLD(int curNode, int cost, int prev) {
	if(chainHead[chainNo] == -1) {
		chainHead[chainNo] = curNode; // Assign chain head
	}
	chainInd[curNode] = chainNo;
	posInBase[curNode] = ptr; // Position of this node in baseArray which we will use in Segtree
	baseArray[ptr++] = cost;

	int sc = -1, ncost;
	// Loop to find special child
	for(int i=0; i<adj[curNode].size(); i++) if(adj[curNode][i] != prev) {
		if(sc == -1 || subsize[sc] < subsize[adj[curNode][i]]) {
			sc = adj[curNode][i];
			ncost = costs[curNode][i];
		}
	}

	if(sc != -1) {
		// Expand the chain
		HLD(sc, ncost, curNode);
	}

	for(int i=0; i<adj[curNode].size(); i++) if(adj[curNode][i] != prev) {
		if(sc != adj[curNode][i]) {
			// New chains at each normal node
			chainNo++;
			HLD(adj[curNode][i], costs[curNode][i], curNode);
		}
	}
}

/*
 * dfs used to set parent of a node, depth of a node, subtree size of a node
 */
void dfs(int cur, int prev, int _depth=0) {
	pa[0][cur] = prev;
	depth[cur] = _depth;
	subsize[cur] = 1;
	for(int i=0; i<adj[cur].size(); i++)
		if(adj[cur][i] != prev) {
			otherEnd[indexx[cur][i]] = adj[cur][i];
			dfs(adj[cur][i], cur, _depth+1);
			subsize[cur] += subsize[adj[cur][i]];
		}
}

int main() {
	int t;
	scanf("%d ", &t);
	while(t--) {
		ptr = 0;
		int n;
		scanf("%d", &n);
		// Cleaning step, new test case
		for(int i=0; i<n; i++) {
			adj[i].clear();
			costs[i].clear();
			indexx[i].clear();
			chainHead[i] = -1;
			for(int j=0; j<LN; j++) pa[j][i] = -1;
		}
		for(int i=1; i<n; i++) {
			int u, v, c;
			scanf("%d %d %d", &u, &v, &c);
			u--; v--;
			adj[u].push_back(v);
			costs[u].push_back(c);
			indexx[u].push_back(i-1);
			adj[v].push_back(u);
			costs[v].push_back(c);
			indexx[v].push_back(i-1);
		}

		chainNo = 0;
		dfs(root, -1); // We set up subsize, depth and parent for each node
		HLD(root, -1, -1); // We decomposed the tree and created baseArray
		make_tree(1, 0, ptr); // We use baseArray and construct the needed segment tree

		// Below Dynamic programming code is for LCA.
		for(int i=1; i<LN; i++)
			for(int j=0; j<n; j++)
				if(pa[i-1][j] != -1)
					pa[i][j] = pa[i-1][pa[i-1][j]];

		while(1) {
			char s[100];
			scanf("%s", s);
			if(s[0]=='D') {
				break;
			}
			int a, b;
			scanf("%d %d", &a, &b);
			if(s[0]=='Q') {
				query(a-1, b-1);
			} else {
				change(a-1, b);
			}
		}
	}
}

/**problem:

You are given a tree (an acyclic undirected connected graph) with N nodes, and edges numbered 1, 2, 3...N-1.

We will ask you to perfrom some instructions of the following form:

    CHANGE i ti : change the cost of the i-th edge to ti
    or
    QUERY a b : ask for the maximum edge cost on the path from node a to node b

Input

The first line of input contains an integer t, the number of test cases (t <= 20). t test cases follow.

For each test case:

    In the first line there is an integer N (N <= 10000),
    In the next N-1 lines, the i-th line describes the i-th edge: a line with three integers a b c denotes an edge between a, b of cost c (c <= 1000000),
    The next lines contain instructions "CHANGE i ti" or "QUERY a b",
    The end of each test case is signified by the string "DONE".

There is one blank line between successive tests.

*/


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

/**
*  Description : Useful for all same types of operations on each node on a tree
*               cnt[] storing all number of color of a type in sub tree of v
*               res[] storing sum of color which apears maximum numbers of time
*
*  Time complexity : O(nlogn)
*  Source :https://codeforces.com/blog/entry/44351
*/


vector<int>g[MX];

int sz[MX];
void getsz(int v, int p){
    sz[v] = 1;  // every vertex has itself in its subtree
    for(auto u : g[v])
        if(u != p){
            getsz(u, v);
            sz[v] += sz[u]; // add size of child u to its parent(v)
        }
}

ll cnt[MX],col[MX],ans[MX],res[MX],mxn;
bool big[MX];

void add(int v, int p, int x){
    //ans[cnt[col[v]]] -= col[v];
    cnt[col[v]] += x;
    //ans[cnt[col[v]]] += col[v];
    //mxn = max(cnt[col[v]],mxn);
    for(auto u: g[v])
        if(u != p && !big[u])
            add(u, v, x);
}


void dfs(int v, int p, bool keep){
    int mx = -1, bigChild = -1;
    for(auto u : g[v])
       if(u != p && sz[u] > mx)
          mx = sz[u], bigChild = u;
    //run a dfs on small childs and clear them from cnt
    for(auto u : g[v])
        if(u != p && u != bigChild)
            dfs(u, v, 0);  //
    // actual processing of vertex v starts from here
    //mxn = 0;
    if(bigChild != -1)
        dfs(bigChild, v, 1), big[bigChild] = 1;  // bigChild marked as big and not cleared from cnt
    // calculating ans
    add(v, p, 1);
    //res[v] = ans[mxn];
    /** here access the result for each node. if needed then access on add() function
        now cnt[c] is the number of vertices in subtree of vertex v that has color c.
        You can answer the queries easily.
    */
    if(bigChild != -1)
        big[bigChild] = 0;
    if(keep == 0)
        add(v, p, -1);
}



int main()
{
    ll i,j,a,b,ts,cn=0,cas=0,n,m,x,y,sum=0,mn=INT_MAX,u,v;
    //freopen("test.txt","r",stdin);
    cin>>n;
    for(int i = 1; i <= n; i++) {
		cin >> col[i];
	}
    for(i=1;i<n;i++){
        cin>>u>>v;
        g[u].push_back(v);
        g[v].push_back(u);
    }
    getsz(1,0);
    dfs(1,0,1);
    for(int i = 1; i <= n; i++) {
		cout << res[i] << ' ';
	}
    return 0;
}

/**
    Name : 2SAT problem
*/

#include <bits/stdc++.h>
using namespace std;
typedef long long ll;
typedef long double ld;
typedef unsigned long long ull;
typedef pair<int,int> pii;
typedef vector<int> vi;
typedef vector<ll> vll;
#define ntm(i, n) for(int i = 0; i<(n); ++i)
#define rep(i, a, b) for(int i = a; i<(b); ++i)
#define all(x) (x).begin(),(x).end()
#define trav(x, v) for(auto &x : (v))
#define sz(x) (int)(x).size()
#define mp(a, b) make_pair(a, b)
#define X first
#define Y second
#define pb(x) push_back(x)
#define PS(x) cout << " " << (x)
#define P1(x) cout << (x) << endl
#define P2(x,y) cout << (x) << " " << (y) << endl
#define P3(x,y,z) cout << (x) << " " << (y) << " " << (z) << endl
#define SQR(x) ((x)*(x))
#define CUBE(x) ((x)*(x)*(x))
#define PF(x, y) fixed << setprecision(y) << x

using namespace std;

struct TwoSat {
    int N;
    vector<vi> gr;
    vi values; // 0 = false , 1 = true
    TwoSat(int n = 0) : N(n), gr(2*n) {}
    int add_var() {
        gr.emplace_back();
        gr.emplace_back();
        return N++;
    }
    void either(int f, int j) {
        f = max(2*f, -1-2*f);
        j = max(2*j, -1-2*j);
        //P3("either: ", f, j);
        gr[f^1].push_back(j);
        gr[j^1].push_back(f);
    }
    void set_value(int x) { either(x, x); }
    void at_most_one(const vi& li) {
        if (sz(li) <= 1) return;
        int cur = ~li[0];
        rep(i,2,sz(li)) {
            int next = add_var();
            either(cur, ~li[i]);
            either(cur, next);
            either(~li[i], next);
            cur = ~next;
        }
        either(cur, ~li[1]);
    }
    vi val, comp, z; int time = 0;
    int dfs(int i) {
        int low = val[i] = ++time, x; z.push_back(i);
        trav(e, gr[i]) if (!comp[e])
            low = min(low, val[e] ?: dfs(e));
        ++time;
        if (low == val[i]) do {
            x = z.back(); z.pop_back();
            comp[x] = time;
            if (values[x>>1] == -1)
                values[x>>1] = !(x&1);
        } while (x != i);
        return val[i] = low;
    }
    bool solve() {
        values.assign(N, -1);
        val.assign(2*N, 0); comp = val;
        rep(i,0,2*N) if (!comp[i]) dfs(i);
        rep(i,0,N) if (comp[2*i] == comp[2*i+1]) return 0;
        return 1;
    }
};

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

/**
    DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG!
    DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG!
    DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG!
    DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG! DYNAMIC PROGRAMIG!
*/

/**
*   Name : Digit DP
*	Description : generate all number upto the number give into the string st[]
				Given implementation will count how many numbers and those numbers digits sum divisible by n
*   Complexity : O(multiplicatins of parameters)
*	Use : set n = mod number | set ln = length(st)
*/


char st[MX];
int n,ln;
int dp[11][3][101][101];

int bct(int pos,bool lst,int sum,int mult){
    if(pos>=ln){
        if(sum%n==0&&mult%n==0) return 1;
        else return 0;
    }
    if(dp[pos][lst][sum][mult]!=-1)
        return dp[pos][lst][sum][mult];
    int x=0,y;
    for(int i=0; i<=9; i++){
        if(lst){
            y = st[pos]-'0';
            if(i>y) break;
            x += bct(pos+1,i==y,(sum+i)%n,(mult*10+i)%n);
        }
        else{
            x += bct(pos+1,0,(sum+i)%n,(mult*10+i)%n);
        }
    }
    return dp[pos][lst][sum][mult] = x;
}

/**
    FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING!
    FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING!
    FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING!
    FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING! FLOW + MATCHING!
*/

/**
    Name : Bipartite Matching
    Description : Find maximum matching. we can find bipartite graph from left or right array
                  lt = left, rt =right, ed = adjesency matrix
                  left and right graph can contains nodes with same same or number
    Time Complexity : O(EV)
    Some Properties :
        * in a graph if we push both side then number of pair is ceil(Bipartite/2)
        * A graph is bipartite, if it doesn't contain any odd cycle.
        0. V = A u B
        1. Konigs minimax: In a bipartite graph v = t and a = p
        2. Hall's theorem: G has a matching from A to B, iif, for all X (subset of A), |adj(X)| >= |X|
        3. The Marriage Theorem: A perfect match exist iff, it has a matching from A to B and |A| = |B|
        4. Dilworth Theorem: For a finite partially ordered poset S,
        the largest antichain is equal to the minimum chain partition and vice versa.
        Consider a bipartite graph built this way, A = S, B = S and and edge from x(A) to y(B) exists if x <= y.
        Let, c be the minimum chain partition. Then, v + c = n.
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

// clear edge
ed[i].push_back(j) // i may same as j


/**
    Name : Flow Edmonds Karp
    Time Complexity : O(VE^2)
    Use : if u to v have mutiple path dont add path multiple time.
    just add the path cost to capacity matrix
    * reverse edge must be push. if unidirectional graph then push reverse edge
    seting the reverse edge capacity zero
    Source :https://cp-algorithms.com/graph/edmonds_karp.html
*/


int n;
vector<vector<int>> capacity;
vector<vector<int>> adj;

int bfs(int s, int t, vector<int>& parent) {
    fill(parent.begin(), parent.end(), -1);
    parent[s] = -2;
    queue<pair<int, int>> q;
    q.push({s, INF});

    while (!q.empty()) {
        int cur = q.front().first;
        int flow = q.front().second;
        q.pop();

        for (int next : adj[cur]) {
            if (parent[next] == -1 && capacity[cur][next]) {
                parent[next] = cur;
                int new_flow = min(flow, capacity[cur][next]);
                if (next == t)
                    return new_flow;
                q.push({next, new_flow});
            }
        }
    }

    return 0;
}

int maxflow(int s, int t) {
    int flow = 0;
    vector<int> parent(n);
    int new_flow;

    while (new_flow = bfs(s, t, parent)) {
        flow += new_flow;
        int cur = t;
        while (cur != s) {
            int prev = parent[cur];
            capacity[prev][cur] -= new_flow;
            capacity[cur][prev] += new_flow;
            cur = prev;
        }
    }

    return flow;
}

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


/**
    STRING! STRING! STRING! STRING! STRING! STRING! STRING!
    STRING! STRING! STRING! STRING! STRING! STRING! STRING!
    STRING! STRING! STRING! STRING! STRING! STRING! STRING!
    STRING! STRING! STRING! STRING! STRING! STRING! STRING!
*/

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



/// **** Rabin-Karp

#define MD 1000000009
using namespace std;

ll pm=31,powr[MX],hs[MX];

void powCalc(){
    powr[0] = 1;
    for(int i=1;i<MX;i++){
        powr[i] = (powr[i-1]*pm)%MD;
    }
}

vector<int> ans;

void rabin_karp(char txt[],char ptt[]){
    int P = strlen(ptt),T=strlen(txt),i,crr;

    memset(hs,0,sizeof hs);

    for(i=0; i<T; i++){
        hs[i+1] = (hs[i]+(txt[i]-'a'+1)*powr[i])%MD;
    }
    ll hss = 0;
    for(i=0; i<P; i++){
        hss = (hss+(ptt[i]-'a'+1)*powr[i])%MD;
    }

    for(i=0;i+P-1<T;i++){
        crr = (hs[i+P] - hs[i] + MD)%MD;

        if(crr == hss*powr[i]%MD){
            ans.push_back(i+1);
        }
    }
}


/**
    Name : Suffix array
    Description : it gives all suffixes in sorted form on an array
    the array contains starting position of suffixes in sorted order
    Longest common prefix array contains long match with the previous suffix in sorted form
    Complexity : O(n^2Logn)
    Ex :

    0 banana                          5 a
    1 anana     Sort the Suffixes     3 ana
    2 nana      ---------------->     1 anana
    3 ana        alphabetically       0 banana
    4 na                              4 na
    5 a                               2 nana
    So the suffix array for "banana" is {5, 3, 1, 0, 4, 2}

    Uses :

        * number of uniqe substring of a string : n(n+1)/2 − ∑ lcp[i] (summation of all lcp is duplicate string as lcp means longest match with the prev string)
        * longest common substring by sliding window   (add special character and joint all strings then calculate lcp and by sliding window we can determine the longest substring)
        * longest repeated substring : maximum value of lcp is longest repeated substring, number of the maximum number is the longest repeated substring is the longest repeated substring
        * Finding a substring in a string : by binary search on the lcp array
*/

#define MX 10005

int sfa[10009],pos[10009],tmp[10009],lcp[10009],gap=1,ln;
char ss[10009];
bool scmp(int a,int b)
{
    if(pos[a]!=pos[b])return pos[a]<pos[b];
    a+=gap;
    b+=gap;
    return (a<ln&&b<ln)?pos[a]<pos[b]:a>b;

}
void buildsa()
{
    int i,j;
    ln=strlen(ss);
    for(i=0;i<=ln;i++){
        sfa[i]=i;
        pos[i]=ss[i];
    }
    for(gap=1;;gap*=2){
        sort(sfa,sfa+ln,scmp);
        for(i=0;i<ln;i++){
            tmp[i+1]=tmp[i]+scmp(sfa[i],sfa[i+1]);
        }
        for(i=0;i<ln;i++)pos[sfa[i]]=tmp[i];
        if(tmp[i]==ln-1)break;
    }
}
void buildlcp()
{
    int i,j,k;
    for(i=0,k=0;i<ln;i++){
        if(pos[i]==ln-1)continue;
        for(j=sfa[pos[i]+1];ss[i+k]==ss[j+k];)k++;
        lcp[pos[i]]=k;
        if(k)k--;
    }
    lcp[ln-1] = 0;
}


/**
    BACK TRACKING! BACK TRACKING! BACK TRACKING! BACK TRACKING!
    BACK TRACKING! BACK TRACKING! BACK TRACKING! BACK TRACKING!
    BACK TRACKING! BACK TRACKING! BACK TRACKING! BACK TRACKING!
    BACK TRACKING! BACK TRACKING! BACK TRACKING! BACK TRACKING!
*/


/**
*   Name : all permutation generation
*   This code will print all nPr
*   While ara contains n elements
*   a is starring postion of array
*/

int ara[200],mr[200],k;

void permutation(int a,int n){
    int i;
    if(a==n+1){
        for(i=1;i<=n;i++) printf("%c",ara[i]);
        printf("\n");
        return ;
    }
    for(i=1;i<=n;i++) if(!mr[i]){
        mr[i] = 1;
        ara[a] = i;
        permutation(a+1,n);
        mr[i] = 0;
    }
}

/**
*   Name : All combination generation
*   This code will print all nCr
*   While ara is elements of n
*/

int ara[20],sln[20],n,r;

void combinition(int pos,int last){
    int i;
    if(pos>r){
        for(i=1;i<=r;i++){
            if(i==r) printf("%d\n",sln[i]);
            else printf("%d ",sln[i]);
        }
        return;
    }

    for(i=last+1;i<=n;i++){
        sln[pos] = ara[i];
        combinition(pos+1,i);
    }
}

/**
*   Name : m-Coloring problem
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

/**
    COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY!
    COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY!
    COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY!
    COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY! COMPUTATIONAL GEOMETRY!
*/


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



/// *** Jarvis Algorihtm [n*n]


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

double distd(point p1,point p2){
    return sqrt((p1.x-p2.x)*(p1.x-p2.x) + (p1.y-p2.y)*(p1.y-p2.y));
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
double convexHull(int n)
{
	vector<point> hull;
	int l = 0;
	for (int i = 1; i < n; i++)
		if (points[i].x < points[l].x)
			l = i;
	int p = l, q;
	do
	{
		hull.push_back(points[p]);
		q = (p+1)%n;
		for (int i = 0; i < n; i++)
		{
		if (orientation(points[p], points[i], points[q]) == 2)
			q = i;
		}

		p = q;

	} while (p != l);
    double ans = distd(hull[0],hull[hull.size()-1]);
	for (int i = 0; i<hull.size()-1; i++){
        ans += distd(hull[i],hull[i+1]);
        //cout<<hull[i].x<<" "<<hull[i].y<<endl;
	}
	return ans;
}

/**
    BINARY EXPONENTIATION! BINARY EXPONENTIATION!
    BINARY EXPONENTIATION! BINARY EXPONENTIATION!
    BINARY EXPONENTIATION! BINARY EXPONENTIATION!
    BINARY EXPONENTIATION! BINARY EXPONENTIATION!
*/


/**
    Name : Bigmod
    Description : calculate (n^r)%MD
    Complexity : O(log r)
*/

#define MD 1000007

ll bigmod(ll n, ll r){
    if(r==0) return 1;
    if(r==1) return n%MD;
    if(r%2==1) return (bigmod(n,r-1)*n)%MD;
    else{
        ll x = bigmod(n,r/2);
        return (x*x)%MD;
    }
}

/**
    Name : bigmod without recursion
    Description : calculate (n^r)%MD
    Complexity : O(log r)
*/


#define MD 1000007

ll bigmod(ll n, ll r){
    ll x = 1;
    while (r > 0){
        if (r&1) x = (x*n)%MD;
        n = (n*n)%MD;
        r >>= 1;
    }
    return x % MD;
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

/**
    COMBINATORICS! COMBINATORICS! COMBINATORICS!
    COMBINATORICS! COMBINATORICS! COMBINATORICS!
    COMBINATORICS! COMBINATORICS! COMBINATORICS!
    COMBINATORICS! COMBINATORICS! COMBINATORICS!
*/

/// *** ncr [n*r] [but it is dp]

ll cr[902][902];
ll ncr(int n,int r){
    if(n==r)return 1;
    if(r==1)return n;
    if(cr[n][r])return cr[n][r];
    return cr[n][r] = ncr(n-1,r)+ncr(n-1,r-1);
}

/**
    Name : Inclusion-Exclusion
    Description : compute the union of sets
    Time Complexity : O(2^n)
    Space : O(n)
    Source : collected form hacker earth
*/

int n; // the number of sets in the set A
int result = 0; //final result, the cardinality of sum of all subsets of A
for(int b = 0; b < (1 << n); ++b){
     vector<int> indices;
     for(int k = 0; k < n; ++k){
          if(b & (1 << k)){
                // we could work with this values
                indices.push_back(k);
          }
     }
     int cardinality = intersectionCardinality(indices); // intersections
     if(indices.size() % 2 == 1) result += cardinality;
     else result -= cardinality;
}
cout << result << endl; //printing the final result


/**
    Name : Catalan Number
    n’th catalan number : Cn = (1/(n+1))(2n n) = (2n n) - (2n n+1)
                =(2n)!/((2n-n)!*n!*(n+1))
                =(2n)!/((n)!*(n+1)!)

    C(0) = C(1) = 1
    C(4) = C0*C3 + C1*C2 + C2*C1 + C3*C0
    C(n) =C(0)*(n-1) +C(1)*C(n-2) +.......... +C(i)*C(n-i-1)+........+C(n-1)*C(0)

    * the number of ways a polygon with n+2 sides can be cut into triangles
    * the number of ways to use n rectangles to tile a stairstep shape (1, 2, ..., n−1, n).
    * Cn counts the number of expressions containing n pairs of parentheses which are correctly matched
    * Cn is the number of different ways (n + 1) factors can be completely parenthesized  catalan(3) = 5 , ex :  ((ab)c)d,   (a(bc))d ,   (ab)(cd) ,   a((bc)d) ,   a(b(cd))
    * Count the number of possible Binary Search Trees with n keys
    * Total number of full binary tree : catalan(n) * factorial(n)
    * the number of paths of length 2n through an n-by-n grid that do not rise above the main diagonal
    * Cn is the number of ways to form a “mountain ranges” with n upstrokes and n down-strokes that all stay above the original line.The mountain range interpretation is that the mountains will never go below the horizontal.
    * Cn is the number of ways that the vertices of a convex 2n-gon can be paired so that the line segments joining paired vertices do not intersect. This is precisely the condition that guarantees that the paired edges can be identified (sewn together) to form a closed surface of genus zero (a topological 2-sphere).

*/

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

/**
    NUMBER THEORY! NUMBER THEORY! NUMBER THEORY! NUMBER THEORY!
    NUMBER THEORY! NUMBER THEORY! NUMBER THEORY! NUMBER THEORY!
    NUMBER THEORY! NUMBER THEORY! NUMBER THEORY! NUMBER THEORY!
    NUMBER THEORY! NUMBER THEORY! NUMBER THEORY! NUMBER THEORY!
*/

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

/**
    Name : Sieve O(n)
    Description : pr stores all primes and lp stores lowest prime factors
    Complexity : O(n)
*/

const int N = 10000000;
int lp[N+1];
vector<int> pr;

void sieve(){
    int i,j;
    for (i=2; i<=N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            pr.push_back (i);
        }
        for (j=0; j<(int)pr.size() && pr[j]<=lp[i] && i*pr[j]<=N; ++j)
            lp[i*pr[j]] = pr[j];
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



/// *** Is prime [sqrt(n)]

int isprime(int n){
    if(n==2)return 1;
    if(!(n%2)||n<2)return 0;
    int i,sq=sqrt(n)+2;
    for(i=3; i<sq; i+=2)if(!(n%i))return 0;
    return 1;
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

/**
*	Description : normal prime factorization
*	Complexity : O(sqrt(n))
*	Use : store prime up to squreroot(n) by efficient sieve, in = number of prime
*/


void primeFactorization(ll n){
	ll i,cn=0;
	for(i=1; i <= in && (ll)pr[i]*(ll)pr[i] <= n; i++){
        cn = 0;
    	while(n%pr[i] == 0){
        	n /= pr[i];
        	cn++;
        }
        // cn is the prime power
    }
    if(n>1) // then n is a prime
}

/**
    Description : prime factorization
    Complexity : O(log n)
    Limitation : memory complexity O(n) so it works up to 10^7
*/

const int N = 10000000;
int lp[N+1];
vector<int> pr;

void sieve(){
    int i,j;
    for (i=2; i<=N; ++i) {
        if (lp[i] == 0) {
            lp[i] = i;
            pr.push_back (i);
        }
        for (j=0; j<(int)pr.size() && pr[j]<=lp[i] && i*pr[j]<=N; ++j)
            lp[i*pr[j]] = pr[j];
    }
}

vector<int> prime_fact(int n){
    vector<int> pf;
    while (n != 1){
        pf.push_back(lp[n]);
        n = n / lp[n];
    }
    return pf;
}

/**
*   Name : NOD
*   Description : Therefore if the prime factorization of n is p1^e1⋅pi2^e2⋯pk^ek,
*		        where pi are distinct prime numbers, then the number of divisors is: d(n)=(e1+1)⋅(e2+1)⋯(ek+1)
*	Complexity : O(sqrt(n)) , Squre root porjnto jotogula prime ase
*	Use : store prime up to squreroot(n) by efficient sieve, in = number of prime
*/



ll NOD(ll n){
	ll ans=1,x,count;
	for(int i=0; (ll)prime[i] * (ll)prime[i]<=n && i<=in ;i++){
		int( n%prime[i] == 0){
			x=prime[i];
			count =1;
			while(n%x==0) {
				n/=x;
				count++;
			}
			ans*= count;
		}
	}
	if(n>1) ans*=2;
	return ans;
}

/**
*   Name : SOD
*   Description : Find sum of divisor of a given number
*				σ(n)=(p1^(e1+1)-1)/(p1−1) * (p2^(e2+1)−1)/(p2−1)***(pk^(ek+1)−1)(pk−1)
*	Complexity : O(sqrt(n)) , Squre root porjnto jotogula prime ase
*	Use : store prime up to squreroot(n) by efficient sieve, in = number of prime
*/


ll SOD(ll n){
    ll i,cnt,sum=1,tmp,B;
    for(i=1;i<=in&&pr[i]*pr[i]<=n;i++){
            cnt=0;
            while(n%pr[i]==0){
                cnt++;
                n /= pr[i];
            }
            tmp = pr[i];
            if(cnt){
                sum *= ((bgM(tmp,cnt+1)-1+MD)%MD)*(modi(tmp-1))%MD;
            }
            sum = sum%MD;
        }
        tmp = n;
        if(n>1){
            sum *= ((bgM(tmp,2)-1+MD)%MD)*(modi(tmp-1))%MD;
            sum = sum%MD;
        }
    return sum;
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

/**
*   Name : Totient Function
*   Description : Find number of co-prime of a number less than that number
*	Complexity : O(sqrt(n))
*   source : https://cp-algorithms.com/algebra/phi-function.html
*/

int phi(int n) {
    int result = n;
    for (int i = 2; i * i <= n; i++) {
        if(n % i == 0) {
            while(n % i == 0)
                n /= i;
            result -= result / i;
        }
    }
    if(n > 1)
        result -= result / n;
    return result;
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


/**
    Name : Chinese remainder theorem

    Description : lets think about an equation:   x % num[i] = rem[i];
    Given num[] and rem[] list, you have to calculate an x ,
    so that for every index i, corrsponding num[i] and rem[i] should satisfy the equation;
*/


ll num[15], rem[15], n;

pll extnuc(ll a,ll b)
{
    if(b==0)return pll(1,0);
    pll d=extnuc(b,a%b);
    return pll(d.second,d.first-d.second*(a/b));
}

ll modi(ll n, ll M)
{
    pll d=extnuc(n,M);
    return ((d.first%M)+M)%M;
}

ll chinese_remainder(int k)
{
    ll prod = 1;

    for(int i = 0; i < k; ++i)
        prod *= num[i];

    ll res = 0;

    for(int i = 0; i < k; ++i) {
        ll pp = prod / num[i];

        int inv = modi(pp, num[i]);
        if(inv<0) inv += num[i];

        res += rem[i] * pp * inv;
        res %= prod;
    }
    return res;
}

int main()
{
    scanf("%lld", &n);
    for(i=0;i<n;i++)scanf("%lld %lld", &num[i], &rem[i]);
    printf("Case %d: %lld\n", ++tc, chinese_remainder(n));
    return 0;
}

/**
    Name : Stern Brocot tree
    here n = numerator d = denominator
    prn = precissions numerator prd =
    O(log n)
    stern brocot tree can generate any fractional number
*/

void stern_brocot(ll n, ll d,ll prn,ll prd){
    ll ln=0,ld=1,rn=1,rd=0,mn,md;
    while(1){
        mn = ln+rn;
        md = ld+rd;
        if(n*md<d*mn){
            rn = mn;
            rd = md;
            /** condition for precissition and breaking condition
                it could be different for different problem
                here we are finding a/b which is near to c/d and a/b<c/d
                and b is as small as possible
            */
            if(prd*(d*mn-n*md)<=d*md*prn&&mn){
                x = mn;
                y = md;
                break;
            }
        }
        else{
            ln = mn;
            ld = md;
        }
    }
}




