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
