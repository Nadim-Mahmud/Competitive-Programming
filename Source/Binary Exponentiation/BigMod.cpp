/**
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
