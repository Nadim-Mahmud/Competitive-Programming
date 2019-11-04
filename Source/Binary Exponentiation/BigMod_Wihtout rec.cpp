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
